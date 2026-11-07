<!-- SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk> -->
<!-- SPDX-License-Identifier: CC-BY-SA-4.0 -->

# Recovering remote stages

All eight SSH stages use the same runner: `learn_errors`, `denoise`,
`chimera_removal`, `assign_taxonomy`, `phylogeny_align`, `phylogeny_tree`,
`phylogeny_add`, and `phylogeny_place`.

After restarting MetaManifold, run the interrupted workflow again. Reference
trees and placements whose saved status says running without a live local job
are marked interrupted, enabling Build again. The stage reconnects to its
remote command, or downloads its completed outputs, without uploading again.
This is recovery on an explicit run request, not automatic job submission at
server startup. The local queue remains an in-memory queue.

The runner requires Unix locally (including WSL), and Linux on the remote host,
with Bash, OpenSSH, util-linux `flock`/`setsid`, and coreutils `nohup`, `dd` and
`sha256sum`. The configured staging filesystem must implement advisory locks
and atomic rename correctly. It must retain the stage directory across the
interruption. No scheduler or new Julia package is required.

## What is retained

`<stage log>.remote.json` records one workflow's attempt. Keep it with the
workflow when restarting. Its sibling `.lock` is an advisory OS lock: the kernel
releases ownership when a caller dies, and the file itself is deliberately not
unlinked. Another simultaneous caller for the same stage log is refused.
Different workflows receive separate attempt directories even for equal data.

An attempt's identity includes the input SHA-256 checksums and relative names,
the command with a fixed staging placeholder, stage name, target host and root,
thread setting, output names, worker version, and any externally stored database
checksums supplied by the caller. Uploads are checked against those checksums
before launch. A changed identity creates a separate attempt and archives the
old journal as `.superseded-<id>`; it never borrows the old command's outputs.
An already running superseded command may finish independently.
Phylogeny outputs are retained only when the journal identity matches and its
retrieval receipt verifies the current local files. Otherwise they are removed
before contacting SSH, so a failed new attempt cannot expose an older tree.

The remote `.metamanifold` directory contains the command, worker, log, process
identity, locks and exit marker. The worker is detached with `nohup setsid` and
has no SSH standard streams. Repeating a launch request after losing its reply
cannot start another workload while the attempt's execution lock is held.
Successful exit markers also prevent relaunch. The log is read in bounded
chunks on each connection.
Control scripts are uploaded to unique temporary names, checksummed, and renamed
into place. A partial upload is never executed as the installed status worker.

The exit status is published by rename. Outputs are downloaded to unique
temporary local files and checked against remote SHA-256 checksums. All required
downloads must verify before any output is published. A transfer failure keeps
the remote directory and journal. After retrieval a receipt records the local
output hashes, then remote cleanup runs. The receipt survives cleanup so that a
restart before the caller commits its checkpoint can recover the verified files.
Receipt recovery does not repeat preparation checks against staging inputs that
cleanup has removed. A separate callback restores caller metadata when needed.
The publication of multiple files is not a single filesystem transaction; a
restart before the receipt is committed fetches them again from the retained
remote directory.

## Cancellation, failures and retention

Cancelling a local job writes a durable cancellation request, fences remote
launch with a tombstone, and has a watchdog inside the worker's own process
group kill that group. It waits for the execution lock to be released before
removing staging. The controller never signals a PID read from a previous run:
even a process-start-time check followed by a kill would have a PID-reuse race.
The boot ID and process start time are retained only for diagnosis.
Cancellation handling covers recovery preparation too, including input hashing
and command construction. It reads the saved journal to stop the original
attempt even if the new run has not yet reconnected or its settings changed.
If the wrapper dies without recording completion, its watchdog stops the
remaining group, releasing the lock so the retained checkpoint can be retried.
A cancellation that cannot reach the host is reported as unacknowledged
and retained locally. Running that stage again first retries cancellation,
including when settings now select different content or another host.
Network partitions cannot guarantee immediate remote cancellation.
To abandon an interrupted computation, run the workflow again and cancel its
new local job before deleting the workflow and its recovery journals.

A command failure keeps its files. A subsequent explicit retry clears the
failure marker under the execution lock and runs in the same directory. Retry
waits for the terminating wrapper to release that lock after publishing its
exit marker. A
remote reboot likewise permits retry once no command holds that lock. IQ-TREE
remote commands omit `-redo`, allowing its own checkpoint recovery. DADA2 and
other tools without internal checkpointing restart their command after a remote
process failure; they still retain their uploads. A local interruption while the
remote process stays alive does not restart these commands.

Successful retrieval and acknowledged cancellation remove their attempt's
directory. The maintenance sweep removes terminal directories older than seven
days only when their execution lock is free. It never deletes a running command
because its directory is old. Sweeping occurs during later runner activity;
without further activity, administrators can remove abandoned terminal attempts
using the same worker's `cleanup` operation. Superseded journals record the old
identity and directory for diagnosis. Do not delete a live attempt by hand.

## Verification and its boundary

Run the regression suite with the pinned toolchain:

```sh
julia --project=. -t 2 test/runtests.jl test_remote_stages.jl test_phylogeny.jl test_jobs.jl
julia --project=. -t 2 test/remote_protocol.jl
julia --project=. -t 2 test/remote_ssh.jl
bash proofs/tests/axiom-audit.sh
cd proofs/agda && agda MetaManifold/All.agda
```

The remote protocol tests run the actual Bash worker and OS process groups over
a local transport. They cover a real Julia-process SIGKILL and fresh-process
recovery, connection failure, completed work, partial downloads, process-group
cancellation, distinct concurrent workflows, changed content/settings, and tool
checkpoint retry. The standalone `remote_protocol.jl` entry point loads the same
production runner and logger without R or the web server. `remote_ssh.jl` creates
a dedicated loopback OpenSSH server, keys and staging root, then kills the Julia
caller and closes its SSH connection. It checks reattachment to live and completed
commands without uploads, and cancellation through SSH. It requires Linux,
`openssh-server`, and the system privilege-separation directory `/run/sshd`.
CI runs this fixture as well as the full unit and integration suite.

`MetaManifold.RemoteExecution` in the existing Agda proof suite proves properties
of an explicit transition system, by induction over arbitrary finite traces:

| Theorem | Protocol property | Implementation boundary |
|---|---|---|
| `one-live-command` | At most one live workload per attempt. | Exclusive `run.lock`; inherited by the synchronous workload. |
| `trace-key` | An attempt never changes owner/content identity. | Exact journal-identity comparison and a new directory for a changed identity. |
| `different-attempts-stay-different` | Independent attempt identities remain distinct through any traces. | Fresh UUID directory names and distinct workflow journals. |
| `cleanup-needs-no-live-command` | Cleanup cannot enter the removed phase with a live workload. | Cleanup must acquire the execution lock. |
| `fetch-needs-success` | A new fetched state requires successful exit. | Exit status zero and all required output checksums verified. |
| `different-content-cannot-attach` | Unequal keys cannot reattach. | Checked identity equality. |
| `disconnect-preserves-state` | Transport loss does not transition the remote state. | Detached session with redirected standard streams. |
| `duplicate-launch-preserves-count` | A duplicate launch while running adds no command. | Nonblocking execution lock in the worker. |

The model distinguishes a cancellation request from stop acknowledgement. It
does not assume that requesting cancellation has already killed the command.
Repeated start requests and disconnections are stuttering transitions.

These are machine-checked protocol theorems, not a proof of equivalence between
Agda, Julia and Bash. OS locking, process groups, atomic rename, truthful exit
status, synchronous commands, filesystem durability, SSH integrity, SHA-256
collision resistance and fresh UUID allocation are trusted implementation
assumptions. A remote command must wait for its work rather than launch an
untracked daemon. Keep remote executables and their runtime environment unchanged
through an attempt and its recovery; this protocol does not snapshot installed
tools or libraries. Host failure can destroy data or a tool checkpoint; indefinite
network failure prevents liveness. Hash equality is not a mathematical proof of
byte equality for all possible inputs. Tests and these theorems therefore do not
justify a claim of unconditional, 100% end-to-end correctness.
