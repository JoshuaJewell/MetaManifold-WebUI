#!/usr/bin/env bash
# SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>
# SPDX-License-Identifier: MPL-2.0

# Private protocol for remote_exec.jl. The remote host must be Linux, with bash,
# flock, setsid, nohup and sha256sum. No SSH channel owns the workload's lifetime.
set -u
stage_dir=$(cd -- "$(dirname -- "$0")/.." && pwd -P) || exit 1
cd -- "$stage_dir" || exit 1
control=.metamanifold
boot_id() { cat /proc/sys/kernel/random/boot_id; }
birth() { awk '{print $22}' "/proc/$1/stat" 2>/dev/null; }
idle() { flock -n "$control/run.lock" true; }
# Publish the wrapper's exit status atomically. Nothing can be renamed onto a
# full or read-only directory, so a failure here is reported rather than
# ignored: an absent marker would keep the controller polling for a terminal
# state that never arrives, and an empty one would read as an unknown status.
publish_exit() {
    { printf '%s\n' "$1" > "$control/exit.tmp"; } 2>/dev/null &&
        mv -- "$control/exit.tmp" "$control/exit" 2>/dev/null
}
# The retryable failure state: `exit 1` when that marker can still be written,
# and otherwise the `interrupted` state, which a status probe reports when the
# launch marker is present and no lock is held. The launch marker is therefore
# deliberately kept; removing it would let a probe read the attempt as `ready`
# and relaunch a command whose checkpoint is the only record of what ran.
interrupted() {
    publish_exit 1 && return 0
    [ -e "$control/pid" ] || publish_exit 1
    return 1
}

case "${1:-}" in
    run)
        # The lock is inherited by the workload: even if its wrapper dies, a
        # surviving child prevents another command from starting in this attempt.
        exec 9>"$control/run.lock"
        # A status/idle probe takes this lock briefly too. Do not mistake that
        # transient owner for an already running workload and lose the launch.
        flock -w 1 9 || exit 0
        exec 8>"$control/control.lock"
        flock -x 8 || exit 1
        [ ! -e "$control/cancel" ] && [ ! -e "$control/exit" ] &&
            [ -e "$control/ready" ] || exit 0
        printf '%s %s %s\n' "$$" "$(boot_id)" "$(birth "$$")" > "$control/pid.tmp"
        mv -- "$control/pid.tmp" "$control/pid" || exit 1
        flock -u 8
        exec 8>&-
        # A member of this process group observes cancellation and kills its
        # own group. The reconnecting controller never signals a stored PID:
        # checking /proc and then signalling would still have a PID-reuse race.
        worker_pid=$$
        (
            while :; do
                # If the wrapper dies before it can publish an exit marker,
                # treat the attempt as interrupted: stop its children, mark it
                # failed, and record the outcome so the retained lock does not
                # hang the controller. An unrecorded completion would otherwise
                # keep this watcher -- and the execution lock it inherits -- for
                # the life of the process group.
                if [ -e "$control/cancel" ]; then
                    kill -KILL -- "-$worker_pid" 2>/dev/null
                    exit 1
                fi
                if ! kill -0 "$worker_pid" 2>/dev/null; then
                    kill -KILL -- "-$worker_pid" 2>/dev/null
                    if ! publish_exit 1; then
                        interrupted
                    fi
                    exit 1
                fi
                [ ! -e "$control/exit" ] || exit 0
                sleep 0.1
            done
        ) &
        watcher=$!
        bash "$control/command.sh"
        result=$?
        # A full or read-only filesystem must not leave the controller polling a
        # barrier that can never appear. `interrupted` is the retryable state.
        if ! publish_exit "$result"; then
            interrupted
            exit 1
        fi
        wait "$watcher"
        exit "$result"
        ;;
    start)
        [ -e "$control/ready" ] || exit 1
        # Retransmitting this request after losing its reply is safe. Only one
        # worker can acquire run.lock; a terminal attempt cannot launch again.
        # Recorded here, not by the detached child: if its log cannot be opened
        # on a full disk it never runs, so waiting for a marker only it can
        # write would report a failed launch as a success.
        nohup setsid bash "$control/worker.sh" run >> "$control/log" 2>&1 < /dev/null &
        child=$!
        # Acknowledging only the shell fork lets an immediate status probe win
        # run.lock before the worker. Wait for a durable launch marker before
        # the caller polls. Startup failure is bounded and retryable, not an
        # apparently successful start followed by an endless ready state. A
        # child that has already exited without a marker (a missing log, an
        # open failure) is a failed launch, not one still starting up.
        for ((n=0; n<100; n++)); do
            if [ -e "$control/pid" ] || [ -e "$control/exit" ] || [ -e "$control/cancel" ]; then
                exit 0
            fi
            kill -0 "$child" 2>/dev/null || break
            sleep 0.1
        done
        echo 'Remote worker did not acknowledge startup; retry the stage' >&2
        exit 1
        ;;
    status)
        if [ -e "$control/cancel" ]; then
            echo cancelled
        elif [ -f "$control/exit" ]; then
            printf 'exit '; cat "$control/exit"
        elif ! idle; then
            echo running
        elif [ -e "$control/pid" ]; then
            # A reboot or a killed wrapper: do not guess that outputs are valid.
            echo interrupted
        elif [ -e "$control/ready" ]; then
            echo ready
        else
            echo preparing
        fi
        ;;
    retry)
        exec 9>"$control/run.lock"
        # The terminal marker is visible just before the old wrapper and its
        # watcher release their inherited descriptors. Wait for that handoff;
        # observing an exit marker does not imply the lock is already free.
        flock -w 10 9 || exit 1
        exec 8>"$control/control.lock"
        flock -x 8 || exit 1
        [ ! -e "$control/cancel" ] || exit 1
        # Only an explicitly retried failure loses its exit marker. Checkpoints
        # and logs stay in place; a successful attempt is never rerun.
        if [ -f "$control/exit" ]; then
            [ "$(cat "$control/exit")" != 0 ] || exit 1
            rm -- "$control/exit"
        fi
        # The launch marker outlives the wrapper and its watcher by design, so a
        # crash between the wrapper's exit and its retried relaunch would make a
        # later status probe read a stale `interrupted` even though this attempt
        # is running. Clear it under the locks a probe needs.
        rm -f -- "$control/pid"
        ;;
    cancel)
        exec 8>"$control/control.lock"
        flock -x 8 || exit 1
        # A tombstone fences a worker launched just before cancellation. It must
        # hold control.lock and observe this before publishing a PID or starting.
        touch "$control/cancel"
        flock -u 8
        exec 8>&-
        for ((n=0; n<100; n++)); do
            if idle; then exit 0; fi
            sleep 0.1
        done
        echo 'Remote process has not stopped; retaining its directory' >&2
        exit 1
        ;;
    cleanup)
        # Fetch or cancel must have finished first. Never unlink a live lock and
        # create a second inode on which a new worker could acquire a new lock.
        exec 9>"$control/run.lock"
        flock -w 10 9 || exit 1
        exec 8>"$control/control.lock"
        flock -x 8 || exit 1
        [ -e "$control/exit" ] || [ -e "$control/cancel" ] || exit 1
        case "$(basename -- "$stage_dir")" in mm_*) ;; *) exit 1 ;; esac
        cd .. || exit 1
        rm -rf -- "$stage_dir"
        ;;
    log)
        [[ "${2:-}" =~ ^[0-9]+$ ]] || exit 1
        [ -f "$control/log" ] || exit 0
        # Bounded reads let cancellation interrupt even a very noisy command.
        dd if="$control/log" bs=65536 skip="$2" count=65536 iflag=skip_bytes,count_bytes status=none
        ;;
    *) echo 'Unknown remote execution operation' >&2; exit 2 ;;
esac
