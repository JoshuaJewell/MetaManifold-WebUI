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

case "${1:-}" in
    run)
        # The lock is inherited by the workload: even if its wrapper dies, a
        # surviving child prevents another command from starting in this attempt.
        exec 9>"$control/run.lock"
        flock -n 9 || exit 0
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
                # If the wrapper dies before publishing an exit marker, stop
                # its children too. Otherwise this watcher would keep the
                # execution lock forever after an unrecorded completion.
                if [ -e "$control/cancel" ] || ! kill -0 "$worker_pid" 2>/dev/null; then
                    kill -KILL -- "-$worker_pid"
                    exit 1
                fi
                [ ! -e "$control/exit" ] || exit 0
                sleep 0.1
            done
        ) &
        watcher=$!
        bash "$control/command.sh"
        result=$?
        printf '%s\n' "$result" > "$control/exit.tmp"
        mv -- "$control/exit.tmp" "$control/exit"
        wait "$watcher"
        exit "$result"
        ;;
    start)
        [ -e "$control/ready" ] || exit 1
        # Retransmitting this request after losing its reply is safe. Only one
        # worker can acquire run.lock; a terminal attempt cannot launch again.
        nohup setsid bash "$control/worker.sh" run >> "$control/log" 2>&1 < /dev/null &
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
