"""Run per-cluster tree inference as concurrent jobs within one thread budget.

IQ-TREE parallelises over alignment site patterns, and cluster alignments are
short: a 167-taxon cluster of the Nipah binder set has 303 patterns. Giving one
IQ-TREE run many threads leaves each thread a handful of patterns and the run
spends its time synchronising - 30 iterations on that cluster took 50 s at 8
threads but 82 s at 32. Running several clusters at once, each with a modest
thread count, uses the same budget far better.

Each job is a separate ``viral-phylo tree`` process writing to its own log, so
concurrent tool output never interleaves and jobs share no Python state.
"""

import os
import signal
import subprocess
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import Any, Dict, List, Optional, Tuple


def _as_count(value) -> Optional[int]:
    try:
        n = int(value)
    except (TypeError, ValueError):
        return None
    return n if n > 0 else None


def plan_tree_jobs(threads, tree_threads="8", tree_jobs=None) -> Tuple[int, str, Optional[str]]:
    """Split the run's thread cap across concurrent tree jobs.

    Returns ``(jobs, threads_per_job, note)``. ``jobs x threads_per_job`` never
    exceeds the cap. With no cap (``--threads AUTO``) trees run one at a time with
    the uncapped setting, as before. ``note`` explains any adjustment made.
    """
    cap = _as_count(threads)
    if cap is None:
        return 1, str(threads), None
    per_job = min(_as_count(tree_threads) or 8, cap)
    fit = max(1, cap // per_job)
    requested = _as_count(tree_jobs)
    if requested is None:
        return fit, str(per_job), None
    if requested > fit:
        return fit, str(per_job), (f"--tree-jobs {requested} x {per_job} threads exceeds --threads {cap}; "
                                   f"running {fit} job(s) at a time")
    return requested, str(per_job), None


def _tail(path: str, lines: int = 15) -> str:
    try:
        with open(path, "r", encoding="utf-8", errors="replace") as handle:
            return "".join(handle.readlines()[-lines:])
    except OSError:
        return ""


def _signal_group(proc: subprocess.Popen, sig: int) -> None:
    try:
        os.killpg(proc.pid, sig)
    except (ProcessLookupError, PermissionError):
        pass


def run_tree_jobs(jobs: List[Dict[str, Any]], max_parallel: int,
                  env: Optional[Dict[str, str]] = None) -> List[Dict[str, Any]]:
    """Run tree jobs, at most ``max_parallel`` at a time, largest first.

    Each job is ``{"name", "taxa", "cmd", "log"}``. Returns one result per job with
    ``returncode``, ``seconds`` and, on failure, the tail of its log. If this
    process is interrupted (Ctrl-C, SIGTERM via SystemExit), every running job is
    terminated before the exception propagates, so no tree is left running.
    """
    # Longest-first scheduling: the biggest clusters dominate the stage, so start
    # them immediately and let small ones fill the remaining slots.
    ordered = sorted(jobs, key=lambda j: -(j.get("taxa") or 0))
    total = len(ordered)
    lock = threading.Lock()
    running: Dict[str, subprocess.Popen] = {}
    stopping = threading.Event()
    done = [0]

    def _say(msg: str) -> None:
        with lock:
            print(msg, flush=True)

    def _run(job: Dict[str, Any]) -> Dict[str, Any]:
        if stopping.is_set():
            return {**job, "returncode": None, "seconds": 0.0, "error": "not started: run interrupted"}
        os.makedirs(os.path.dirname(job["log"]) or ".", exist_ok=True)
        start = time.time()
        with open(job["log"], "w", encoding="utf-8") as log:
            log.write("$ " + " ".join(job["cmd"]) + "\n\n")
            log.flush()
            # Own process group, so stopping a job also stops the IQ-TREE /
            # VeryFastTree it launched rather than orphaning it.
            proc = subprocess.Popen(job["cmd"], stdout=log, stderr=subprocess.STDOUT, env=env,
                                    start_new_session=True)
            with lock:
                running[job["name"]] = proc
                n_running = len(running)
            _say(f"[Trees] started  {job['name']} ({job['taxa']} taxa) - {n_running} running; log: {job['log']}")
            rc = proc.wait()
            with lock:
                running.pop(job["name"], None)
        seconds = time.time() - start
        result = {**job, "returncode": rc, "seconds": round(seconds, 1)}
        with lock:
            done[0] += 1
            finished = done[0]
        if rc == 0:
            _say(f"[Trees] finished {job['name']} in {seconds / 60:.1f} min ({finished}/{total})")
        else:
            result["error"] = f"exit status {rc}"
            result["log_tail"] = _tail(job["log"])
            _say(f"[Trees] [!] FAILED {job['name']} (exit {rc}) after {seconds / 60:.1f} min "
                 f"({finished}/{total}); see {job['log']}")
        return result

    results: List[Dict[str, Any]] = []
    pool = ThreadPoolExecutor(max_workers=max(1, max_parallel))
    try:
        futures = [pool.submit(_run, job) for job in ordered]
        for future in as_completed(futures):
            results.append(future.result())
    except BaseException:
        stopping.set()
        with lock:
            live = list(running.values())
        for proc in live:
            _signal_group(proc, signal.SIGTERM)
        for proc in live:
            try:
                proc.wait(timeout=30)
            except subprocess.TimeoutExpired:
                _signal_group(proc, signal.SIGKILL)
        raise
    finally:
        pool.shutdown(wait=True, cancel_futures=True)
    order = {job["name"]: i for i, job in enumerate(ordered)}
    return sorted(results, key=lambda r: order[r["name"]])
