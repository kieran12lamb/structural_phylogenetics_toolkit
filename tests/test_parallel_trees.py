"""Concurrent per-cluster tree jobs within the pipeline's thread cap."""

import os
import signal
import sys
import tempfile
import textwrap
import time
import unittest

from viral_phylo.parallel_trees import plan_tree_jobs, run_tree_jobs


class TestPlanTreeJobs(unittest.TestCase):
    def test_cap_is_divided_into_jobs(self):
        self.assertEqual(plan_tree_jobs("32", "8"), (4, "8", None))

    def test_per_job_threads_never_exceed_the_cap(self):
        self.assertEqual(plan_tree_jobs("4", "8"), (1, "4", None))

    def test_requested_jobs_are_clamped_to_what_fits(self):
        jobs, per_job, note = plan_tree_jobs("32", "8", 10)
        self.assertEqual((jobs, per_job), (4, "8"))
        self.assertIn("exceeds --threads 32", note)

    def test_fewer_jobs_may_be_requested(self):
        self.assertEqual(plan_tree_jobs("32", "8", 2), (2, "8", None))

    def test_auto_runs_one_uncapped_job(self):
        self.assertEqual(plan_tree_jobs("AUTO", "8"), (1, "AUTO", None))


def _job(tmp, name, taxa, code):
    return {"name": name, "taxa": taxa, "log": os.path.join(tmp, name, "tree_job.log"),
            "cmd": [sys.executable, "-c", textwrap.dedent(code)]}


class TestRunTreeJobs(unittest.TestCase):
    def test_jobs_run_concurrently_largest_first(self):
        with tempfile.TemporaryDirectory() as tmp:
            stamp = os.path.join(tmp, "starts")
            code = f"import time; open({stamp!r}, 'a').write('{{}} %f\\n' % time.time()); time.sleep(1.0)"
            jobs = [_job(tmp, f"c{n}", n, code.replace("{}", f"c{n}")) for n in (5, 50, 20, 10)]
            t0 = time.time()
            results = run_tree_jobs(jobs, max_parallel=2)
            elapsed = time.time() - t0
            self.assertTrue(all(r["returncode"] == 0 for r in results))
            # 4 one-second jobs, 2 at a time: about 2 s, not 4.
            self.assertLess(elapsed, 3.5)
            starts = [line.split()[0] for line in open(stamp)]
            self.assertEqual(set(starts[:2]), {"c50", "c20"}, "the two largest clusters must start first")

    def test_failure_is_reported_with_its_log(self):
        with tempfile.TemporaryDirectory() as tmp:
            jobs = [_job(tmp, "ok", 10, "print('fine')"),
                    _job(tmp, "bad", 20, "import sys; print('IQ-TREE exploded'); sys.exit(2)")]
            results = {r["name"]: r for r in run_tree_jobs(jobs, max_parallel=2)}
            self.assertEqual(results["ok"]["returncode"], 0)
            self.assertEqual(results["bad"]["returncode"], 2)
            self.assertIn("IQ-TREE exploded", results["bad"]["log_tail"])
            self.assertTrue(open(results["bad"]["log"]).read().startswith("$ "), "log starts with the command")

    @unittest.skipUnless(hasattr(signal, "setitimer"), "needs POSIX timers")
    def test_interrupt_stops_jobs_and_their_children(self):
        """A killed pipeline must not leave IQ-TREE running: the job's own child dies too."""
        with tempfile.TemporaryDirectory() as tmp:
            pidfile = os.path.join(tmp, "grandchild.pid")
            code = f"""
                import subprocess, sys, time
                p = subprocess.Popen([sys.executable, "-c", "import time; time.sleep(60)"])
                open({pidfile!r}, "w").write(str(p.pid))
                time.sleep(60)
            """
            jobs = [_job(tmp, "long", 100, code)]

            def _interrupt(signum, frame):
                raise KeyboardInterrupt

            previous = signal.signal(signal.SIGALRM, _interrupt)
            try:
                signal.setitimer(signal.ITIMER_REAL, 1.5)
                with self.assertRaises(KeyboardInterrupt):
                    run_tree_jobs(jobs, max_parallel=1)
            finally:
                signal.setitimer(signal.ITIMER_REAL, 0)
                signal.signal(signal.SIGALRM, previous)
            grandchild = int(open(pidfile).read())
            for _ in range(50):
                try:
                    os.kill(grandchild, 0)
                except ProcessLookupError:
                    break
                try:
                    os.waitpid(grandchild, os.WNOHANG)
                except ChildProcessError:
                    pass
                time.sleep(0.1)
            else:
                os.kill(grandchild, signal.SIGKILL)
                self.fail("the job's child process survived the interrupt")


if __name__ == "__main__":
    unittest.main()
