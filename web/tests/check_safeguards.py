#!/usr/bin/env python3
"""Do the tests notice when a safeguard is taken out?

    python tests/check_safeguards.py            (from web/; about a minute)

For each case below, one safeguard is removed in a throw-away copy of this folder and
the test that should guard it is run there. The test must fail. A case whose text is
no longer in the source counts as a failure of this check, so the list cannot go stale
unnoticed. Before the cases, the guarding tests are run on an unchanged copy, where
they must pass (otherwise "the test failed" would prove nothing).

This is not part of the normal test run. Run it after changing a safeguard or its
test, and add a case when you add a safeguard.

Each case: (file in synvoy_web/, text to replace, replacement, test that must fail).
"""

import os
import secrets
import shutil
import signal
import subprocess
import sys
import tempfile
import time
from pathlib import Path

WEB = Path(__file__).resolve().parents[1]

# Everything this check starts carries this in its environment. With a safeguard gone,
# a stand-in search can survive its test (that is what some cases are about), and one
# kind never ends by itself. They are found by this mark and ended after every case.
MARK = "SAFEGUARD_CHECK_RUN=" + secrets.token_hex(8)

T_API = "tests/test_api.py::"
T_RUN = "tests/test_runner.py::"
T_ENV = "tests/test_environment.py::"
T_INP = "tests/test_inputs.py::"
T_RES = "tests/test_results.py::"
T_GEN = "tests/test_genomes.py::"

CASES = [
    # ── files and requests ──
    ("api.py", "            resolved.relative_to(self.site_dir)\n", "            pass\n",
     T_API + "test_links_out_of_the_pages_folder_are_not_followed"),
    ("api.py", "        if relative not in listed or not path.is_file():",
     "        if not path.is_file():",
     T_API + "test_only_listed_files_of_a_result_are_served"),
    ("api.py", '        if content_type != "application/json":', "        if False:",
     T_API + "test_request_bodies_are_checked"),
    ("api.py", '        headers["Content-Security-Policy"] = FIGURE_PAGE_CSP', "        pass",
     T_API + "test_example_figures_get_the_sandbox_policy"),
    ("net.py", '    if parts.scheme != "https" or parts.hostname not in ALLOWED_HOSTS:',
     "    if False:",
     T_INP + "test_requests_go_to_the_two_databases_only"),
    # ── limits, ownership, privacy ──
    ("jobs.py", "        if active >= settings.max_active_per_client:", "        if False:",
     T_API + "test_limits_per_visitor_and_on_the_queue"),
    ("jobs.py", '        if remaining is None:\n            raise QueryError("not_yours"',
     '        if False:\n            raise QueryError("not_yours"',
     T_RUN + "test_cancel_needs_the_token_of_the_submitter"),
    ("jobs.py", '        text = self._secret + "|" + (address or "unknown")\n'
                '        return hashlib.sha256(text.encode("utf-8")).hexdigest()[:32]',
     '        return address or "unknown"',
     T_API + "test_visitor_addresses_are_not_stored"),
    ("runner.py", '            self.store.forget_visitors(job["id"])', "            pass",
     T_RUN + "test_nothing_links_a_finished_search_to_a_visitor"),
    ("jobs.py", '        "code": version["fingerprint"],\n', "",
     T_INP + "test_cache_key_follows_what_decides_the_result"),
    # ── processes and files of a search ──
    ("procs.py", "        _signal(pids, signal.SIGKILL)\n", "        pass\n",
     T_RUN + "test_time_limit_ends_a_search_that_ignores_requests_to_stop"),
    ("fsutil.py", "    if resolved == resolved_root or resolved_root not in resolved.parents:",
     "    if False:",
     T_RUN + "test_remove_tree_never_leaves_its_root"),
    ("runner.py",
     "        self.stop(job, handle)                    # nothing should be left; make sure\n", "",
     T_RUN + "test_process_left_behind_by_a_finished_search_is_stopped"),
    ("runner.py", '        elif job["attempt"] < self.settings.max_attempts:', "        elif False:",
     T_RUN + "test_search_that_died_with_the_machine_is_resumed"),
    ("runner.py", '        if not unchanged:\n'
                  '            write_atomic(query_file, request["query"]["fasta"])\n'
                  '        os.utime(str(query_file), (job["created_at"], job["created_at"]))\n',
     '        write_atomic(query_file, request["query"]["fasta"])\n',
     T_RUN + "test_search_that_died_with_the_machine_is_resumed"),
    ("runner.py", "            if resume:\n                # An earlier attempt ran with the older code.",
     "            if False:\n                # An earlier attempt ran with the older code.",
     T_RUN + "test_search_started_again_after_a_code_update_begins_from_scratch"),
    ("runner.py", '        if free < need:\n            return "low_disk"',
     '        if False:\n            return "low_disk"',
     T_RUN + "test_search_waits_when_the_disk_is_short"),
    ("runner.py",
     '        env = {k: v for k, v in os.environ.items() if not k.startswith("SYNVOY_WEB_")}',
     "        env = dict(os.environ)",
     T_RUN + "test_service_settings_are_not_passed_to_the_pipeline"),
    # ── the environment of a search ──
    ("run_job.sh", '    [ -f "$entry/../pyvenv.cfg" ] && continue\n', "",
     T_ENV + "test_search_gets_the_conda_python_whatever_is_in_front"),
    ("run_job.sh",
     '    [ -n "${VIRTUAL_ENV:-}" ] && [ "$entry" = "${VIRTUAL_ENV}/bin" ] && continue\n', "",
     T_ENV + "test_virtualenv_named_only_by_the_variable_is_dropped_too"),
    ("run_job.sh", '    PATH="$CONDA_PREFIX/bin${rest:+:$rest}"\n', "    :\n",
     T_ENV + "test_search_gets_the_conda_python_whatever_is_in_front"),
    ("run_job.sh", "unset VIRTUAL_ENV\n", "",
     T_ENV + "test_search_gets_the_conda_python_whatever_is_in_front"),
    ("runner.py", "                problem = self.probe_task_environment()",
     "                problem = None",
     T_ENV + "test_start_up_check_of_the_worker_looks_at_the_tasks_python"),
    ("runner.py", "        if not used.startswith(str(newest)):", "        if False:",
     T_ENV + "test_start_up_check_of_the_worker_looks_at_the_tasks_python"),
    ("runner.py",
     '        if error_code == "server_setup":\n            self.pause_after_setup_failure(job["id"])',
     '        if False:\n            self.pause_after_setup_failure(job["id"])',
     T_ENV + "test_broken_software_pauses_the_queue"),
    ("runner.py",
     '            if self.ctx.now() < int(paused.get("until") or 0):\n                return',
     "            if False:\n                return",
     T_ENV + "test_broken_software_pauses_the_queue"),
    ("runner.py", 'r"parasail is not importable|ModuleNotFoundError|No module named"',
     'r"NEVER_IN_A_LOG"',
     T_ENV + "test_broken_software_pauses_the_queue"),
    ("runner.py", '        if paused.get("code") == "operator":\n            return',
     "        if False:\n            return",
     T_ENV + "test_operator_can_pause_and_resume"),
    # ── genomes fetched from NCBI ──
    ("ncbi.py", "    return value if ACCESSION_RE.match(value) else None", "    return value",
     T_GEN + "test_only_assembly_accessions_with_a_version_are_taken"),
    ("ncbi.py", "            if file_md5(part) != md5:", "            if False:",
     T_GEN + "test_damaged_file_is_fetched_again_and_refused_the_second_time"),
    ("ncbi.py", "                if size == 0 and block[:1] != first_byte:", "                if False:",
     T_GEN + "test_file_that_is_not_what_its_name_says_is_refused"),
    ("genomes.py", '        candidates = [r for r in rows if r["state"] == READY and not r["pinned"]\n',
     '        candidates = [r for r in rows if r["state"] == READY\n',
     T_GEN + "test_pinned_and_needed_genomes_are_never_removed"),
    ("genomes.py", '                      and r["accession"] not in protect]', "                      ]",
     T_GEN + "test_pinned_and_needed_genomes_are_never_removed"),
    ("genomes.py", "            if need > self.budget - locked:", "            if False:",
     T_GEN + "test_pinned_and_needed_genomes_are_never_removed"),
    ("genomes.py", "            locked = self._locked_bytes(rows, protect, accession) + stray\n",
     "            locked = self._locked_bytes(rows, protect, accession)\n",
     T_GEN + "test_files_in_no_record_count_against_the_size_of_the_store"),
    ("genomes.py", '            used = stray + sum(r["disk_bytes"] for r in rows if r["accession"] != keep)\n',
     '            used = sum(r["disk_bytes"] for r in rows if r["accession"] != keep)\n',
     T_GEN + "test_files_in_no_record_count_against_the_size_of_the_store"),
    ("genomes.py", '            if kept.is_file() and suffix == "fna" and row["ready_at"]:',
     '            if kept.is_file() and suffix == "fna":',
     T_GEN + "test_file_the_store_did_not_fetch_is_not_taken_for_the_genome"),
    ("genomes.py", "        skip.update(self._operators_own())\n", "",
     T_GEN + "test_genomes_the_operator_registered_are_never_taken_for_strays"),
    ("genomes.py", "        if not self.own_files()[1]:", "        if False:",
     T_GEN + "test_genomes_the_operator_registered_are_never_taken_for_strays"),
    ("genomes.py", '            if origin == "visitor" and new_request:', "            if False:",
     T_GEN + "test_one_visitor_cannot_have_many_genomes_fetched"),
    ("genomes.py", '        if size_mb > self.settings.max_genome_mb and origin != "operator":',
     "        if False:",
     T_GEN + "test_genome_larger_than_the_limits_is_refused_before_any_download"),
    ("genomes.py", '            if row is None or row["state"] != WANTED:\n                return row',
     "            if row is None:\n                return row",
     T_GEN + "test_download_cut_off_by_a_restart_is_wanted_again"),
    ("runner.py", '            if state == "waiting":\n                continue',
     "            if False:\n                continue",
     T_GEN + "test_search_waits_for_its_genome_and_then_runs"),
    ("runner.py", '            if state == "failed":\n                self.finalise(job, db.FAILED, "genome_download_failed", detail)\n                continue',
     '            if False:\n                continue',
     T_GEN + "test_search_fails_with_the_reason_when_its_genome_cannot_be_fetched"),
    # ── results ──
    ("runner.py", "            if self.code_changed_since_launch(job):", "            if False:",
     T_RUN + "test_code_updated_under_a_running_search_is_flagged"),
    ("results.py", '"GOI_placeholder" not in newick.read_text()[:2000]', "True",
     T_RES + "test_stand_in_tree_is_not_offered"),
    ("results.py", '    reason = _SET_ASIDE.get(record.get("goi_class") or "")', "    reason = None",
     T_RES + "test_labels_and_counts"),
    ("results.py", '    label = re.sub(r"[^A-Za-z0-9-]+", "_", short_name(query))',
     "    label = short_name(query)",
     T_RES + "test_archive_name_is_made_from_the_query_name_safely"),
    ("gallery.py", "    if hits:\n        shutil.rmtree(str(bundle))", "    if False:\n        shutil.rmtree(str(bundle))",
     T_RES + "test_example_is_published_under_its_own_name_not_under_the_runs_label"),
    ("jobs.py", '    release = version != "unknown" and tag in (version, "v" + version) and status == ""',
     '    release = version != "unknown"',
     "tests/test_inputs.py::test_version_number_is_shown_for_a_tagged_release_only"),
]


def copy_to(folder):
    shutil.copytree(str(WEB), folder, ignore=shutil.ignore_patterns(
        "var", "data", "__pycache__", ".pytest_cache"))
    return folder


def end_leftovers():
    """End every process started under this check that is still there. Returns how many."""
    mark = MARK.encode("utf-8")
    ended = 0
    for _ in range(3):                                  # a dying process may start another
        found = 0
        for name in os.listdir("/proc"):
            if not name.isdigit() or int(name) == os.getpid():
                continue
            try:
                with open("/proc/%s/environ" % name, "rb") as handle:
                    if mark not in handle.read().split(b"\0"):
                        continue
                os.kill(int(name), signal.SIGKILL)
                found += 1
            except OSError:                             # gone, or not ours to read
                continue
        ended += found
        if not found:
            break
        time.sleep(0.3)
    return ended


def run_tests(folder, tests):
    env = {k: v for k, v in os.environ.items()
           if k != "VIRTUAL_ENV" and not k.startswith("SYNVOY_WEB_")}
    name, value = MARK.split("=", 1)
    env[name] = value
    try:
        return subprocess.run([sys.executable, "-m", "pytest", "-q", "-x", "-p", "no:cacheprovider"]
                              + list(tests), cwd=folder, env=env, capture_output=True, text=True,
                              timeout=900)
    finally:
        end_leftovers()


def main():
    scratch = tempfile.mkdtemp(prefix="synvoy_web_safeguards_")
    not_caught = 0
    try:
        clean = copy_to(os.path.join(scratch, "unchanged"))
        done = run_tests(clean, sorted({case[3] for case in CASES}))
        print("unchanged copy: " + (done.stdout.strip().splitlines() or ["no output"])[-1])
        if done.returncode != 0:
            print(done.stdout[-3000:] + done.stderr[-1000:])
            print("The guarding tests do not pass on an unchanged copy; nothing was checked.")
            return 2
        for number, (name, old, new, test) in enumerate(CASES):
            folder = copy_to(os.path.join(scratch, "case%d" % number))
            path = os.path.join(folder, "synvoy_web", name)
            with open(path) as handle:
                text = handle.read()
            if text.count(old) != 1:
                print("CASE DOES NOT APPLY (%d places)  %s  %r" % (text.count(old), name, old[:60]))
                not_caught += 1
            else:
                with open(path, "w") as handle:
                    handle.write(text.replace(old, new, 1))
                done = run_tests(folder, [test])
                caught = done.returncode == 1            # 1 = a test failed; else the run broke
                print("%-8s %-11s %s" % ("caught" if caught else "MISSED", name, test.split("::")[-1]))
                if not caught:
                    print("    " + "\n    ".join(done.stdout.strip().splitlines()[-6:]))
                    not_caught += 1
            shutil.rmtree(folder, ignore_errors=True)
    finally:
        end_leftovers()
        shutil.rmtree(scratch, ignore_errors=True)
    print("%d cases, %d not caught" % (len(CASES), not_caught))
    return 1 if not_caught else 0


if __name__ == "__main__":
    sys.exit(main())
