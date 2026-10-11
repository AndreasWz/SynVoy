"""The HTTP interface: what it answers, what it refuses, and that it never hands out
anything it should not."""

import json
import time
import urllib.error
import urllib.request

import pytest

from conftest import PROTEIN, REPO_DIR, call, call_json, run_until
from synvoy_web import db
from synvoy_web.api import App
from synvoy_web.server import serve_in_thread

MISSING = "Aa0000000000000000000a"          # looks like an id, belongs to nothing


@pytest.fixture
def app(env):
    env.beat()
    return App(env.ctx, log=lambda message: None)


def start(app, env, name="ok", remote="10.0.0.1", **kwargs):
    status, data = call_json(app, "POST", "/api/v1/searches", env.payload(name=name, **kwargs),
                             remote=remote)
    assert status in (200, 201), data
    return data


def finish(env, job_id):
    worker = env.worker()
    run_until(worker, lambda: env.store.get_job(job_id)["state"] in db.FINAL_STATES)
    return env.store.get_job(job_id)


# ── reading ─────────────────────────────────────────────────────────────────

def test_status_tells_limits_and_version(app, env):
    status, data = call_json(app, "GET", "/api/v1/status")
    assert status == 200 and data["accepting"] is True and data["api"] == "v1"
    assert data["limits"]["max_targets"] == env.settings.max_targets
    # No version number is claimed for a checkout that is not a tagged release: the
    # number in the VERSION file is shown only when the checkout is at that tag.
    assert set(data["pipeline"]) == {"label", "release", "commit", "modified"}
    number = (REPO_DIR / "VERSION").read_text().strip()
    if data["pipeline"]["release"]:
        assert data["pipeline"]["label"] == number
    else:
        assert data["pipeline"]["label"] == "pre-release" and number not in json.dumps(data)
    assert data["queue"] == {"waiting": 0, "running": 0}
    assert data["features"] == {"email": False}


def test_library_lists_genomes_without_paths(app, env):
    status, data = call_json(app, "GET", "/api/v1/library")
    assert status == 200 and len(data["genomes"]) == 6
    text = json.dumps(data)
    assert str(env.genomes_dir) not in text and "checksum" not in text and "fasta" not in text
    home = [g for g in data["genomes"] if g["id"] == "homus_primus"][0]
    assert home["can_be_home"] is True and home["size_mb"] == 12.0


def test_library_orders_targets_by_relatedness(app):
    status, data = call_json(app, "GET", "/api/v1/library", query="home=homus_primus")
    assert status == 200
    # shared taxonomy levels with the home genome: 3, 3, 2, 1, 0
    assert data["targets"] == ["alpha_beta", "big_one", "gamma_delta", "epsilon_zeta", "eta_theta"]
    assert call_json(app, "GET", "/api/v1/library", query="home=nobody")[0] == 404


def test_lookup_accession_finds_the_home_genome(app, env):
    status, data = call_json(app, "POST", "/api/v1/lookup", {"accession": " p01501 "})
    assert status == 200, data
    assert data["name"] == "Melittin" and data["gene"] == "MELT" and data["length"] == len(PROTEIN)
    assert data["home"]["id"] == "homus_primus"            # "(strain T1)" did not get in the way
    assert data["warnings"] == ["short_protein"] and data["acceptable"] is True
    assert "sequence" not in data
    calls = len(env.fetch.calls)
    call_json(app, "POST", "/api/v1/lookup", {"accession": "P01501"})
    assert len(env.fetch.calls) == calls                   # the second answer came from the cache


def test_lookup_of_a_pasted_sequence(app):
    status, data = call_json(app, "POST", "/api/v1/lookup", {"sequence": ">my prot\n" + PROTEIN})
    assert status == 200 and data["length"] == len(PROTEIN) and data["name"] == "my"
    status, data = call_json(app, "POST", "/api/v1/lookup", {"sequence": "ACGT" * 30})
    assert status == 400 and data["error"]["code"] == "sequence_is_dna"


def test_lookup_limits_requests_that_leave_the_server(make_env):
    env = make_env(max_lookups_per_hour=2)
    env.beat()
    app = App(env.ctx, log=lambda message: None)
    for accession in ("Q00001", "Q00002"):
        assert call_json(app, "POST", "/api/v1/lookup", {"accession": accession})[0] == 404
    status, data = call_json(app, "POST", "/api/v1/lookup", {"accession": "Q00003"})
    assert status == 429 and data["error"]["code"] == "too_many_lookups"
    # Another visitor is not affected, and neither is an answer that is already stored.
    assert call_json(app, "POST", "/api/v1/lookup", {"accession": "Q00004"}, remote="10.0.0.2")[0] == 404


# ── starting a search ───────────────────────────────────────────────────────

def test_submit_then_follow_to_the_result(app, env):
    data = start(app, env)
    assert data["state"] == "queued" and data["reused"] is False and len(data["token"]) >= 20
    assert data["queue"] == {"ahead": 0, "running": 0, "genomes": []}
    assert "fasta" not in json.dumps(data) and "client" not in json.dumps(data)
    job_id = data["id"]
    assert call_json(app, "GET", "/api/v1/searches/%s/result" % job_id)[0] == 409
    assert finish(env, job_id)["state"] == db.DONE

    status, view = call_json(app, "GET", "/api/v1/searches/" + job_id)
    assert status == 200 and view["state"] == "done" and view["runtime_s"] is not None
    assert "token" not in view
    status, result = call_json(app, "GET", "/api/v1/searches/%s/result" % job_id)
    assert status == 200 and result["summary"]["found"] == 2

    name = result["loci"][0]["files"]["grid"].split("/")[-1]
    status, headers, body = call(app, "GET", "/api/v1/searches/%s/files/%s" % (job_id, name))
    assert status == 200 and headers["Content-Type"] == "image/svg+xml" and body.startswith(b"<svg")
    assert "sandbox" in headers["Content-Security-Policy"]
    assert "Content-Disposition" not in headers
    status, headers, _ = call(app, "GET", "/api/v1/searches/%s/files/gene_models.zip" % job_id)
    assert status == 200 and headers["Content-Disposition"].startswith("attachment")
    page = result["loci"][0]["files"]["ribbons_page"].split("/")[-1]
    _, headers, _ = call(app, "GET", "/api/v1/searches/%s/files/%s" % (job_id, page))
    assert "allow-scripts" in headers["Content-Security-Policy"]


def test_same_search_twice_gives_the_same_search(app, env):
    first = start(app, env)
    status, second = call_json(app, "POST", "/api/v1/searches", env.payload(), remote="10.0.0.2")
    assert status == 200 and second["id"] == first["id"] and second["reused"] is True
    assert second["token"] != first["token"]
    finish(env, first["id"])
    status, third = call_json(app, "POST", "/api/v1/searches", env.payload(), remote="10.0.0.3")
    assert status == 200 and third["id"] == first["id"] and third["state"] == "done"
    assert len(env.job_ids()) == 1


def test_a_different_question_is_a_different_search(app, env):
    ids = {start(app, env)["id"],
           start(app, env, name="ok_other")["id"],
           start(app, env, targets=["alpha_beta"])["id"],
           start(app, env, targets=["alpha_beta", "eta_theta", "epsilon_zeta", "gamma_delta"])["id"]}
    assert len(ids) == 4
    # The order in which the genomes were ticked does not make a new search.
    assert start(app, env, targets=["eta_theta", "epsilon_zeta", "alpha_beta"])["id"] in ids


def test_submit_by_accession_uses_the_database_record(app, env):
    status, data = call_json(app, "POST", "/api/v1/searches",
                             {"query": {"accession": "P01501"}, "targets": ["alpha_beta"]})
    assert status == 201, data
    assert data["home"]["id"] == "homus_primus" and data["query"]["label"] == "MELT"
    stored = env.store.get_job(data["id"])["request"]["query"]["fasta"]
    assert stored.startswith(">sp|P01501|MEL_TEST Melittin") and PROTEIN[:60] in stored


@pytest.mark.parametrize("change, code", [
    ({"targets": []}, "no_targets"),
    ({"targets": ["alpha_beta", "alpha_beta"]}, "bad_request"),
    ({"targets": ["nobody"]}, "target_unknown"),
    ({"targets": ["homus_primus"]}, "target_is_home"),
    ({"targets": ["big_one"]}, "genome_too_large"),
    ({"targets": "alpha_beta"}, "bad_request"),
    ({"targets": [7]}, "bad_request"),
    ({"extra": 1}, "bad_request"),
    ({"email": "someone@example.org"}, "email_off"),
    ({"query": {"sequence": PROTEIN}}, "home_missing"),
    ({"query": {"sequence": PROTEIN, "home": "nobody"}}, "home_unknown"),
    ({"query": {"sequence": PROTEIN, "home": "alpha_beta"}}, "home_without_annotation"),
    ({"query": {"sequence": PROTEIN[:10], "home": "homus_primus"}}, "sequence_too_short"),
    ({"query": {"sequence": PROTEIN + "*" + PROTEIN, "home": "homus_primus"}}, "sequence_invalid"),
    ({"query": {"sequence": ">a\n" + PROTEIN + "\n>b\n" + PROTEIN, "home": "homus_primus"}},
     "several_sequences"),
    ({"query": {"sequence": PROTEIN, "accession": "P01501"}}, "bad_request"),
    ({"query": {"accession": "P01501-2"}}, "accession_isoform"),
    ({"query": {"accession": "not an accession"}}, "accession_unknown"),
    ({"query": {"accession": "Q99999"}}, "accession_not_found"),
    ({"query": "P01501"}, "bad_request"),
])
def test_requests_that_are_refused(app, env, change, code):
    body = env.payload()
    body.update(change)
    status, data = call_json(app, "POST", "/api/v1/searches", body)
    assert 400 <= status < 500 and data["error"]["code"] == code, data
    assert data["error"]["message"]
    assert env.job_ids() == []


def test_limits_on_number_and_size_of_genomes(make_env):
    env = make_env(max_targets=2, max_total_mb=1000)
    env.beat()
    app = App(env.ctx, log=lambda message: None)
    status, data = call_json(app, "POST", "/api/v1/searches",
                             env.payload(targets=["alpha_beta", "epsilon_zeta", "eta_theta"]))
    assert status == 400 and data["error"]["code"] == "too_many_targets"
    status, data = call_json(app, "POST", "/api/v1/searches",
                             env.payload(targets=["alpha_beta", "gamma_delta"]))
    assert status == 400 and data["error"]["code"] == "selection_too_large"


def test_limits_per_visitor_and_on_the_queue(make_env):
    env = make_env(max_active_per_client=1, max_queue=2, max_submissions_per_day=3)
    env.beat()
    app = App(env.ctx, log=lambda message: None)
    start(app, env, "ok_1", remote="10.0.0.1")
    status, data = call_json(app, "POST", "/api/v1/searches", env.payload(name="ok_2"),
                             remote="10.0.0.1")
    assert status == 429 and data["error"]["code"] == "too_many_active"
    start(app, env, "ok_3", remote="10.0.0.2")
    status, data = call_json(app, "POST", "/api/v1/searches", env.payload(name="ok_4"),
                             remote="10.0.0.3")
    assert status == 503 and data["error"]["code"] == "queue_full"
    # A search that exists already is still handed out when the queue is full.
    assert call_json(app, "POST", "/api/v1/searches", env.payload(name="ok_1"),
                     remote="10.0.0.3")[0] == 200


def test_limit_per_day(make_env):
    env = make_env(max_submissions_per_day=2)
    env.beat()
    app = App(env.ctx, log=lambda message: None)
    start(app, env, "ok_1")
    start(app, env, "ok_2")
    status, data = call_json(app, "POST", "/api/v1/searches", env.payload(name="ok_3"))
    assert status == 429 and data["error"]["code"] == "too_many_submissions"
    assert call_json(app, "POST", "/api/v1/searches", env.payload(name="ok_3"),
                     remote="10.0.0.2")[0] == 201


def test_no_search_is_taken_when_none_could_be_run(make_env):
    env = make_env()
    app = App(env.ctx, log=lambda message: None)
    status, data = call_json(app, "POST", "/api/v1/searches", env.payload())
    assert status == 503 and data["error"]["code"] == "worker_down"      # no worker yet
    env.beat()
    env.store.set_meta("worker_heartbeat", int(time.time()) - 10000)
    assert call_json(app, "POST", "/api/v1/searches", env.payload())[1]["error"]["code"] == "worker_down"
    env.beat()
    env.store.set_meta("worker_preflight", {"ok": False, "problems": ["no nextflow"]})
    status, data = call_json(app, "GET", "/api/v1/status")
    assert data["accepting"] is False and data["reason"] == "worker_not_ready"
    assert "nextflow" not in json.dumps(data)                # details stay with the operator
    assert env.job_ids() == []


def test_operator_can_close_the_service(make_env):
    env = make_env(accepting=False, closed_message="Back on Monday.")
    env.beat()
    app = App(env.ctx, log=lambda message: None)
    status, data = call_json(app, "POST", "/api/v1/searches", env.payload())
    assert status == 503 and data["error"]["message"] == "Back on Monday."


# ── cancelling ──────────────────────────────────────────────────────────────

def test_cancel_over_http(app, env):
    data = start(app, env)
    path = "/api/v1/searches/%s/cancel" % data["id"]
    assert call_json(app, "POST", path, {"token": "wrong"})[0] == 403
    assert call_json(app, "POST", path, {})[0] == 403
    assert call_json(app, "POST", path, ["x"])[0] == 403
    status, view = call_json(app, "POST", path, {"token": data["token"]})
    assert status == 200 and view["cancelled"] is True and view["state"] == "cancelled"
    assert call_json(app, "POST", path, {"token": data["token"]})[0] == 409
    assert call_json(app, "POST", "/api/v1/searches/%s/cancel" % MISSING, {"token": "x"})[0] == 404


# ── what is refused at the door ─────────────────────────────────────────────

def test_request_bodies_are_checked(app, env):
    path = "/api/v1/searches"
    assert call(app, "POST", path, b"{not json")[0] == 400
    assert call(app, "POST", path, env.payload(), headers={"Content-Type": "text/plain"})[0] == 415
    assert call(app, "POST", path, env.payload(),
                headers={"Content-Type": "application/x-www-form-urlencoded"})[0] == 415
    assert call(app, "POST", path, env.payload(), headers={"Content-Length": None})[0] == 411
    assert call(app, "POST", path, env.payload(), headers={"Content-Length": "abc"})[0] == 411
    assert call(app, "POST", path, b"x" * (200 * 1024))[0] == 413
    assert call(app, "POST", path, b'{"a":1}', headers={"Content-Length": "9999"})[0] == 400
    assert call(app, "POST", path, b"\xff\xfe")[0] == 400
    assert env.job_ids() == []


def test_unknown_addresses_and_methods(app):
    assert call_json(app, "GET", "/api/v1/nothing")[0] == 404
    assert call_json(app, "GET", "/api/v2/status")[0] == 404
    status, headers, _ = call(app, "DELETE", "/api/v1/status")
    assert status == 405 and headers["Allow"] == "GET"
    status, headers, _ = call(app, "GET", "/api/v1/searches")
    assert status == 405 and headers["Allow"] == "POST"
    assert call(app, "POST", "/index.html", {})[0] == 405
    for bad in ("short", "A" * 23, "../../etc/passwd", "-" + "a" * 21, "a" * 21 + "!"):
        assert call_json(app, "GET", "/api/v1/searches/" + bad)[0] == 404
    assert call_json(app, "GET", "/api/v1/searches/" + MISSING)[0] == 404
    assert call_json(app, "GET", "/api/v1/searches/%s/result" % MISSING)[0] == 404


def test_only_listed_files_of_a_result_are_served(app, env):
    job_id = start(app, env)["id"]
    finish(env, job_id)
    base = "/api/v1/searches/%s/files/" % job_id
    for name in ("result.json", "not_listed.txt", "..", "../result.json", "%2e%2e", ".hidden",
                 "..%2f..%2fstate", "synvoy_report.json/", "a/b"):
        assert call(app, "GET", base + name)[0] == 404, name
    # A file that exists in the bundle folder but is not in the result is not served either.
    extra = env.settings.jobs_dir / job_id / "result" / "files" / "extra.txt"
    extra.write_text("x")
    assert call(app, "GET", base + "extra.txt")[0] == 404
    assert call(app, "GET", base + "synvoy_report.json")[0] == 200


def test_expired_result_says_so(app, env):
    job_id = start(app, env)["id"]
    finish(env, job_id)
    env.store.update_job(job_id, expires_at=int(time.time()) - 1)
    env.worker().collect_garbage()
    status, data = call_json(app, "GET", "/api/v1/searches/%s/result" % job_id)
    assert status == 410 and data["error"]["code"] == "expired"
    assert call_json(app, "GET", "/api/v1/searches/" + job_id)[1]["state"] == "expired"


def test_failed_search_shows_its_reason(app, env):
    job_id = start(app, env, "nogene")["id"]
    finish(env, job_id)
    status, view = call_json(app, "GET", "/api/v1/searches/" + job_id)
    assert view["state"] == "failed" and view["error"]["code"] == "gene_not_in_home"
    assert "Homus primus" in view["error"]["message"]
    assert call_json(app, "GET", "/api/v1/searches/%s/result" % job_id)[0] == 409


def test_an_error_inside_the_server_gives_no_details(env, monkeypatch):
    messages = []
    app = App(env.ctx, log=messages.append)
    monkeypatch.setattr("synvoy_web.jobs.status_view", lambda ctx: 1 / 0)
    status, data = call_json(app, "GET", "/api/v1/status")
    assert status == 500 and data["error"]["code"] == "internal_error"
    assert "ZeroDivision" not in json.dumps(data) and "Traceback" not in json.dumps(data)
    assert any("ZeroDivisionError" in m for m in messages)     # the log has it


# ── pages and headers ───────────────────────────────────────────────────────

def test_pages_are_served_with_strict_headers(app):
    status, headers, body = call(app, "GET", "/")
    assert status == 200 and b"<title>home</title>" in body
    assert headers["Content-Type"].startswith("text/html")
    policy = headers["Content-Security-Policy"]
    assert "script-src 'self'" in policy and "unsafe-inline" not in policy
    assert "frame-ancestors 'none'" in policy
    assert headers["X-Content-Type-Options"] == "nosniff"
    assert headers["Referrer-Policy"] == "no-referrer"       # a search id must not leak in a link
    assert call(app, "GET", "/assets/site.css")[1]["Content-Type"].startswith("text/css")
    status, headers, body = call(app, "HEAD", "/index.html")
    assert status == 200 and body == b"" and int(headers["Content-Length"]) > 0


def test_pages_are_revalidated_not_resent(app):
    _, headers, _ = call(app, "GET", "/index.html")
    status, _, body = call(app, "GET", "/index.html", headers={"If-None-Match": headers["ETag"]})
    assert status == 304 and body == b""


@pytest.mark.parametrize("path", [
    "/../var/secret", "/assets/../../var/secret", "/.git/config", "/assets/.hidden",
    "/..", "/assets/..", "/nothing.html", "/assets", "/assets/", "//etc/passwd",
])
def test_nothing_outside_the_pages_folder_is_served(app, env, path):
    (env.data_dir / "secret").write_text("secret")
    (env.site_dir / "assets" / ".hidden").write_text("hidden")
    status, _, body = call(app, "GET", path)
    assert status == 404 and b"secret" not in body and b"hidden" not in body


def test_links_out_of_the_pages_folder_are_not_followed(app, env):
    (env.data_dir / "secret").write_text("secret")
    (env.site_dir / "leak.txt").symlink_to(env.data_dir / "secret")
    assert call(app, "GET", "/leak.txt")[0] == 404


def test_strange_paths_are_refused(app):
    for path in ("/caf\xe9", "/a\x00b", "/a\\b", "/\x7f"):
        assert call(app, "GET", path)[0] == 400


def test_example_figures_get_the_sandbox_policy(app, env):
    folder = env.site_dir / "data" / "examples" / "demo" / "files"
    folder.mkdir(parents=True)
    (folder / "plot.html").write_text("<!doctype html><script>1</script>")
    (env.site_dir / "data" / "gallery.json").write_text("{}")
    _, headers, _ = call(app, "GET", "/data/examples/demo/files/plot.html")
    assert headers["Content-Security-Policy"].startswith("sandbox allow-scripts")
    _, headers, _ = call(app, "GET", "/data/gallery.json")
    assert "script-src 'self'" in headers["Content-Security-Policy"]


def test_site_without_examples_has_an_empty_list(app, env):
    """A new installation has no examples built. The start page asks for the list:
    it gets an empty one, not an error."""
    assert not (env.site_dir / "data").exists()
    status, headers, body = call(app, "GET", "/data/gallery.json")
    assert status == 200 and json.loads(body) == {"examples": []}
    assert headers["Content-Type"].startswith("application/json")
    assert call(app, "GET", "/data/examples/demo/result.json")[0] == 404    # only the list
    assert call(app, "GET", "/data/other.json")[0] == 404


def test_api_answers_are_never_stored_by_caches(app, env):
    for path in ("/api/v1/status", "/api/v1/searches/" + MISSING):
        _, headers, _ = call(app, "GET", path)
        assert headers["Cache-Control"] == "no-store"
        assert headers["X-Robots-Tag"] == "noindex, nofollow"


def test_other_sites_may_call_only_when_listed(make_env):
    env = make_env(cors_origins=["https://pages.example.org"])
    env.beat()
    app = App(env.ctx, log=lambda message: None)
    _, headers, _ = call(app, "GET", "/api/v1/status", headers={"Origin": "https://pages.example.org"})
    assert headers["Access-Control-Allow-Origin"] == "https://pages.example.org"
    _, headers, _ = call(app, "GET", "/api/v1/status", headers={"Origin": "https://evil.example"})
    assert "Access-Control-Allow-Origin" not in headers
    status, headers, _ = call(app, "OPTIONS", "/api/v1/searches",
                              headers={"Origin": "https://pages.example.org"})
    assert status == 204 and "POST" in headers["Access-Control-Allow-Methods"]
    status, headers, _ = call(app, "OPTIONS", "/api/v1/searches",
                              headers={"Origin": "https://evil.example"})
    assert "Access-Control-Allow-Origin" not in headers
    assert "https://pages.example.org" in call(app, "GET", "/")[1]["Content-Security-Policy"]


def test_visitor_address_behind_a_proxy(make_env):
    direct = make_env()
    proxied = make_env(trust_proxy=True)
    environ = {"REMOTE_ADDR": "127.0.0.1", "HTTP_X_FORWARDED_FOR": "6.6.6.6, 203.0.113.7"}
    assert App(direct.ctx).client_address(environ) == "127.0.0.1"     # header is not believed
    assert App(proxied.ctx).client_address(environ) == "203.0.113.7"  # what the proxy saw
    assert App(proxied.ctx).client_address({"REMOTE_ADDR": "127.0.0.1"}) == "127.0.0.1"


def test_visitor_addresses_are_not_stored(app, env):
    start(app, env, remote="203.0.113.77")
    dump = "\n".join(env.store._conn().iterdump())
    assert "203.0.113.77" not in dump


def test_log_shortens_search_ids(env):
    messages = []
    app = App(env.ctx, log=messages.append)
    call(app, "GET", "/api/v1/searches/" + MISSING)
    assert MISSING not in messages[-1] and MISSING[:6] in messages[-1]


def test_health_check(app):
    status, _, body = call(app, "GET", "/healthz")
    assert status == 200 and body == b"ok\n"


def test_over_a_real_socket(app, env):
    server, _ = serve_in_thread(app)
    base = "http://127.0.0.1:%d" % server.server_address[1]
    try:
        with urllib.request.urlopen(base + "/api/v1/status", timeout=10) as response:
            assert json.loads(response.read())["service"] == "synvoy-web"
        request = urllib.request.Request(
            base + "/api/v1/searches", data=json.dumps(env.payload()).encode("utf-8"),
            headers={"Content-Type": "application/json"})
        with urllib.request.urlopen(request, timeout=10) as response:
            assert response.status == 201 and json.loads(response.read())["state"] == "queued"
        with pytest.raises(urllib.error.HTTPError) as caught:
            urllib.request.urlopen(base + "/api/v1/searches/" + MISSING, timeout=10)
        assert caught.value.code == 404
        assert json.loads(caught.value.read())["error"]["code"] == "not_found"
        with urllib.request.urlopen(base + "/", timeout=10) as response:
            assert b"<title>home</title>" in response.read()
    finally:
        server.shutdown()
        server.server_close()
