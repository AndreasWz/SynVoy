# SynVoy web

A website for SynVoy: a visitor names a protein, picks species from the genomes kept
on the server (or names one that is not there, and the server fetches it from NCBI),
and gets a link to a page that shows the progress and then the result.

It has two parts that can be used apart:

- **The pages** (`site/`). Plain HTML, CSS and JavaScript, no build step. Without a
  server behind them they show the examples of finished searches.
- **The service** (`synvoy_web/`). A small HTTP API, a queue in SQLite, and a worker
  that runs the pipeline. Python standard library only.

The pipeline is not changed and not imported. The worker starts
`nextflow run main.nf` exactly as the command line does and reads the files the
pipeline publishes (`docs/OUTPUT.md`).

## Try it on this machine

You need what a normal SynVoy run needs (the conda environment `synvoy_env`) and the
demo genomes (`demo/fetch_data.sh`).

```bash
cd web

# 1. Settings: where conda is. Everything else has a default.
mkdir -p var
echo '{"conda_base": "'"$(conda info --base)"'"}' > var/settings.json
export SYNVOY_WEB_CONFIG=var/settings.json

# 2. Genomes the site may search. With --gff a genome can also be the starting genome.
D=../local_data/demo
python -m synvoy_web library add --species "Saccharomyces cerevisiae" --group yeasts \
    --group-name "Yeasts" --fasta $D/home/Saccharomyces_cerevisiae.fna \
    --gff $D/home/Saccharomyces_cerevisiae.gff --assembly GCF_000146045.2
for s in Naumovozyma_castellii Kluyveromyces_lactis Lachancea_thermotolerans; do
    python -m synvoy_web library add --species "${s/_/ }" --group yeasts --fasta $D/genomes/$s.fna
done

# 3. Is everything in place?  Then start pages, API and worker in one process.
python -m synvoy_web check
python -m synvoy_web run          # http://127.0.0.1:8080/
```

On the page, enter `D6VTK4` (yeast STE2). On a 12-core laptop that search took
5 min 14 s.

The first search builds the pipeline's tool environment (several GB, several
minutes, needs the network). It is kept in `var/conda_cache` and reused by every
later search. To reuse one you already have, set `"conda_cache_dir"` to it.

## How it is built

```
 browser ──► reverse proxy (TLS) ──► API process ──► SQLite file ◄── worker process
                                        │   (queue, limits, cache)        │
                                        │                                 ├─ starts: run_job.sh → nextflow run main.nf
                                        └─ serves pages + result files    ├─ watches: /proc, exit code file, time limit
                                                                          └─ builds: result.json + figures  (results.py)
```

| File | What it does |
|---|---|
| `synvoy_web/config.py` | Settings: defaults, JSON file, environment. An unknown key is an error. |
| `synvoy_web/db.py` | The SQLite store. Searches, who asked, rate-limit events, look-up cache. |
| `synvoy_web/library.py` | The genome library file and "closest relatives first". |
| `synvoy_web/genomes.py` | The genome store: genomes fetched from NCBI, kept within a size, least used removed first. |
| `synvoy_web/ncbi.py` | What NCBI has for a species or accession; resumable, checksum-verified download. |
| `synvoy_web/query.py` | Checks a pasted sequence; looks an accession up at UniProt / NCBI. |
| `synvoy_web/jobs.py` | Accepts a search: checks, limits, result cache. Starts nothing. |
| `synvoy_web/runner.py` | The worker. Starts, watches, stops, resumes, cleans up; fetches genomes that searches wait for. |
| `synvoy_web/run_job.sh` | Sets up the environment of one search and records its exit code. |
| `synvoy_web/procs.py` | Finds and stops all processes of one search. |
| `synvoy_web/progress.py` | Reads how far a running search is from the pipeline's output. |
| `synvoy_web/results.py` | Turns the pipeline's output folder into the result page's data. |
| `synvoy_web/api.py` | The HTTP API and the pages, with their security headers. |
| `synvoy_web/gallery.py` | Builds the examples on the start page from finished runs. |
| `synvoy_web/cli.py` | `python -m synvoy_web ...` |
| `site/` | The pages. `site/config.js` says where the service is. |
| `deploy/` | Examples: settings, service units, reverse proxy, pipeline settings for a server. |
| `tests/` | See "Tests" below. |

A search is in one state at a time: `queued` → `running` → `done`, or `failed`,
`cancelled`, and later `expired`. Its folder is `<data_dir>/jobs/<id>/`. When it ends,
the working files are deleted; the result (a few hundred kB) and the logs stay until
the result expires.

### What a result page shows

- **Every genome that was searched has a row**, also one without a call. The row says
  how many of the gene's neighbours lie side by side in that genome (`neighbours` in
  `result.json`: counted from the gene models the search wrote, among the neighbours
  that were searched as proteins) and, where no call was made, how many weak matches
  of the gene itself were seen (`leads`). A genome without the gene is a result.
  The figures draw these genomes too; that is the pipeline's default since 2026-10-09
  (`hide_goi_absent_tracks`, `docs/OUTPUT.md`).
- **Everything in one file**: `files/synvoy_<name>.zip` holds the figures, the gene
  models, the report, the table and a read-me. It is built with the result, so the
  examples have it as well.
- **What made the result.** A version number is shown only when the checkout stands
  exactly on the git tag of its `VERSION` file (`v1.0.0`) and has no local changes.
  Anything else reads "pre-release" with the commit. `version_label` in the settings
  replaces this for an installation without git history.

### The same question is answered once

A search is identified by the query sequence, the genomes, the pipeline's code (a hash
over `bin/`, `modules/`, `main.nf`, `nextflow.config`, ...) and the run settings. A
second request with the same four gets the first one's result page at once, or joins
it while it is still running. Changing the pipeline's code changes the hash, so results
of older code are never handed out for new code.

## What keeps it running

Each line names the test that fails when the safeguard is taken out
(`tests/test_runner.py` unless another file is named). `tests/check_safeguards.py`
does exactly that for 31 of them: it removes one safeguard at a time in a copy and
checks that the test notices.

| What can go wrong | What the service does | Test |
|---|---|---|
| The worker or the whole service is restarted during a search | The search keeps running. The new worker finds it through `/proc` and its exit-code file and carries on. | `test_restarted_worker_takes_over_a_running_search`; browser run-through |
| The machine reboots or a search is killed | The search is queued again and started with `-resume`; steps that had finished are reused. After `max_attempts` (2) it fails with a reason. | `test_search_that_died_with_the_machine_is_resumed`; `test_end_to_end.py` |
| The pipeline is updated while a search runs | The result page carries a warning and the result is not given to anyone else. A search that is started again after the update begins from scratch, so nothing computed with the older code goes into its result. | `test_code_updated_under_a_running_search_is_flagged`, `test_search_started_again_after_a_code_update_begins_from_scratch` |
| A search never ends | It is stopped after `max_runtime_hours` (12). Every process it started is ended: asked first, then forced. | `test_time_limit_ends_a_search_that_ignores_requests_to_stop` |
| A visitor cancels | Same stop. Only the browser that started a search can cancel it; a search two visitors wait for stops when both have left. | `test_cancel_*`, `test_search_shared_by_two_visitors_*` |
| The disk fills up | A search starts only if the free space covers its estimated need plus `min_free_gb`. Otherwise it waits and the site stops taking searches. | `test_search_waits_when_the_disk_is_short` |
| The pipeline's software is broken (a missing package, conda gone) | The search fails with "the server could not prepare its software", the queue pauses for 30 minutes or until `resume`, waiting searches keep their place. | `tests/test_environment.py` |
| Something else's Python is in front of the pipeline's in `PATH` | `run_job.sh` removes virtualenvs from `PATH` and puts the conda environment first; `check` and the worker test what a task would really run. | `tests/test_environment.py` |
| One search makes the worker raise an error | After three errors in a row that search is failed; the others go on. | `test_one_bad_search_does_not_stop_the_others` |
| Two workers on one data folder | The second refuses to start. | `test_only_one_worker_per_data_folder` |
| A step of the pipeline fails for single genomes (it carries on) | The result page says how many steps failed. | `test_search_runs_to_a_result` |
| The connection to NCBI breaks during a genome download | The download goes on from the byte it stopped at, also after a restart of the worker. NCBI being down is tried again after 4 and 16 minutes, then the waiting search fails with that reason. | `tests/test_genomes.py`: `test_broken_connection_*`, `test_stopping_the_worker_*`, `test_ncbi_being_down_*` |
| A genome arrives damaged, or is not what its name says | It is compared with NCBI's checksum before it is unpacked, fetched once more, and refused the second time. Nothing of it stays in the store. | `test_damaged_file_*`, `test_file_that_is_not_what_its_name_says_*` |
| The genome store is full | The genomes used least in the last `store_score_days` go first. A genome that is pinned, or that a waiting or running search needs, is never removed; if nothing can go, the request is refused at once with the reason. | `test_least_used_genome_*`, `test_pinned_and_needed_*` |
| A genome is removed while a search that needs it waits | It is fetched again before the search starts. | `test_genome_removed_while_a_search_waited_*` |
| One visitor asks for many or huge genomes | Size is checked at NCBI before anything is downloaded; a limit per visitor and day; only assembly accessions with a version reach a download address. | `test_genome_larger_than_the_limits_*`, `test_one_visitor_*`, `test_what_a_request_may_and_may_not_name` |
| Two processes fetch the same genome (worker and `library fetch`) | The second finds it taken and waits. | `test_download_cut_off_by_a_restart_*` |
| The store's folder holds files that are in no record (the database was lost or replaced; files were copied in by hand) | They count against the store's size, so the folder as a whole stays within the space that was set. The store never deletes them by itself and never takes one for a genome: a genome is fetched and checked even when a file of its name lies there. `library status` and `check` name them; `library sweep --yes` deletes them. Files of genomes registered with `library add` are never among them. | `test_files_in_no_record_*`, `test_file_the_store_did_not_fetch_*`, `test_genomes_the_operator_registered_*` |
| Too many or too large requests | Limits per visitor (running searches, searches per day, look-ups per hour), on the queue, on genome number and size, on the request size. | `tests/test_api.py` |
| Someone asks for a file that is not theirs to have | Only files listed in a result are served; nothing outside `site/` is; search ids are 128 random bits and are not logged in full. | `tests/test_api.py` |
| A name or sequence from a visitor ends up in a page | Pages insert text as text only, load scripts from this site only, and figures made by the pipeline are shown in a sandbox. | `tests/test_api.py`; browser run-through |
| UniProt or NCBI do not answer | The form says so and offers pasting the sequence. Answers are kept for 30 days. | `tests/test_inputs.py` |
| The pipeline changes a log line or a report field the site reads | A test compares the site's expectations with `main.nf`, `bin/` and `generate_report.py`. Progress display degrades; a search never fails because of it. | `tests/test_results.py` ("what the web layer expects of the pipeline") |
| The page loses the connection | It keeps the last state, says so, and retries with growing pauses. | browser run-through |

The counts on a result page are computed from the report and compared with the
report's own headline numbers (`self_check` in `result.json`). If they ever differ,
the page shows a warning instead of a wrong number.

## The API

All under `/api/v1/`. Errors are `{"error": {"code", "message", "field"}}`; `message`
is written for the visitor and is shown as it is.

| Call | Purpose |
|---|---|
| `GET /status` | Is the service taking searches; limits; pipeline version; queue length. |
| `GET /library[?home=<id>]` | Genomes on the server; with `home`, the others ordered by relatedness (`?taxon=<NCBI id>` orders them for a species whose genome is not there). |
| `POST /lookup` `{"accession"}` or `{"sequence"}` | What the protein is and which genome on the server it belongs to. |
| `POST /genomes/find` `{"term", "purpose"?}` | What NCBI has for a species name or assembly accession: which is on the server already, which could be fetched. |
| `POST /searches` `{"query": {"accession"} \| {"sequence","name","home"}, "targets": [...], "email"?}` | Start a search. A genome (in `targets`, or as `home`) is its id on the server or `{"accession": "GCF_..."}`. Returns the search's id, state and a `token` for cancelling. |
| `GET /searches/<id>` | State, queue position, the genomes it waits for, progress, error. |
| `POST /searches/<id>/cancel` `{"token"}` | Withdraw. |
| `GET /searches/<id>/result` | The result (`result.json`) when done. |
| `GET /searches/<id>/files/<name>` | One figure or download listed in the result. |
| `GET /healthz` (not under `/api/v1/`) | `ok` when the database answers. For monitoring. |

The pages can be hosted elsewhere than the service: put the service's address in
`site/config.js` and the pages' address under `cors_origins` in the settings. With
`apiBase: null` the pages work from plain files and show the examples only.

## The genome library

`python -m synvoy_web library add|list|check|remove`. `add` counts the genome, takes a
checksum and asks NCBI Taxonomy for the lineage (that is what orders "closest relatives
first"; `--no-network` or `--lineage` skip the question). The files stay where they are;
the library file (`<data_dir>/library.json`) points to them. Give the service read
access only.

A genome with `--gff` can be a starting genome ("home"). Target genomes need no
annotation, and none is passed to the pipeline unless `use_target_annotations` is set.

The library is re-read when the file changes; no restart is needed.

## Genomes that are not on the server: the genome store

With `genome_store_gb` above 0, the server takes that much space for genomes it fetches
from NCBI itself (`genome_store_dir`, default `<data_dir>/genomes`). They come in two ways:

- **The operator fills part of it.** `python -m synvoy_web library fetch "Apis mellifera"
  GCA_928718305.1 --annotation` takes species names (NCBI's reference genome is used)
  or assembly accessions; `--from-file list.txt --fill-to 50` stops when the store is
  half full; `--pin` keeps a genome whatever happens; `--annotation` also gets the
  gene annotation, without which a genome cannot be a starting genome.
- **Visitors name what is missing.** The search form has "A species that is not
  listed". The search then waits ("Getting a genome first"), the worker downloads the
  genome, compares it with NCBI's checksum, unpacks it, and starts the search. The
  genome stays for later searches.

Every search that starts counts as one use of each of its genomes. When a new genome
needs room, the genomes with the fewest uses in the last `store_score_days` (90) days
are removed first, then the ones unused for longest. Never removed: a pinned genome,
and a genome that a waiting or running search needs. A genome of the store is known by
its accession with version, which names one fixed set of sequences, so it can be
removed and fetched again without changing a result; a result that used it stays valid.

`python -m synvoy_web library status` shows the space, every genome with its uses, and
which would go first. `library pin|unpin|evict ACCESSION` are the manual controls.
The genomes registered with `library add` are not part of the store and are never
removed. Limits a visitor meets: `max_genome_mb` per genome, `max_targets` genomes per
search, `max_fetches_per_day` new genomes per visitor.

What the store knows about its genomes is in the database. If the database is lost or
replaced while the folder stays, the files in it are in no record any more: they still
count against the store's size, `library status` lists them, and `library sweep --yes`
deletes them (a genome that is asked for again is fetched again, not taken from such a
file). Give the store a folder of its own; your own genomes can lie in it, but a
separate folder is clearer.

## Examples on the start page

```bash
python -m synvoy_web gallery build --spec my_examples.json
```

writes `site/data/`. An example is a pipeline output folder, or a search that was run
on the server (`"job": "<id>"`). See `deploy/gallery.example.json`. Nothing a visitor
does puts a search there; the operator chooses. The command warns when a file it
publishes contains a path of the machine the run was made on.

A run that was started under an internal label can be published under another name:
`"rename": {"label": "name"}` replaces it in every text file, and `"must_not_contain":
["label"]` makes the build stop, naming the file, if the word is still somewhere (also
inside the archives). The query is always published under the example's own name.

## Running it on a server

1. A user for the service, a checkout of this repository, the conda environment
   (`./install.sh`), the genomes.
2. `deploy/settings.example.json` → `/etc/synvoy-web/settings.json`. At least
   `data_dir`, `conda_base`, `trust_proxy: true`, `public_base_url`,
   `operator_name`, `contact_email`.
3. `deploy/server.config` → pipeline memory and cores for this machine. **Its numbers
   are examples that were not measured.** Then run a real search with the largest
   genomes you will offer (`tests/test_end_to_end.py`, below) and adjust.
4. `python -m synvoy_web library add ...` for every genome of your own, and, with
   `genome_store_gb` set, `python -m synvoy_web library fetch --from-file genomes.txt
   --fill-to 50` for the ones NCBI has. Then `python -m synvoy_web check`. It must end
   without `[FAIL]`.
5. The two units in `deploy/` (`synvoy-web-api`, `synvoy-web-worker`) and the
   `deploy/Caddyfile` (TLS, request size, logs without visitor addresses).
6. Fill in `site/imprint.html` and `site/privacy.html`. Both carry a notice for the
   operator that says what is missing. A public site run from Germany needs both.

The built-in HTTP server is enough behind a reverse proxy for a site of this size. Any
WSGI server works as well: `gunicorn --workers 2 --threads 4 synvoy_web.wsgi:application`.

### Day to day

```bash
python -m synvoy_web jobs list                 # searches and their states
python -m synvoy_web jobs show <id>            # one search, and where its folder is
python -m synvoy_web jobs cancel <id>
python -m synvoy_web pause "Back at 14:00."    # no new searches, none started; running ones finish
python -m synvoy_web resume
python -m synvoy_web library status            # the genome store: space, uses, what goes first
python -m synvoy_web library sweep             # files in the store's folder that are in no record
```

The worker writes what it does to its standard error (the journal, under systemd):
every start, end, pause and clean-up. Lines starting with `PAUSED` or `NOT READY` need
a person. The logs of one search are in `<data_dir>/jobs/<id>/logs/`.

Updating the pipeline: `pause`, wait until no search is running
(`jobs list --state running`), update the checkout, `resume`. A search that runs
while the code changes uses the workflow it started with and whatever scripts are on
disk when each step starts. The worker notices this when the search ends: the result
page then carries a warning and the result is not reused for anyone else. A search
that was interrupted and is started again after the update begins from scratch.
Results of older code are never handed out for new code.

Back-up: the settings file and `library.json`. Everything else can be made again; after
a loss of the database, see `library sweep` above for the genome store's folder.

## Tests

```bash
cd web
python -m pytest tests                          # about 300 tests, 20 s, no network, no pipeline
python tests/check_safeguards.py                # do the tests notice a removed safeguard? 1 min
SYNVOY_WEB_BROWSER=1 python -m pytest tests/test_browser.py -s      # real Chrome, about 1 min
SYNVOY_WEB_E2E_CONFIG=var/settings.json \
SYNVOY_WEB_E2E_QUERY=../local_data/demo/query/PGK1.faa \
SYNVOY_WEB_E2E_HOME=saccharomyces_cerevisiae \
    python -m pytest tests/test_end_to_end.py -s                    # the real pipeline, about 10 min
```

- The normal run uses a stand-in for the pipeline (`tests/fake_pipeline.py`) that
  behaves like it from the outside and can be told to fail, hang, crash or leave
  processes behind.
- The browser run-through (`tests/browser/run.mjs`, needs Node 22 and Chrome; name them
  with `SYNVOY_WEB_NODE` and `SYNVOY_WEB_CHROME` when they are not on the `PATH`) drives a
  real browser through the paths a visitor takes: search by accession to the result,
  the same search again, paste a sequence and cancel it, what the form refuses, broken
  links, the service going away and coming back, a species that has to be fetched
  from (a stand-in for) NCBI, a small screen, and the pages without a service.
  `SYNVOY_WEB_SHOTS=<folder>` also saves pictures of the pages.
- NCBI is a stand-in in all of these (`conftest.FakeNcbi`): it can be down, drop a
  connection part-way, ignore the request to resume, or send a damaged file.
- The end-to-end test runs the real pipeline on real genomes, cancels one search and
  kills another outright. Run it while nobody edits the pipeline: a search that is
  started again after a code change rightly reuses nothing, and the test then skips.
  On a small machine run it under a memory limit and never next to another pipeline:
  `systemd-run --user --unit=synvoy-web-e2e --collect -p MemoryMax=10G ... python -m
  pytest tests/test_end_to_end.py -s`.
- Nothing here needs `var/` or built examples: the tests and the browser run-through
  pass on a fresh checkout (the start page then says that there are no examples yet).

## What has been checked, and what has not

**With the real pipeline, on the code before the genome store was added** (laptop,
12 cores, 15 GB; the four yeast genomes of the demo; yeast STE2 as the query):

- A search through the site gives the same report as `./run_synvoy.sh` with the same
  inputs (the two `synvoy_report.json` are identical).
- Restarting the service during a search did not interrupt it (14 processes kept
  running; the new worker took the search over).
- Cancelling a running search ended all of its processes (16) within a second of the
  request and removed its working files.
- A search whose processes were all killed (16, with 13 steps finished) was started
  again by the worker and came to the right result. All 13 finished steps were reused.
- These runs were started from a shell with a Python virtualenv first in `PATH`, the
  state in which searches failed before `run_job.sh` cleaned the `PATH`.
- The pipeline's code was being edited while three of these searches ran. The worker
  noticed each time (its log; that the result is then flagged and kept out of the
  cache is checked with the stand-in pipeline only).
- The service under gunicorn (2 workers, 4 threads each), 40 requests at once.

Two faults were found only by these runs, not by the stand-in pipeline: the virtualenv
in `PATH`, and a query file that was written again for the second attempt, after which
Nextflow reused 5 of 13 finished steps instead of all. After a change to `runner.py`
or `run_job.sh`, run `tests/test_end_to_end.py`.

**On the present code** (with the genome store, the new result page and the archive):

- One search with the real pipeline through the site: yeast SYP1 against three yeasts.
  All three genomes are in every figure and in the table, also the one without a call.
- One genome fetched from the real NCBI by `library fetch`: the *S. cerevisiae*
  reference with its annotation, 12 MB + 12 MB, checksums compared, 9 seconds.
- The suite on Python 3.10 and 3.12; `tests/check_safeguards.py` (50 safeguards
  removed one at a time, each noticed by a test); the browser run-through (63 checks)
  in Chrome; all of these also in a copy that holds only the files under version
  control.

**Not repeated on the present code:** cancelling a running search and killing one
outright, with the real pipeline (`tests/test_end_to_end.py`). Both pass with the
stand-in pipeline. The worker has changed since the runs above (it settles the genomes
of a search before it starts, and fetches genomes in a thread of its own), so run that
test once on a machine where no other pipeline runs and nobody edits the pipeline,
before the service goes public.

Not checked:

- **A search on a genome that the store fetched, with the real pipeline.** The two
  halves are checked (the fetch against NCBI; the search on stored genomes with the
  stand-in pipeline), not the two together.
- **Genomes larger than yeast.** The limits (`max_genome_mb`, `max_total_mb`), the
  disk estimate (`disk_factor`), the space a fetched genome is counted with and
  `deploy/server.config` are starting values. A download of a 1 Gb genome, or the
  store deleting genomes because it is full, has only happened with small stand-in
  files.
- **NCBI under load.** How NCBI answers many metadata requests or several large
  downloads in a row from one server (its rate limits) was not tried.
- **More than one search at a time** (`max_running` above 1) with the real pipeline.
- **E-mail** against a real mail server. The message and the deletion of the address
  are tested; the sending is not.
- **Load.** Nobody has sent this service many requests at once from many addresses.
- The service units and the Caddyfile were written, not run.

Known limits:

- The gene tree is switched off by the `low_mem` profile, which is the default. The
  result page then has no "Gene tree" tab.
- The pipeline does not publish the proteins of the gene models it finds in a search
  (only those of the rescue steps), so the result page offers the gene models (GFF)
  and not a protein FASTA file.
- Progress inside the genome search ("3 of 10 genomes") moves in the steps the
  pipeline works in (waves of genomes), not genome by genome.
- A visitor is told apart by network address. Many people behind one address (a
  university network) share the limits.
