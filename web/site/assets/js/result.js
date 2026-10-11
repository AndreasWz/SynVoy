// result.html — one address for the whole life of a search:
//   ?id=<search id>     waiting, running (updates by itself), then the result
//   ?example=<name>     a finished example kept with the pages

import { ApiError, cancelSearch, getResult, getSearch, getStatic, hasService, searchUrl } from "./api.js";
import { byId, copyText, el, icon, params, replace, species } from "./dom.js";
import { forget, remember, tokenFor } from "./recent.js";
import { renderResult } from "./render.js";
import { initShell } from "./shell.js";
import * as words from "./words.js";

const root = byId("view");
const ID_RE = /^[A-Za-z0-9][A-Za-z0-9_-]{21}$/;
const SLUG_RE = /^[a-z0-9][a-z0-9_-]{0,40}$/;
const POLL_MS = 5000;
const MAX_WAIT_MS = 60000;

function message(title, text, extra) {
  document.title = title + " · SynVoy";
  replace(root, el("section", { class: "stack" },
    el("h1", { class: "page-title", text: title }),
    el("p", { class: "lead", text: text }),
    extra || el("div", { class: "row" },
      hasService ? el("a", { class: "button", href: "search.html", text: "Start a search" }) : null,
      el("a", { class: "button secondary", href: "index.html#examples", text: "See the examples" }))));
}

// ── examples ───────────────────────────────────────────────────────────────

async function showExample(slug) {
  if (!SLUG_RE.test(slug)) {
    message("Example not found", "There is no example with this name.");
    return;
  }
  const base = "data/examples/" + slug + "/";
  const result = await getStatic(base + "result.json");
  if (!result || !result.summary) {
    message("Example not found", "This example is not on the site, or it could not be loaded. Reload the page to try again.");
    return;
  }
  renderResult(root, result, (file) => base + file, { example: true, hasService });
}

// ── a search on the server ─────────────────────────────────────────────────

function linkBox() {
  const url = window.location.href;
  const done = el("span", { class: "small muted", role: "status" });
  return el("div", { class: "card linkbox" },
    el("div", { class: "stack-s" },
      el("div", { class: "smaller muted", text: "Your link. Keep it to come back later; anyone who has it can see the result." }),
      el("div", { class: "url", text: url })),
    el("div", { class: "row" }, done,
      el("button", { type: "button", class: "button secondary", id: "copy-link", text: "Copy link", onclick: async () => {
        done.textContent = (await copyText(url)) ? "Copied." : "Select the link and copy it by hand.";
      } })));
}

function elapsed(job) {
  if (!job.started_at) return "";
  const seconds = (Date.now() - new Date(job.started_at).getTime()) / 1000;
  return seconds > 0 ? "Running for " + words.duration(seconds) + "." : "";
}

function targetsLine(job) {
  const out = ["Searching in: "];
  job.targets.forEach((t, i) => { if (i) out.push(", "); out.push(species(t.species)); });
  out.push(". Home genome: ", species(job.home.species || ""), ".");
  return el("p", { class: "small muted" }, out);
}

function stepsList(job) {
  const steps = (job.progress && job.progress.steps) || [];
  return el("ol", { class: "steps" }, steps.map((step) => {
    const text = words.stepText(step, job);
    const detail = step.detail || {};
    const ticks = step.key === "search" && step.state === "active" && detail.genomes_total
      ? el("div", { class: "ticks", role: "img", "aria-label": detail.genomes_done + " of " + detail.genomes_total + " genomes searched" },
        Array.from({ length: detail.genomes_total }, (_, i) => el("span", { class: i < detail.genomes_done ? "on" : "" })))
      : null;
    return el("li", { class: step.state },
      icon(step.state),
      el("div", { class: "stack-s" },
        el("div", { class: "title", text: text.title }),
        text.detail ? el("div", { class: "detail", text: text.detail }) : null,
        ticks));
  }));
}

function cancelControl(job, refresh) {
  const token = tokenFor(job.id);
  if (!token) return null;
  const note = el("span", { class: "small", role: "status" });
  const button = el("button", { type: "button", class: "button danger", id: "cancel-search", text: "Cancel this search" });
  button.addEventListener("click", async () => {
    if (!window.confirm("Stop this search? What it has computed so far is discarded.")) return;
    button.disabled = true;
    note.textContent = "Cancelling…";
    try {
      const outcome = await cancelSearch(job.id, token);
      if (!outcome.cancelled) {
        note.textContent = "Someone else is waiting for the same search, so it keeps running.";
      }
      refresh();
    } catch (error) {
      note.textContent = error.message;
      button.disabled = false;
    }
  });
  return el("div", { class: "row" }, button, note);
}

// What is on the screen for a waiting or running search. The page asks the server
// every few seconds; it is drawn again only when something a visitor can see has
// changed. Drawing it again each time would take the keyboard focus off the button
// someone is about to press, and make a screen reader start over.
let shown = { signature: null, clock: null };

function runningNote(job) {
  return [elapsed(job), job.attempt > 1 ? "The server was restarted during this search; it was picked up again where it stopped." : ""].filter(Boolean).join(" ");
}

function showActive(job, refresh, connection) {
  const queued = job.state === "queued";
  const name = words.title(job.query, job.targets.length, "");
  const signature = JSON.stringify([job.state, job.queue, job.progress, job.attempt, connection]);
  if (signature === shown.signature && shown.clock && document.body.contains(shown.clock)) {
    shown.clock.textContent = runningNote(job);
    return;
  }
  const focused = document.activeElement ? document.activeElement.id : "";
  const clock = el("p", { class: "small muted", text: runningNote(job) });
  shown = { signature, clock };
  document.title = (queued ? "Waiting: " : "Running: ") + name + " · SynVoy";
  let body;
  if (queued) {
    const ahead = (job.queue ? job.queue.ahead + job.queue.running : 0);
    const coming = (job.queue && job.queue.genomes) || [];
    body = el("section", { class: "card stack" },
      el("h2", { text: coming.length ? "Getting " + words.plural(coming.length, "a genome", "genomes") + " first" : "Waiting for its turn" }),
      coming.length ? el("div", { class: "stack" },
        el("p", { text: "This search needs " + words.plural(coming.length, "a genome that is", "genomes that are") + " not on the server yet. " + words.plural(coming.length, "It is", "They are") + " fetched from NCBI and checked; the search starts by itself after that." }),
        el("ul", { class: "offers", id: "coming-genomes" }, coming.map((genome) => el("li", { class: "offer" },
          el("div", { class: "what" }, species(genome.species), " ", el("span", { class: "facts-line", text: genome.accession }),
            el("div", { class: "fetching" },
              el("span", { class: "small muted", text: words.fetchText(genome) }),
              genome.state === "fetching" && genome.total_bytes > 0 ? el("div", { class: "meter", role: "img", "aria-label": words.fetchText(genome) },
                el("span", { css: { width: Math.min(100, Math.round(100 * genome.done_bytes / genome.total_bytes)) + "%" } })) : null))))))
        : el("p", { text: ahead === 0 ? "Your search is next and starts within a few seconds."
          : ahead + " " + words.plural(ahead, "search is", "searches are") + " ahead of yours. Searches are run in turn." }),
      el("p", { class: "small muted", text: "A search takes from a few minutes to several hours, depending on the number and size of the genomes." }));
  } else {
    body = el("section", { class: "card stack" },
      el("h2", { text: "Progress" }),
      stepsList(job),
      clock);
  }
  replace(root,
    el("section", { class: "stack" },
      el("div", { class: "row" },
        el("h1", { class: "page-title", text: name }),
        el("span", { class: "badge " + job.state, text: queued ? "Waiting" : "Running" })),
      linkBox(),
      el("p", { class: "muted", text: "You can close this page. The search keeps going, and this page updates by itself." }),
      connection ? el("div", { class: "notice warn", role: "status", text: connection }) : null),
    el("div", { class: "columns" },
      el("div", { class: "main stack" }, body),
      el("div", { class: "side stack" },
        el("div", { class: "card tight stack-s" },
          el("h2", { text: "This search" }),
          el("p", { class: "small" }, el("strong", { text: words.queryName(job.query) }),
            job.query.length ? ", " + job.query.length + " residues" : "",
            job.query.accession ? " (" + job.query.accession + ")" : ""),
          targetsLine(job),
          cancelControl(job, refresh)))));
  if (focused) {
    const again = document.getElementById(focused);
    if (again) again.focus();
  }
}

function showEnded(job) {
  const name = words.title(job.query, job.targets.length, "");
  const again = el("a", { class: "button", href: "search.html" + (job.query.accession ? "?acc=" + encodeURIComponent(job.query.accession) : ""), text: "Start a new search" });
  if (job.state === "failed") {
    const reference = el("p", { class: "small muted", text: "If you ask the operator about it, give this reference: " + job.id.slice(0, 8) });
    message("The search could not be completed", job.error ? job.error.message : "The search stopped with an error.",
      el("div", { class: "stack" }, el("p", { class: "muted", text: name + "." }), reference, el("div", { class: "row" }, again)));
  } else if (job.state === "cancelled") {
    forget(job.id);
    message("This search was cancelled", name + ". Nothing was kept.", el("div", { class: "row" }, again));
  } else if (job.state === "expired") {
    forget(job.id);
    message("This result is no longer kept", name + ". Results are deleted after a while to make room. The same search gives the same result again.", el("div", { class: "row" }, again));
  }
}

async function followSearch(id) {
  if (!ID_RE.test(id)) {
    message("Search not found", "This link is not complete. Check that you copied all of it.");
    return;
  }
  if (!hasService) {
    message("Search not found", "This copy of the site has no search service, so it has no searches.");
    return;
  }
  let wait = POLL_MS;
  let timer = null;
  let last = null;
  let stopped = false;
  let asked = false;

  const schedule = (ms) => {
    clearTimeout(timer);
    if (!stopped) timer = setTimeout(tick, ms);
  };
  // After something the visitor did (cancel), look at once, whatever the tab's state.
  const refresh = () => { asked = true; schedule(0); };

  async function tick() {
    if (document.hidden && !asked) {  // do not load the server from a tab nobody looks at
      schedule(POLL_MS);
      return;
    }
    asked = false;
    let job;
    try {
      job = await getSearch(id);
    } catch (error) {
      if (error instanceof ApiError && error.status === 404) {
        stopped = true;
        forget(id);
        message("Search not found", "There is no search with this link. It may have been deleted, or the link is not complete.");
        return;
      }
      wait = Math.min(wait * 2, MAX_WAIT_MS);
      const note = (error.message || "The server did not answer.") + " Trying again in " + Math.round(wait / 1000) + " s.";
      if (last) showActive(last, refresh, note);
      else message("Loading the search…", note, el("div", {}));
      schedule(wait);
      return;
    }
    wait = POLL_MS;
    if (job.state === "queued" || job.state === "running") {
      last = job;
      showActive(job, refresh, "");
      schedule(POLL_MS);
      return;
    }
    stopped = true;
    if (job.state !== "done") {
      showEnded(job);
      return;
    }
    try {
      const result = await getResult(id);
      remember({ id, token: tokenFor(id), title: words.title(result.query, result.summary.genomes, result.scope) });
      renderResult(root, result, (file) => searchUrl(id, "/" + file), { hasService, expiresAt: job.expires_at });
    } catch (error) {
      if (error instanceof ApiError && error.temporary) {
        stopped = false;
        message("Loading the result…", error.message, el("div", {}));
        schedule(POLL_MS);
      } else {
        message("The result cannot be shown", error.message);
      }
    }
  }

  document.addEventListener("visibilitychange", () => { if (!document.hidden && !stopped) schedule(0); });
  tick();
}

async function main() {
  initShell();
  const query = params();
  if (query.get("example")) await showExample(query.get("example"));
  else if (query.get("id")) await followSearch(query.get("id"));
  else message("No search given", "This page shows one search. Open it from the link you got when you started the search, or pick an example.");
}

main().catch((error) => {
  message("Something went wrong on this page", "Reload the page to try again.");
  console.error(error);
});
