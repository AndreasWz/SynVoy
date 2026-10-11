// search.html — the form that starts a search.
//
// The server checks everything again; what is checked here only saves a round trip.
// Nothing is sent until "Start search" is pressed, except the look-up of an accession
// and, where the server fetches genomes, the question to NCBI for a species name.

import { findGenomes, getLibrary, lookup, startSearch } from "./api.js";
import { byId, el, params, replace, species } from "./dom.js";
import { remember } from "./recent.js";
import { initShell } from "./shell.js";
import * as words from "./words.js";

const state = {
  mode: "accession",
  found: null,          // answer of the look-up for the accession in the box
  library: null,        // {groups, genomes}
  byId: new Map(),
  order: [],            // target ids, closest relatives of the home genome first
  related: null,        // ids near enough to the home genome to be ticked by default
  selected: new Set(),  // ids of genomes on the server
  extra: new Map(),     // accession -> genome at NCBI that is to be fetched and searched
  group: "all",
  home: "",             // id of the starting genome on the server
  homeExtra: null,      // or: a genome at NCBI that is to be fetched as starting genome
  limits: null,
  downloads: false,     // does this server fetch genomes from NCBI?
  email: false,
  busy: false,
};

const $ = byId;
const LETTERS = /[A-Za-z]/g;

function setError(id, text) {
  const node = $(id);
  if (!node) return;
  node.textContent = text || "";
  node.hidden = !text;
  const field = { "acc-error": "acc", "seq-error": "seq", "home-error": "home", "email-error": "email",
    "ncbi-error": "ncbi-term", "home-ncbi-error": "home-term" }[id];
  if (field && $(field)) $(field).setAttribute("aria-invalid", text ? "true" : "false");
}

function clearErrors() {
  for (const id of ["acc-error", "seq-error", "home-error", "targets-error", "email-error", "ncbi-error", "home-ncbi-error"]) setError(id, "");
  $("form-error").hidden = true;
}

function plainOrganism(found) {
  return ((found && found.organism) || "").replace(/\s*\([^()]*\)/g, "").trim();
}

// ── 1. the gene ────────────────────────────────────────────────────────────

function setMode(mode) {
  state.mode = mode;
  $("mode-accession").setAttribute("aria-pressed", String(mode === "accession"));
  $("mode-sequence").setAttribute("aria-pressed", String(mode === "sequence"));
  $("accession-pane").hidden = mode !== "accession";
  $("sequence-pane").hidden = mode !== "sequence";
  clearErrors();
  refreshHomePane();
}

function homeOptions() {
  const select = $("home");
  const homes = state.library ? state.library.genomes.filter((g) => g.can_be_home) : [];
  replace(select, el("option", { value: "", text: homes.length ? "Choose a species…" : "No genome with an annotation on this server" }),
    homes.map((g) => el("option", { value: g.id, text: g.species + (g.common_name ? " (" + g.common_name + ")" : "") })));
  select.value = state.home && homes.some((g) => g.id === state.home) ? state.home : "";
}

function refreshHomePane() {
  // The starting genome is asked for when it cannot be told from the accession.
  const needed = state.mode === "sequence" || (state.found && !state.found.home);
  $("home-pane").hidden = !needed;
  $("home-ncbi").hidden = !(needed && state.downloads);
  $("home-hint").textContent = state.mode === "sequence"
    ? "The gene is first located in this genome, and its neighbours there are the landmarks."
    : "The genome of this species is not on the server. Choose the closest relative that is: the protein is located there by similarity."
      + (state.downloads ? " Or have the species’ own genome fetched, below." : "");
  showHomeExtra();
}

function showFound() {
  const found = state.found;
  const box = $("lookup-result");
  if (!found) {
    replace(box);
    return;
  }
  const from = found.organism
    ? [", from ", species(plainOrganism(found)), found.organism_common ? " (" + found.organism_common.toLowerCase() + ")" : ""]
    : [];
  const parts = [el("div", { class: "notice" },
    el("div", {}, el("strong", { text: found.gene ? found.gene + " · " + (found.name || found.accession) : (found.name || found.accession) }), from),
    el("div", { class: "small", text: found.length + " residues. " + (found.home
      ? "Its genome and gene annotation are on this server."
      : "Its genome is not on this server.") }))];
  if (!found.acceptable) {
    parts.push(el("div", { class: "notice error", text: found.problem }));
  } else if ((found.warnings || []).includes("short_protein")) {
    parts.push(el("div", { class: "notice warn", text: "This is a short protein. In our tests, proteins under 100 residues gave few confident results." }));
  }
  replace(box, parts);
}

async function doLookup() {
  const text = $("acc").value.trim();
  setError("acc-error", "");
  state.found = null;
  showFound();
  if (!text) {
    setError("acc-error", "Enter a protein accession.");
    return false;
  }
  $("lookup").disabled = true;
  replace($("lookup-result"), el("div", { class: "notice plain", text: "Looking up " + text + "…" }));
  try {
    const found = await lookup({ accession: text });
    if ($("acc").value.trim() !== text) return false;      // the box was changed meanwhile
    state.found = found;
    showFound();
    refreshHomePane();
    if (found.home) {
      state.homeExtra = null;
      await setHome(found.home.id);
    } else {
      if (!state.home && !state.homeExtra) showSpecies();
      if (state.downloads && plainOrganism(found)) {         // offer the species' own genome
        $("home-term").value = plainOrganism(found);
        find("home");
      }
    }
    return true;
  } catch (error) {
    replace($("lookup-result"));
    setError("acc-error", error.message);
    return false;
  } finally {
    $("lookup").disabled = false;
  }
}

function sequenceHint() {
  const text = $("seq").value;
  const body = text.split(/\r?\n/).filter((line) => !line.startsWith(">")).join("");
  const count = (body.match(LETTERS) || []).length;
  const limits = state.limits;
  let hint = count ? count + " residues." : "";
  if (limits && count && count < limits.min_query_length) hint += " The shortest accepted is " + limits.min_query_length + ".";
  if (limits && count > limits.max_query_length) hint += " The longest accepted is " + limits.max_query_length + ".";
  if (count && count < 100 && (!limits || count >= limits.min_query_length)) hint += " Short proteins gave few confident results in our tests.";
  $("seq-count").textContent = hint;
}

// ── genomes at NCBI ────────────────────────────────────────────────────────

const FIND = {
  target: { term: "ncbi-term", button: "ncbi-find", error: "ncbi-error", list: "ncbi-offers" },
  home: { term: "home-term", button: "home-find", error: "home-ncbi-error", list: "home-offers" },
};

function offerFacts(c) {
  return [c.accession, c.assembly, c.level, words.megabases(c.size_mb),
    c.reference ? "reference genome" : "", c.annotated ? "with gene annotation" : "no gene annotation"]
    .filter(Boolean).join(" · ");
}

function offerAction(kind, c) {
  const again = () => showOffers(kind, shownOffers[kind]);
  if (kind === "home") {
    if (c.on_server) {
      return el("button", { type: "button", class: "button secondary small", text: "Use it (on this server)",
        onclick: async () => { state.homeExtra = null; await setHome(c.on_server); again(); } });
    }
    if (!c.usable) return el("span", { class: "small", text: c.problem });
    if (state.homeExtra && state.homeExtra.accession === c.accession) return el("span", { class: "found small", text: "Starting genome" });
    return el("button", { type: "button", class: "button secondary small", text: "Use as starting genome",
      onclick: async () => { await setHomeExtra(c); again(); } });
  }
  if (c.on_server) {
    if (state.selected.has(c.on_server)) return el("span", { class: "found small", text: "Selected" });
    if ((state.home && c.on_server === state.home)) return el("span", { class: "small", text: "This is the starting genome" });
    return el("button", { type: "button", class: "button secondary small", text: "Select it (on this server)",
      onclick: () => { addListed(c.on_server); again(); } });
  }
  if (!c.usable) return el("span", { class: "small", text: c.problem });
  if (state.extra.has(c.accession)) return el("span", { class: "found small", text: "Added" });
  if (state.homeExtra && state.homeExtra.accession === c.accession) return el("span", { class: "small", text: "This is the starting genome" });
  return el("button", { type: "button", class: "button secondary small", text: "Add to the search",
    onclick: () => { addExtra(c); again(); } });
}

const shownOffers = { target: [], home: [] };

function showOffers(kind, candidates) {
  shownOffers[kind] = candidates;
  replace($(FIND[kind].list), candidates.map((c) => el("li", { class: "offer" + (c.usable || c.on_server ? "" : " unusable") },
    el("div", { class: "what" },
      el("div", {}, species(c.species), c.common_name ? " (" + c.common_name + ")" : ""),
      el("div", { class: "facts-line", text: offerFacts(c) })),
    offerAction(kind, c))));
}

async function find(kind) {
  const ids = FIND[kind];
  const term = $(ids.term).value.trim();
  setError(ids.error, "");
  if (!term) {
    setError(ids.error, "Type a species name or an assembly accession.");
    return;
  }
  $(ids.button).disabled = true;
  replace($(ids.list), el("li", { class: "offer unusable", text: "Asking NCBI for " + term + "…" }));
  try {
    const data = await findGenomes({ term, purpose: kind });
    if ($(ids.term).value.trim() !== term) return;            // changed meanwhile
    if (!data.candidates.length) {
      replace($(ids.list));
      setError(ids.error, "NCBI has no genome under this name. Check the spelling, or give the assembly accession (GCF_… or GCA_…).");
      return;
    }
    showOffers(kind, data.candidates);
  } catch (error) {
    replace($(ids.list));
    setError(ids.error, error.message);
  } finally {
    $(ids.button).disabled = false;
  }
}

function atLimit() {
  const limits = state.limits || { max_targets: 10 };
  return state.selected.size + state.extra.size >= limits.max_targets;
}

function addExtra(candidate) {
  if (atLimit()) {
    setError("ncbi-error", "The limit of genomes per search is reached. Take one off first.");
    return;
  }
  state.extra.set(candidate.accession, candidate);
  showSpecies();
}

function addListed(id) {
  if (!state.byId.has(id)) return;                          // listed since the page was loaded
  if (atLimit() && !state.selected.has(id)) {
    setError("ncbi-error", "The limit of genomes per search is reached. Take one off first.");
    return;
  }
  state.selected.add(id);
  showSpecies();
}

function showHomeExtra() {
  const chosen = state.homeExtra;
  const note = $("home-chosen");
  if (!note) return;
  note.hidden = !chosen;
  if (!chosen) return;
  replace(note, el("strong", { text: "Starting genome: " }), species(chosen.species), " (" + chosen.accession + "). ",
    "It is fetched from NCBI with its gene annotation before the search starts. ",
    el("button", { type: "button", class: "linklike", text: "Take another", onclick: () => {
      state.homeExtra = null; state.order = []; showHomeExtra(); showOffers("home", shownOffers.home); showSpecies();
    } }));
}

async function setHomeExtra(candidate) {
  state.homeExtra = candidate;
  state.home = "";
  $("home").value = "";
  state.extra.delete(candidate.accession);
  setError("home-error", "");
  showHomeExtra();
  try {
    const data = await getLibrary(null, candidate.taxon_id || "");
    if (state.homeExtra !== candidate) return;
    state.order = data.targets || state.library.genomes.filter((g) => g.available).map((g) => g.id);
    state.related = data.related ? new Set(data.related) : null;
    state.group = "all";
    selectClosest();
  } catch (error) {
    state.order = state.library.genomes.filter((g) => g.available).map((g) => g.id);
    state.related = new Set();                       // the order is not known: tick nothing
  }
  showSpecies();
}

// ── 2. the species ─────────────────────────────────────────────────────────

async function setHome(id) {
  state.home = id;
  $("home").value = id || "";
  setError("home-error", "");
  if (id) {
    state.homeExtra = null;
    showHomeExtra();
  }
  if (!id) {
    if (!state.homeExtra) {
      state.order = [];
      state.selected.clear();
    }
    showSpecies();
    return;
  }
  try {
    const data = await getLibrary(id);
    if (state.home !== id) return;
    state.order = data.targets || [];
    state.related = data.related ? new Set(data.related) : null;
    const home = state.byId.get(id);
    state.group = home && !fromNcbi(home) && operatorGroups().length > 1 ? home.group : "all";
    selectClosest();
  } catch (error) {
    setError("home-error", error.message);
    state.order = [];
  }
  showSpecies();
}

function selectable(genome) {
  return genome.available && (!state.limits || genome.size_mb <= state.limits.max_genome_mb);
}

function inGroup(genome) {
  return state.group === "all" || genome.group === state.group;
}

// Genomes fetched from NCBI belong to no group of the operator's. They are listed
// whatever group is chosen, at their place by relatedness.
const fromNcbi = (genome) => genome.origin === "ncbi";
const operatorGroups = () => state.library.groups.filter((g) => g.id !== "ncbi");

function selectClosest() {
  const limits = state.limits || { max_targets: 10, max_total_mb: Infinity };
  state.selected.clear();
  let total = 0;
  for (const c of state.extra.values()) total += c.size_mb;
  for (const id of state.order) {
    const genome = state.byId.get(id);
    if (!genome || !selectable(genome) || !(inGroup(genome) || fromNcbi(genome))) continue;
    // The operator's genomes were chosen and grouped by hand. Those fetched from NCBI
    // are whatever visitors asked for: tick only the ones near the starting genome.
    if (fromNcbi(genome) && state.related && !state.related.has(id)) continue;
    if (state.selected.size + state.extra.size >= limits.max_targets || total + genome.size_mb > limits.max_total_mb) continue;
    state.selected.add(id);
    total += genome.size_mb;
  }
}

function totals() {
  let total = 0;
  for (const id of state.selected) total += (state.byId.get(id) || { size_mb: 0 }).size_mb;
  for (const c of state.extra.values()) total += c.size_mb;
  return total;
}

function showSpecies() {
  const hasHome = !!state.home || !!state.homeExtra;
  const ready = hasHome && (state.order.length > 0 || state.extra.size > 0 || state.downloads);
  $("species-empty").hidden = ready;
  $("species-pane").hidden = !ready;
  $("ncbi-pane").hidden = !(ready && state.downloads);
  if (!ready) {
    $("species-empty").textContent = hasHome
      ? "There is no other genome on this server to search."
      : "Choose the gene first. The species are then listed with the closest relatives on top.";
    return;
  }
  const limits = state.limits || { max_targets: 10, max_total_mb: Infinity, max_genome_mb: Infinity };
  const groups = operatorGroups();
  const sizes = (list) => {
    if (!list.length) return "";
    const low = Math.min(...list.map((g) => g.size_mb)), high = Math.max(...list.map((g) => g.size_mb));
    return low === high ? words.megabases(low) : words.megabases(low) + " to " + words.megabases(high);
  };
  const targets = state.order.map((id) => state.byId.get(id)).filter(Boolean);
  $("group-pane").hidden = groups.length < 2;
  if (groups.length > 1) {
    const choice = (value, name, list) => el("label", { class: "choice" },
      el("input", { type: "radio", name: "group", value, checked: state.group === value, onchange: () => { state.group = value; selectClosest(); showSpecies(); } }),
      el("span", {}, el("span", { class: "name", text: name }),
        el("span", { class: "note", text: list.length + " " + words.plural(list.length, "genome", "genomes") + (list.length ? ", " + sizes(list) : "") })));
    replace($("groups"), choice("all", "All groups", targets),
      groups.map((g) => choice(g.id, g.name, targets.filter((t) => t.group === g.id))));
  }
  const full = atLimit();
  const shown = targets.filter((genome) => inGroup(genome) || fromNcbi(genome));
  // Genomes that are to be fetched come first: they are part of this search only.
  const fetched = [...state.extra.values()].map((c) => el("label", { class: "check" },
    el("input", { type: "checkbox", value: c.accession, checked: true, onchange: () => {
      state.extra.delete(c.accession);
      showSpecies();
      showOffers("target", shownOffers.target);
    } }),
    el("span", {}, species(c.species), " ", el("span", { class: "tag", text: "from NCBI" }), " ",
      el("span", { class: "meta", text: words.megabases(c.size_mb) + ", fetched first" }))));
  replace($("species-list"), fetched, shown.map((genome) => {
    const checked = state.selected.has(genome.id);
    const usable = selectable(genome);
    const disabled = !usable || (!checked && full);
    const why = !genome.available ? "not available at the moment"
      : !usable ? "above the limit of " + words.megabases(limits.max_genome_mb) + " per genome" : "";
    return el("label", { class: "check" + (disabled ? " disabled" : ""), title: why || null },
      el("input", { type: "checkbox", value: genome.id, checked, disabled, onchange: (event) => {
        if (event.target.checked) state.selected.add(genome.id); else state.selected.delete(genome.id);
        showSpecies();
      } }),
      el("span", {}, species(genome.species), " ",
        fromNcbi(genome) ? [el("span", { class: "tag", text: "from NCBI" }), " "] : null,
        el("span", { class: "meta", text: words.megabases(genome.size_mb) + (why ? ", " + why : "") })));
  }));
  const total = totals();
  const count = state.selected.size + state.extra.size;
  const over = total > limits.max_total_mb;
  $("count").textContent = count + " of " + limits.max_targets + " genomes selected";
  $("size").textContent = words.megabases(total) + " in total" + (Number.isFinite(limits.max_total_mb) ? " (limit " + words.megabases(limits.max_total_mb) + ")" : "");
  $("meter").className = "meter" + (over ? " over" : "");
  $("meter-fill").style.width = Math.min(100, 100 * count / limits.max_targets) + "%";
  setError("targets-error", over ? "The chosen genomes add up to more than the limit per search. Take one off." : "");
  $("limit-note").hidden = !full;
}

// ── 3. sending ─────────────────────────────────────────────────────────────

function fieldFor(error) {
  return { accession: "acc-error", sequence: "seq-error", name: "seq-error", home: "home-error",
    targets: "targets-error", email: "email-error" }[error.field] || null;
}

function homeItem() {
  return state.homeExtra ? { accession: state.homeExtra.accession } : state.home;
}

async function submit(event) {
  event.preventDefault();
  if (state.busy) return;
  clearErrors();
  const query = {};
  if (state.mode === "accession") {
    const text = $("acc").value.trim();
    if (!state.found || state.found.accession.toUpperCase().split(".")[0] !== text.toUpperCase().split(".")[0]) {
      if (!(await doLookup())) { $("acc").focus(); return; }
    }
    if (!state.found.acceptable) { $("acc").focus(); return; }
    query.accession = state.found.accession;
    if (!state.found.home) {
      if (!homeItem()) { setError("home-error", "Choose the starting genome."); $("home").focus(); return; }
      query.home = homeItem();
    }
  } else {
    if (!$("seq").value.trim()) { setError("seq-error", "Paste a protein sequence."); $("seq").focus(); return; }
    if (!homeItem()) { setError("home-error", "Choose the species the protein comes from."); $("home").focus(); return; }
    query.sequence = $("seq").value;
    query.home = homeItem();
    if ($("seq-name").value.trim()) query.name = $("seq-name").value.trim();
  }
  const targets = state.order.filter((id) => state.selected.has(id))
    .concat([...state.extra.keys()].map((accession) => ({ accession })));
  if (!targets.length) {
    setError("targets-error", "Choose at least one genome to search.");
    $("species-pane").scrollIntoView({ block: "center" });
    return;
  }
  const body = { query, targets };
  const mail = state.email ? $("email").value.trim() : "";
  if (mail) body.email = mail;

  state.busy = true;
  const button = $("submit");
  button.disabled = true;
  button.textContent = "Starting…";
  try {
    const job = await startSearch(body);
    remember({ id: job.id, token: job.token, title: words.title(job.query, job.targets.length, "") });
    window.location.assign("result.html?id=" + encodeURIComponent(job.id));
  } catch (error) {
    const target = fieldFor(error);
    if (target) {
      setError(target, error.message);
      const node = $(target);
      if (node && node.scrollIntoView) node.scrollIntoView({ block: "center" });
    } else {
      $("form-error").textContent = error.message;
      $("form-error").hidden = false;
    }
    state.busy = false;
    button.disabled = false;
    button.textContent = "Start search";
  }
}

// ── start ──────────────────────────────────────────────────────────────────

function showLimits(status) {
  const l = status.limits;
  replace($("limits"),
    el("li", { text: "Up to " + l.max_targets + " genomes per search" }),
    el("li", { text: "Genomes up to " + words.megabases(l.max_genome_mb) + " each, " + words.megabases(l.max_total_mb) + " together" }),
    el("li", { text: "Proteins of " + l.min_query_length + " to " + words.number(l.max_query_length) + " residues" }),
    el("li", { text: l.max_active_per_client + " " + words.plural(l.max_active_per_client, "search", "searches") + " at a time per visitor, " + l.max_submissions_per_day + " per day" }),
    el("li", { text: "A search is stopped after " + l.max_runtime_hours + " hours" }),
    el("li", { text: l.downloads
      ? "Genomes that are not here are fetched from NCBI, up to " + l.max_fetches_per_day + " per visitor and day"
      : "Only genomes kept on this server" }));
  $("retention").textContent = "Results are kept for " + l.retention_days + " days. The same search started again gives the same result.";
  const ahead = status.queue.waiting + status.queue.running;
  $("queue-note").textContent = ahead === 0 ? "No search is running; yours starts at once."
    : ahead + " " + words.plural(ahead, "search is", "searches are") + " ahead of yours.";
}

function unavailable(text) {
  $("search-form").hidden = true;
  const box = $("unavailable");
  replace(box, el("strong", { text: "Searches cannot be started right now. " }), text);
  box.hidden = false;
}

async function main() {
  const service = await initShell();
  if (!service.live) {
    unavailable(service.message + " The examples on the start page are real results.");
    return;
  }
  const status = service.data;
  state.limits = status.limits;
  state.downloads = !!status.limits.downloads;
  state.email = !!(status.features && status.features.email);
  $("email-pane").hidden = !state.email;
  showLimits(status);
  if (!status.accepting) {
    // The button stays usable: a search that has been run before is still answered.
    $("queue-note").textContent = (status.message || "No new searches at the moment.")
      + " A search that was run before is still answered.";
  }
  try {
    state.library = await getLibrary();
  } catch (error) {
    unavailable(error.message);
    return;
  }
  state.byId = new Map(state.library.genomes.map((g) => [g.id, g]));
  if (!state.downloads && !state.library.genomes.some((g) => g.can_be_home)) {
    unavailable("This server has no genomes to search yet.");
    return;
  }
  homeOptions();
  $("mode-accession").addEventListener("click", () => setMode("accession"));
  $("mode-sequence").addEventListener("click", () => setMode("sequence"));
  $("lookup").addEventListener("click", doLookup);
  $("acc").addEventListener("keydown", (event) => {
    if (event.key === "Enter") { event.preventDefault(); doLookup(); }
  });
  $("acc").addEventListener("input", () => {
    if (state.found) { state.found = null; showFound(); refreshHomePane(); }
    setError("acc-error", "");
  });
  $("seq").addEventListener("input", sequenceHint);
  $("home").addEventListener("change", (event) => setHome(event.target.value));
  for (const kind of ["target", "home"]) {
    $(FIND[kind].button).addEventListener("click", () => find(kind));
    $(FIND[kind].term).addEventListener("keydown", (event) => {
      if (event.key === "Enter") { event.preventDefault(); find(kind); }
    });
    $(FIND[kind].term).addEventListener("input", () => setError(FIND[kind].error, ""));
  }
  $("pick-closest").addEventListener("click", () => { selectClosest(); showSpecies(); });
  $("pick-none").addEventListener("click", () => {
    state.selected.clear(); state.extra.clear(); showSpecies(); showOffers("target", shownOffers.target);
  });
  $("search-form").addEventListener("submit", submit);
  showSpecies();
  const given = (params().get("acc") || "").trim().slice(0, 20);
  if (given) {
    $("acc").value = given;
    doLookup();
  }
}

main().catch((error) => {
  unavailable("This page could not be set up. Reload the page to try again.");
  console.error(error);
});
