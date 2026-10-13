// Draws a finished result. Used for a search run on the server and for the examples;
// the only difference between the two is where the files are (``fileUrl``).

import { copyText, el, replace, species } from "./dom.js";
import * as words from "./words.js";

function where(call) {
  if (!call.chrom) return "";
  return call.chrom + ":" + words.number(call.start) + "-" + words.number(call.end);
}

function speciesList(genomes) {
  const out = [];
  genomes.forEach((genome, index) => {
    if (index) out.push(", ");
    out.push(species(genome.species));
    if (genome.low_quality) out.push("†");
  });
  return out;
}

function download(result, key) {
  return (result.downloads || []).find((item) => item.key === key) || null;
}

// ── top ────────────────────────────────────────────────────────────────────

function intro(result, fileUrl, options) {
  const s = result.summary;
  const q = result.query || {};
  const bits = [options.example ? "Example" : "Finished", words.software(result.pipeline)];
  if (result.timing && result.timing.runtime_s) bits.push("took " + words.duration(result.timing.runtime_s));
  const about = [words.queryLongName(q), q.about, q.accession,
    q.length ? q.length + " residues" : ""].filter(Boolean);
  const all = download(result, "all");
  const copied = el("span", { class: "small muted", role: "status" });
  return el("section", { class: "result-head" },
    el("div", { class: "stack" },
      el("div", { class: "eyebrow", text: bits.join(" · ") }),
      el("h1", { class: "page-title", text: words.title(q, s.genomes, result.scope) }),
      about.length || q.organism ? el("p", { class: "muted" }, about.join(" · "),
        q.organism ? [about.length ? " · from " : "From ", species(q.organism.replace(/\s*\([^()]*\)/g, ""))] : null) : null,
      el("p", { class: "lead-strong", text: words.headline(result) })),
    el("div", { class: "actions" },
      all ? el("a", { class: "button", id: "download-all", href: fileUrl(all.file), download: all.file.split("/").pop(),
        text: "Download everything (" + words.bytes(all.bytes) + ")" }) : null,
      el("button", { type: "button", class: "button secondary", id: "copy-result-link", text: "Copy link", onclick: async () => {
        copied.textContent = (await copyText(window.location.href)) ? "Copied." : "Copy the address from the address bar.";
      } }),
      copied));
}

// Links to the parts of the page. The bar stays under the header while the page scrolls.
function subnav(parts) {
  const links = parts.map(([id, text]) => el("a", { href: "#" + id, text }));
  const bar = el("div", { class: "subnav" }, el("div", { class: "inner" },
    el("nav", { "aria-label": "On this page" }, links)));
  if ("IntersectionObserver" in window) {
    const seen = new Map();
    const observer = new IntersectionObserver((entries) => {
      for (const entry of entries) seen.set(entry.target.id, entry.isIntersecting);
      const current = parts.map(([id]) => id).find((id) => seen.get(id));
      for (const link of links) {
        if (link.getAttribute("href") === "#" + current) link.setAttribute("aria-current", "true");
        else link.removeAttribute("aria-current");
      }
    }, { rootMargin: "-120px 0px -55% 0px" });
    // The sections exist once the page has been put together.
    setTimeout(() => parts.forEach(([id]) => {
      const node = document.getElementById(id);
      if (node) observer.observe(node);
    }), 0);
  }
  return bar;
}

function notes(result) {
  const s = result.summary;
  const out = [];
  const low = result.genomes.filter((g) => g.low_quality);
  if (low.length) {
    out.push(el("div", { class: "notice warn" },
      el("strong", { text: "† Fragmented assembly: " }), speciesList(low.map((g) => ({ species: g.species }))),
      ". These genomes failed the quality check and were searched anyway. A gene can be missing there only because the assembly is broken at that place."));
  }
  if (s.failed_tasks > 0) {
    out.push(el("div", { class: "notice warn" },
      el("strong", { text: s.failed_tasks + " " + words.plural(s.failed_tasks, "checking step", "checking steps") + " failed. " }),
      "The search went on without " + words.plural(s.failed_tasks, "it", "them") + ", so the result for single genomes can be incomplete."));
  }
  if (result.self_check && result.self_check.agrees === false) {
    out.push(el("div", { class: "notice error" },
      el("strong", { text: "The numbers on this page do not match the report file. " }),
      "Use the report file (under “Take it with you”) and tell the operator of this site."));
  }
  for (const note of result.notes || []) {
    if (note && note.text) out.push(el("div", { class: "notice warn", text: note.text }));
  }
  return out;
}

// ── figures ────────────────────────────────────────────────────────────────

// A figure at its own size in a window over the page, for figures with many genomes.
function viewer() {
  const image = el("img", { alt: "" });
  const title = el("strong", {});
  let natural = 0, zoom = 1;
  const apply = () => { if (natural) image.style.width = Math.round(natural * zoom) + "px"; };
  const dialog = el("dialog", { class: "viewer", "aria-label": "Figure, enlarged" },
    el("div", { class: "bar-top" }, title,
      el("div", { class: "row" },
        el("button", { type: "button", class: "button secondary small", text: "Smaller", onclick: () => { zoom = Math.max(0.4, zoom / 1.25); apply(); } }),
        el("button", { type: "button", class: "button secondary small", text: "Larger", onclick: () => { zoom = Math.min(5, zoom * 1.25); apply(); } }),
        el("button", { type: "button", class: "button small", text: "Close", onclick: () => dialog.close() }))),
    el("div", { class: "scroll" }, image));
  image.addEventListener("load", () => { natural = image.naturalWidth || 1200; apply(); });
  dialog.addEventListener("click", (event) => { if (event.target === dialog) dialog.close(); });
  return {
    node: dialog,
    open(src, alt, name) {
      zoom = 1; natural = 0;
      image.style.width = "";
      image.alt = alt;
      image.src = src;
      title.textContent = name;
      if (dialog.showModal) dialog.showModal(); else window.open(src, "_blank", "noopener");
    },
  };
}

function figures(result, fileUrl) {
  const loci = (result.loci || []).filter((l) => l.files && (l.files.grid || l.files.ribbons || l.files.positions || l.files.tree_page));
  if (!loci.length) return null;
  const state = { locus: 0, tab: null };
  const tabsFor = (locus) => {
    const f = locus.files;
    return [["grid", f.grid], ["ribbons", f.ribbons], ["positions", f.positions], ["tree", f.tree_page]]
      .filter((pair) => pair[1]).map((pair) => pair[0]);
  };
  const tabBar = el("div", { class: "segmented", role: "group", "aria-label": "Figure" });
  const locusBar = el("div", { class: "segmented", role: "group", "aria-label": "Place in the home genome" });
  const frame = el("div", { class: "figure-frame" });
  const below = el("div", { class: "stack-s" });
  const large = viewer();
  const enlarge = el("button", { type: "button", class: "button secondary small", id: "enlarge-figure", text: "Enlarge" });

  function draw() {
    const locus = loci[state.locus];
    const tabs = tabsFor(locus);
    if (!tabs.includes(state.tab)) state.tab = tabs[0];
    replace(tabBar, tabs.map((key) => el("button", {
      type: "button", "aria-pressed": key === state.tab ? "true" : "false", text: words.FIGURES[key].tab,
      onclick: () => { state.tab = key; draw(); },
    })));
    replace(locusBar, loci.length < 2 ? [] : loci.map((l, index) => el("button", {
      type: "button", "aria-pressed": index === state.locus ? "true" : "false", text: "Place " + l.n,
      onclick: () => { state.locus = index; draw(); },
    })));
    locusBar.hidden = loci.length < 2;
    const info = words.FIGURES[state.tab];
    enlarge.hidden = state.tab === "tree";
    if (state.tab === "tree") {
      replace(frame, el("iframe", { src: fileUrl(locus.files.tree_page), title: info.alt, sandbox: "", loading: "lazy" }));
    } else {
      const src = fileUrl(locus.files[state.tab]);
      replace(frame, el("img", { src, alt: info.alt, loading: "lazy" }));
      enlarge.onclick = () => large.open(src, info.alt, info.tab);
    }
    const interactive = state.tab === "grid" ? locus.files.grid_page : state.tab === "ribbons" ? locus.files.ribbons_page : null;
    replace(below,
      el("p", { class: "small", text: info.how }),
      loci.length > 1 ? el("p", { class: "small muted", text: "The gene has " + loci.length + " places in its home genome. Each has its own figures." + (locus.chrom ? " Shown: " + locus.chrom + ":" + words.number(locus.start) + "-" + words.number(locus.end) + "." : "") }) : null,
      interactive ? el("p", { class: "small" }, el("a", { href: fileUrl(interactive), target: "_blank", rel: "noopener noreferrer", text: "Open the interactive version in a new tab" }), " (hover for details, zoom).") : null);
  }
  draw();
  return el("section", { class: "card stack", id: "figures" },
    el("div", { class: "row-between" },
      el("h2", { text: "The neighbourhood, species by species" }),
      el("div", { class: "row" }, locusBar, tabBar, enlarge)),
    frame, below, large.node);
}

// ── table ──────────────────────────────────────────────────────────────────

const RANK = { found: 0, probable: 1, candidate: 2, set_aside: 3 };

function counted(genome) {
  return genome.calls.filter((c) => c.label !== "set_aside");
}

// Genomes for the table: the clearest result first, then by similarity; genomes
// without a call follow, those with the most of the neighbourhood first.
function ordered(genomes) {
  const key = (genome) => {
    const calls = counted(genome);
    if (genome.status === "not_searched") return [9, 0];
    if (!calls.length) return [8, -((genome.neighbours || {}).together || 0)];
    const best = calls.reduce((a, b) => (RANK[b.label] < RANK[a.label]
      || (RANK[b.label] === RANK[a.label] && (b.identity || 0) > (a.identity || 0)) ? b : a));
    return [RANK[best.label], -(best.identity || 0)];
  };
  return genomes.map((genome, index) => ({ genome, index, key: key(genome) }))
    .sort((a, b) => a.key[0] - b.key[0] || a.key[1] - b.key[1] || a.index - b.index)
    .map((entry) => entry.genome);
}

function pill(label) {
  return el("span", { class: "pill p-" + label, text: words.LABEL[label] || label });
}

function neighbourCell(genome, rows) {
  const n = genome.neighbours;
  if (!n || !n.total) return el("td", { class: "neigh", rowspan: rows > 1 ? rows : null });
  const share = Math.max(0, Math.min(1, n.together / n.total));
  const more = words.neighbourMore(genome);
  return el("td", { class: "neigh", rowspan: rows > 1 ? rows : null },
    el("div", { class: "meter", role: "img", "aria-label": n.together + " of " + n.total + " neighbouring genes found side by side" },
      el("span", { css: { width: Math.round(share * 100) + "%" } })),
    el("div", { class: "count", text: words.neighbourCount(genome) }),
    more ? el("div", { class: "more", text: more }) : null);
}

function speciesCell(genome, rows) {
  const name = [species(genome.species)];
  if (genome.low_quality) name.push("†");
  return el("td", { rowspan: rows > 1 ? rows : null }, name,
    rows > 1 ? el("div", { class: "copy-n", text: rows + " copies" }) : null);
}

function genomeTable(result, fileUrl) {
  const s = result.summary;
  const rows = [];
  for (const genome of ordered(result.genomes)) {
    const calls = counted(genome);
    if (!calls.length) {
      const skipped = genome.status === "not_searched";
      rows.push(el("tr", { class: "none" },
        speciesCell(genome, 1),
        el("td", {}, pill(skipped ? "not_searched" : "none")),
        el("td", { class: "num" }), el("td", { class: "num" }),
        neighbourCell(genome, 1),
        el("td", { class: "why", text: skipped ? "The assembly is too fragmented to search." : words.notFoundText(genome) })));
      continue;
    }
    calls.forEach((call, index) => {
      rows.push(el("tr", {},
        index === 0 ? speciesCell(genome, calls.length) : null,
        el("td", {}, pill(call.label)),
        el("td", { class: "num", text: words.percent(call.identity) }),
        el("td", { class: "num", text: call.coverage === null ? "" : words.percent(call.coverage * 100) }),
        index === 0 ? neighbourCell(genome, calls.length) : null,
        el("td", { class: "where", text: where(call) })));
    });
  }
  const aside = [];
  for (const genome of result.genomes) {
    for (const call of genome.calls.filter((c) => c.label === "set_aside")) {
      const reason = call.reason === "related_gene"
        ? "matches " + (call.related_gene ? "the related gene " + call.related_gene : "a related gene") + " better than the query"
        : "matches only a short piece of the query";
      aside.push(el("li", {}, species(genome.species), ", ", where(call), ": ", reason,
        call.identity === null ? "" : " (" + words.percent(call.identity) + " identical)", "."));
    }
  }
  const chip = (n, text) => el("span", { class: "chip" }, el("b", { text: String(n) }), text);
  const table = download(result, "table");
  return el("section", { class: "card stack", id: "genomes" },
    el("div", { class: "row-between" },
      el("h2", { text: "Genome by genome" }),
      el("div", { class: "chips", "aria-label": "Tally" },
        chip(s.found, "found"), chip(s.candidate, words.plural(s.candidate, "candidate", "candidates")),
        chip(s.none, "not found"), s.not_searched ? chip(s.not_searched, "not searched") : null)),
    el("p", { class: "small muted", text: "Every genome that was searched is listed, also where the gene was not found: the bar shows how many of the gene’s neighbours lie side by side there." }),
    el("div", { class: "table-scroll" },
      el("table", { class: "calls" },
        el("thead", {}, el("tr", {},
          el("th", { scope: "col", text: "Species" }), el("th", { scope: "col", text: "Result" }),
          el("th", { scope: "col", class: "num", text: "Identical" }),
          el("th", { scope: "col", class: "num", text: "Covered" }),
          el("th", { scope: "col", text: "Neighbouring genes found" }),
          el("th", { scope: "col", text: "Where in the genome, or what was seen" }))),
        el("tbody", {}, rows))),
    aside.length ? el("div", { class: "stack-s" },
      el("h3", { text: "Set aside: " + aside.length + " " + words.plural(aside.length, "call that is not the gene", "calls that are not the gene") }),
      el("ul", { class: "small" }, aside)) : null,
    table ? el("p", { class: "small" }, el("a", { href: fileUrl(table.file), download: "genome_by_genome.tsv", text: "This table as a file" }), " (TSV, opens in a spreadsheet).") : null);
}

// ── explanations and files ─────────────────────────────────────────────────

function meanings(result) {
  const home = (result.home && result.home.species) || "the home species";
  return el("div", { class: "card stack-s", id: "words" },
    el("h2", { text: "What the words mean" }),
    el("div", { class: "words small" },
      el("div", {}, pill("found"), " ", pill("probable"), " A gene with matching sequence, inside a neighbourhood that matches ", species(home), ". Only these count as the gene."),
      el("div", {}, pill("candidate"), " The neighbourhood matches. What sits in it could be the gene, but its own sequence does not show it. A lost gene leaves the same neighbourhood behind."),
      el("div", {}, pill("none"), " Nothing passed the checks. The gene may be absent, too different, or in a part of the genome that was not assembled. A neighbourhood that is there without the gene is a result too."),
      el("div", {}, el("strong", { text: "Neighbouring genes found " }), "counts the genes next to the query in ", species(home), " that lie side by side in this genome. The more of them, the surer it is that the right place was searched."),
      el("div", {}, el("strong", { text: "Identical " }), "is the share of matching residues; ", el("strong", { text: "covered " }), "is how much of the query protein the copy spans. A high identity over a small part says little.")));
}

function downloads(result, fileUrl, options) {
  const items = [];
  const add = (name, format, file, size) => items.push(el("li", {},
    el("a", { href: fileUrl(file), download: file.split("/").pop() },
      el("span", { class: "name", text: name }),
      el("span", { class: "format", text: format + (size ? " · " + words.bytes(size) : "") }))));
  for (const item of result.downloads || []) {
    const w = words.DOWNLOADS[item.key];
    if (w && item.key !== "all") add(w[0], w[1], item.file, item.bytes);
  }
  const many = (result.loci || []).length > 1;
  for (const locus of result.loci || []) {
    const tag = many ? " (place " + locus.n + ")" : "";
    const f = locus.files || {};
    if (f.print) add("Figure at print size" + tag, "SVG", f.print);
    if (f.print_legend) add("Legend text for that figure" + tag, "TXT", f.print_legend);
    if (f.print_values) add("Values behind the figure" + tag, "TSV", f.print_values);
    if (f.gene_names) add("Gene names used in the figures" + tag, "TSV", f.gene_names);
    if (f.tree_newick) add("Gene tree" + tag, "Newick", f.tree_newick);
    if (f.grid_fragments) add("Grid with fragments drawn" + tag, "SVG", f.grid_fragments);
  }
  const accession = result.query && result.query.accession;
  const all = download(result, "all");
  return el("div", { class: "card stack", id: "files" },
    el("h2", { text: "Take it with you" }),
    all ? el("div", { class: "download-all" },
      el("div", { class: "stack-s" },
        el("strong", { text: "Everything in one file" }),
        el("span", { class: "small muted", text: "Figures, gene models, the table, the report and a read-me." })),
      el("a", { class: "button", href: fileUrl(all.file), download: all.file.split("/").pop(), text: "Download ZIP · " + words.bytes(all.bytes) })) : null,
    el("div", { class: "stack-s" },
      all ? el("h3", { text: "Single files" }) : null,
      el("ul", { class: "downloads" }, items)),
    options.hasService ? el("div", { class: "row" },
      el("a", { class: "button secondary", href: "search.html" + (accession ? "?acc=" + encodeURIComponent(accession) : ""), text: accession ? "Search other species" : "Start a new search" })) : null);
}

function specialists(result) {
  const s = result.summary;
  const p = result.pipeline || {};
  const q = result.query || {};
  const h = result.home || {};
  const fact = (name, value) => (value ? [el("dt", { text: name }), el("dd", { text: value })] : []);
  const loci = (result.loci || []).filter((l) => l.chrom).map((l) =>
    l.chrom + ":" + words.number(l.start) + "-" + words.number(l.end)
    + (l.neighbours ? " (" + l.neighbours + " neighbouring genes)" : "")).join("; ");
  return el("details", { class: "card", id: "details" },
    el("summary", { text: "For specialists: settings, software and checks" }),
    el("dl", { class: "facts" },
      fact("Software", words.softwareLine(p)),
      fact("Settings", [p.profile ? "profile " + p.profile : "", p.preset || ""].filter(Boolean).join(", ")),
      fact("Query", [q.accession, q.gene, q.length ? q.length + " residues" : "", q.organism].filter(Boolean).join(", ")),
      fact("Home genome", [h.species, h.assembly].filter(Boolean).join(", ")),
      fact(words.plural((result.loci || []).length, "Place in the home genome", "Places in the home genome"), loci),
      fact("Calls counted as the gene", s.copies + " (" + s.copies_clear + " high, " + s.copies_probable + " medium confidence)"),
      fact("Candidates, set aside", s.candidates + " ambiguous, " + s.set_aside + " set aside"),
      fact("Weak matches not listed as calls", String(s.weak_leads) + " low-confidence " + words.plural(s.weak_leads, "model", "models") + " (in the report file)"),
      fact("Steps that failed", String(s.failed_tasks || 0)),
      fact("Finished", result.timing && result.timing.finished_at ? result.timing.finished_at.replace("T", " ").replace("Z", " UTC") : "")),
    el("p", { class: "small muted" }, "The labels of the report file are HIGH (found), MEDIUM (probable) and AMBIGUOUS (candidate). ",
      el("a", { href: "about.html#words", text: "How a call gets its label" }), "."));
}

export function renderResult(root, result, fileUrl, options) {
  options = options || {};
  document.title = words.title(result.query, result.summary.genomes, result.scope) + " · SynVoy";
  const kept = options.expiresAt ? "This result is kept until " + words.day(options.expiresAt) + ". " : "";
  const figure = figures(result, fileUrl);
  const parts = [figure ? ["figures", "Figures"] : null, ["genomes", "Genome by genome"],
    ["words", "What the words mean"], ["files", "Files"], ["details", "Details"]].filter(Boolean);
  replace(root,
    intro(result, fileUrl, options),
    subnav(parts),
    notes(result),
    figure,
    genomeTable(result, fileUrl),
    el("section", { class: "pair" }, meanings(result), downloads(result, fileUrl, options)),
    specialists(result),
    el("p", { class: "smaller muted" }, kept, "The same gene, genomes and software give the same result every time. ",
      el("a", { href: "about.html#cite", text: "How to cite" }), "."));
}
