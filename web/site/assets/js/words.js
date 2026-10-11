// Every sentence the pages build from numbers lives here, so the wording is in one
// place and the same on the start page, while a search runs, and on the result.

export const LABEL = { found: "Found", probable: "Probable", candidate: "Candidate", set_aside: "Set aside",
  none: "Not found", not_searched: "Not searched" };

export function plural(n, one, many) {
  return n === 1 ? one : many;
}

export function number(n) {
  return Number(n).toLocaleString("en-US");
}

export function percent(value, digits) {
  if (value === null || value === undefined || Number.isNaN(Number(value))) return "";
  return Number(value).toFixed(digits || 0) + " %";
}

export function megabases(mb) {
  if (mb === null || mb === undefined) return "";
  if (mb >= 1000) return (mb / 1000).toFixed(1) + " Gb";
  return mb < 0.5 ? "under 1 Mb" : Math.round(mb) + " Mb";
}

export function bytes(n) {
  if (n >= 1048576) return (n / 1048576).toFixed(1) + " MB";
  if (n >= 1024) return Math.round(n / 1024) + " kB";
  return n + " B";
}

// 45 s, 8 min, 1 h 08 min
export function duration(seconds) {
  if (seconds === null || seconds === undefined) return "";
  const s = Math.max(0, Math.round(Number(seconds)));
  if (s < 90) return s + " s";
  const minutes = Math.round(s / 60);
  if (minutes < 60) return minutes + " min";
  const rest = minutes % 60;
  return Math.floor(minutes / 60) + " h " + String(rest).padStart(2, "0") + " min";
}

export function day(iso) {
  const date = new Date(iso);
  if (Number.isNaN(date.getTime())) return "";
  return date.toLocaleDateString("en-GB", { day: "numeric", month: "long", year: "numeric" });
}

// The short name a search goes by: the gene symbol when there is one ("SYP1"), else
// the name that was given or looked up.
export function queryName(query) {
  return (query && (query.gene || query.name || query.label || query.accession)) || "The gene";
}

// The long name, when it says more than the short one.
export function queryLongName(query) {
  const long = (query && query.name) || "";
  return long && long !== queryName(query) ? long : "";
}

// "SynVoy 1.0.0" for a release. Anything else is "SynVoy pre-release": a number is
// shown only when the server runs exactly that tagged release.
export function software(pipeline) {
  return "SynVoy " + ((pipeline && pipeline.label) || "pre-release");
}

export function softwareLine(pipeline) {
  let text = software(pipeline);
  if (pipeline && pipeline.commit) {
    text += ", commit " + String(pipeline.commit).slice(0, 10) + (pipeline.modified ? " (with local changes)" : "");
  }
  return text;
}

// "Apyrase in 33 ants and wasps" / "STE2 in 3 genomes"
export function title(query, genomes, scope) {
  const where = scope ? genomes + " " + scope : genomes + " " + plural(genomes, "genome", "genomes");
  return queryName(query) + " in " + where;
}

function ofThe(n, total) {
  if (total === 1) return "the genome";
  if (n === total) return n === 2 ? "both genomes" : "all " + total + " genomes";
  return n + " of the " + total + " genomes";
}

// The paragraph under the title of a result.
export function headline(result) {
  const s = result.summary;
  const name = queryName(result.query);
  const searched = s.genomes - (s.not_searched || 0);
  const parts = [];
  if (searched === 0) {
    parts.push("None of the genomes could be searched.");
  } else if (s.found === 0 && s.candidate === 0) {
    parts.push(name + " was not found in " + (searched === 1 ? "the genome" : "any of the " + searched + " genomes") + ".");
  } else {
    if (s.found > 0) {
      parts.push((s.clear > 0 ? name + " was found in " : "A probable copy of " + name + " was found in ")
        + ofThe(s.found, searched) + ".");
    }
    if (s.candidate > 0) {
      parts.push((s.found > 0 ? "In " + s.candidate + " more" : "In " + ofThe(s.candidate, searched))
        + ", the neighbourhood is there and holds a candidate that could not be confirmed.");
    }
    if (s.none > 0) {
      parts.push((s.none === 1 ? "In the other one" : "In the other " + s.none) + " it was not found; the table says how much of its neighbourhood is there.");
    }
  }
  if (s.not_searched > 0) {
    parts.push(s.not_searched + " " + plural(s.not_searched, "genome was", "genomes were")
      + " left out because the assembly is too fragmented.");
  }
  return parts.join(" ");
}

// One line for an example card on the start page.
export function cardSentence(card) {
  const parts = [];
  if (card.found > 0) {
    parts.push((card.clear > 0 ? "Found in " : "Probable in ")
      + (card.found === card.genomes ? (card.genomes === 2 ? "both" : "all " + card.genomes)
        : card.found + " of " + card.genomes) + " genomes");
  } else {
    parts.push("Not found in " + card.genomes + " " + plural(card.genomes, "genome", "genomes"));
  }
  if (card.candidate > 0) parts.push("a candidate in " + card.candidate + (card.found > 0 ? " more" : ""));
  if (card.none > 0 && (card.found > 0 || card.candidate > 0)) parts.push("not found in " + card.none);
  let text = parts.join(", ") + ".";
  if (card.identity_min !== null && card.identity_min !== undefined && card.found > 0) {
    const low = Math.round(card.identity_min), high = Math.round(card.identity_max);
    text += low === high ? " " + low + " % identical." : " " + low + " to " + high + " % identical.";
  }
  return text;
}

// "10 of 20" neighbouring genes: how many of the home gene's neighbours lie side by
// side in this genome, and what else was seen of the neighbourhood.
export function neighbourCount(genome) {
  const n = genome.neighbours;
  if (!n || !n.total) return "";
  return n.together + " of " + n.total;
}

export function neighbourMore(genome) {
  const n = genome.neighbours;
  if (!n || !n.total) return "";
  if (n.found === 0) return "none found";
  if (n.found > n.together) return n.found + " in all, in " + n.places + " " + plural(n.places, "place", "places");
  return "";
}

// What a genome without a call still tells.
export function notFoundText(genome) {
  const n = genome.neighbours;
  const leads = genome.leads;
  const parts = [];
  if (n && n.total) {
    // Five neighbours side by side is what the search itself takes as a firm neighbourhood.
    if (n.found === 0) parts.push("None of the " + n.total + " neighbouring genes was found either.");
    else if (n.together >= 5) parts.push("Its neighbourhood is there; the gene was not found in it.");
    else if (n.together >= 2) parts.push("Part of its neighbourhood is there (" + n.together + " genes side by side); the gene was not found in it.");
    else parts.push("Only single neighbouring genes were found, too few to tell where the gene would be.");
  }
  if (leads && leads.n) {
    let text = leads.n + " weak " + plural(leads.n, "match", "matches") + " of the gene";
    if (leads.identity !== null && leads.identity !== undefined) {
      text += ", the best " + percent(leads.identity) + " identical";
      if (leads.coverage !== null && leads.coverage !== undefined) text += " over " + percent(leads.coverage * 100) + " of the protein";
    }
    parts.push(text + ": not enough for a call.");
  }
  return parts.join(" ");
}

// A genome a waiting search needs first.
export function fetchText(genome) {
  if (genome.state === "fetching" && genome.total_bytes > 0) {
    return "Fetching from NCBI: " + Math.min(100, Math.round(100 * genome.done_bytes / genome.total_bytes)) + " %";
  }
  if (genome.state === "fetching") return "Fetching from NCBI";
  return "Waiting to be fetched from NCBI";
}

export const STEP_TITLES = {
  locate: ["Locating the gene in its home genome", "Gene located in its home genome"],
  neighbours: ["Reading the neighbouring genes", "Neighbouring genes read"],
  search: ["Searching the genomes", "Genomes searched"],
  check: ["Checking each call against related genes", "Calls checked against related genes"],
  report: ["Writing figures and summary", "Figures and summary written"],
};

export function stepText(step, job) {
  const titles = STEP_TITLES[step.key] || [step.key, step.key];
  const detail = step.detail || {};
  const home = (job.home && job.home.species) || "the home genome";
  let text = "";
  if (step.key === "locate") {
    if (step.state === "done" && detail.loci) {
      text = detail.loci === 1 ? "One place found in " + home + "."
        : detail.loci + " places found in " + home + ". Each is searched for separately.";
    } else if (step.state === "active") {
      text = "Looking for the protein in the genome of " + home + ".";
    }
  } else if (step.key === "neighbours") {
    text = step.state === "waiting" ? "" : "The genes next to it are the landmarks for the search.";
  } else if (step.key === "search") {
    const bits = [];
    if (detail.loci_total) bits.push("Place " + Math.min(detail.loci_done + 1, detail.loci_total) + " of " + detail.loci_total + ".");
    if (step.state === "active" && detail.genomes_total) {
      bits.push(detail.genomes_done + " of " + detail.genomes_total + " genomes done.");
    } else if (step.state === "active") {
      bits.push("This is the long step.");
    }
    text = bits.join(" ");
  } else if (step.key === "check") {
    text = "A call that matches another gene of " + home + " better is set aside.";
  } else if (step.key === "report") {
    text = step.state === "waiting" ? "Waiting." : "";
  }
  return { title: step.state === "done" ? titles[1] : titles[0], detail: text };
}

export const DOWNLOADS = {
  all: ["Everything in one file", "ZIP"],
  table: ["The table, genome by genome", "TSV"],
  report: ["Summary of all calls", "JSON"],
  gene_models: ["Gene models per genome", "GFF, zipped"],
  query: ["The query sequence", "FASTA"],
};

export const FIGURES = {
  grid: {
    tab: "Grid",
    alt: "Grid figure: one row per genome, one column per neighbouring gene, with an arrow where that gene was placed.",
    how: "Each row is a genome, each column one of the genes next to the query gene in its home genome. An arrow means that gene was found in the same neighbourhood. The paler the arrow, the more the protein has changed. The column of the query gene shows one small arrow per copy.",
  },
  ribbons: {
    tab: "Ribbons",
    alt: "Ribbon figure: the neighbourhood drawn to scale in each genome, with ribbons joining matching genes.",
    how: "The neighbourhood drawn to scale, one line per genome. Ribbons join the same gene in neighbouring lines, so rearrangements show as crossing ribbons.",
  },
  positions: {
    tab: "Positions",
    alt: "Position figure: where each gene of the neighbourhood lies in each genome.",
    how: "Where each gene lies on its chromosome or scaffold, genome by genome.",
  },
  anchors: {
    tab: "Landmarks",
    alt: "Landmark figure: the neighbouring genes found in each genome, at their positions.",
    how: "The neighbouring genes that were found in each genome, at their positions. Copies of the query gene are drawn as notched arrows.",
  },
  tree: {
    tab: "Gene tree",
    alt: "Tree of the recovered proteins.",
    how: "A tree of the proteins that were recovered. It is a rough guide: short or partial proteins place badly.",
  },
};
