// index.html — the examples, and the searches started from this browser.

import { getSearch, getStatic, hasService } from "./api.js";
import { byId, el, replace } from "./dom.js";
import { forget, recentSearches } from "./recent.js";
import { initShell } from "./shell.js";
import * as words from "./words.js";

function bar(card) {
  const part = (cls, n) => (n > 0 ? el("span", { class: cls, css: { flexGrow: String(n) } }) : null);
  return el("div", { class: "bar", role: "img",
    "aria-label": card.found + " found, " + card.candidate + " candidate, " + card.none + " not found" },
    part("b-found", card.found), part("b-candidate", card.candidate), part("b-none", card.none));
}

function exampleCard(card) {
  const sub = [card.about, card.length ? card.length + " residues" : ""].filter(Boolean).join(", ");
  const foot = [card.genomes + " " + (card.scope || words.plural(card.genomes, "genome", "genomes")),
    card.runtime_s ? words.duration(card.runtime_s) : ""].filter(Boolean).join(" · ");
  return el("a", { class: "card", href: "result.html?example=" + encodeURIComponent(card.slug) },
    el("div", { class: "stack-s" },
      el("div", { class: "title", text: card.name }),
      sub ? el("div", { class: "smaller muted", text: sub }) : null),
    bar(card),
    el("div", { class: "small", text: words.cardSentence(card) }),
    el("div", { class: "mono smaller muted", text: foot }));
}

async function showExamples() {
  const box = byId("example-cards");
  const listing = await getStatic("data/gallery.json");
  const cards = listing && Array.isArray(listing.examples) ? listing.examples : [];
  if (!cards.length) {
    replace(box, el("p", { class: "muted", text: "There are no examples on this site yet." }));
    byId("example-legend").hidden = true;
    return;
  }
  replace(box, cards.map(exampleCard));
}

const STATE_WORD = { queued: "Waiting", running: "Running", done: "Finished", failed: "Stopped with an error" };

async function showRecent() {
  const list = recentSearches().slice(0, 6);
  if (!hasService || !list.length) return;
  const section = byId("recent");
  const items = byId("recent-list");
  section.hidden = false;
  const rows = list.map((entry) => {
    const badge = el("span", { class: "badge", text: "…" });
    const row = el("li", { class: "row-between" },
      el("a", { href: "result.html?id=" + encodeURIComponent(entry.id), text: entry.title || "Search" }),
      el("span", { class: "row" },
        el("span", { class: "smaller muted", text: words.day(entry.created) }), badge));
    getSearch(entry.id).then((job) => {
      if (job.state === "cancelled" || job.state === "expired") {
        forget(entry.id);
        row.remove();
        return;
      }
      badge.textContent = STATE_WORD[job.state] || job.state;
      badge.className = "badge " + job.state;
    }, (error) => {
      if (error.status === 404) { forget(entry.id); row.remove(); } else badge.textContent = "?";
    }).finally(() => { if (!items.children.length) section.hidden = true; });
    return row;
  });
  replace(items, rows);
}

initShell();
showExamples();
showRecent();
