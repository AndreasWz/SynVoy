// What every page shares: the state of the search service (asked once per page), the
// notice shown when it takes no searches, and the version line in the footer.

import { getStatus, hasService } from "./api.js";
import { el } from "./dom.js";
import { software } from "./words.js";

let statusPromise = null;

// Resolves to {live: true, data} or {live: false, message}. Never rejects.
export function serviceStatus() {
  if (!statusPromise) {
    statusPromise = !hasService
      ? Promise.resolve({ live: false, none: true, message: "Searches cannot be started on this copy of the site. The examples are real results." })
      : getStatus().then(
        (data) => ({ live: true, data }),
        (error) => ({ live: false, none: false, message: error.message }));
  }
  return statusPromise;
}

function versionText(pipeline) {
  if (!pipeline) return "SynVoy";
  let text = software(pipeline);
  if (pipeline.commit) text += " (" + String(pipeline.commit).slice(0, 7) + (pipeline.modified ? ", modified" : "") + ")";
  return text;
}

// The header stays at the top of the window. It gets a shadow once the page has
// moved under it, and on a narrow screen its links fold into a menu.
function initHeader() {
  document.documentElement.classList.add("js");
  const header = document.querySelector(".site-header");
  if (!header) return;
  const onScroll = () => header.classList.toggle("is-stuck", window.scrollY > 4);
  window.addEventListener("scroll", onScroll, { passive: true });
  onScroll();
  const toggle = header.querySelector(".nav-toggle");
  if (!toggle) return;
  const set = (open) => {
    if (open) header.setAttribute("data-open", ""); else header.removeAttribute("data-open");
    toggle.setAttribute("aria-expanded", String(open));
  };
  toggle.addEventListener("click", () => set(!header.hasAttribute("data-open")));
  document.addEventListener("keydown", (event) => { if (event.key === "Escape") set(false); });
  for (const link of header.querySelectorAll(".site-nav a")) link.addEventListener("click", () => set(false));
}

export async function initShell() {
  initHeader();
  const here = (window.location.pathname.split("/").pop() || "index.html");
  for (const link of document.querySelectorAll(".site-nav a[href]")) {
    if (link.getAttribute("href") === here) link.setAttribute("aria-current", "page");
  }
  const state = await serviceStatus();
  const banner = document.getElementById("service-banner");
  if (state.live) {
    const data = state.data;
    for (const node of document.querySelectorAll("[data-version]")) node.textContent = versionText(data.pipeline);
    const operator = data.operator || {};
    for (const node of document.querySelectorAll("[data-operator]")) {
      if (operator.name) node.textContent = operator.name;
    }
    if (!data.accepting && banner) {
      banner.replaceChildren(el("div", { class: "wrap" },
        el("strong", { text: "No new searches at the moment. " }),
        data.message || "", " Results that exist can still be opened."));
      banner.hidden = false;
    }
  } else if (banner && !state.none && document.body.dataset.needsService !== undefined) {
    banner.replaceChildren(el("div", { class: "wrap" },
      el("strong", { text: "The search service does not answer. " }), state.message));
    banner.hidden = false;
  }
  // Texts that hold only where the server fetches genomes from NCBI on request.
  const downloads = !!(state.live && state.data.limits && state.data.limits.downloads);
  for (const node of document.querySelectorAll("[data-if-downloads]")) node.hidden = !downloads;
  for (const node of document.querySelectorAll("[data-if-service]")) node.hidden = !!state.none;
  for (const node of document.querySelectorAll("[data-if-no-service]")) node.hidden = !state.none;
  return state;
}
