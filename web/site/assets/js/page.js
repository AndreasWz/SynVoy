// Pages with fixed text (about, privacy, imprint): only the shared parts are filled in.

import { initShell } from "./shell.js";

initShell().then((state) => {
  if (!state.live) return;
  const limits = state.data.limits;
  for (const node of document.querySelectorAll("[data-limit]")) {
    const value = limits[node.dataset.limit];
    if (value !== undefined) node.textContent = String(value);
  }
  const operator = state.data.operator || {};
  for (const node of document.querySelectorAll("[data-contact]")) {
    if (operator.contact) node.textContent = operator.contact;
  }
});
