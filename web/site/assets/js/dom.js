// Small helpers for building the page. Text always goes in as text, never as markup,
// so nothing that comes from a search (names, messages) can turn into code.

export function el(tag, attrs, ...children) {
  const node = document.createElement(tag);
  for (const [key, value] of Object.entries(attrs || {})) {
    if (value === null || value === undefined || value === false) continue;
    if (key === "class") node.className = value;
    else if (key === "text") node.textContent = value;
    else if (key === "css") Object.assign(node.style, value);   // not a style attribute: the page's policy forbids those
    else if (key.startsWith("on") && typeof value === "function") node.addEventListener(key.slice(2), value);
    else if (value === true) node.setAttribute(key, "");
    else node.setAttribute(key, String(value));
  }
  append(node, children);
  return node;
}

export function append(node, children) {
  for (const child of children.flat(Infinity)) {
    if (child === null || child === undefined || child === false) continue;
    node.append(child instanceof Node ? child : document.createTextNode(String(child)));
  }
  return node;
}

export function replace(node, ...children) {
  node.replaceChildren();
  return append(node, children);
}

export const byId = (id) => document.getElementById(id);

export function species(name) {
  return el("em", { class: "species", text: name });
}

const SVG = "http://www.w3.org/2000/svg";
export function icon(kind) {
  // done: ticked circle, active: clock, waiting: empty circle, failed: crossed circle
  const colors = { done: "#0B6B57", active: "#1F4E8C", waiting: "#8A9A9F", failed: "#8E1B1B" };
  const svg = document.createElementNS(SVG, "svg");
  svg.setAttribute("viewBox", "0 0 24 24");
  svg.setAttribute("width", "24");
  svg.setAttribute("height", "24");
  svg.setAttribute("fill", "none");
  svg.setAttribute("stroke", colors[kind] || colors.waiting);
  svg.setAttribute("stroke-width", "2.5");
  svg.setAttribute("stroke-linecap", "round");
  svg.setAttribute("stroke-linejoin", "round");
  svg.setAttribute("aria-hidden", "true");
  svg.setAttribute("class", "mark");
  const circle = document.createElementNS(SVG, "circle");
  circle.setAttribute("cx", "12");
  circle.setAttribute("cy", "12");
  circle.setAttribute("r", "10");
  svg.append(circle);
  const paths = { done: "M7.5 12.5l3 3 6-6.5", active: "M12 7v5l3.5 2", failed: "M8.5 8.5l7 7M15.5 8.5l-7 7" };
  if (paths[kind]) {
    const path = document.createElementNS(SVG, "path");
    path.setAttribute("d", paths[kind]);
    svg.append(path);
  }
  return svg;
}

export function params() {
  return new URLSearchParams(window.location.search);
}

// Copy text; resolves to true when it worked.
export async function copyText(text) {
  try {
    await navigator.clipboard.writeText(text);
    return true;
  } catch (error) {
    try {
      const area = el("textarea", { readonly: true });
      area.value = text;
      area.style.position = "fixed";
      area.style.opacity = "0";
      document.body.append(area);
      area.select();
      const ok = document.execCommand("copy");
      area.remove();
      return ok;
    } catch (error2) {
      return false;
    }
  }
}
