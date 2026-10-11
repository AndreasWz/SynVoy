// Talking to the search service. Every call ends in one of two ways: the decoded
// answer, or an ApiError with a message that can be shown to the visitor. There is
// no third way (a hanging request is ended by the time limit below).

const settings = window.SYNVOY || {};
export const hasService = settings.apiBase !== null && settings.apiBase !== undefined;
const BASE = (hasService ? String(settings.apiBase) : "").replace(/\/+$/, "") + "/api/v1";

export class ApiError extends Error {
  constructor(status, code, message, field) {
    super(message);
    this.status = status;       // 0 = no answer at all
    this.code = code;
    this.field = field || null;
  }
  get offline() { return this.status === 0; }
  get temporary() { return this.status === 0 || this.status === 429 || this.status >= 500; }
}

export function searchUrl(id, tail) {
  return BASE + "/searches/" + encodeURIComponent(id) + (tail || "");
}

export async function request(method, path, body, timeoutMs) {
  if (!hasService) {
    throw new ApiError(0, "no_service", "Searches cannot be started on this copy of the site.");
  }
  const controller = new AbortController();
  const timer = setTimeout(() => controller.abort(), timeoutMs || 20000);
  let response;
  try {
    response = await fetch(BASE + path, {
      method,
      headers: body === undefined ? {} : { "Content-Type": "application/json" },
      body: body === undefined ? undefined : JSON.stringify(body),
      signal: controller.signal,
      cache: "no-store",
      credentials: "omit",
      referrerPolicy: "no-referrer",
    });
  } catch (error) {
    throw new ApiError(0, controller.signal.aborted ? "timeout" : "network",
      "The server did not answer. Check your connection; this page keeps trying.");
  } finally {
    clearTimeout(timer);
  }
  let data = null;
  try {
    data = await response.json();
  } catch (error) {
    data = null;
  }
  if (response.ok && data !== null) return data;
  if (data && data.error && typeof data.error.message === "string") {
    throw new ApiError(response.status, data.error.code || "error", data.error.message, data.error.field);
  }
  // Not an answer of the service: a proxy error page, or pages hosted without a service.
  const message = response.status === 404
    ? "The search service was not found at this address."
    : response.status === 429
      ? "Too many requests. Wait a moment and try again."
      : "The search service is not available at the moment (" + response.status + ").";
  throw new ApiError(response.status, "no_api", message);
}

export const getStatus = () => request("GET", "/status");
// The genomes on the server. With ``home`` (an id) or ``taxon`` (an NCBI taxon id) the
// answer also orders them, closest relatives first.
export const getLibrary = (home, taxon) => request("GET", "/library"
  + (home ? "?home=" + encodeURIComponent(home) : taxon ? "?taxon=" + encodeURIComponent(taxon) : ""));
export const findGenomes = (body) => request("POST", "/genomes/find", body, 45000);
export const lookup = (body) => request("POST", "/lookup", body, 40000);
export const startSearch = (body) => request("POST", "/searches", body, 40000);
export const getSearch = (id) => request("GET", "/searches/" + encodeURIComponent(id));
export const getResult = (id) => request("GET", "/searches/" + encodeURIComponent(id) + "/result");
export const cancelSearch = (id, token) => request("POST", "/searches/" + encodeURIComponent(id) + "/cancel", { token });

// A file kept next to the pages (examples). Null when it is not there.
export async function getStatic(path) {
  try {
    const response = await fetch(path, { cache: "no-cache", credentials: "omit" });
    if (!response.ok) return null;
    return await response.json();
  } catch (error) {
    return null;
  }
}
