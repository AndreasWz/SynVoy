// The searches started from this browser: their links and the tokens that allow
// cancelling them. Kept in the browser only. If the browser keeps no data (private
// window, blocked storage) everything still works; the visitor then has to keep the
// link, and cancelling is not offered.

const KEY = "synvoy.searches.v1";
const MAX = 20;

function read() {
  try {
    const list = JSON.parse(window.localStorage.getItem(KEY) || "[]");
    return Array.isArray(list) ? list.filter((e) => e && typeof e.id === "string") : [];
  } catch (error) {
    return [];
  }
}

function write(list) {
  try {
    window.localStorage.setItem(KEY, JSON.stringify(list.slice(0, MAX)));
    return true;
  } catch (error) {
    return false;
  }
}

export function recentSearches() {
  return read();
}

export function remember(entry) {
  const list = read().filter((e) => e.id !== entry.id);
  list.unshift({ id: entry.id, token: entry.token || "", title: entry.title || "", created: entry.created || new Date().toISOString() });
  return write(list);
}

export function tokenFor(id) {
  const found = read().find((e) => e.id === id);
  return found && found.token ? found.token : "";
}

export function forget(id) {
  write(read().filter((e) => e.id !== id));
}
