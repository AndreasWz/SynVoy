// Browser run-through of the site: a real Chrome, the real pages, the stand-in pipeline.
//
//   node run.mjs            (Node 22 or newer, google-chrome on the PATH)
//
// It starts its own throw-away service (service.py), clicks through the paths a visitor
// takes, and fails if a step does not happen or the browser reports any script error
// or blocked request. No packages are needed: Chrome is driven over its debugging port.

import { spawn, spawnSync } from "node:child_process";
import { cpSync, existsSync, mkdirSync, mkdtempSync, rmSync, writeFileSync } from "node:fs";
import net from "node:net";
import os from "node:os";
import path from "node:path";
import { fileURLToPath } from "node:url";

const HERE = path.dirname(fileURLToPath(import.meta.url));
const WEB = path.resolve(HERE, "..", "..");
const PYTHON = process.env.SYNVOY_WEB_PYTHON || "python3";
const CHROME = process.env.SYNVOY_WEB_CHROME || "google-chrome";
// A checkout whose examples have not been built is a valid site without examples.
const HAS_EXAMPLES = existsSync(path.join(WEB, "site", "data", "gallery.json"));
// SYNVOY_WEB_SHOTS=<folder>: also save pictures of the pages at some points.
const SHOTS = process.env.SYNVOY_WEB_SHOTS || "";
const PROTEIN = "MKFLVNVALVFMVVYISYIYAAPEPEPAPEPEAEADAEADPEAGIGAVLKVLTTGLPALISWIKRKRQQG";

const sleep = (ms) => new Promise((resolve) => setTimeout(resolve, ms));

function freePort() {
  return new Promise((resolve, reject) => {
    const server = net.createServer();
    server.once("error", reject);
    server.listen(0, "127.0.0.1", () => {
      const { port } = server.address();
      server.close(() => resolve(port));
    });
  });
}

function waitForLine(child, pattern, stream, timeoutMs) {
  return new Promise((resolve, reject) => {
    let seen = "";
    const timer = setTimeout(() => reject(new Error("no '" + pattern + "' within " + timeoutMs + " ms; got: " + seen.slice(-400))), timeoutMs);
    child[stream].on("data", (chunk) => {
      seen += chunk;
      const match = seen.match(pattern);
      if (match) { clearTimeout(timer); resolve(match); }
    });
    child.once("exit", (code) => { clearTimeout(timer); reject(new Error("process ended (" + code + "): " + seen.slice(-400))); });
  });
}

class Page {
  constructor(socket, sessionId) {
    this.socket = socket;
    this.sessionId = sessionId;
    this.problems = [];
  }

  send(method, params) { return this.socket.call(method, params || {}, this.sessionId); }

  async goto(url) {
    await this.send("Page.navigate", { url });
    await this.waitFor("document.readyState === 'complete'", 15000, "page load of " + url);
  }

  async evaluate(expression) {
    const answer = await this.send("Runtime.evaluate", { expression, awaitPromise: true, returnByValue: true });
    if (answer.exceptionDetails) {
      throw new Error("in page: " + (answer.exceptionDetails.exception?.description || answer.exceptionDetails.text));
    }
    return answer.result.value;
  }

  async waitFor(expression, timeoutMs, what) {
    const end = Date.now() + (timeoutMs || 20000);
    let last;
    while (Date.now() < end) {
      try { last = await this.evaluate("(() => { try { return (" + expression + "); } catch (e) { return false; } })()"); } catch (e) { last = false; }
      if (last) return last;
      await sleep(150);
    }
    const body = await this.evaluate("document.body ? document.body.innerText.slice(0, 600) : ''").catch(() => "");
    throw new Error("timed out waiting for " + (what || expression) + "\n--- page text ---\n" + body);
  }

  click(selector) {
    return this.evaluate("(() => { const n = document.querySelector(" + JSON.stringify(selector) + "); if (!n) throw new Error('no element " + selector + "'); n.click(); return true; })()");
  }

  fill(selector, text) {
    return this.evaluate("(() => { const n = document.querySelector(" + JSON.stringify(selector) + "); if (!n) throw new Error('no element " + selector + "'); n.value = " + JSON.stringify(text)
      + "; n.dispatchEvent(new Event('input', { bubbles: true })); n.dispatchEvent(new Event('change', { bubbles: true })); return true; })()");
  }

  text(selector) {
    return this.evaluate("(() => { const n = document.querySelector(" + JSON.stringify(selector) + "); return n ? n.innerText : null; })()");
  }

  // A picture of the whole page, when pictures were asked for. Taken from the top of the
  // page: the header stays in view while scrolling and would lie across the picture.
  async shot(name) {
    if (!SHOTS) return;
    mkdirSync(SHOTS, { recursive: true });
    // A smooth scroll that is still under way carries on after a jump to the top, so
    // the page is put back until it has stayed there for a few frames.
    await this.evaluate("new Promise((done) => { let still = 0, frames = 0; const step = () => { frames += 1; if (window.scrollY !== 0) { window.scrollTo({ top: 0, behavior: 'instant' }); still = 0; } else { still += 1; } if (still >= 6 || frames > 400) done(); else requestAnimationFrame(step); }; step(); })");
    const { data } = await this.send("Page.captureScreenshot", { format: "png", captureBeyondViewport: true });
    writeFileSync(path.join(SHOTS, name + ".png"), Buffer.from(data, "base64"));
  }
}

class Chrome {
  static async launch() {
    const profile = mkdtempSync(path.join(os.tmpdir(), "synvoy-chrome-"));
    const child = spawn(CHROME, ["--headless=new", "--no-sandbox", "--disable-gpu", "--remote-debugging-port=0",
      "--user-data-dir=" + profile, "--window-size=1280,1600", "about:blank"], { stdio: ["ignore", "ignore", "pipe"] });
    child.stderr.setEncoding("utf8");
    const match = await waitForLine(child, /DevTools listening on (ws:\/\/\S+)/, "stderr", 30000);
    const chrome = new Chrome(child, profile);
    await chrome.connect(match[1]);
    return chrome;
  }

  constructor(child, profile) {
    this.child = child;
    this.profile = profile;
    this.nextId = 1;
    this.pending = new Map();
    this.pages = new Map();
  }

  connect(url) {
    return new Promise((resolve, reject) => {
      this.ws = new WebSocket(url);
      this.ws.addEventListener("open", () => resolve());
      this.ws.addEventListener("error", (event) => reject(new Error("cannot reach Chrome: " + (event.message || ""))));
      this.ws.addEventListener("message", (event) => this.onMessage(JSON.parse(event.data)));
    });
  }

  onMessage(message) {
    if (message.id && this.pending.has(message.id)) {
      const { resolve, reject } = this.pending.get(message.id);
      this.pending.delete(message.id);
      if (message.error) reject(new Error(message.error.message)); else resolve(message.result);
      return;
    }
    const page = this.pages.get(message.sessionId);
    if (!page) return;
    const p = message.params || {};
    if (message.method === "Runtime.exceptionThrown") {
      page.problems.push("script error: " + (p.exceptionDetails.exception?.description || p.exceptionDetails.text));
    } else if (message.method === "Runtime.consoleAPICalled" && (p.type === "error" || p.type === "assert")) {
      page.problems.push("console error: " + p.args.map((a) => a.value ?? a.description ?? "").join(" "));
    } else if (message.method === "Log.entryAdded" && p.entry.level === "error" && !page.expectFailures) {
      // The service saying "no" (4xx) to a request of the page is an answer the page
      // handles; the run-through provokes several on purpose. Anything else counts:
      // a missing file, a blocked script or style, a refused connection.
      const refusal = /status of 4\d\d/.test(p.entry.text) && (p.entry.url || "").includes("/api/v1/");
      // Plain file hosting has no list of examples to hand out when none were built.
      const noExamples = !HAS_EXAMPLES && /status of 404/.test(p.entry.text)
        && (p.entry.url || "").endsWith("/data/gallery.json");
      if (!refusal && !noExamples) page.problems.push("browser: " + p.entry.text + (p.entry.url ? " (" + p.entry.url + ")" : ""));
    }
  }

  call(method, params, sessionId) {
    const id = this.nextId++;
    return new Promise((resolve, reject) => {
      this.pending.set(id, { resolve, reject });
      this.ws.send(JSON.stringify(sessionId ? { id, method, params, sessionId } : { id, method, params }));
    });
  }

  async newPage() {
    const { targetId } = await this.call("Target.createTarget", { url: "about:blank" });
    const { sessionId } = await this.call("Target.attachToTarget", { targetId, flatten: true });
    const page = new Page(this, sessionId);
    this.pages.set(sessionId, page);
    await page.send("Page.enable");
    await page.send("Runtime.enable");
    await page.send("Log.enable");
    return page;
  }

  close() {
    try { this.ws.close(); } catch (e) { /* already closed */ }
    this.child.kill("SIGKILL");
    try { rmSync(this.profile, { recursive: true, force: true }); } catch (e) { /* best effort */ }
  }
}

// ── the run-through ──────────────────────────────────────────────────────────

const steps = [];
function ok(text) { steps.push(text); console.log("  ok  " + text); }
function check(condition, text) {
  if (!condition) throw new Error("check failed: " + text);
  ok(text);
}

async function api(base, pathname) {
  const response = await fetch(base + pathname);
  return { status: response.status, body: await response.json().catch(() => null) };
}

async function withService(root, port) {
  const child = spawn(PYTHON, [path.join(HERE, "service.py"), root, String(port)], {
    stdio: ["ignore", "pipe", "pipe"], env: { ...process.env, VIRTUAL_ENV: "", SYNVOY_WEB_CONFIG: "" } });
  child.stdout.setEncoding("utf8");
  child.stderr.setEncoding("utf8");
  let errors = "";
  child.stderr.on("data", (chunk) => { errors += chunk; });
  await waitForLine(child, /ready/, "stdout", 30000).catch((error) => { throw new Error(error.message + "\n" + errors.slice(-800)); });
  return child;
}

async function main() {
  const scratch = mkdtempSync(path.join(os.tmpdir(), "synvoy-browser-"));
  const port = await freePort();
  const base = "http://127.0.0.1:" + port;
  let service = await withService(path.join(scratch, "service"), port);
  const chrome = await Chrome.launch();
  let staticServer = null;
  try {
    const page = await chrome.newPage();

    console.log("Start page");
    await page.goto(base + "/");
    await page.waitFor("document.querySelector('[data-version]').innerText.includes('SynVoy pre-release')", 15000, "what runs, in the footer");
    check(!/SynVoy \d/.test(await page.text("[data-version]")), "no version number is claimed for a checkout that is not a release");
    const cards = await page.evaluate("document.querySelectorAll('#example-cards a.card').length");
    ok("start page loads; " + cards + " example card(s)");
    await page.evaluate("window.scrollTo(0, 700)");
    await page.waitFor("document.querySelector('.site-header').classList.contains('is-stuck')", 5000, "the header after scrolling");
    check(await page.evaluate("Math.abs(document.querySelector('.site-header').getBoundingClientRect().top) < 1"), "the header stays at the top while the page scrolls");
    await page.evaluate("window.scrollTo(0, 0)");
    await page.shot("start");
    check(await page.evaluate("document.querySelector('#recent').hidden"), "no list of own searches before any was started");

    console.log("Search by accession, through to the result");
    await page.goto(base + "/search.html?acc=P01501");
    await page.waitFor("document.querySelector('#lookup-result').innerText.includes('Melittin')", 15000, "the look-up answer");
    check((await page.text("#lookup-result")).includes("Homus primus"), "look-up names the protein and its species");
    await page.waitFor("document.querySelectorAll('#species-list input[type=checkbox]').length >= 4", 15000, "the species list");
    const order = await page.evaluate("[...document.querySelectorAll('#species-list em')].map((n) => n.innerText)");
    check(order[0] === "Alpha beta" && order[order.length - 1] === "Eta theta", "species are listed closest relatives first: " + order.join(", "));
    check(await page.evaluate("[...document.querySelectorAll('#species-list label')].find((l) => l.innerText.includes('Magnus maximus')).querySelector('input').disabled"),
      "a genome above the size limit cannot be ticked");
    check((await page.text("#count")).startsWith("4 of 10"), "the closest usable genomes are ticked by default");
    await page.click("#submit");
    await page.waitFor("location.pathname.endsWith('/result.html') && new URLSearchParams(location.search).get('id')", 15000, "the move to the search page");
    const id = await page.evaluate("new URLSearchParams(location.search).get('id')");
    await page.waitFor("document.querySelectorAll('table.calls tbody tr').length === 4", 30000, "the result");
    check((await page.text("h1")).startsWith("MELT in 4 genomes"), "result title");
    check((await page.text(".result-head")).includes("Melittin"), "the long name of the protein is given under the title");
    const tally = await page.evaluate("[...document.querySelectorAll('.chips .chip')].map((n) => n.innerText.replace(/\\s+/g, ' ').trim())");
    check(tally.join(" / ") === "2 found / 1 candidate / 1 not found", "tally of the result: " + tally.join(" / "));
    const labels = await page.evaluate("[...document.querySelectorAll('table.calls tbody tr .pill')].map((n) => n.innerText)");
    check(labels.join(",") === "Found,Probable,Candidate,Not found", "every genome has a row, in the order Found, Probable, Candidate, Not found");
    const rows = await page.evaluate("[...document.querySelectorAll('table.calls tbody tr')].map((r) => r.innerText.replace(/\\s+/g, ' '))");
    check(rows[0].includes("6 of 6") && rows[3].includes("2 of 6"), "each row says how many neighbouring genes lie side by side there");
    check(rows[3].includes("neighbourhood is there") && rows[3].includes("the gene was not found in it") && rows[3].includes("2 weak matches"),
      "the genome without a call says what was found there: " + rows[3].slice(0, 160));
    check(await page.evaluate("document.querySelectorAll('.stat').length") === 0, "the three number tiles are gone");
    await page.waitFor("(() => { const i = document.querySelector('.figure-frame img'); return i && i.complete && i.naturalWidth > 0; })()", 15000, "the figure");
    const firstFigure = await page.evaluate("document.querySelector('.figure-frame img').getAttribute('src')");
    check(firstFigure.includes("/api/v1/searches/" + id + "/files/") && firstFigure.endsWith("anchor_grid.svg"), "the grid figure comes from the search's own files");
    await page.evaluate("[...document.querySelectorAll('.segmented button')].find((b) => b.innerText === 'Ribbons').click()");
    await page.waitFor("document.querySelector('.figure-frame img').getAttribute('src').endsWith('synteny_plot_view.svg')", 5000, "the ribbon figure after the click");
    ok("the figure tabs switch the figure");
    await page.click("#enlarge-figure");
    await page.waitFor("(() => { const d = document.querySelector('dialog.viewer'); const i = d && d.querySelector('img'); return d.open && i.complete && i.naturalWidth > 0; })()", 10000, "the enlarged figure");
    await page.evaluate("[...document.querySelectorAll('dialog.viewer button')].find((b) => b.innerText === 'Close').click()");
    check(await page.evaluate("!document.querySelector('dialog.viewer').open"), "a figure can be enlarged and closed again");
    const downloads = await page.evaluate("Promise.all([...document.querySelectorAll('.downloads a')].map((a) => fetch(a.href).then((r) => r.status)))");
    check(downloads.length >= 4 && downloads.every((s) => s === 200), "all " + downloads.length + " download links answer");
    const everything = await page.evaluate("fetch(document.querySelector('#download-all').href).then(async (r) => [r.status, r.headers.get('content-type'), (await r.arrayBuffer()).byteLength, document.querySelector('#download-all').getAttribute('download')])");
    check(everything[0] === 200 && everything[1] === "application/zip" && everything[2] > 1000 && everything[3] === "synvoy_MELT.zip",
      "everything comes in one file: " + everything[3] + ", " + everything[2] + " bytes");
    check(await page.evaluate("document.title").then((t) => t.startsWith("MELT in 4 genomes")), "the browser tab is named after the search");
    await page.evaluate("document.querySelector('#genomes').scrollIntoView()");
    await page.waitFor("(() => { const b = document.querySelector('.subnav').getBoundingClientRect(); const h = document.querySelector('.site-header').getBoundingClientRect(); return Math.abs(h.top) < 1 && Math.abs(b.top - h.bottom) < 2; })()", 5000, "the jump bar under the header");
    ok("header and jump bar stay in place on a long result");
    await page.shot("result");

    await page.goto(base + "/");
    await page.waitFor("!document.querySelector('#recent').hidden && document.querySelector('#recent-list .badge').innerText === 'Finished'", 15000, "the list of own searches");
    ok("the start page lists the search as finished");

    console.log("The same search again");
    await page.goto(base + "/search.html?acc=P01501");
    await page.waitFor("document.querySelectorAll('#species-list input:checked').length === 4", 15000, "the form");
    await page.click("#submit");
    await page.waitFor("document.querySelectorAll('table.calls tbody tr').length === 4", 20000, "the stored result");
    check(await page.evaluate("new URLSearchParams(location.search).get('id')") === id, "an identical search leads to the same result, without a new run");

    console.log("Paste a sequence, watch it run, cancel it");
    await page.goto(base + "/search.html");
    await page.waitFor("document.querySelector('#home').options.length > 1", 15000, "the form");
    await page.click("#mode-sequence");
    await page.fill("#seq", ">slow_one\n" + PROTEIN);
    check((await page.text("#seq-count")).startsWith("70 residues"), "the residues are counted while typing");
    await page.click("#submit");
    await page.waitFor("!document.querySelector('#home-error').hidden", 5000, "the request to choose a starting genome");
    ok("without a starting genome the form asks for one");
    await page.fill("#home", "homus_primus");
    await page.waitFor("document.querySelectorAll('#species-list input:checked').length === 4", 15000, "the species list");
    await page.click("#submit");
    await page.waitFor("document.querySelector('.badge.running')", 30000, "the running search");
    await page.waitFor("document.querySelector('.steps') && document.querySelector('.steps').innerText.includes('1 of 4 genomes done')", 30000, "the progress");
    ok("a running search shows its steps and the genome count");
    const running = await page.evaluate("new URLSearchParams(location.search).get('id')");
    check(await page.evaluate("!!document.querySelector('.button.danger')"), "the browser that started it may cancel it");
    // The page asks the server every 5 s. Focus on a button must survive that.
    await page.evaluate("document.querySelector('#cancel-search').focus(); window.__button = document.querySelector('#cancel-search'); true");
    await sleep(11000);
    check(await page.evaluate("document.activeElement && document.activeElement.id === 'cancel-search'"), "the keyboard focus stays on a button while the page updates");
    check(await page.evaluate("window.__button === document.querySelector('#cancel-search')"), "an unchanged view is not drawn again");
    check(/Running for \d+ s/.test(await page.text(".steps ~ p")), "the running time keeps counting");

    const other = await chrome.newPage();          // someone else with the link
    await other.goto(base + "/result.html?id=" + running);
    await other.waitFor("document.querySelector('.badge.running')", 20000, "the search in a second browser tab");
    ok("the link works for anyone who has it");

    await page.send("Page.bringToFront");
    await page.evaluate("window.confirm = () => true");
    await page.click(".button.danger");
    await page.waitFor("document.querySelector('h1').innerText === 'This search was cancelled'", 30000, "the cancelled view");
    check((await api(base, "/api/v1/searches/" + running)).body.state === "cancelled", "the server has cancelled it");
    await other.send("Page.bringToFront");
    await other.waitFor("document.querySelector('h1').innerText === 'This search was cancelled'", 30000, "the cancelled view in the other tab");
    ok("the other tab notices by itself");

    console.log("What the form refuses");
    await page.send("Page.bringToFront");
    await page.goto(base + "/search.html");
    await page.waitFor("document.querySelector('#home').options.length > 1", 15000, "the form");
    await page.fill("#acc", "NOPE");
    await page.click("#lookup");
    await page.waitFor("!document.querySelector('#acc-error').hidden", 10000, "the message for a wrong accession");
    check((await page.text("#acc-error")).includes("not a UniProt or NCBI"), "a wrong accession is explained");
    await page.click("#mode-sequence");
    await page.fill("#seq", "ACGTTGCAAC".repeat(8));
    await page.fill("#home", "homus_primus");
    await page.waitFor("document.querySelectorAll('#species-list input:checked').length === 4", 15000, "the species list");
    await page.click("#submit");
    await page.waitFor("!document.querySelector('#seq-error').hidden", 10000, "the message for DNA");
    check((await page.text("#seq-error")).includes("DNA"), "DNA is refused with a reason, next to the field");
    check(await page.evaluate("location.pathname.endsWith('/search.html')"), "and the form stays where it is");
    await page.fill("#seq", PROTEIN);
    await page.click("#pick-none");
    await page.click("#submit");
    await page.waitFor("!document.querySelector('#targets-error').hidden", 5000, "the message for no genome");
    ok("no genome ticked: asked to choose one");
    check(await page.evaluate("document.querySelector('#submit').disabled === false"), "the button can be used again after an error");

    console.log("Links that lead nowhere");
    await page.goto(base + "/result.html?id=Aa0000000000000000000a");
    await page.waitFor("document.querySelector('h1').innerText === 'Search not found'", 15000, "the not-found view");
    ok("an unknown search link says so");
    await page.goto(base + "/result.html?id=%3Cscript%3Ealert(1)%3C/script%3E");
    await page.waitFor("document.querySelector('h1').innerText === 'Search not found'", 15000, "the not-found view");
    ok("a malformed link is not sent to the server");
    await page.goto(base + "/result.html?example=../../etc/passwd");
    await page.waitFor("document.querySelector('h1').innerText === 'Example not found'", 15000, "the missing example view");
    ok("an example name with a path in it is refused");

    if (cards > 0) {
      console.log("An example");
      await page.goto(base + "/");
      await page.waitFor("document.querySelector('#example-cards a.card')", 15000, "the example cards");
      await page.click("#example-cards a.card");
      await page.waitFor("document.querySelectorAll('table.calls tbody tr').length > 0", 20000, "the example");
      await page.waitFor("(() => { const i = document.querySelector('.figure-frame img'); return i && i.complete && i.naturalWidth > 0; })()", 15000, "the example's figure");
      ok("an example opens with its figure: " + (await page.text("h1")));
    }

    console.log("The server goes away while a search is watched");
    await page.goto(base + "/search.html");
    await page.waitFor("document.querySelector('#home').options.length > 1", 15000, "the form");
    await page.click("#mode-sequence");
    await page.fill("#seq", ">slow_two\n" + PROTEIN);
    await page.fill("#home", "homus_primus");
    await page.waitFor("document.querySelectorAll('#species-list input:checked').length === 4", 15000, "the species list");
    await page.click("#submit");
    await page.waitFor("document.querySelector('.badge.running')", 30000, "the running search");
    const watched = await page.evaluate("new URLSearchParams(location.search).get('id')");
    page.expectFailures = true;                    // failed requests are the point of this part
    service.kill("SIGKILL");
    await page.waitFor("document.querySelector('.notice.warn') && document.querySelector('.notice.warn').innerText.includes('Trying again')", 30000, "the notice that the server does not answer");
    check(await page.evaluate("!!document.querySelector('.steps')"), "the last known progress stays on the page");
    service = await withService(path.join(scratch, "service"), port);
    await page.waitFor("document.querySelector('.badge.running') && !document.querySelector('.notice.warn')", 90000, "the page to recover");
    ok("when the server is back the page carries on by itself");
    check((await api(base, "/api/v1/searches/" + watched)).body.attempt === 1, "and the search itself was never interrupted");
    page.expectFailures = false;

    console.log("A species that is not on the server");
    await page.send("Page.bringToFront");
    await page.goto(base + "/search.html?acc=P01501");
    await page.waitFor("document.querySelectorAll('#species-list input:checked').length === 4", 15000, "the form");
    check(!(await page.evaluate("document.querySelector('#ncbi-pane').hidden")), "the form offers to fetch a species that is not listed");
    await page.fill("#ncbi-term", "Gigas maximus");
    await page.click("#ncbi-find");
    await page.waitFor("document.querySelector('#ncbi-offers').innerText.includes('Gigas maximus') && !document.querySelector('#ncbi-find').disabled", 15000, "what NCBI has");
    check((await page.text("#ncbi-offers")).includes("above the limit") && !(await page.evaluate("!!document.querySelector('#ncbi-offers button')")),
      "a genome above the size limit is shown and cannot be added: " + (await page.text("#ncbi-offers")).replace(/\s+/g, " ").slice(0, 140));
    await page.fill("#ncbi-term", "Nomen nescio");
    await page.click("#ncbi-find");
    await page.waitFor("!document.querySelector('#ncbi-error').hidden", 15000, "the answer for an unknown name");
    ok("a name NCBI does not know is said to be unknown");
    await page.fill("#ncbi-term", "Formica nova");
    await page.click("#ncbi-find");
    await page.waitFor("document.querySelector('#ncbi-offers button')", 15000, "the genome NCBI offers");
    check((await page.text("#ncbi-offers")).includes("GCA_000000101.1"), "NCBI's genome for the species is offered with its accession");
    await page.click("#ncbi-offers button");
    await page.waitFor("document.querySelector('#species-list .tag')", 5000, "the added genome in the list");
    check((await page.text("#count")).startsWith("5 of 10"), "the added genome counts towards the limit");
    await page.shot("form-with-ncbi-genome");
    await page.click("#submit");
    await page.waitFor("location.pathname.endsWith('/result.html')", 15000, "the move to the search page");
    await page.waitFor("document.querySelector('#coming-genomes') && document.querySelector('#coming-genomes').innerText.includes('Formica nova')", 20000, "the note about the genome that is fetched");
    check((await page.text("h2")).startsWith("Getting a genome first"), "a search that waits for a genome says which one");
    await page.shot("waiting-for-a-genome");
    await page.waitFor("document.querySelectorAll('table.calls tbody tr').length === 5", 90000, "the result with the fetched genome");
    const fiveRows = await page.evaluate("[...document.querySelectorAll('table.calls tbody tr')].map((r) => r.innerText.replace(/\\s+/g, ' '))");
    check(fiveRows.some((t) => t.includes("Formica nova")), "the fetched genome is in the result");
    check(fiveRows.some((t) => t.includes("None of the 6 neighbouring genes was found either")),
      "a genome in which nothing was found has its row and says so");
    await page.goto(base + "/search.html?acc=P01501");
    await page.waitFor("[...document.querySelectorAll('#species-list label')].some((l) => l.innerText.includes('Formica nova') && l.innerText.includes('from NCBI') && l.querySelector('input').checked)", 15000, "the fetched genome in the list");
    ok("a genome fetched once is on the server for the next search, listed with the others");

    console.log("A starting genome from NCBI");
    const chosenHome = "!document.querySelector('#home-chosen').hidden && document.querySelector('#home-chosen').innerText.includes('Apis nova')";
    const noHome = "document.querySelector('#home-chosen').hidden && document.querySelector('#home').value === '' && [...document.querySelector('#home').options].every((o) => !o.value.startsWith('ncbi:'))";
    const useIt = "[...document.querySelectorAll('#home-offers button')].find((b) => b.innerText === 'Use as starting genome').click()";
    await page.goto(base + "/search.html");
    await page.waitFor("document.querySelector('#home').options.length > 1", 15000, "the form");
    await page.click("#mode-sequence");
    await page.fill("#seq", ">home_from_ncbi\n" + PROTEIN);
    await page.waitFor("!document.querySelector('#home-ncbi').hidden", 5000, "the field for a starting genome from NCBI");
    await page.fill("#home-term", "Apis nova");
    await page.click("#home-find");
    await page.waitFor("document.querySelector('#home-offers button')", 15000, "the genome NCBI offers as starting genome");
    await page.evaluate(useIt);
    await page.waitFor(chosenHome, 5000, "the chosen starting genome");
    ok("a genome at NCBI can be made the starting genome");
    check(await page.evaluate("document.querySelector('#home').selectedOptions[0].text.startsWith('Apis nova (fetched from NCBI')"),
      "the list of starting genomes says which one is chosen");
    await page.waitFor("document.querySelectorAll('#species-list input').length > 0", 10000, "the species to search");
    await page.shot("home-from-ncbi");
    await page.click("#home-remove");
    await page.waitFor(noHome, 5000, "the starting genome to be taken back");
    check(await page.evaluate("!document.querySelector('#species-empty').hidden && document.querySelectorAll('#species-list input:checked').length === 0"),
      "Remove takes the starting genome back: none is chosen, no species are ticked");
    check(await page.evaluate("[...document.querySelectorAll('#home-offers button')].map((b) => b.innerText).join('|')") === "Use as starting genome",
      "it can be chosen again");
    await page.evaluate(useIt);
    await page.waitFor(chosenHome, 5000, "the starting genome chosen again");
    check(await page.evaluate("[...document.querySelectorAll('#home-offers button')].map((b) => b.innerText).join('|')") === "Remove",
      "the row of the chosen genome says so at once and offers to remove it");
    await page.evaluate("[...document.querySelectorAll('#home-offers button')].find((b) => b.innerText === 'Remove').click()");
    await page.waitFor(noHome, 5000, "the starting genome to be taken back in its row");
    ok("that Remove button takes it back as well");
    await page.evaluate(useIt);
    await page.waitFor(chosenHome, 5000, "the starting genome chosen a third time");
    await page.fill("#home", "");
    await page.waitFor(noHome, 5000, "the starting genome to be taken back through the list");
    ok("choosing no species in the list takes it back too");
    await page.evaluate(useIt);
    await page.waitFor(chosenHome, 5000, "the starting genome chosen a fourth time");
    await page.fill("#home", "homus_primus");
    await page.waitFor("document.querySelector('#home-chosen').hidden && document.querySelector('#home').value === 'homus_primus' && document.querySelectorAll('#species-list input:checked').length >= 4", 15000, "the genome of the server as starting genome");
    ok("a genome of the server chosen in the list replaces it");
    // A choice made for one protein must not stay when another protein is looked up.
    await page.goto(base + "/search.html?acc=Q00001");
    await page.waitFor("[...document.querySelectorAll('#home-offers button')].some((b) => b.innerText === 'Use as starting genome')", 20000, "the genome of the protein's own species, offered");
    await page.evaluate(useIt);
    await page.waitFor(chosenHome, 5000, "the starting genome for the first protein");
    await page.fill("#acc", "Q00002");
    await page.click("#lookup");
    await page.waitFor("document.querySelector('#lookup-result').innerText.includes('Formica nova')", 15000, "the second protein");
    await page.waitFor(noHome, 5000, "the starting genome of the first protein to be dropped");
    ok("a starting genome chosen for one protein is not kept for a protein of another species");

    console.log("A small screen");
    await page.send("Emulation.setDeviceMetricsOverride", { width: 400, height: 800, deviceScaleFactor: 1, mobile: true });
    await page.goto(base + "/");
    await page.waitFor("getComputedStyle(document.querySelector('.nav-toggle')).display !== 'none'", 10000, "the menu button");
    check(await page.evaluate("getComputedStyle(document.querySelector('.site-nav')).display === 'none'"), "the links fold into a menu");
    await page.click(".nav-toggle");
    check(await page.evaluate("getComputedStyle(document.querySelector('.site-nav')).display !== 'none' && document.querySelector('.nav-toggle').getAttribute('aria-expanded') === 'true'"), "the menu opens");
    check(await page.evaluate("document.documentElement.scrollWidth <= window.innerWidth + 1"), "the page does not scroll sideways");
    await page.shot("small-screen");
    await page.send("Emulation.clearDeviceMetricsOverride");

    console.log("Pages without a service");
    const copy = path.join(scratch, "pages");
    cpSync(path.join(WEB, "site"), copy, { recursive: true });
    writeFileSync(path.join(copy, "config.js"), "window.SYNVOY = { apiBase: null };\n");
    const staticPort = await freePort();
    staticServer = spawn(PYTHON, ["-m", "http.server", String(staticPort), "--bind", "127.0.0.1", "--directory", copy], { stdio: "ignore" });
    await sleep(1200);
    const plain = "http://127.0.0.1:" + staticPort;
    const solo = await chrome.newPage();
    await solo.goto(plain + "/index.html");
    await solo.waitFor("document.querySelector('[data-if-no-service]') && !document.querySelector('[data-if-no-service]').hidden", 15000, "the note that this copy has no service");
    check(await solo.evaluate("[...document.querySelectorAll('[data-if-service]')].every((n) => n.hidden)"), "everything that needs the service is hidden");
    if (cards > 0) {
      await solo.waitFor("document.querySelector('#example-cards a.card')", 15000, "the example cards");
      await solo.click("#example-cards a.card");
      await solo.waitFor("document.querySelectorAll('table.calls tbody tr').length > 0", 20000, "the example");
      ok("the examples work from plain files");
    }
    await solo.goto(plain + "/search.html");
    await solo.waitFor("!document.querySelector('#unavailable').hidden && document.querySelector('#search-form').hidden", 15000, "the notice on the search page");
    ok("the search page says that searches cannot be started there");
    check(solo.problems.length === 0, "no errors in the browser on the pages without a service");

    const problems = [...page.problems, ...other.problems];
    if (problems.length) throw new Error("the browser reported:\n  " + problems.join("\n  "));
    ok("no script errors and no blocked requests in the browser");
    console.log("\n" + steps.length + " checks passed");
  } finally {
    chrome.close();
    service.kill("SIGKILL");
    if (staticServer) staticServer.kill("SIGKILL");
    // Stand-in pipelines of the searches left running: stop them by their marker.
    spawnSync(PYTHON, ["-c", "import os,signal\nfor n in os.listdir('/proc'):\n    if not n.isdigit(): continue\n    try: data = open('/proc/%s/environ' % n, 'rb').read()\n    except OSError: continue\n"
      + "    if b'SYNVOY_WEB_JOB=' in data and " + JSON.stringify(scratch) + ".encode() in open('/proc/%s/cmdline' % n, 'rb').read():\n        os.kill(int(n), signal.SIGKILL)\n"]);
    try { rmSync(scratch, { recursive: true, force: true }); } catch (e) { /* best effort */ }
  }
}

main().then(() => process.exit(0), (error) => {
  console.error("\nFAILED: " + error.message);
  process.exit(1);
});
