/* ===========================================================================
   Log-normal Outlier Explorer

   Four views over one pre-exported simulation. The whole site is static: the
   two JSON files below are the entire index, and each replicate's per-
   observation columns arrive as a ~23 KB binary fetched on demand.

   Regenerate data/ with `Rscript export_data.R`; meta.json carries the byte
   layout so nothing here hard-codes an offset.

   Adapted from the linear-regression explorer. The design here is
   (pct, delta, sim), the parameters are (mu, sigma2), there is no OLS and no
   x column, and the likelihood grid has infeasible cells that are shaded.
   =========================================================================== */

const d3 = window.d3;

/* ---- palette --------------------------------------------------------------
   Clay is spoken for: it marks outliers and true parameter values, the two
   things every view is ultimately about. The two fits get cool hues so a
   fitted rule never reads as an outlier. */
const PAPER = "#fdfcf9";
const PANEL = "#f2eee5";
const INK   = "#1f1d1a";
const MUTED = "#1f1d1a";
const RULE  = "#1f1d1a";
const AXIS  = "#1f1d1a";
const CLAY  = "#b8452f";
const CALM  = "#1f1d1a";   // ordinary observations
const SHADE = "#dedad0";   // infeasible likelihood cells

const METHODS = [
  { key: "gel_bias", label: "ETEL, biased start", colour: "#1a8f6a" },
  { key: "gel",      label: "ETEL, true start",   colour: "#7b3fa0" }
];

/* Filled in at boot once meta is known: truth comes from the export. */
const PARAMS = [
  { key: "mu",     label: "μ̂",  tex: "$\\hat{\\mu}$",      truth: 0 },
  { key: "sigma2", label: "σ̂²", tex: "$\\hat{\\sigma}^2$", truth: 0 }
];

const VIEWS = [
  { key: "sample",     label: "Simulated data" },
  { key: "estimates",  label: "Density of the estimates" },
  { key: "weights",    label: "Inner products and implied probabilities" },
  { key: "likelihood", label: "Likelihood surface" }
];

const f2 = d3.format(".2f");
const f3 = d3.format(".3f");
const fg = d3.format("~g");
const f1 = d3.format(".1f");
const fe = d3.format(".3e");

/* ---- state ---------------------------------------------------------------- */
const state = {
  view: "sample",
  iPct: 0,
  iDelta: 0,
  iSim: 0,
  deltaLo: 0,
  deltaHi: 0,
  methods: new Set(METHODS.map(m => m.key)),
  gel: "gel_bias",
  trim: true,
  logY: false,
  logWeights: true,
  llWindow: 400,
  hover: null
};

let meta = null;
let coefs = null;
let nPct, nDelta, nSim, nObs, nLL, nMu, nS2, dMu, dS2;

/* Handlers a view installs so hovering one panel can highlight the same
   observation in its siblings without re-rendering a thousand circles. */
let hoverHandlers = [];
function setHover(i) {
  if (state.hover === i) return;
  state.hover = i;
  hoverHandlers.forEach(fn => fn(i));
}

/* ===========================================================================
   data access
   =========================================================================== */

const slot = (iPct, iDelta, iSim) => (iPct * nDelta + iDelta) * nSim + iSim;

const repCache = new Map();

function repUrl(iPct, iDelta, iSim) {
  return `data/rep/${meta.pct[iPct]}-${meta.delta_slug[iDelta]}-${meta.sim[iSim]}.bin`;
}

/* Cached by URL, and the *promise* is cached rather than its value, so two
   overlapping renders of the same replicate share one request. */
function loadRep(iPct, iDelta, iSim) {
  const url = repUrl(iPct, iDelta, iSim);
  let p = repCache.get(url);
  if (!p) {
    p = fetch(url)
      .then(r => {
        if (!r.ok) throw new Error(`${url}: ${r.status}`);
        return r.arrayBuffer();
      })
      .then(decodeRep);
    repCache.set(url, p);
  }
  return p;
}

function decodeRep(buf) {
  const L = meta.layout;
  return {
    y: new Float32Array(buf, L.y, nObs),
    gtRaw: {
      gel_bias: new Float32Array(buf, L.gt_bias, nObs),
      gel:      new Float32Array(buf, L.gt, nObs)
    },
    ll:      new Float32Array(buf, L.ll, nLL),
    hull:    new Uint8Array(buf, L.hull, nLL),
    outlier: new Uint8Array(buf, L.outlier, nObs),
    _pt: {}, _lpt: {}
  };
}

/* log pi_i = (lambda' g_i) - logsumexp(lambda' g), exactly as the R fitter
   computes it, so only the inner products are shipped. The log version is
   kept separately: a collapsed fit puts weights of 1e-300 and below on most
   points, which exp() turns into an honest zero, and the log-scale panels
   want the number before that happens. */
function logImplied(rec, key) {
  if (rec._lpt[key]) return rec._lpt[key];
  const g = rec.gtRaw[key];
  const mx = d3.max(g);
  let s = 0;
  for (let i = 0; i < nObs; i++) s += Math.exp(g[i] - mx);
  const lse = mx + Math.log(s);
  const lp = new Float64Array(nObs);
  for (let i = 0; i < nObs; i++) lp[i] = g[i] - lse;
  rec._lpt[key] = lp;
  return lp;
}

function implied(rec, key) {
  if (rec._pt[key]) return rec._pt[key];
  const lp = logImplied(rec, key);
  const w = Float64Array.from(lp, v => Math.exp(v));
  rec._pt[key] = w;
  return w;
}

/* Deviation on the log scale from the *true* log-mean, not from any fit: it
   is the only measure of outlyingness that does not already depend on the
   estimator being judged, and the contamination is a shift on this scale. */
function logDev(rec) {
  if (rec._dev) return rec._dev;
  const mu = meta.theta_true.mu;
  const r = new Float64Array(nObs);
  for (let i = 0; i < nObs; i++) r[i] = Math.log(rec.y[i]) - mu;
  rec._dev = r;
  return r;
}

/* The replicate estimates at one (pct, delta) for one method and parameter,
   with any failed fit dropped. */
function coefRun(mkey, param, iPct, iDelta) {
  const a = coefs[mkey][param];
  const from = slot(iPct, iDelta, 0);
  return a.slice(from, from + nSim).filter(Number.isFinite);
}

const trueMean = () => Math.exp(meta.theta_true.mu + meta.theta_true.sigma2 / 2);

function fitAt(mkey, s) {
  const mu = coefs[mkey].mu[s], s2 = coefs[mkey].sigma2[s];
  if (!Number.isFinite(mu) || !Number.isFinite(s2)) return null;
  return { mu, s2, ll: coefs[mkey].ll[s], conv: coefs[mkey].conv[s] };
}

const fitMean = f => (f ? Math.exp(f.mu + f.s2 / 2) : null);

const activeMethods = () => METHODS.filter(m => state.methods.has(m.key));

/* ===========================================================================
   chart plumbing
   =========================================================================== */

function frame(sel, w, h) {
  const host = d3.select(sel);
  host.selectAll("*").remove();
  const svg = host.append("svg")
    .attr("viewBox", `0 0 ${w} ${h}`)
    .attr("width", "100%")
    .style("height", "auto")
    .style("overflow", "visible");
  addDownloadButtons(host, svg, sel.replace("#", ""));
  return svg;
}

/* ---- figure download ------------------------------------------------------
   Axis labels are MathJax SVG, which survives export as-is. Tick labels and
   any foreignObject fallbacks are swapped for plain <text> with a Unicode
   rendering of the same label, so the file opens anywhere. */

const TEX_SYM = {
  "\\mu": "μ", "\\sigma": "σ", "\\delta": "δ", "\\pi": "π", "\\lambda": "λ",
  "\\epsilon": "ε", "\\top": "⊤", "\\log": "log", "\\exp": "exp"
};
const SUBDIG = { 0: "₀", 1: "₁", 2: "₂", 3: "₃", 4: "₄",
                 5: "₅", 6: "₆", 7: "₇", 8: "₈", 9: "₉" };
const SUPDIG = { 0: "⁰", 1: "¹", 2: "²", 3: "³", 4: "⁴",
                 5: "⁵", 6: "⁶", 7: "⁷", 8: "⁸", 9: "⁹" };

function texToUnicode(src) {
  let s = src;
  s = s.replace(/\\text\{([^}]*)\}/g, "$1");
  s = s.replace(/\\mathrm\{([^}]*)\}/g, "$1");
  s = s.replace(/\\hat\{([^}]*)\}/g, "$1̂");
  s = s.replace(/\\[;,!]/g, " ");
  Object.keys(TEX_SYM).forEach(k => { s = s.split(k).join(TEX_SYM[k]); });
  let prev;
  do { prev = s; s = s.replace(/_\{([^{}]*)\}/, "_$1"); } while (s !== prev);
  do { prev = s; s = s.replace(/\^\{([^{}]*)\}/, "^$1"); } while (s !== prev);
  s = s.replace(/_([0-9])/g, (_, d) => SUBDIG[d]);
  s = s.replace(/\^([0-9])/g, (_, d) => SUPDIG[d]);
  return s.replace(/[{}$\\]/g, "").replace(/\s+/g, " ").trim();
}

function svgString(node, opts = {}) {
  const NS = "http://www.w3.org/2000/svg";
  const clone = node.cloneNode(true);

  if (!opts.keepMath) {
    clone.querySelectorAll("foreignObject[data-tex]").forEach(fo => {
      const t = document.createElementNS(NS, "text");
      t.setAttribute("transform", fo.getAttribute("transform") || "");
      t.setAttribute("text-anchor", fo.getAttribute("data-anchor") || "middle");
      t.setAttribute("dominant-baseline", "middle");
      t.setAttribute("font-size", 11);
      t.setAttribute("font-family", "Georgia, serif");
      t.setAttribute("fill", MUTED);
      t.textContent = texToUnicode(fo.getAttribute("data-tex"));
      fo.parentNode.replaceChild(t, fo);
    });
  }

  const vb = (clone.getAttribute("viewBox") || "0 0 640 470").split(/\s+/).map(Number);
  clone.setAttribute("xmlns", NS);
  clone.setAttribute("width", vb[2]);
  clone.setAttribute("height", vb[3]);

  // Transparent PNGs are a nuisance to place in a document.
  const bg = document.createElementNS(NS, "rect");
  bg.setAttribute("x", vb[0]); bg.setAttribute("y", vb[1]);
  bg.setAttribute("width", vb[2]); bg.setAttribute("height", vb[3]);
  bg.setAttribute("fill", "#ffffff");
  const BGS = new Set([PAPER, PANEL, "#f5f1e7", "#f6f3ec", "#fdfcf9"]
    .map(c => c.toLowerCase()));
  clone.querySelectorAll("rect").forEach(el => {
    const f = (el.getAttribute("fill") || "").toLowerCase();
    if (BGS.has(f)) el.setAttribute("fill", "#ffffff");
  });
  clone.insertBefore(bg, clone.firstChild);

  return { text: new XMLSerializer().serializeToString(clone), w: vb[2], h: vb[3] };
}

function saveBlob(blob, name) {
  const url = URL.createObjectURL(blob);
  const a = document.createElement("a");
  a.href = url; a.download = name;
  document.body.appendChild(a); a.click(); a.remove();
  setTimeout(() => URL.revokeObjectURL(url), 1000);
}

/* PDF by printing rather than converting: the browser already renders this
   figure correctly, and its PDF output is vector. */
function printPdf(svg, id) {
  const { text, w, h } = svgString(svg.node(), { keepMath: true });
  const win = window.open("", "_blank");
  if (!win) {
    alert("Allow pop-ups for this site to export PDF.");
    return;
  }
  win.document.write(
    `<!doctype html><html><head><meta charset="utf-8">` +
    `<title>${figureName(id)}</title>` +
    `<style>` +
    `@page { size: ${(w / 96).toFixed(3)}in ${(h / 96).toFixed(3)}in; margin: 0; }` +
    `html, body { margin: 0; padding: 0; background: #fff; ` +
    `-webkit-print-color-adjust: exact; print-color-adjust: exact; }` +
    `svg { display: block; width: ${w}px; height: ${h}px; }` +
    `</style></head><body>${text}</body></html>`
  );
  win.document.close();
  const go = () => { win.focus(); win.print(); };
  if (win.document.fonts && win.document.fonts.ready) {
    win.document.fonts.ready.then(() => setTimeout(go, 250));
  } else {
    setTimeout(go, 800);
  }
}

/* Name files after the cell they came from, so a folder of exports is still
   readable in a month. */
function figureName(id) {
  const pct = meta.pct[state.iPct];
  const dl = meta.delta_slug[state.iDelta];
  return `${id}-${pct}pct-delta${dl}-rep${meta.sim[state.iSim]}`;
}

function addDownloadButtons(host, svg, id) {
  const bar = host.append("div").attr("class", "dl");

  bar.append("button").attr("type", "button").text("SVG")
    .on("click", () => {
      const { text } = svgString(svg.node());
      saveBlob(new Blob([text], { type: "image/svg+xml;charset=utf-8" }),
               `${figureName(id)}.svg`);
    });

  bar.append("button").attr("type", "button").text("PNG")
    .on("click", () => {
      const { text, w, h } = svgString(svg.node());
      const scale = 3;                       // 3× so it holds up in print
      const img = new Image();
      img.onload = () => {
        const c = document.createElement("canvas");
        c.width = w * scale; c.height = h * scale;
        const ctx = c.getContext("2d");
        ctx.scale(scale, scale);
        ctx.drawImage(img, 0, 0);
        c.toBlob(b => saveBlob(b, `${figureName(id)}.png`), "image/png");
      };
      img.src = "data:image/svg+xml;base64," +
                btoa(unescape(encodeURIComponent(text)));
    });

  bar.append("button").attr("type", "button").text("PDF")
    .on("click", () => printPdf(svg, id));
}

/* LaTeX to native SVG paths. Nothing here depends on a font being installed,
   so the same output works on screen and in a downloaded SVG. EX converts
   MathJax's ex units to our pixel sizes. */
function mathSvg(tex, px, colour) {
  if (!window.MathJax || !MathJax.tex2svg) return null;
  const out = MathJax.tex2svg(tex, { display: false });
  const node = out.querySelector("svg");
  if (!node) return null;

  const EX = px * 0.45;
  const w = parseFloat(node.getAttribute("width")) * EX;
  const h = parseFloat(node.getAttribute("height")) * EX;
  node.setAttribute("width", w);
  node.setAttribute("height", h);
  node.removeAttribute("style");
  node.setAttribute("fill", colour);
  node.querySelectorAll('[fill="currentColor"]').forEach(el =>
    el.setAttribute("fill", colour));
  return { node, w, h };
}

/* Target text sizes, in real screen pixels. */
const TICK_PX = 13;
const LABEL_PX = 20;

function unitOf(node) {
  const el = node.node ? node.node() : node;
  const svg = el.ownerSVGElement || el;
  return +svg.getAttribute("data-unit") || 1;
}

function styleAxes(svg) {
  const px = TICK_PX * unitOf(svg);
  svg.selectAll(".domain, .tick line").attr("stroke", AXIS);
  svg.selectAll(".tick text")
    .attr("fill", MUTED)
    .attr("font-size", px);
  mathTicks(svg, px);
  return svg;
}

/* Tick labels through the same path, so digits match the labels and survive
   export. Position comes from the text's own box rather than per-axis rules. */
function mathTicks(svg, px = 12) {
  svg.selectAll(".tick text").each(function () {
    const raw = this.textContent.replace(/,/g, "").replace(/−/g, "-").trim();
    if (!raw) return;
    let box;
    try { box = this.getBBox(); } catch (e) { return; }
    const m = mathSvg(raw, px, MUTED);
    if (!m) return;
    m.node.setAttribute("x", box.x + box.width / 2 - m.w / 2);
    m.node.setAttribute("y", box.y + box.height / 2 - m.h / 2);
    this.parentNode.appendChild(m.node);
    this.remove();
  });
}

function axisLabel(g, x, y, text, anchor = "middle", rotate = 0) {
  const px = LABEL_PX * unitOf(g);
  const isTex = text.startsWith("$") && text.endsWith("$") && text.length > 2;
  const m = isTex ? mathSvg(text.slice(1, -1), px, MUTED) : null;

  if (!m) {
    g.append("text")
      .attr("transform", `translate(${x},${y}) rotate(${rotate})`)
      .attr("text-anchor", anchor)
      .attr("font-size", px)
      .attr("fill", MUTED)
      .text(isTex ? texToUnicode(text.slice(1, -1)) : text);
    return;
  }

  const dx = anchor === "end" ? -m.w : anchor === "start" ? 0 : -m.w / 2;
  m.node.setAttribute("x", dx);
  m.node.setAttribute("y", -m.h / 2);
  g.append("g")
    .attr("transform", `translate(${x},${y}) rotate(${rotate})`)
    .node().appendChild(m.node);
}

function pad(extent, frac = 0.05) {
  const [lo, hi] = extent;
  const s = (hi - lo) || Math.abs(hi) || 1;
  return [lo - s * frac, hi + s * frac];
}

/* Gaussian kernel density on a fixed grid. Silverman's rule with the robust
   spread, which keeps a run of 50 estimates from being oversmoothed by the one
   replicate that blew up. */
function kde(values, grid) {
  const n = values.length;
  if (!n) return new Float64Array(grid.length);
  const sorted = Float64Array.from(values).sort(d3.ascending);
  const sd = d3.deviation(values) || 0;
  const iqr = d3.quantileSorted(sorted, 0.75) - d3.quantileSorted(sorted, 0.25);
  const spread = Math.min(sd || Infinity, (iqr / 1.349) || Infinity);
  const range = (grid[grid.length - 1] - grid[0]) || 1;
  const bw = Math.max(0.9 * (isFinite(spread) ? spread : range / 6) * Math.pow(n, -0.2),
                      range / 400);

  const out = new Float64Array(grid.length);
  const c = 1 / (n * bw * Math.sqrt(2 * Math.PI));
  for (let k = 0; k < n; k++) {
    const v = values[k];
    for (let i = 0; i < grid.length; i++) {
      const z = (grid[i] - v) / bw;
      if (z > -5 && z < 5) out[i] += Math.exp(-0.5 * z * z);
    }
  }
  for (let i = 0; i < grid.length; i++) out[i] *= c;
  return out;
}

function quantileOf(values, p) {
  const s = Float64Array.from(values).sort(d3.ascending);
  return d3.quantileSorted(s, p);
}

const dnorm = (x, mu, sd) =>
  Math.exp(-0.5 * ((x - mu) / sd) ** 2) / (sd * Math.sqrt(2 * Math.PI));

/* A hover marker that lives above the point cloud and is moved rather than
   redrawn. Returns the update function a view registers as a hover handler. */
function hoverRing(g, colour = INK) {
  const ring = g.append("circle")
    .attr("r", 7)
    .attr("fill", "none")
    .attr("stroke", colour)
    .attr("stroke-width", 1.4)
    .style("display", "none");
  return (i, cx, cy) => {
    if (i == null || cx == null || !isFinite(cx) || !isFinite(cy)) {
      ring.style("display", "none"); return;
    }
    ring.attr("cx", cx).attr("cy", cy).style("display", null);
  };
}

/* One line of fit summary, used by several readouts. */
function fitSummary(mm, f) {
  if (!f) return `<b style="color:${mm.colour}">${mm.label}</b> <span class="dim">fit failed</span>`;
  return `<b style="color:${mm.colour}">${mm.label}</b>` +
    ` <span class="dim">μ̂</span> ${f3(f.mu)}` +
    ` <span class="dim">σ̂²</span> ${f3(f.s2)}` +
    ` <span class="dim">ll</span> ${f1(f.ll)}` +
    (f.conv === 0 ? "" : ` <b style="color:${CLAY}">not converged</b>`);
}

/* ===========================================================================
   view: one sample
   =========================================================================== */

function renderSample(rec) {
  drawSampleScatter(rec);
  drawSampleCoefs(rec);
  drawSampleDens();
  sampleReadout(null, rec);
}

function drawSampleScatter(rec) {
  const W = 640, H = 470, m = { l: 62, r: 16, t: 10, b: 44 };
  const svg = frame("#sampleScatter", W, H);
  const xs = d3.scaleLinear([1, nObs], [m.l, W - m.r]);
  const ext = d3.extent(rec.y);
  const ys = state.logY
    ? d3.scaleLog([ext[0] / 1.2, ext[1] * 1.2], [H - m.b, m.t])
    : d3.scaleLinear([0, ext[1] * 1.04], [H - m.b, m.t]).nice();

  svg.append("rect")
    .attr("x", m.l).attr("y", m.t)
    .attr("width", W - m.r - m.l).attr("height", H - m.b - m.t)
    .attr("fill", PAPER).attr("stroke", RULE);

  svg.append("g").attr("transform", `translate(0,${H - m.b})`)
    .call(d3.axisBottom(xs).ticks(8).tickSize(4));
  svg.append("g").attr("transform", `translate(${m.l},0)`)
    .call(state.logY
      ? d3.axisLeft(ys).ticks(6, "~g").tickSize(4)
      : d3.axisLeft(ys).ticks(7).tickSize(4));
  axisLabel(svg, (m.l + W - m.r) / 2, H - 6, "Observation index");
  axisLabel(svg, 15, (m.t + H - m.b) / 2, state.logY ? "$y$ (log axis)" : "$y$", "middle", -90);

  const clip = "clip-sample";
  svg.append("clipPath").attr("id", clip).append("rect")
    .attr("x", m.l).attr("y", m.t)
    .attr("width", W - m.r - m.l).attr("height", H - m.b - m.t);
  const plot = svg.append("g").attr("clip-path", `url(#${clip})`);

  const hline = (v, colour, dash, label) => {
    if (v == null || !isFinite(ys(v))) return;
    plot.append("line")
      .attr("x1", m.l).attr("x2", W - m.r)
      .attr("y1", ys(v)).attr("y2", ys(v))
      .attr("stroke", colour)
      .attr("stroke-width", dash ? 1.4 : 1.8)
      .attr("stroke-dasharray", dash || null);
    if (label) {
      plot.append("text")
        .attr("x", W - m.r - 5).attr("y", ys(v) - 5)
        .attr("text-anchor", "end").attr("font-size", 10).attr("fill", colour)
        .attr("font-family", "ui-monospace, Menlo, monospace")
        .attr("paint-order", "stroke").attr("stroke", PAPER).attr("stroke-width", 3)
        .text(label);
    }
  };

  // Ordinary points first, outliers over them: with 10% contamination the clay
  // would otherwise be buried under the cloud it is supposed to explain.
  const idx = d3.range(nObs);
  plot.append("g").selectAll("circle")
    .data(idx.filter(i => !rec.outlier[i])).join("circle")
    .attr("cx", i => xs(i + 1)).attr("cy", i => ys(rec.y[i]))
    .attr("r", 2.4).attr("fill", CALM).attr("fill-opacity", 0.8);
  plot.append("g").selectAll("circle")
    .data(idx.filter(i => rec.outlier[i])).join("circle")
    .attr("cx", i => xs(i + 1)).attr("cy", i => ys(rec.y[i]))
    .attr("r", 3.8).attr("fill", CLAY).attr("fill-opacity", 0.9);

  const s = slot(state.iPct, state.iDelta, state.iSim);
  hline(trueMean(), INK, "6 4", "true mean");
  activeMethods().forEach(mm => {
    hline(fitMean(fitAt(mm.key, s)), mm.colour, null, null);
  });

  const ring = hoverRing(plot);
  const finder = d3.Delaunay.from(idx, i => xs(i + 1), i => ys(rec.y[i]));
  svg.append("rect")
    .attr("x", m.l).attr("y", m.t)
    .attr("width", W - m.r - m.l).attr("height", H - m.b - m.t)
    .attr("fill", "transparent")
    .on("pointermove", ev => {
      const [px, py] = d3.pointer(ev);
      setHover(finder.find(px, py));
    })
    .on("pointerleave", () => setHover(null));

  hoverHandlers.push(i => {
    ring(i, i == null ? null : xs(i + 1), i == null ? null : ys(rec.y[i]));
    sampleReadout(i, rec);
  });

  styleAxes(svg);
}

/* Strip plot: where this replicate's estimate falls among the 50 drawn at the
   same setting. One lane per method, two rows for the two parameters. */
function drawSampleCoefs(rec) {
  const W = 430, rowH = 104, gap = 16, m = { l: 10, r: 12, t: 28, b: 22 };
  const svg = frame("#sampleCoefs", W, PARAMS.length * rowH + gap);
  const shown = activeMethods();
  const s = slot(state.iPct, state.iDelta, state.iSim);

  PARAMS.forEach((p, r) => {
    const g = svg.append("g").attr("transform", `translate(0,${r * (rowH + gap)})`);
    const runs = shown.map(mm => coefRun(mm.key, p.key, state.iPct, state.iDelta));
    const all = runs.flat().concat([p.truth]);
    const dom = state.trim && all.length > 10
      ? pad([Math.min(p.truth, quantileOf(all, 0.02)), Math.max(p.truth, quantileOf(all, 0.98))], 0.08)
      : pad(d3.extent(all), 0.06);
    const xs = d3.scaleLinear(dom, [m.l + 4, W - m.r]).nice();

    g.append("g").attr("transform", `translate(0,${rowH - m.b})`)
      .call(d3.axisBottom(xs).ticks(5).tickSize(3));

    const clip = `clip-coefs-${r}`;
    g.append("clipPath").attr("id", clip).append("rect")
      .attr("x", m.l).attr("y", 0).attr("width", W - m.r - m.l).attr("height", rowH - m.b);
    const plot = g.append("g").attr("clip-path", `url(#${clip})`);

    plot.append("line")
      .attr("x1", xs(p.truth)).attr("x2", xs(p.truth))
      .attr("y1", m.t - 6).attr("y2", rowH - m.b)
      .attr("stroke", CLAY).attr("stroke-dasharray", "4 3");

    const laneH = (rowH - m.b - m.t) / Math.max(shown.length, 1);
    shown.forEach((mm, k) => {
      const cy = m.t + laneH * (k + 0.5);
      plot.append("g").selectAll("line").data(runs[k]).join("line")
        .attr("x1", v => xs(v)).attr("x2", v => xs(v))
        .attr("y1", cy - 5).attr("y2", cy + 5)
        .attr("stroke", mm.colour).attr("stroke-opacity", 0.38);
      const f = fitAt(mm.key, s);
      if (f) {
        const v = p.key === "mu" ? f.mu : f.s2;
        plot.append("circle")
          .attr("cx", xs(v)).attr("cy", cy)
          .attr("r", 3.8)
          .attr("fill", f.conv === 0 ? mm.colour : PAPER)
          .attr("stroke", mm.colour).attr("stroke-width", 1.6);
      }
      g.append("text")
        .attr("x", m.l).attr("y", cy - 8)
        .attr("font-size", 9.5).attr("fill", MUTED)
        .attr("font-family", "ui-monospace, Menlo, monospace")
        .text(mm.label);
    });

    g.append("text")
      .attr("x", m.l).attr("y", 11)
      .attr("font-size", 12).attr("fill", INK)
      .text(`${p.label}   true ${p.truth}`);
  });

  styleAxes(svg);
}

/* The two components of the mixture the sample is drawn from, on the raw
   scale: clean LN(mu, sigma2) against contaminated LN(mu + delta, sigma2).
   Mirrors plot_lognormal_dens() in 00_plot_lognormal_density.R, including
   its x range of 0 to exp(max(mu, mu + delta) + 3 sigma). */
const dlnorm = (x, mu, sd) =>
  x <= 0 ? 0 : Math.exp(-0.5 * ((Math.log(x) - mu) / sd) ** 2) / (x * sd * Math.sqrt(2 * Math.PI));

function drawSampleDens() {
  const W = 430, H = 210, m = { l: 48, r: 12, t: 12, b: 34 };
  const svg = frame("#sampleDens", W, H);
  const tt = meta.theta_true;
  const sd = Math.sqrt(tt.sigma2);
  const delta = meta.delta[state.iDelta];
  const xMax = Math.exp(Math.max(tt.mu, tt.mu + delta) + 3 * sd);
  const grid = d3.range(500).map(k => xMax * k / 499);

  const comps = [
    { label: "Clean",   colour: "#36648b", v: grid.map(x => dlnorm(x, tt.mu, sd)) },
    { label: "Outlier", colour: CLAY,      v: grid.map(x => dlnorm(x, tt.mu + delta, sd)) }
  ];

  const xs = d3.scaleLinear([0, xMax], [m.l, W - m.r]).nice();
  const ys = d3.scaleLinear([0, d3.max(comps, c => d3.max(c.v)) || 1], [H - m.b, m.t]).nice();

  svg.append("rect")
    .attr("x", m.l).attr("y", m.t)
    .attr("width", W - m.r - m.l).attr("height", H - m.b - m.t)
    .attr("fill", PAPER).attr("stroke", RULE);

  svg.append("g").attr("transform", `translate(0,${H - m.b})`)
    .call(d3.axisBottom(xs).ticks(6).tickSize(4));
  svg.append("g").attr("transform", `translate(${m.l},0)`)
    .call(d3.axisLeft(ys).ticks(4).tickSize(4).tickFormat(d3.format(".2~g")));
  axisLabel(svg, (m.l + W - m.r) / 2, H - 4, "$y$");
  axisLabel(svg, 14, (m.t + H - m.b) / 2, "Density", "middle", -90);

  const clip = "clip-dens";
  svg.append("clipPath").attr("id", clip).append("rect")
    .attr("x", m.l).attr("y", m.t - 2)
    .attr("width", W - m.r - m.l).attr("height", H - m.b - m.t + 2);
  const plot = svg.append("g").attr("clip-path", `url(#${clip})`);

  const line = d3.line().x((_, k) => xs(grid[k])).y(v => ys(v));
  comps.forEach(c => {
    plot.append("path").datum(c.v)
      .attr("d", line)
      .attr("fill", "none").attr("stroke", c.colour).attr("stroke-width", 1.8);
  });

  // Legend inside the panel, top right, in the same idiom as the R plot.
  const lg = svg.append("g").attr("transform", `translate(${W - m.r - 8},${m.t + 12})`);
  comps.forEach((c, k) => {
    const row = lg.append("g").attr("transform", `translate(0,${k * 15})`);
    row.append("line").attr("x1", -60).attr("x2", -44).attr("y1", 0).attr("y2", 0)
      .attr("stroke", c.colour).attr("stroke-width", 1.8);
    row.append("text").attr("x", -40).attr("y", 3.5)
      .attr("font-size", 10).attr("fill", INK)
      .attr("font-family", "ui-monospace, Menlo, monospace")
      .text(c.label);
  });

  styleAxes(svg);
}

function sampleReadout(i, rec) {
  const el = d3.select("#sampleReadout");
  const nOut = d3.sum(rec.outlier);
  const s = slot(state.iPct, state.iDelta, state.iSim);
  if (i == null) {
    const bits = [
      `<span class="dim">n</span> ${nObs}` +
      ` <span class="dim">outliers</span> <b style="color:${CLAY}">${nOut}</b>` +
      ` <span class="dim">(${f1(100 * nOut / nObs)}%)</span>`
    ];
    activeMethods().forEach(mm => bits.push(fitSummary(mm, fitAt(mm.key, s))));
    el.html(bits.join(`<span class="sep">|</span>`));
    return;
  }
  const d = logDev(rec)[i];
  el.html(
    `<span class="dim">obs</span> ${i + 1}` +
    `<span class="sep">|</span><span class="dim">y</span> ${f3(rec.y[i])}` +
    `<span class="sep">|</span><span class="dim">log y</span> ${f3(Math.log(rec.y[i]))}` +
    `<span class="sep">|</span><span class="dim">log y − μ</span> ${f3(d)}` +
    `<span class="sep">|</span>` +
    (rec.outlier[i]
      ? `<b style="color:${CLAY}">contaminated (shift ${fg(meta.delta[state.iDelta])})</b>`
      : `<span class="dim">ordinary</span>`)
  );
}

/* ===========================================================================
   view: sampling distributions
   =========================================================================== */

function renderEstimates() {
  drawRidges();
  drawEstPanels();
}

function drawRidges() {
  const W = 850, H = 580, gutter = 34;
  const facetW = (W - gutter) / 2;
  const svg = frame("#ridges", W, H);
  const shown = activeMethods();
  if (!shown.length) {
    svg.append("text").attr("x", W / 2).attr("y", H / 2)
      .attr("text-anchor", "middle").attr("fill", MUTED).attr("font-size", 14)
      .text("Select at least one method.");
    return;
  }

  const eLo = state.deltaLo, eHi = state.deltaHi;
  const rows = d3.range(eLo, eHi + 1);
  const nRows = rows.length;

  /* Both facets register their paths under the same (delta, method) key, so
     hovering a ridge in one panel lights up the same run in the other. */
  const ridgePaths = new Map();

  function applyHover(key) {
    ridgePaths.forEach((arr, k) => {
      const on = key == null || k === key;
      arr.forEach(o => {
        o.area.attr("fill-opacity", key == null ? 0.20 : (on ? 0.45 : 0.05));
        o.outline
          .attr("stroke-width", key != null && on ? 2.4 : 1.2)
          .attr("stroke-opacity", key == null ? 1 : (on ? 1 : 0.2));
      });
    });
  }

  function ridgeReadout(e, mkey) {
    const el = d3.select("#estReadout");
    if (e == null) { el.html(""); return; }
    const mm = METHODS.find(x => x.key === mkey);
    const bits = [
      `<b style="color:${mm.colour}">${mm.label}</b>`,
      `<span class="dim">δ</span> ${fg(meta.delta[e])}`
    ];
    PARAMS.forEach(p => {
      const v = coefRun(mkey, p.key, state.iPct, e);
      if (!v.length) { bits.push(`<span class="dim">${p.label}</span> no fits`); return; }
      const mean = d3.mean(v);
      const rmse = Math.sqrt(d3.mean(v, d => (d - p.truth) ** 2));
      bits.push(
        `<span class="dim">${p.label}</span> ${f3(mean)}` +
        ` <span class="dim">bias</span> ${f3(mean - p.truth)}` +
        ` <span class="dim">sd</span> ${f3(d3.deviation(v) || 0)}` +
        ` <span class="dim">rmse</span> ${f3(rmse)}`
      );
    });
    el.html(bits.join(`<span class="sep">|</span>`));
  }

  PARAMS.forEach((p, fi) => {
    const m = { l: fi === 0 ? 46 : 30, r: 12, t: 26, b: 58 };
    const x0 = fi * (facetW + gutter);
    const g = svg.append("g").attr("transform", `translate(${x0},0)`);

    const runs = [];
    rows.forEach(e => shown.forEach(mm => runs.push(coefRun(mm.key, p.key, state.iPct, e))));
    const pooled = runs.flat();
    if (!pooled.length) return;

    // A biased start sends a few fits a long way out; trimming keeps the bulk
    // legible, and the toggle exists because the excursions are the point at
    // the far end of the sweep.
    const dom = state.trim
      ? pad([Math.min(p.truth, quantileOf(pooled, 0.01)), Math.max(p.truth, quantileOf(pooled, 0.99))], 0.06)
      : pad(d3.extent(pooled.concat([p.truth])), 0.03);
    const xs = d3.scaleLinear(dom, [m.l, facetW - m.r]);
    const grid = d3.range(200).map(i => dom[0] + (dom[1] - dom[0]) * i / 199);

    const RIDGE = 1.8;                       // ridge height, in row heights
    const rowH = (H - m.t - m.b) / (nRows + RIDGE - 0.5);
    const rowY = e => H - m.b - ((e - eLo) + 0.5) * rowH;

    // One density per (delta, method), scaled against the tallest in the facet
    // so the ridges share a vertical unit -- ggridges' default.
    const dens = [];
    let top = 0;
    rows.forEach(e => {
      shown.forEach(mm => {
        const z = kde(coefRun(mm.key, p.key, state.iPct, e), grid);
        top = Math.max(top, d3.max(z));
        dens.push({ e, mm, z });
      });
    });

    const ridge = rowH * RIDGE;
    const ys = z => (top > 0 ? (z / top) * ridge : 0);

    g.append("rect")
      .attr("x", m.l).attr("y", m.t)
      .attr("width", facetW - m.r - m.l).attr("height", H - m.b - m.t)
      .attr("fill", PAPER).attr("stroke", RULE);

    const clip = `clip-ridge-${fi}`;
    g.append("clipPath").attr("id", clip).append("rect")
      .attr("x", m.l).attr("y", m.t)
      .attr("width", facetW - m.r - m.l).attr("height", H - m.b - m.t);
    const plot = g.append("g").attr("clip-path", `url(#${clip})`);

    rows.forEach(e => {
      plot.append("line")
        .attr("x1", m.l).attr("x2", facetW - m.r)
        .attr("y1", rowY(e)).attr("y2", rowY(e))
        .attr("stroke", RULE).attr("stroke-width", 0.7);
    });

    plot.append("line")
      .attr("x1", xs(p.truth)).attr("x2", xs(p.truth))
      .attr("y1", m.t).attr("y2", H - m.b)
      .attr("stroke", CLAY).attr("stroke-dasharray", "5 4").attr("stroke-width", 1.1);

    const area = d3.area()
      .x((_, i) => xs(grid[i]))
      .y1(v => -ys(v))
      .y0(0)
      .curve(d3.curveBasis);
    const outline = d3.line()
      .x((_, i) => xs(grid[i]))
      .y(v => -ys(v))
      .curve(d3.curveBasis);

    // Drawn from the top row down, so nearer (lower) ridges overlap the ones
    // behind them rather than being hidden by them.
    dens.sort((a, b) => d3.descending(a.e, b.e));
    dens.forEach(({ e, mm, z }) => {
      const row = plot.append("g").attr("transform", `translate(0,${rowY(e)})`);
      const aPath = row.append("path").datum(z)
        .attr("d", area)
        .attr("fill", mm.colour).attr("fill-opacity", 0.20)
        .attr("stroke", "none");
      const oPath = row.append("path").datum(z)
        .attr("d", outline)
        .attr("fill", "none")
        .attr("stroke", mm.colour).attr("stroke-width", 1.2);

      const key = `${e}|${mm.key}`;
      if (!ridgePaths.has(key)) ridgePaths.set(key, []);
      ridgePaths.get(key).push({ area: aPath, outline: oPath });
    });

    g.append("g").attr("transform", `translate(0,${H - m.b})`)
      .call(d3.axisBottom(xs).ticks(5).tickSize(4));

    axisLabel(g, (m.l + facetW - m.r) / 2, H - 14, p.tex);

    if (fi === 0) {
      g.append("line")
        .attr("x1", m.l).attr("x2", m.l)
        .attr("y1", m.t).attr("y2", H - m.b)
        .attr("stroke", AXIS);
      g.append("g").selectAll("text").data(rows).join("text")
        .attr("x", m.l - 8).attr("y", e => rowY(e) + 3.5)
        .attr("text-anchor", "end")
        .attr("font-size", TICK_PX * unitOf(svg))
        .attr("font-family", "ui-monospace, Menlo, monospace")
        .attr("fill", MUTED)
        .text(e => meta.delta[e]);
      axisLabel(g, 13, H / 2, "$\\delta$", "middle", -90);
    }
    g.append("text")
      .attr("x", m.l).attr("y", 14)
      .attr("font-size", 12).attr("fill", INK)
      .text(`${p.label}   true ${p.truth}`);

    /* Pick the ridge the pointer is actually inside, front row first, so the
       one drawn on top is the one that responds. Falling back to the nearest
       baseline means the gaps between ridges still select something. */
    g.append("rect")
      .attr("x", m.l).attr("y", m.t)
      .attr("width", facetW - m.r - m.l).attr("height", H - m.b - m.t)
      .attr("fill", "transparent")
      .on("pointermove", function (ev) {
        const [px, py] = d3.pointer(ev, g.node());
        let gi = Math.round((xs.invert(px) - dom[0]) / (dom[1] - dom[0]) * 199);
        gi = Math.max(0, Math.min(199, gi));

        const inside = dens.filter(d => {
          const base = rowY(d.e);
          return py <= base && py >= base - ys(d.z[gi]);
        });
        let best;
        if (inside.length) {
          const front = d3.max(inside, d => rowY(d.e));
          best = inside.filter(d => rowY(d.e) === front)
            .reduce((a, b) =>
              Math.abs(py - (rowY(a.e) - ys(a.z[gi]))) <=
              Math.abs(py - (rowY(b.e) - ys(b.z[gi]))) ? a : b);
        } else {
          best = dens.reduce((a, b) =>
            Math.abs(rowY(a.e) - py) <= Math.abs(rowY(b.e) - py) ? a : b);
        }
        applyHover(`${best.e}|${best.mm.key}`);
        ridgeReadout(best.e, best.mm.key);
      })
      .on("pointerleave", () => { applyHover(null); ridgeReadout(null); });
  });

  styleAxes(svg);
}

/* Bias and RMSE across the delta sweep, one small panel each per parameter.
   The delta grid skips 0, so the line is broken there rather than drawn
   across a setting that was never simulated. */
function drawEstPanels() {
  const W = 420, H = 126, gap = 18, m = { l: 58, r: 10, t: 16, b: 22 };
  const shown = activeMethods();
  const panels = [
    { stat: "bias", param: PARAMS[0], label: "bias  μ̂" },
    { stat: "bias", param: PARAMS[1], label: "bias  σ̂²" },
    { stat: "rmse", param: PARAMS[0], label: "RMSE  μ̂" },
    { stat: "rmse", param: PARAMS[1], label: "RMSE  σ̂²" }
  ];
  const svg = frame("#estPanels", W, panels.length * H + (panels.length - 1) * gap);
  if (!shown.length) return;

  const xs = d3.scaleLinear(pad(d3.extent(meta.delta), 0.04), [m.l, W - m.r]);
  const step = nDelta > 1 ? d3.min(d3.pairs(meta.delta), ([a, b]) => b - a) : 1;

  panels.forEach((p, r) => {
    const truth = p.param.truth;
    const series = shown.map(mm => {
      const pts = d3.range(nDelta).map(e => {
        const run = coefRun(mm.key, p.param.key, state.iPct, e);
        const y = !run.length ? null
          : p.stat === "bias" ? d3.mean(run) - truth
          : Math.sqrt(d3.mean(run, v => (v - truth) ** 2));
        return { e, x: meta.delta[e], y: Number.isFinite(y) ? y : null };
      });
      const withGaps = [];
      pts.forEach((d, k) => {
        if (k > 0 && d.x - pts[k - 1].x > 1.5 * step) withGaps.push({ e: -1, x: NaN, y: null });
        withGaps.push(d);
      });
      return { mm, pts, v: withGaps };
    });

    const g = svg.append("g").attr("transform", `translate(0,${r * (H + gap)})`);
    const vals = series.flatMap(s => s.pts.map(d => d.y)).filter(v => v != null);
    const ext = d3.extent(vals);
    const ys = d3.scaleLinear(
      p.stat === "bias" ? [Math.min(ext[0], 0), Math.max(ext[1], 0)] : [0, ext[1] || 1],
      [H - m.b, m.t]
    ).nice();

    g.append("rect")
      .attr("x", m.l).attr("y", m.t)
      .attr("width", W - m.r - m.l).attr("height", H - m.b - m.t)
      .attr("fill", PAPER).attr("stroke", RULE);

    if (p.stat === "bias" && ys.domain()[0] < 0) {
      g.append("line")
        .attr("x1", m.l).attr("x2", W - m.r)
        .attr("y1", ys(0)).attr("y2", ys(0))
        .attr("stroke", CLAY).attr("stroke-dasharray", "4 3");
    }

    g.append("g").attr("transform", `translate(0,${H - m.b})`)
      .call(d3.axisBottom(xs).tickValues(meta.delta).tickSize(3));
    g.append("g").attr("transform", `translate(${m.l},0)`)
      .call(d3.axisLeft(ys).ticks(3).tickSize(3).tickFormat(d3.format(".2~g")));

    const line = d3.line().defined(d => d.y != null).x(d => xs(d.x)).y(d => ys(d.y));
    series.forEach(s => {
      g.append("path").datum(s.v)
        .attr("fill", "none").attr("stroke", s.mm.colour)
        .attr("stroke-width", 1.5).attr("d", line);
      g.append("g").selectAll("circle").data(s.pts.filter(d => d.y != null)).join("circle")
        .attr("cx", d => xs(d.x)).attr("cy", d => ys(d.y))
        .attr("r", 2).attr("fill", s.mm.colour);
      const here = s.pts[state.iDelta];
      if (here && here.y != null) {
        g.append("circle")
          .attr("cx", xs(here.x)).attr("cy", ys(here.y))
          .attr("r", 3.4).attr("fill", s.mm.colour);
      }
    });

    g.append("line")
      .attr("x1", xs(meta.delta[state.iDelta])).attr("x2", xs(meta.delta[state.iDelta]))
      .attr("y1", m.t).attr("y2", H - m.b)
      .attr("stroke", MUTED).attr("stroke-dasharray", "2 3").attr("stroke-opacity", 0.7);

    g.append("text").attr("x", m.l).attr("y", 10)
      .attr("font-size", 11).attr("fill", MUTED).text(p.label);
  });

  styleAxes(svg);
}

/* ===========================================================================
   view: weights
   =========================================================================== */

function renderWeights(rec) {
  const gt  = rec.gtRaw[state.gel];
  const pt  = implied(rec, state.gel);
  const lpt = logImplied(rec, state.gel);
  const dev = logDev(rec);
  const marks = [];

  const shown = state.logWeights ? lpt : pt;
  const rule = state.logWeights ? Math.log(1 / nObs) : 1 / nObs;
  const ruleLabel = state.logWeights ? "log(1/N)" : "1/N";
  const yLabel = state.logWeights ? "$\\log \\pi_i$" : "$\\pi_i$";

  marks.push(drawIndexPanel("#gtPlot", rec, gt, {
    rule: 0, ruleLabel: "0", yLabel: "$\\lambda^{\\top} g_i$"
  }));
  marks.push(drawIndexPanel("#ptPlot", rec, shown, { rule, ruleLabel, yLabel }));
  marks.push(drawWeightDev(rec, shown, dev, { rule, yLabel }));
  marks.push(drawWeightLorenz(rec, pt));

  hoverHandlers.push(i => {
    marks.forEach(fn => fn(i));
    weightsReadout(i, rec, gt, pt, lpt, dev);
  });
  weightsReadout(null, rec, gt, pt, lpt, dev);
}

function drawIndexPanel(sel, rec, values, opt) {
  const { rule, ruleLabel, yLabel } = opt;
  const W = 640, H = 470, m = { l: 70, r: 14, t: 12, b: 54 };
  const svg = frame(sel, W, H);
  const xs = d3.scaleLinear([1, nObs], [m.l, W - m.r]);
  const ys = d3.scaleLinear(pad(d3.extent(values), 0.05), [H - m.b, m.t]).nice();

  svg.append("rect")
    .attr("x", m.l).attr("y", m.t)
    .attr("width", W - m.r - m.l).attr("height", H - m.b - m.t)
    .attr("fill", PAPER).attr("stroke", RULE);

  svg.append("g").attr("transform", `translate(0,${H - m.b})`)
    .call(d3.axisBottom(xs).ticks(8).tickSize(4));
  svg.append("g").attr("transform", `translate(${m.l},0)`)
    .call(d3.axisLeft(ys).ticks(6).tickSize(4).tickFormat(d3.format(".3~g")));
  axisLabel(svg, (m.l + W - m.r) / 2, H - 4, "Observation index");
  axisLabel(svg, 16, (m.t + H - m.b) / 2, yLabel, "middle", -90);

  const clip = `clip${sel.replace("#", "")}`;
  svg.append("clipPath").attr("id", clip).append("rect")
    .attr("x", m.l).attr("y", m.t)
    .attr("width", W - m.r - m.l).attr("height", H - m.b - m.t);
  const plot = svg.append("g").attr("clip-path", `url(#${clip})`);

  if (rule != null && rule >= ys.domain()[0] && rule <= ys.domain()[1]) {
    plot.append("line")
      .attr("x1", m.l).attr("x2", W - m.r)
      .attr("y1", ys(rule)).attr("y2", ys(rule))
      .attr("stroke", MUTED).attr("stroke-dasharray", "4 3").attr("stroke-width", 2);
    plot.append("text")
      .attr("x", W - m.r - 5).attr("y", ys(rule) - 5)
      .attr("text-anchor", "end").attr("font-size", 10).attr("fill", MUTED)
      .attr("font-family", "ui-monospace, Menlo, monospace")
      .text(ruleLabel);
  }

  const idx = d3.range(nObs);
  plot.append("g").selectAll("circle")
    .data(idx.filter(i => !rec.outlier[i])).join("circle")
    .attr("cx", i => xs(i + 1)).attr("cy", i => ys(values[i]))
    .attr("r", 2.6).attr("fill", CALM).attr("fill-opacity", 0.82);
  plot.append("g").selectAll("circle")
    .data(idx.filter(i => rec.outlier[i])).join("circle")
    .attr("cx", i => xs(i + 1)).attr("cy", i => ys(values[i]))
    .attr("r", 4).attr("fill", CLAY).attr("fill-opacity", 1);

  const ring = hoverRing(plot);
  const finder = d3.Delaunay.from(idx, i => xs(i + 1), i => ys(values[i]));
  svg.append("rect")
    .attr("x", m.l).attr("y", m.t)
    .attr("width", W - m.r - m.l).attr("height", H - m.b - m.t)
    .attr("fill", "transparent")
    .on("pointermove", ev => {
      const [px, py] = d3.pointer(ev);
      setHover(finder.find(px, py));
    })
    .on("pointerleave", () => setHover(null));

  styleAxes(svg);
  return i => ring(i, i == null ? null : xs(i + 1), i == null ? null : ys(values[i]));
}

function drawWeightDev(rec, values, dev, opt) {
  const W = 640, H = 470, m = { l: 70, r: 14, t: 12, b: 48 };
  const svg = frame("#ptDev", W, H);
  const xs = d3.scaleLinear(pad(d3.extent(dev), 0.04), [m.l, W - m.r]).nice();
  const ys = d3.scaleLinear(pad(d3.extent(values), 0.05), [H - m.b, m.t]).nice();

  svg.append("rect")
    .attr("x", m.l).attr("y", m.t)
    .attr("width", W - m.r - m.l).attr("height", H - m.b - m.t)
    .attr("fill", PAPER).attr("stroke", RULE);

  svg.append("g").attr("transform", `translate(0,${H - m.b})`)
    .call(d3.axisBottom(xs).ticks(6).tickSize(4));
  svg.append("g").attr("transform", `translate(${m.l},0)`)
    .call(d3.axisLeft(ys).ticks(5).tickSize(4).tickFormat(d3.format(".3~g")));
  axisLabel(svg, (m.l + W - m.r) / 2, H - 4, "$\\log y_i - \\mu$");
  axisLabel(svg, 16, (m.t + H - m.b) / 2, opt.yLabel, "middle", -90);

  const clip = "clip-ptdev";
  svg.append("clipPath").attr("id", clip).append("rect")
    .attr("x", m.l).attr("y", m.t)
    .attr("width", W - m.r - m.l).attr("height", H - m.b - m.t);
  const plot = svg.append("g").attr("clip-path", `url(#${clip})`);

  if (opt.rule >= ys.domain()[0] && opt.rule <= ys.domain()[1]) {
    plot.append("line")
      .attr("x1", m.l).attr("x2", W - m.r)
      .attr("y1", ys(opt.rule)).attr("y2", ys(opt.rule))
      .attr("stroke", MUTED).attr("stroke-dasharray", "4 3");
  }
  if (0 > xs.domain()[0] && 0 < xs.domain()[1]) {
    plot.append("line")
      .attr("x1", xs(0)).attr("x2", xs(0))
      .attr("y1", m.t).attr("y2", H - m.b)
      .attr("stroke", MUTED).attr("stroke-dasharray", "2 3").attr("stroke-opacity", 0.6);
  }

  const idx = d3.range(nObs);
  plot.append("g").selectAll("circle")
    .data(idx.filter(i => !rec.outlier[i])).join("circle")
    .attr("cx", i => xs(dev[i])).attr("cy", i => ys(values[i]))
    .attr("r", 2.2).attr("fill", CALM).attr("fill-opacity", 0.82);
  plot.append("g").selectAll("circle")
    .data(idx.filter(i => rec.outlier[i])).join("circle")
    .attr("cx", i => xs(dev[i])).attr("cy", i => ys(values[i]))
    .attr("r", 3.6).attr("fill", CLAY).attr("fill-opacity", 1);

  const ring = hoverRing(plot);
  const finder = d3.Delaunay.from(idx, i => xs(dev[i]), i => ys(values[i]));
  svg.append("rect")
    .attr("x", m.l).attr("y", m.t)
    .attr("width", W - m.r - m.l).attr("height", H - m.b - m.t)
    .attr("fill", "transparent")
    .on("pointermove", ev => {
      const [px, py] = d3.pointer(ev);
      setHover(finder.find(px, py));
    })
    .on("pointerleave", () => setHover(null));

  styleAxes(svg);
  return i => ring(i, i == null ? null : xs(dev[i]), i == null ? null : ys(values[i]));
}

/* Lorenz curve for the implied probabilities: observations sorted by weight,
   largest first, against cumulative share of the total. Uniform weighting is
   the diagonal, so the gap between curve and diagonal is the departure from
   equal weighting. */
function drawWeightLorenz(rec, pt) {
  const W = 640, H = 470, m = { l: 70, r: 14, t: 12, b: 48 };
  const svg = frame("#ptLorenz", W, H);

  const order = d3.range(nObs).sort((a, b) => pt[b] - pt[a]);
  const total = d3.sum(pt) || 1;
  const cum = new Float64Array(nObs);
  let run = 0;
  order.forEach((i, k) => { run += pt[i] / total; cum[k] = run; });

  const rank = new Int32Array(nObs);
  order.forEach((i, k) => { rank[i] = k; });

  const xs = d3.scaleLinear([0, 1], [m.l, W - m.r]);
  const ys = d3.scaleLinear([0, 1], [H - m.b, m.t]);
  const frac = k => (k + 1) / nObs;

  svg.append("rect")
    .attr("x", m.l).attr("y", m.t)
    .attr("width", W - m.r - m.l).attr("height", H - m.b - m.t)
    .attr("fill", PAPER).attr("stroke", RULE);

  svg.append("g").attr("transform", `translate(0,${H - m.b})`)
    .call(d3.axisBottom(xs).ticks(6).tickSize(4).tickFormat(d3.format(".0%")));
  svg.append("g").attr("transform", `translate(${m.l},0)`)
    .call(d3.axisLeft(ys).ticks(6).tickSize(4).tickFormat(d3.format(".0%")));
  axisLabel(svg, (m.l + W - m.r) / 2, H - 6, "share of the sample, heaviest first");
  axisLabel(svg, 16, (m.t + H - m.b) / 2, "share of total weight", "middle", -90);

  const clip = "clip-ptlorenz";
  svg.append("clipPath").attr("id", clip).append("rect")
    .attr("x", m.l).attr("y", m.t)
    .attr("width", W - m.r - m.l).attr("height", H - m.b - m.t);
  const plot = svg.append("g").attr("clip-path", `url(#${clip})`);

  plot.append("line")
    .attr("x1", xs(0)).attr("y1", ys(0))
    .attr("x2", xs(1)).attr("y2", ys(1))
    .attr("stroke", MUTED).attr("stroke-dasharray", "4 3").attr("stroke-width", 1.6);

  const line = d3.line().x((_, k) => xs(frac(k))).y((_, k) => ys(cum[k]));
  plot.append("path").datum(d3.range(nObs))
    .attr("d", line)
    .attr("fill", "none")
    .attr("stroke", INK)
    .attr("stroke-width", 1.8);

  const nOut = d3.sum(rec.outlier);
  if (nOut > 0 && nOut < nObs) {
    const heldByOut = order.slice(0, nOut).reduce((s, i) => s + pt[i] / total, 0);
    plot.append("line")
      .attr("x1", xs(nOut / nObs)).attr("x2", xs(nOut / nObs))
      .attr("y1", ys(0)).attr("y2", ys(1))
      .attr("stroke", CLAY).attr("stroke-opacity", 0.45)
      .attr("stroke-dasharray", "3 3");
    plot.append("text")
      .attr("x", xs(nOut / nObs) + 6).attr("y", ys(0.04))
      .attr("font-size", 10).attr("fill", CLAY)
      .attr("paint-order", "stroke").attr("stroke", PAPER).attr("stroke-width", 3)
      .text(`top ${nOut} hold ${d3.format(".1%")(heldByOut)}`);
  }

  const ring = hoverRing(plot);
  const finder = d3.Delaunay.from(d3.range(nObs), k => xs(frac(k)), k => ys(cum[k]));
  svg.append("rect")
    .attr("x", m.l).attr("y", m.t)
    .attr("width", W - m.r - m.l).attr("height", H - m.b - m.t)
    .attr("fill", "transparent")
    .on("pointermove", ev => {
      const [px, py] = d3.pointer(ev);
      setHover(order[finder.find(px, py)]);
    })
    .on("pointerleave", () => setHover(null));

  styleAxes(svg);
  return i => {
    if (i == null) { ring(null, null, null); return; }
    const k = rank[i];
    ring(i, xs(frac(k)), ys(cum[k]));
  };
}

function weightsReadout(i, rec, gt, pt, lpt, dev) {
  const el = d3.select("#weightsReadout");
  const mm = METHODS.find(x => x.key === state.gel);
  const s = slot(state.iPct, state.iDelta, state.iSim);

  if (i == null) {
    let outW = 0, nOut = 0;
    for (let k = 0; k < nObs; k++) if (rec.outlier[k]) { outW += pt[k]; nOut++; }
    const share = nOut ? outW / (nOut / nObs) : 0;   // 1 means "no downweighting"
    el.html(
      fitSummary(mm, fitAt(state.gel, s)) +
      `<span class="sep">|</span>` +
      `<span class="dim">outliers hold</span> <b style="color:${CLAY}">${f2(100 * outW)}%</b>` +
      ` <span class="dim">of the weight, against</span> ${f1(100 * nOut / nObs)}%` +
      ` <span class="dim">of the sample &mdash; a factor of</span> <b>${f2(share)}</b>`
    );
    return;
  }
  el.html(
    `<span class="dim">obs</span> ${i + 1}` +
    `<span class="sep">|</span><span class="dim">λ⊤gᵢ</span> ${f3(gt[i])}` +
    `<span class="sep">|</span><span class="dim">πᵢ</span> ${fe(pt[i])}` +
    ` <span class="dim">(${d3.format(".3~g")(pt[i] * nObs)}× uniform)</span>` +
    `<span class="sep">|</span><span class="dim">log πᵢ</span> ${f2(lpt[i])}` +
    `<span class="sep">|</span><span class="dim">log y − μ</span> ${f3(dev[i])}` +
    (rec.outlier[i] ? `<span class="sep">|</span><b style="color:${CLAY}">outlier</b>` : "")
  );
}

/* ===========================================================================
   view: likelihood surface
   =========================================================================== */

/* d3.contours treats grid value k as sitting at the centre of cell k, so
   grid coordinate g maps to index g - 0.5. */
const HALF = 0.5;
const gxToMu = gx => meta.mu[0] + (gx - HALF) * dMu;
const gyToS2 = gy => meta.sigma2[0] + (gy - HALF) * dS2;

function renderLikelihood(rec) {
  // The export is row-major over mu; d3.contours wants the first axis varying
  // fastest, so transpose once here and mu becomes the horizontal axis.
  const z = new Float64Array(nLL);
  const hull = new Uint8Array(nLL);
  for (let i = 0; i < nMu; i++)
    for (let j = 0; j < nS2; j++) {
      z[j * nMu + i] = rec.ll[i * nS2 + j];
      hull[j * nMu + i] = rec.hull[i * nS2 + j];
    }

  let iMax = -1;
  for (let k = 0; k < nLL; k++) {
    if (!Number.isFinite(z[k])) continue;
    if (iMax < 0 || z[k] > z[iMax]) iMax = k;
  }
  if (iMax < 0) {
    d3.select("#llSurface").html(`<p class="caption">No feasible grid cell for this replicate.</p>`);
    d3.select("#llProfiles").html("");
    d3.select("#llReadout").html("");
    return;
  }
  const maxLL = z[iMax];
  const maxI = iMax % nMu, maxJ = (iMax - maxI) / nMu;

  const floor = maxLL - state.llWindow;
  const zc = Float64Array.from(z, v => (Number.isFinite(v) ? Math.max(v, floor) : floor));

  const info = { maxLL, floor, maxI, maxJ, hull };
  drawLLSurface(zc, z, info);
  drawLLProfiles(z, info);

  const tt = meta.theta_true;
  const llTrue = nodeValue(z, hull, tt.mu, tt.sigma2);
  const s = slot(state.iPct, state.iDelta, state.iSim);
  const bits = [
    `<span class="dim">grid max</span> <b>${f1(maxLL)}</b>` +
    ` <span class="dim">at</span> (${fg(meta.mu[maxI])}, ${fg(meta.sigma2[maxJ])})`,
    `<span class="dim">at the truth</span> ` +
    `<b style="color:${CLAY}">${llTrue == null ? "infeasible" : f1(llTrue)}</b>` +
    (llTrue == null ? "" : ` <span class="dim">gap</span> ${f1(maxLL - llTrue)}`)
  ];
  activeMethods().forEach(mm => bits.push(fitSummary(mm, fitAt(mm.key, s))));
  d3.select("#llReadout").html(bits.join(`<span class="sep">|</span>`));
}

/* Value at the grid node nearest (mu, s2), or null when that node is
   infeasible or off the grid. The truth sits exactly on a node. */
function nodeValue(z, hull, t1, t2) {
  const i = Math.round((t1 - meta.mu[0]) / dMu), j = Math.round((t2 - meta.sigma2[0]) / dS2);
  if (i < 0 || j < 0 || i >= nMu || j >= nS2) return null;
  const k = j * nMu + i;
  return hull[k] && Number.isFinite(z[k]) ? z[k] : null;
}

/* Bilinear read of the surface at an arbitrary (mu, s2), for the hover
   readout; the grid is evenly spaced in both. */
function bilinear(z, t1, t2) {
  const u = (t1 - meta.mu[0]) / dMu, v = (t2 - meta.sigma2[0]) / dS2;
  if (u < 0 || v < 0 || u > nMu - 1 || v > nS2 - 1) return null;
  const i = Math.min(Math.floor(u), nMu - 2), j = Math.min(Math.floor(v), nS2 - 2);
  const a = u - i, b = v - j;
  const at = (p, q) => z[q * nMu + p];
  const out = at(i, j) * (1 - a) * (1 - b) + at(i + 1, j) * a * (1 - b) +
              at(i, j + 1) * (1 - a) * b + at(i + 1, j + 1) * a * b;
  return Number.isFinite(out) ? out : null;
}

/* Bounding box of every contour ring, in data units. */
function contourBBox(contours) {
  let x0 = Infinity, x1 = -Infinity, y0 = Infinity, y1 = -Infinity;
  contours.forEach(c => c.coordinates.forEach(poly => poly.forEach(ring =>
    ring.forEach(([gx, gy]) => {
      const x = gxToMu(gx), y = gyToS2(gy);
      if (x < x0) x0 = x; if (x > x1) x1 = x;
      if (y < y0) y0 = y; if (y > y1) y1 = y;
    })
  )));
  return isFinite(x0) ? [x0, x1, y0, y1] : null;
}

function drawLLSurface(zc, z, info) {
  const W = 640, H = 470, m = { l: 66, r: 16, t: 12, b: 44 };
  const svg = frame("#llSurface", W, H);
  const tt = meta.theta_true;
  const s = slot(state.iPct, state.iDelta, state.iSim);

  const levels = d3.ticks(info.floor, info.maxLL, 11).filter(v => v > info.floor);
  const contours = d3.contours().size([nMu, nS2]).thresholds(levels)(zc);

  // Zoom to where there is something to see: the contours, the truth, the
  // grid maximum and the fits, with the grid edges as a hard limit.
  let lo1 = Infinity, hi1 = -Infinity, lo2 = Infinity, hi2 = -Infinity;
  const inc = (a, b) => {
    if (!Number.isFinite(a) || !Number.isFinite(b)) return;
    lo1 = Math.min(lo1, a); hi1 = Math.max(hi1, a);
    lo2 = Math.min(lo2, b); hi2 = Math.max(hi2, b);
  };
  const bbox = contourBBox(contours);
  if (bbox) { inc(bbox[0], bbox[2]); inc(bbox[1], bbox[3]); }
  inc(tt.mu, tt.sigma2);
  inc(meta.mu[info.maxI], meta.sigma2[info.maxJ]);
  const fits = activeMethods().map(mm => ({ mm, f: fitAt(mm.key, s) })).filter(d => d.f);
  fits.forEach(({ mm, f }) => {
    inc(f.mu, f.s2);
    if (mm.key === "gel_bias") inc(coefs.gel_bias.start_mu[s], coefs.gel_bias.start_sigma2[s]);
  });
  const gridLo1 = meta.mu[0] - dMu / 2, gridHi1 = meta.mu[nMu - 1] + dMu / 2;
  const gridLo2 = meta.sigma2[0] - dS2 / 2, gridHi2 = meta.sigma2[nS2 - 1] + dS2 / 2;
  let dom1 = [gridLo1, gridHi1], dom2 = [gridLo2, gridHi2];
  if (isFinite(lo1)) {
    const p1 = Math.max((hi1 - lo1) * 0.12, 3 * dMu);
    const p2 = Math.max((hi2 - lo2) * 0.12, 3 * dS2);
    dom1 = [Math.max(gridLo1, lo1 - p1), Math.min(gridHi1, hi1 + p1)];
    dom2 = [Math.max(gridLo2, lo2 - p2), Math.min(gridHi2, hi2 + p2)];
  }
  const xs = d3.scaleLinear(dom1, [m.l, W - m.r]);
  const ys = d3.scaleLinear(dom2, [H - m.b, m.t]);

  // Contour coordinates arrive in grid units.
  const path = d3.geoPath(d3.geoTransform({
    point(gx, gy) { this.stream.point(xs(gxToMu(gx)), ys(gyToS2(gy))); }
  }));

  const fill = d3.scaleLinear()
    .domain([info.floor, info.floor + (info.maxLL - info.floor) * 0.55, info.maxLL])
    .range(["#eef4f8", "#cfdce4", "#6f97ac"])
    .interpolate(d3.interpolateLab)
    .clamp(true);
  svg.append("rect")
    .attr("x", m.l).attr("y", m.t)
    .attr("width", W - m.r - m.l).attr("height", H - m.b - m.t)
    .attr("fill", PAPER).attr("stroke", RULE);

  const clip = "clip-ll";
  svg.append("clipPath").attr("id", clip).append("rect")
    .attr("x", m.l).attr("y", m.t)
    .attr("width", W - m.r - m.l).attr("height", H - m.b - m.t);
  const plot = svg.append("g").attr("clip-path", `url(#${clip})`);

  // Infeasible cells first, so the contours sit over them.
  const cells = [];
  for (let j = 0; j < nS2; j++) for (let i = 0; i < nMu; i++) {
    if (info.hull[j * nMu + i]) continue;
    const mu = meta.mu[i], s2 = meta.sigma2[j];
    if (mu + dMu / 2 < dom1[0] || mu - dMu / 2 > dom1[1]) continue;
    if (s2 + dS2 / 2 < dom2[0] || s2 - dS2 / 2 > dom2[1]) continue;
    cells.push([mu, s2]);
  }
  plot.append("g").selectAll("rect").data(cells).join("rect")
    .attr("x", d => xs(d[0] - dMu / 2))
    .attr("y", d => ys(d[1] + dS2 / 2))
    .attr("width", xs(dMu) - xs(0) + 0.5)
    .attr("height", ys(0) - ys(dS2) + 0.5)
    .attr("fill", SHADE);

  plot.append("g").selectAll("path").data(contours).join("path")
    .attr("d", path)
    .attr("fill", c => fill(c.value))
    .attr("stroke", INK)
    .attr("stroke-opacity", 0.4)
    .attr("stroke-width", 0.95);

  const anchorX = xs(meta.mu[info.maxI]), anchorY = ys(meta.sigma2[info.maxJ]);
  labelContours(plot, contours, xs, ys, anchorX, anchorY);

  const mark = (cx, cy, colour, text, dy, hollow) => {
    plot.append("circle").attr("cx", cx).attr("cy", cy)
      .attr("r", 4.6)
      .attr("fill", hollow ? PAPER : colour)
      .attr("stroke", hollow ? colour : PAPER).attr("stroke-width", hollow ? 1.6 : 1.2);
    if (text) plot.append("text")
      .attr("x", Math.max(m.l + 26, Math.min(W - m.r - 26, cx)))
      .attr("y", cy + dy)
      .attr("text-anchor", "middle").attr("font-size", 11).attr("fill", colour)
      .attr("paint-order", "stroke").attr("stroke", PAPER).attr("stroke-width", 3.5)
      .text(text);
  };

  // Biased start and the path the optimiser took, drawn under the marks.
  fits.forEach(({ mm, f }) => {
    if (mm.key !== "gel_bias") return;
    const sm = coefs.gel_bias.start_mu[s], ss = coefs.gel_bias.start_sigma2[s];
    if (!Number.isFinite(sm) || !Number.isFinite(ss)) return;
    plot.append("line")
      .attr("x1", xs(sm)).attr("y1", ys(ss)).attr("x2", xs(f.mu)).attr("y2", ys(f.s2))
      .attr("stroke", mm.colour).attr("stroke-width", 1).attr("stroke-dasharray", "3 2");
    mark(xs(sm), ys(ss), mm.colour, "start", -10, true);
  });
  // The estimates carry no text: their colours match the legend below, and
  // when a fit lands on the truth a label would only cover the clay mark.
  fits.forEach(({ mm, f }) => mark(xs(f.mu), ys(f.s2), mm.colour, null, 0, false));
  mark(anchorX, anchorY, INK, "grid max", -11, false);
  mark(xs(tt.mu), ys(tt.sigma2), CLAY, "true (μ, σ²)", 18, false);

  svg.append("g").attr("transform", `translate(0,${H - m.b})`)
    .call(d3.axisBottom(xs).ticks(8).tickSize(4));
  svg.append("g").attr("transform", `translate(${m.l},0)`)
    .call(d3.axisLeft(ys).ticks(8).tickSize(4));
  axisLabel(svg, (m.l + W - m.r) / 2, H - 8, "$\\mu$");
  axisLabel(svg, 15, (m.t + H - m.b) / 2, "$\\sigma^2$", "middle", -90);

  const cross = plot.append("g").style("display", "none");
  const cx = cross.append("circle").attr("r", 3).attr("fill", "none")
    .attr("stroke", INK).attr("stroke-width", 1.2);
  svg.append("rect")
    .attr("x", m.l).attr("y", m.t)
    .attr("width", W - m.r - m.l).attr("height", H - m.b - m.t)
    .attr("fill", "transparent")
    .on("pointermove", ev => {
      const [px, py] = d3.pointer(ev);
      const t1 = xs.invert(px), t2 = ys.invert(py);
      const feasible = nodeValue(z, info.hull, t1, t2) != null;
      const v = feasible ? bilinear(z, t1, t2) : null;
      cross.style("display", null);
      cx.attr("cx", px).attr("cy", py);
      d3.select("#llReadout").select(".probe").remove();
      d3.select("#llReadout").append("span").attr("class", "probe")
        .html(`<span class="sep">|</span><span class="dim">at</span> ` +
              `(${f2(t1)}, ${f2(t2)}) <span class="dim">ll</span> ` +
              `${!feasible ? "infeasible" : v == null ? "—" : f1(v)}`);
    })
    .on("pointerleave", () => {
      cross.style("display", "none");
      d3.select("#llReadout").select(".probe").remove();
    });

  styleAxes(svg);
}

/* Labelled contours: take the first crossing of the horizontal line through
   the maximum, which puts the labels in a readable row because nested
   contours cross that line at increasing distance. */
function labelContours(g, contours, xs, ys, x0, y0) {
  const layer = g.append("g");
  const px = ([gx, gy]) => [xs(gxToMu(gx)), ys(gyToS2(gy))];
  const lo = xs.range()[0], hi = xs.range()[1];
  const dir = (hi - x0) >= (x0 - lo) ? 1 : -1;

  const placed = [];
  contours.forEach(c => {
    let best = null;
    c.coordinates.forEach(poly => poly.forEach(ring => {
      for (let i = 1; i < ring.length; i++) {
        const a = px(ring[i - 1]), b = px(ring[i]);
        if ((a[1] - y0) * (b[1] - y0) > 0) continue;      // no crossing
        const dy = b[1] - a[1];
        const x = dy === 0 ? a[0] : a[0] + ((y0 - a[1]) / dy) * (b[0] - a[0]);
        if (dir * (x - x0) < 8) continue;                 // on the wrong side
        if (best == null || dir * (x - best) < 0) best = x;
      }
    }));
    if (best == null || best < lo + 14 || best > hi - 14) return;
    if (placed.some(x => Math.abs(x - best) < 30)) return;
    placed.push(best);

    layer.append("text")
      .attr("x", best).attr("y", y0)
      .attr("text-anchor", "middle").attr("dy", 3.5)
      .attr("font-size", 9.5)
      .attr("font-family", "ui-monospace, Menlo, monospace")
      .attr("fill", INK)
      .attr("paint-order", "stroke")
      .attr("stroke", "#f5f1e7").attr("stroke-width", 3.5)
      .text(d3.format(".0f")(c.value));
  });
}

function drawLLProfiles(z, info) {
  const W = 420, H = 190, gap = 26, m = { l: 62, r: 12, t: 18, b: 30 };
  const svg = frame("#llProfiles", W, 2 * H + gap);
  const tt = meta.theta_true;

  const panels = [
    {
      label: `μ at σ² = ${fg(meta.sigma2[info.maxJ])}`,
      axis: meta.mu,
      raw: d3.range(nMu).map(i => z[info.maxJ * nMu + i]),
      truth: tt.mu,
      at: info.maxI
    },
    {
      label: `σ² at μ = ${fg(meta.mu[info.maxI])}`,
      axis: meta.sigma2,
      raw: d3.range(nS2).map(j => z[j * nMu + info.maxI]),
      truth: tt.sigma2,
      at: info.maxJ
    }
  ];

  panels.forEach((p, r) => {
    const g = svg.append("g").attr("transform", `translate(0,${r * (H + gap)})`);
    const xs = d3.scaleLinear(d3.extent(p.axis), [m.l, W - m.r]);
    // Clipped at the same window as the contours, so the two panels agree
    // about what counts as "near the peak".
    const ys = d3.scaleLinear([info.floor, info.maxLL], [H - m.b, m.t]).nice();

    g.append("rect")
      .attr("x", m.l).attr("y", m.t)
      .attr("width", W - m.r - m.l).attr("height", H - m.b - m.t)
      .attr("fill", PAPER).attr("stroke", RULE);

    const clip = `clip-prof-${r}`;
    g.append("clipPath").attr("id", clip).append("rect")
      .attr("x", m.l).attr("y", m.t)
      .attr("width", W - m.r - m.l).attr("height", H - m.b - m.t);
    const plot = g.append("g").attr("clip-path", `url(#${clip})`);

    // Infeasible stretches of the profile are shaded, and the line stops there.
    const step = p.axis[1] - p.axis[0];
    plot.append("g").selectAll("rect")
      .data(p.raw.map((v, i) => ({ v, i })).filter(d => !Number.isFinite(d.v))).join("rect")
      .attr("x", d => xs(p.axis[d.i] - step / 2))
      .attr("width", xs(step) - xs(0) + 0.5)
      .attr("y", m.t).attr("height", H - m.b - m.t)
      .attr("fill", SHADE);

    plot.append("line")
      .attr("x1", xs(p.truth)).attr("x2", xs(p.truth))
      .attr("y1", m.t).attr("y2", H - m.b)
      .attr("stroke", CLAY).attr("stroke-dasharray", "4 3");

    plot.append("path")
      .datum(p.raw.map((v, i) => ({ x: p.axis[i], v })))
      .attr("fill", "none").attr("stroke", INK).attr("stroke-width", 1.5)
      .attr("d", d3.line()
        .defined(d => Number.isFinite(d.v))
        .x(d => xs(d.x)).y(d => ys(Math.max(d.v, info.floor))));

    plot.append("circle")
      .attr("cx", xs(p.axis[p.at])).attr("cy", ys(info.maxLL))
      .attr("r", 3.4).attr("fill", INK);

    g.append("g").attr("transform", `translate(0,${H - m.b})`)
      .call(d3.axisBottom(xs).ticks(6).tickSize(3));
    g.append("g").attr("transform", `translate(${m.l},0)`)
      .call(d3.axisLeft(ys).ticks(4).tickSize(3).tickFormat(d3.format(".0f")));
    g.append("text").attr("x", m.l).attr("y", 11)
      .attr("font-size", 11).attr("fill", MUTED).text(p.label);
  });

  styleAxes(svg);
}

/* ===========================================================================
   controls
   =========================================================================== */

function tabGroup(sel, items, isOn, onPick) {
  const host = d3.select(sel);
  host.selectAll("*").remove();
  host.selectAll("button").data(items).join("button")
    .attr("class", d => "tab" + (isOn(d) ? " on" : ""))
    .attr("type", "button")
    .text(d => d.label)
    .on("click", (_, d) => onPick(d));
}

function methodLegend(sel) {
  const host = d3.select(sel);
  host.selectAll("*").remove();
  const g = host.append("span").attr("class", "legend");
  g.selectAll("button").data(METHODS).join("button")
    .attr("class", d => "item" + (state.methods.has(d.key) ? " on" : ""))
    .attr("type", "button")
    .style("color", d => state.methods.has(d.key) ? d.colour : null)
    .html(d => `<span class="swatch"></span>${d.label}`)
    .on("click", (_, d) => {
      // Never let the last method be turned off: an empty panel is a worse
      // answer to "show me nothing" than simply refusing the click.
      if (state.methods.has(d.key)) {
        if (state.methods.size > 1) state.methods.delete(d.key);
      } else {
        state.methods.add(d.key);
      }
      render();
    });
}

function checkbox(host, label, get, set) {
  const l = host.append("label").attr("class", "ctl");
  const input = l.append("input").attr("type", "checkbox").property("checked", get());
  l.append("span").text(label);
  input.on("change", function () { set(this.checked); render(); });
}

function buildKnobs() {
  tabGroup("#views", VIEWS, v => v.key === state.view, v => {
    state.view = v.key;
    render();
  });

  tabGroup("#pctTabs",
    meta.pct.map((p, i) => ({ label: `${+p}%`, i })),
    d => d.i === state.iPct,
    d => { state.iPct = d.i; render(); });

  d3.select("#deltaRange")
    .attr("max", nDelta - 1).property("value", state.iDelta)
    .on("input", function () { state.iDelta = +this.value; render(); });
  d3.select("#simRange")
    .attr("max", nSim - 1).property("value", state.iSim)
    .on("input", function () { state.iSim = +this.value; render(); });

  /* Two ends over the delta grid. Each pushes the other rather than crossing
     it, so the range is always valid without a second guard at draw time. */
  const setEnd = (which, v) => {
    if (which === "lo") {
      state.deltaLo = v;
      if (state.deltaHi < v) state.deltaHi = v;
    } else {
      state.deltaHi = v;
      if (state.deltaLo > v) state.deltaLo = v;
    }
    render();
  };
  d3.select("#deltaLoRange").attr("max", nDelta - 1).property("value", state.deltaLo)
    .on("input", function () { setEnd("lo", +this.value); });
  d3.select("#deltaHiRange").attr("max", nDelta - 1).property("value", state.deltaHi)
    .on("input", function () { setEnd("hi", +this.value); });

  // Per-view knobs, rebuilt on each render so their state always matches.
  const sk = d3.select("#sampleKnobs"); sk.selectAll("*").remove();
  sk.append("span").attr("class", "ctl").append("span").attr("class", "lbl").text("fits shown");
  sk.append("span").attr("id", "sampleLegend");
  methodLegend("#sampleLegend");
  checkbox(sk, "log axis for y", () => state.logY, v => state.logY = v);

  const ek = d3.select("#estKnobs"); ek.selectAll("*").remove();
  ek.append("span").attr("class", "ctl").append("span").attr("class", "lbl").text("methods");
  ek.append("span").attr("id", "estLegend");
  methodLegend("#estLegend");
  checkbox(ek, "trim the outer 1% of estimates", () => state.trim, v => state.trim = v);

  const wk = d3.select("#weightsKnobs"); wk.selectAll("*").remove();
  const wc = wk.append("span").attr("class", "ctl");
  wc.append("span").attr("class", "lbl").text("fit");
  wc.append("span").attr("id", "gelTabs").attr("class", "tabs");
  tabGroup("#gelTabs", METHODS.map(m => ({ label: m.label, key: m.key })),
    d => d.key === state.gel,
    d => { state.gel = d.key; render(); });
  checkbox(wk, "log scale on the weights (log πᵢ)", () => state.logWeights, v => state.logWeights = v);

  const lk = d3.select("#llKnobs"); lk.selectAll("*").remove();
  const lc = lk.append("label").attr("class", "ctl");
  lc.append("span").attr("class", "lbl").text("contour window");
  // Same range as the Shiny app's ll_window slider: 50 to 5000 by 50.
  lc.append("input").attr("type", "range").attr("min", 50).attr("max", 5000)
    .attr("step", 50).property("value", state.llWindow)
    .on("input", function () { state.llWindow = +this.value; render(); });
  lc.append("span").attr("class", "val").attr("id", "llWindowVal").text(state.llWindow);
  lk.append("span").attr("class", "ctl").append("span").attr("class", "lbl").text("fits shown");
  lk.append("span").attr("id", "llLegend");
  methodLegend("#llLegend");
}

function syncKnobs() {
  d3.selectAll("#views .tab").classed("on", d => d.key === state.view);
  d3.selectAll("#pctTabs .tab").classed("on", d => d.i === state.iPct);
  d3.select("#deltaRange").property("value", state.iDelta);
  d3.select("#simRange").property("value", state.iSim);
  const dl = meta.delta[state.iDelta];
  d3.select("#deltaVal").text(dl > 0 ? `+${fg(dl)}` : fg(dl));
  d3.select("#simVal").text(`${meta.sim[state.iSim]} / ${nSim}`);

  /* The density view pools replicates and spans a range of delta, so it gets
     the range ends and neither single-value slider. */
  const onEst = state.view === "estimates";
  d3.select("#simCtl").style("display", onEst ? "none" : "");
  d3.select("#deltaCtl").style("display", onEst ? "none" : "");
  d3.select("#deltaRangeCtl").style("display", onEst ? "" : "none");

  d3.select("#deltaLoRange").property("value", state.deltaLo);
  d3.select("#deltaHiRange").property("value", state.deltaHi);

  const pc = i => (nDelta > 1 ? (i / (nDelta - 1)) * 100 : 0);
  d3.select("#deltaFill")
    .style("left", `${pc(state.deltaLo)}%`)
    .style("width", `${pc(state.deltaHi) - pc(state.deltaLo)}%`);
  d3.select("#deltaLoRange").style("z-index", state.deltaLo > (nDelta - 1) / 2 ? 2 : 1);
  d3.select("#deltaHiRange").style("z-index", state.deltaLo > (nDelta - 1) / 2 ? 1 : 2);

  d3.select("#deltaRangeVal")
    .text(`${fg(meta.delta[state.deltaLo])} … ${fg(meta.delta[state.deltaHi])}`);
  d3.select("#llWindowVal").text(state.llWindow);
  d3.selectAll(".legend .item")
    .classed("on", d => state.methods.has(d.key))
    .style("color", d => state.methods.has(d.key) ? d.colour : null);
  d3.selectAll("#gelTabs .tab").classed("on", d => d.key === state.gel);
  d3.selectAll("main .view").classed("on", function () {
    return this.dataset.view === state.view;
  });
}

/* ---- deep links -----------------------------------------------------------
   The three global knobs plus the view go in the hash, so a particular
   replicate can be sent to someone: #weights/05/-3/17 */
function writeHash() {
  const h = `#${state.view}/${meta.pct[state.iPct]}/${meta.delta_slug[state.iDelta]}/${meta.sim[state.iSim]}`;
  if (location.hash !== h) history.replaceState(null, "", h);
}

function readHash() {
  const parts = location.hash.replace(/^#/, "").split("/");
  if (parts.length < 4) return;
  const [v, p, e, s] = parts;
  if (VIEWS.some(x => x.key === v)) state.view = v;
  const ip = meta.pct.indexOf(p);          if (ip >= 0) state.iPct = ip;
  const ie = meta.delta_slug.indexOf(e);   if (ie >= 0) state.iDelta = ie;
  const is = meta.sim.indexOf(+s);         if (is >= 0) state.iSim = is;
}

/* ===========================================================================
   render
   =========================================================================== */

let renderToken = 0;

async function render() {
  const token = ++renderToken;
  hoverHandlers = [];
  state.hover = null;
  syncKnobs();
  writeHash();

  if (state.view === "estimates") {
    renderEstimates();
    prefetch();
    return;
  }

  let rec;
  try {
    rec = await loadRep(state.iPct, state.iDelta, state.iSim);
  } catch (err) {
    if (token !== renderToken) return;
    d3.select(`main .view[data-view="${state.view}"] .figure`)
      .html(`<p class="caption">Could not load ${repUrl(state.iPct, state.iDelta, state.iSim)}.</p>`);
    return;
  }
  if (token !== renderToken) return;   // a knob moved while the fetch was in flight

  if (state.view === "sample") renderSample(rec);
  else if (state.view === "weights") renderWeights(rec);
  else if (state.view === "likelihood") renderLikelihood(rec);

  prefetch();
}

/* Warm the neighbours of the current cell, so dragging either slider one step
   is already resolved by the time the pointer stops. */
function prefetch() {
  const near = [
    [state.iPct, state.iDelta, state.iSim + 1],
    [state.iPct, state.iDelta, state.iSim - 1],
    [state.iPct, state.iDelta + 1, state.iSim],
    [state.iPct, state.iDelta - 1, state.iSim]
  ];
  for (const [p, e, s] of near) {
    if (e < 0 || e >= nDelta || s < 0 || s >= nSim) continue;
    loadRep(p, e, s).catch(() => {});
  }
}

/* ===========================================================================
   boot
   =========================================================================== */

(async function boot() {
  [meta, coefs] = await Promise.all([
    fetch("data/meta.json").then(r => r.json()),
    fetch("data/coefs.json").then(r => r.json())
  ]);

  nPct = meta.pct.length;
  nDelta = meta.delta.length;
  nSim = meta.sim.length;
  nObs = meta.n_obs;
  nMu = meta.mu.length;
  nS2 = meta.sigma2.length;
  nLL = nMu * nS2;
  dMu = (meta.mu[nMu - 1] - meta.mu[0]) / (nMu - 1);
  dS2 = (meta.sigma2[nS2 - 1] - meta.sigma2[0]) / (nS2 - 1);
  PARAMS.forEach(p => { p.truth = meta.theta_true[p.key]; });

  // Open on a mild shift rather than an extreme one: the smallest positive
  // delta on the grid, otherwise the middle of the sweep.
  const mild = meta.delta.findIndex(d => d > 0);
  state.iDelta = mild >= 0 ? mild : Math.floor(nDelta / 2);
  state.deltaLo = 0;
  state.deltaHi = nDelta - 1;

  readHash();
  buildKnobs();
  window.addEventListener("hashchange", () => { readHash(); buildKnobs(); render(); });
  render();

  if (window.renderMathInElement) {
    renderMathInElement(document.body, {
      delimiters: [
        { left: "$$", right: "$$", display: true },
        { left: "$",  right: "$",  display: false }
      ],
      throwOnError: false
    });
  }
})();
