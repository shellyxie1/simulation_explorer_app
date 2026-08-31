/* ===========================================================================
   Outlier Simulation Explorer

   Four views over one pre-exported simulation. The whole site is static: the
   two JSON files below are the entire index, and each replicate's per-
   observation columns arrive as a ~19 KB binary fetched on demand.

   Regenerate data/ with `Rscript export_data.R`; meta.json carries the byte
   layout so nothing here hard-codes an offset.
   =========================================================================== */

const d3 = window.d3;

/* ---- palette --------------------------------------------------------------
   Clay is spoken for: it marks outliers and true parameter values, the two
   things every view is ultimately about. Estimators are given cool hues so a
   fitted line never reads as an outlier. */
const PAPER = "#fdfcf9";
const PANEL = "#f2eee5";
const INK   = "#1f1d1a";
const MUTED = "#1f1d1a";
const RULE  = "#1f1d1a";
const AXIS  = "#1f1d1a";
const CLAY  = "#b8452f";
const CALM  = "#1f1d1a";   // ordinary observations

const METHODS = [
  { key: "OLS",      label: "OLS",       colour: "#1f6fb4" },
  { key: "gel_bias", label: "ETEL_bias", colour: "#1a8f6a" },
  { key: "gel",      label: "ETEL",      colour: "#7b3fa0" }
];
const GEL_METHODS = METHODS.filter(m => m.key !== "OLS");

const VIEWS = [
  { key: "sample",     label: "Simulated data" },
  { key: "estimates",  label: "Density of the estimates" },
  { key: "weights",    label: "Inner products and Implied probabilities" },
  { key: "likelihood", label: "Likelihood surface" }
];

const f2 = d3.format(".2f");
const f3 = d3.format(".3f");
const f4g = d3.format(".4g");
const fg  = d3.format("~g");
const f1 = d3.format(".1f");

/* ---- state ---------------------------------------------------------------- */
const state = {
  view: "sample",
  iPct: 0,
  iEps: 0,
  iSim: 0,
  methods: new Set(METHODS.map(m => m.key)),
  gel: "gel_bias",
  trim: true,
  logWeights: false,
  llWindow: 340,
  hover: null
};

let meta = null;
let coefs = null;
let nPct, nEps, nSim, nObs, nLL, n1, n2, d1, d2;

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

const slot = (iPct, iEps, iSim) => (iPct * nEps + iEps) * nSim + iSim;

const repCache = new Map();

function repUrl(iPct, iEps, iSim) {
  return `data/rep/${meta.pct[iPct]}-${meta.eps_slug[iEps]}-${meta.sim[iSim]}.bin`;
}

/* Cached by URL, and the *promise* is cached rather than its value, so two
   overlapping renders of the same replicate share one request. */
function loadRep(iPct, iEps, iSim) {
  const url = repUrl(iPct, iEps, iSim);
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
    x: new Float32Array(buf, L.x, nObs),
    y: new Float32Array(buf, L.y, nObs),
    gtRaw: {
      gel_bias: new Float32Array(buf, L.gt_bias, nObs),
      gel:      new Float32Array(buf, L.gt, nObs)
    },
    ll: new Float32Array(buf, L.ll, nLL),
    outlier: new Uint8Array(buf, L.outlier, nObs),
    _pt: {}
  };
}

/* The implied probabilities are exactly softmax(gt * lambda) -- verified
   against the fitted values to 1e-15 -- so only the tilt is shipped. Shifting
   by the maximum before exponentiating keeps the sum finite when the tilt runs
   to -20 or beyond, which it does under heavy contamination. */
function implied(rec, key) {
  if (rec._pt[key]) return rec._pt[key];
  const g = rec.gtRaw[key];
  const mx = d3.max(g);
  const w = new Float64Array(nObs);
  let s = 0;
  for (let i = 0; i < nObs; i++) { w[i] = Math.exp(g[i] - mx); s += w[i]; }
  for (let i = 0; i < nObs; i++) w[i] /= s;
  rec._pt[key] = w;
  return w;
}

/* Residual from the *true* line, not from any fit: it is the only measure of
   outlyingness that does not already depend on the estimator being judged. */
function residuals(rec) {
  if (rec._resid) return rec._resid;
  const b0 = meta.beta_true.b0, b1 = meta.beta_true.b1;
  const r = new Float64Array(nObs);
  for (let i = 0; i < nObs; i++) r[i] = rec.y[i] - (b0 + b1 * rec.x[i]);
  rec._resid = r;
  return r;
}

/* The 50 replicate estimates at one (pct, eps) for one method and parameter. */
function coefRun(mkey, param, iPct, iEps) {
  const a = coefs[mkey][param];
  const from = slot(iPct, iEps, 0);
  return a.slice(from, from + nSim);
}

const activeMethods = () => METHODS.filter(m => state.methods.has(m.key));

/* ===========================================================================
   chart plumbing
   =========================================================================== */

function frame(sel, w, h) {
  const host = d3.select(sel);
  host.selectAll("*").remove();
  return host.append("svg")
    .attr("viewBox", `0 0 ${w} ${h}`)
    .attr("width", "100%")
    .style("height", "auto")
    .style("overflow", "visible");
}

/* Every axis in the site is the same weight and colour; doing it once at the
   end of a draw beats threading the same four attributes through each call. */
function styleAxes(svg) {
  svg.selectAll(".domain, .tick line").attr("stroke", AXIS);
  svg.selectAll(".tick text")
    .attr("fill", MUTED)
    .attr("font-size", 10)
    .attr("font-family", "KaTeX_Main, Georgia, serif");
  return svg;
}

function axisLabel(g, x, y, text, anchor = "middle", rotate = 0) {
  const isTex = text.startsWith("$") && text.endsWith("$") && text.length > 2;

  if (!isTex || typeof katex === "undefined") {
    g.append("text")
      .attr("transform", `translate(${x},${y}) rotate(${rotate})`)
      .attr("text-anchor", anchor)
      .attr("font-size", 11)
      .attr("fill", MUTED)
      .text(text);
    return;
  }

  /* SVG <text> can't hold markup, so KaTeX output goes in a foreignObject:
     a window of ordinary HTML inside the drawing. Box is oversized and the
     content centred, so the label lands on (x, y) whatever its width. */
  const BW = 260, BH = 28;
  const justify = anchor === "end" ? "flex-end"
                : anchor === "start" ? "flex-start" : "center";

  const div = g.append("foreignObject")
    .attr("width", BW).attr("height", BH)
    .attr("x", -BW / 2).attr("y", -BH / 2)
    .attr("transform", `translate(${x},${y}) rotate(${rotate})`)
    .style("overflow", "visible")
    .style("pointer-events", "none")
    .append("xhtml:div")
      .style("width", `${BW}px`).style("height", `${BH}px`)
      .style("display", "flex")
      .style("align-items", "center")
      .style("justify-content", justify)
      .style("font-size", "11px")
      .style("color", MUTED);

  katex.render(text.slice(1, -1), div.node(), { throwOnError: false });
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
    if (i == null || cx == null || !isFinite(cx)) { ring.style("display", "none"); return; }
    ring.attr("cx", cx).attr("cy", cy).style("display", null);
  };
}

/* ===========================================================================
   view: one sample
   =========================================================================== */

function renderSample(rec) {
  drawSampleScatter(rec);
  drawSampleCoefs(rec);
  drawSampleResid(rec);
  sampleReadout(null, rec);
}

function drawSampleScatter(rec) {
  const W = 640, H = 470, m = { l: 54, r: 16, t: 10, b: 44 };
  const svg = frame("#sampleScatter", W, H);
  const xs = d3.scaleLinear(pad(d3.extent(rec.x)), [m.l, W - m.r]);
  const ys = d3.scaleLinear(pad(d3.extent(rec.y)), [H - m.b, m.t]);

  svg.append("rect")
    .attr("x", m.l).attr("y", m.t)
    .attr("width", W - m.r - m.l).attr("height", H - m.b - m.t)
    .attr("fill", PAPER).attr("stroke", RULE);

  svg.append("g").attr("transform", `translate(0,${H - m.b})`)
    .call(d3.axisBottom(xs).ticks(7).tickSize(4));
  svg.append("g").attr("transform", `translate(${m.l},0)`)
    .call(d3.axisLeft(ys).ticks(7).tickSize(4));
  axisLabel(svg, (m.l + W - m.r) / 2, H - 8, "$x$");
  axisLabel(svg, 15, (m.t + H - m.b) / 2, "$y$", "middle", -90);

  const clip = "clip-sample";
  svg.append("clipPath").attr("id", clip).append("rect")
    .attr("x", m.l).attr("y", m.t)
    .attr("width", W - m.r - m.l).attr("height", H - m.b - m.t);
  const plot = svg.append("g").attr("clip-path", `url(#${clip})`);

  const line = (b0, b1, colour, dash) => {
    const [x0, x1] = xs.domain();
    plot.append("line")
      .attr("x1", xs(x0)).attr("y1", ys(b0 + b1 * x0))
      .attr("x2", xs(x1)).attr("y2", ys(b0 + b1 * x1))
      .attr("stroke", colour)
      .attr("stroke-width", dash ? 1.4 : 1.8)
      .attr("stroke-dasharray", dash || null);
  };

  // Ordinary points first, outliers over them: with 10% contamination the clay
  // would otherwise be buried under the cloud it is supposed to explain.
  const idx = d3.range(nObs);
  plot.append("g").selectAll("circle")
    .data(idx.filter(i => !rec.outlier[i])).join("circle")
    .attr("cx", i => xs(rec.x[i])).attr("cy", i => ys(rec.y[i]))
    .attr("r", 2.7).attr("fill", CALM).attr("fill-opacity", 0.85);
  plot.append("g").selectAll("circle")
    .data(idx.filter(i => rec.outlier[i])).join("circle")
    .attr("cx", i => xs(rec.x[i])).attr("cy", i => ys(rec.y[i]))
    .attr("r", 4).attr("fill", CLAY).attr("fill-opacity", 0.9);

  line(meta.beta_true.b0, meta.beta_true.b1, INK, "6 4");
  const s = slot(state.iPct, state.iEps, state.iSim);
  activeMethods().forEach(mm => {
    line(coefs[mm.key].b0[s], coefs[mm.key].b1[s], mm.colour, null);
  });

  const ring = hoverRing(plot);
  const finder = d3.Delaunay.from(idx, i => xs(rec.x[i]), i => ys(rec.y[i]));
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
    ring(i, i == null ? null : xs(rec.x[i]), i == null ? null : ys(rec.y[i]));
    sampleReadout(i, rec);
  });

  styleAxes(svg);
}

/* Strip plot: where this replicate's estimate falls among the 50 drawn at the
   same setting. One lane per method, two rows for the two coefficients. */
function drawSampleCoefs(rec) {
  const W = 430, rowH = 104, gap = 16, m = { l: 10, r: 12, t: 28, b: 22 };
  const params = [
    { key: "b0", label: "β̂₀", truth: meta.beta_true.b0 },
    { key: "b1", label: "β̂₁", truth: meta.beta_true.b1 }
  ];
  const svg = frame("#sampleCoefs", W, params.length * rowH + gap);
  const shown = activeMethods();
  const s = slot(state.iPct, state.iEps, state.iSim);

  params.forEach((p, r) => {
    const g = svg.append("g").attr("transform", `translate(0,${r * (rowH + gap)})`);
    const runs = shown.map(mm => coefRun(mm.key, p.key, state.iPct, state.iEps));
    const all = runs.flat().concat([p.truth]);
    if (!all.length) return;
    const xs = d3.scaleLinear(pad(d3.extent(all), 0.06), [m.l + 4, W - m.r]).nice();

    g.append("g").attr("transform", `translate(0,${rowH - m.b})`)
      .call(d3.axisBottom(xs).ticks(5).tickSize(3));

    g.append("line")
      .attr("x1", xs(p.truth)).attr("x2", xs(p.truth))
      .attr("y1", m.t - 6).attr("y2", rowH - m.b)
      .attr("stroke", CLAY).attr("stroke-dasharray", "4 3");

    const laneH = (rowH - m.b - m.t) / Math.max(shown.length, 1);
    shown.forEach((mm, k) => {
      const cy = m.t + laneH * (k + 0.5);
      g.append("g").selectAll("line").data(runs[k]).join("line")
        .attr("x1", v => xs(v)).attr("x2", v => xs(v))
        .attr("y1", cy - 5).attr("y2", cy + 5)
        .attr("stroke", mm.colour).attr("stroke-opacity", 0.38);
      g.append("circle")
        .attr("cx", xs(coefs[mm.key][p.key][s])).attr("cy", cy)
        .attr("r", 3.6).attr("fill", mm.colour);
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

function drawSampleResid(rec) {
  const W = 430, H = 190, m = { l: 44, r: 12, t: 12, b: 32 };
  const svg = frame("#sampleResid", W, H);
  const r = residuals(rec);
  const xs = d3.scaleLinear(pad(d3.extent(r), 0.03), [m.l, W - m.r]).nice();
  const bins = d3.bin().domain(xs.domain()).thresholds(36);

  const idx = d3.range(nObs);
  const normal = bins(idx.filter(i => !rec.outlier[i]).map(i => r[i]));
  const outs   = bins(idx.filter(i => rec.outlier[i]).map(i => r[i]));
  const top = d3.max(normal, b => b.length) || 1;
  const ys = d3.scaleLinear([0, top], [H - m.b, m.t]).nice();

  svg.append("g").attr("transform", `translate(0,${H - m.b})`)
    .call(d3.axisBottom(xs).ticks(6).tickSize(4));
  svg.append("g").attr("transform", `translate(${m.l},0)`)
    .call(d3.axisLeft(ys).ticks(4).tickSize(4));
  axisLabel(svg, (m.l + W - m.r) / 2, H - 6, "$y - (\\beta_0 + \\beta_1 x)$");

  const bars = (data, fill, opacity) =>
    svg.append("g").selectAll("rect").data(data.filter(b => b.length)).join("rect")
      .attr("x", b => xs(b.x0) + 0.5)
      .attr("width", b => Math.max(1, xs(b.x1) - xs(b.x0) - 1))
      // Counts above the axis maximum are clipped rather than allowed to rescale
      // it: a single tall bin at zero would flatten everything else.
      .attr("y", b => ys(Math.min(b.length, top)))
      .attr("height", b => ys(0) - ys(Math.min(b.length, top)))
      .attr("fill", fill).attr("fill-opacity", opacity);

  bars(normal, CALM, 0.92);
  bars(outs, CLAY, 1);
  styleAxes(svg);
}

function sampleReadout(i, rec) {
  const el = d3.select("#sampleReadout");
  const nOut = d3.sum(rec.outlier);
  if (i == null) {
    el.html(
      `<span class="dim">n</span> ${nObs}` +
      `<span class="sep">|</span><span class="dim">outliers</span> ` +
      `<b style="color:${CLAY}">${nOut}</b> <span class="dim">(${f1(100 * nOut / nObs)}%)</span>` +
      `<span class="sep">|</span><span class="dim">hover a point</span>`
    );
    return;
  }
  const r = residuals(rec)[i];
  el.html(
    `<span class="dim">obs</span> ${i + 1}` +
    `<span class="sep">|</span><span class="dim">x</span> ${f3(rec.x[i])}` +
    `<span class="sep">|</span><span class="dim">y</span> ${f3(rec.y[i])}` +
    `<span class="sep">|</span><span class="dim">error</span> ${f3(r)}` +
    `<span class="sep">|</span>` +
    (rec.outlier[i]
      ? `<b style="color:${CLAY}">injected outlier</b>`
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
  const W = 780, H = 580, gutter = 34;
  const facetW = (W - gutter) / 2;
  const svg = frame("#ridges", W, H);
  const shown = activeMethods();
  if (!shown.length) {
    svg.append("text").attr("x", W / 2).attr("y", H / 2)
      .attr("text-anchor", "middle").attr("fill", MUTED).attr("font-size", 14)
      .text("Select at least one method.");
    return;
  }

  const params = [
    { key: "b0", label: "β̂₀", truth: meta.beta_true.b0 },
    { key: "b1", label: "β̂₁", truth: meta.beta_true.b1 }
  ];
  const s = slot(state.iPct, state.iEps, state.iSim);

  params.forEach((p, fi) => {
    const m = { l: fi === 0 ? 46 : 30, r: 12, t: 26, b: 38 };
    const x0 = fi * (facetW + gutter);
    const g = svg.append("g").attr("transform", `translate(${x0},0)`);

    const runs = [];
    for (let e = 0; e < nEps; e++)
      shown.forEach(mm => runs.push(coefRun(mm.key, p.key, state.iPct, e)));
    const pooled = runs.flat();

    // Contamination sends a few OLS fits a long way out; trimming keeps the
    // bulk legible, and the toggle exists because the excursions are the point
    // at the other end of the sweep.
    const dom = state.trim
      ? pad([quantileOf(pooled, 0.01), quantileOf(pooled, 0.99)], 0.06)
      : pad(d3.extent(pooled), 0.03);
    const xs = d3.scaleLinear(dom, [m.l, facetW - m.r]);
    const grid = d3.range(200).map(i => dom[0] + (dom[1] - dom[0]) * i / 199);

    const RIDGE = 1.8;                       // ridge height, in row heights
    const rowH = (H - m.t - m.b) / (nEps + RIDGE - 0.5);
    const rowY = e => H - m.b - (e + 0.5) * rowH;

    // One density per (eps, method), scaled against the tallest in the facet so
    // the ridges share a vertical unit -- ggridges' default, and the reason a
    // flat ridge reads as "spread out" rather than "different scale".
    const dens = [];
    let top = 0;
    for (let e = 0; e < nEps; e++) {
      shown.forEach(mm => {
        const z = kde(coefRun(mm.key, p.key, state.iPct, e), grid);
        top = Math.max(top, d3.max(z));
        dens.push({ e, mm, z });
      });
    }
    const ridge = rowH * RIDGE;
    const ys = z => (top > 0 ? (z / top) * ridge : 0);

    g.append("rect")
      .attr("x", m.l).attr("y", m.t)
      .attr("width", facetW - m.r - m.l).attr("height", H - m.b - m.t)
      .attr("fill", PAPER).attr("stroke", RULE);

    // Band the selected row so the sampling distribution the other views are
    // drawn from is identifiable at a glance.
    g.append("rect")
      .attr("x", m.l).attr("y", rowY(state.iEps) - rowH / 2)
      .attr("width", facetW - m.r - m.l).attr("height", rowH)
      .attr("fill", PANEL);

    const clip = `clip-ridge-${fi}`;
    g.append("clipPath").attr("id", clip).append("rect")
      .attr("x", m.l).attr("y", m.t)
      .attr("width", facetW - m.r - m.l).attr("height", H - m.b - m.t);
    const plot = g.append("g").attr("clip-path", `url(#${clip})`);

    for (let e = 0; e < nEps; e++) {
      plot.append("line")
        .attr("x1", m.l).attr("x2", facetW - m.r)
        .attr("y1", rowY(e)).attr("y2", rowY(e))
        .attr("stroke", RULE).attr("stroke-width", 0.7);
    }

    plot.append("line")
      .attr("x1", xs(p.truth)).attr("x2", xs(p.truth))
      .attr("y1", m.t).attr("y2", H - m.b)
      .attr("stroke", CLAY).attr("stroke-dasharray", "5 4").attr("stroke-width", 1.1);

    // Drawn from the top row down, so nearer (lower) ridges overlap the ones
    // behind them rather than being hidden by them.
    const area = d3.area()
      .x((_, i) => xs(grid[i]))
      .y1(v => -ys(v))
      .y0(0)
      .curve(d3.curveBasis);
    const outline = d3.line()
      .x((_, i) => xs(grid[i]))
      .y(v => -ys(v))
      .curve(d3.curveBasis);

    dens.sort((a, b) => d3.descending(a.e, b.e));
    dens.forEach(({ e, mm, z }) => {
      const row = plot.append("g").attr("transform", `translate(0,${rowY(e)})`);
      row.append("path").datum(z)
        .attr("d", area)
        .attr("fill", mm.colour).attr("fill-opacity", 0.20)
        .attr("stroke", "none");
      row.append("path").datum(z)
        .attr("d", outline)
        .attr("fill", "none")
        .attr("stroke", mm.colour).attr("stroke-width", 1.2);
    });

    // The selected replicate, on the selected row: the same number the strip
    // plot in the sample view shows as a filled dot.
    shown.forEach(mm => {
      const v = coefs[mm.key][p.key][s];
      plot.append("line")
        .attr("x1", xs(v)).attr("x2", xs(v))
        .attr("y1", rowY(state.iEps) - rowH * 0.42)
        .attr("y2", rowY(state.iEps) + rowH * 0.42)
        .attr("stroke", mm.colour).attr("stroke-width", 2);
    });

    g.append("g").attr("transform", `translate(0,${H - m.b})`)
      .call(d3.axisBottom(xs).ticks(5).tickSize(4));
    if (fi === 0) {
      g.append("line")
        .attr("x1", m.l).attr("x2", m.l)
        .attr("y1", m.t).attr("y2", H - m.b)
        .attr("stroke", AXIS);
      g.append("g").selectAll("text").data(d3.range(nEps)).join("text")
        .attr("x", m.l - 8).attr("y", e => rowY(e) + 3.5)
        .attr("text-anchor", "end")
        .attr("font-size", 10).attr("font-family", "ui-monospace, Menlo, monospace")
        .attr("fill", e => e === state.iEps ? INK : MUTED)
        .text(e => meta.eps[e]);
      axisLabel(g, 13, H / 2, "$\\mu_{\\epsilon_{\\text{out}}}$", "middle", -90);
    }
    g.append("text")
      .attr("x", m.l).attr("y", 14)
      .attr("font-size", 12).attr("fill", INK)
      .text(`${p.label}   true ${p.truth}`);
  });

  styleAxes(svg);
}

/* Bias and RMSE across the sweep, one small panel each per coefficient, in the
   same stacked idiom as the fit diagnostics in ripr-vis. */
function drawEstPanels() {
  const W = 420, H = 126, gap = 18, m = { l: 54, r: 10, t: 16, b: 22 };
  const shown = activeMethods();
  const panels = [
    { stat: "bias", param: "b0", label: "bias  β₀" },
    { stat: "bias", param: "b1", label: "bias  β₁" },
    { stat: "rmse", param: "b0", label: "RMSE  β₀" },
    { stat: "rmse", param: "b1", label: "RMSE  β₁" }
  ];
  const svg = frame("#estPanels", W, panels.length * H + (panels.length - 1) * gap);
  if (!shown.length) return;

  const xs = d3.scaleLinear(d3.extent(meta.eps), [m.l, W - m.r]);

  panels.forEach((p, r) => {
    const truth = meta.beta_true[p.param];
    const series = shown.map(mm => ({
      mm,
      v: d3.range(nEps).map(e => {
        const run = coefRun(mm.key, p.param, state.iPct, e);
        return {
          e,
          y: p.stat === "bias"
            ? d3.mean(run) - truth
            : Math.sqrt(d3.mean(run, v => (v - truth) ** 2))
        };
      })
    }));

    const g = svg.append("g").attr("transform", `translate(0,${r * (H + gap)})`);
    const vals = series.flatMap(s => s.v.map(d => d.y));
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
      .call(d3.axisBottom(xs).tickValues(meta.eps.filter((_, i) => i % 2 === 0)).tickSize(3));
    g.append("g").attr("transform", `translate(${m.l},0)`)
      .call(d3.axisLeft(ys).ticks(3).tickSize(3).tickFormat(d3.format(".2g")));

    const line = d3.line().x(d => xs(meta.eps[d.e])).y(d => ys(d.y));
    series.forEach(s => {
      g.append("path").datum(s.v)
        .attr("fill", "none").attr("stroke", s.mm.colour)
        .attr("stroke-width", 1.5).attr("d", line);
      const here = s.v[state.iEps];
      g.append("circle")
        .attr("cx", xs(meta.eps[here.e])).attr("cy", ys(here.y))
        .attr("r", 3).attr("fill", s.mm.colour);
    });

    g.append("line")
      .attr("x1", xs(meta.eps[state.iEps])).attr("x2", xs(meta.eps[state.iEps]))
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
  const gt = rec.gtRaw[state.gel];
  const pt = implied(rec, state.gel);
  const res = residuals(rec);
  const marks = [];

  marks.push(drawIndexPanel("#gtPlot", rec, gt, {
    W: 640, H: 470, rule: 0, ruleLabel: "0", log: false, big: true
  }));
  marks.push(drawIndexPanel("#ptPlot", rec, pt, {
    W: 640, H: 470, rule: 1 / nObs, ruleLabel: "1/N",
    log: state.logWeights, big: true
  }));
  marks.push(drawWeightResid(rec, pt, res));
  marks.push(drawWeightLorenz(rec, pt));

  hoverHandlers.push(i => {
    marks.forEach(fn => fn(i));
    weightsReadout(i, rec, gt, pt, res);
  });
  weightsReadout(null, rec, gt, pt, res);
}

function drawIndexPanel(sel, rec, values, opt) {
  const { W, H, rule, ruleLabel, log, big } = opt;
  const m = { l: big ? 60 : 56, r: 14, t: 12, b: big ? 42 : 32 };
  const svg = frame(sel, W, H);
  const xs = d3.scaleLinear([1, nObs], [m.l, W - m.r]);

  const positive = Array.from(values).filter(v => v > 0);
  const ys = log && positive.length
    ? d3.scaleLog([d3.min(positive), d3.max(values)], [H - m.b, m.t])
    : d3.scaleLinear(pad(d3.extent(values), 0.05), [H - m.b, m.t]).nice();

  svg.append("rect")
    .attr("x", m.l).attr("y", m.t)
    .attr("width", W - m.r - m.l).attr("height", H - m.b - m.t)
    .attr("fill", PAPER).attr("stroke", RULE);

  svg.append("g").attr("transform", `translate(0,${H - m.b})`)
    .call(d3.axisBottom(xs).ticks(big ? 8 : 5).tickSize(4));
  svg.append("g").attr("transform", `translate(${m.l},0)`)
    .call(log
      ? d3.axisLeft(ys).ticks(4, "0.0e").tickSize(4)
      : d3.axisLeft(ys).ticks(big ? 6 : 4).tickSize(4).tickFormat(d3.format(".3~g")));
  if (big) axisLabel(svg, (m.l + W - m.r) / 2, H - 8, "Observation index");

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
  const r = big ? 2.6 : 1.9;
  plot.append("g").selectAll("circle")
    .data(idx.filter(i => !rec.outlier[i])).join("circle")
    .attr("cx", i => xs(i + 1)).attr("cy", i => ys(values[i]))
    .attr("r", r).attr("fill", CALM).attr("fill-opacity", 0.82);
  plot.append("g").selectAll("circle")
    .data(idx.filter(i => rec.outlier[i])).join("circle")
    .attr("cx", i => xs(i + 1)).attr("cy", i => ys(values[i]))
    .attr("r", r + 1.4).attr("fill", CLAY).attr("fill-opacity", 1);

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

function drawWeightResid(rec, pt, res) {
  const W = 640, H = 470, m = { l: 60, r: 14, t: 12, b: 42 };
  const svg = frame("#ptResid", W, H);
  const xs = d3.scaleLinear(pad(d3.extent(res), 0.04), [m.l, W - m.r]).nice();
  const positive = Array.from(pt).filter(v => v > 0);
  const ys = state.logWeights && positive.length
    ? d3.scaleLog([d3.min(positive), d3.max(pt)], [H - m.b, m.t])
    : d3.scaleLinear(pad(d3.extent(pt), 0.05), [H - m.b, m.t]).nice();

  svg.append("rect")
    .attr("x", m.l).attr("y", m.t)
    .attr("width", W - m.r - m.l).attr("height", H - m.b - m.t)
    .attr("fill", PAPER).attr("stroke", RULE);

  svg.append("g").attr("transform", `translate(0,${H - m.b})`)
    .call(d3.axisBottom(xs).ticks(6).tickSize(4));
  svg.append("g").attr("transform", `translate(${m.l},0)`)
    .call(state.logWeights
      ? d3.axisLeft(ys).ticks(4, "0.0e").tickSize(4)
      : d3.axisLeft(ys).ticks(4).tickSize(4).tickFormat(d3.format(".3~g")));
  axisLabel(svg, (m.l + W - m.r) / 2, H - 6, "$y - (\\beta_0 + \\beta_1 x)$");

  const clip = "clip-ptresid";
  svg.append("clipPath").attr("id", clip).append("rect")
    .attr("x", m.l).attr("y", m.t)
    .attr("width", W - m.r - m.l).attr("height", H - m.b - m.t);
  const plot = svg.append("g").attr("clip-path", `url(#${clip})`);

  const uni = 1 / nObs;
  if (uni >= ys.domain()[0] && uni <= ys.domain()[1]) {
    plot.append("line")
      .attr("x1", m.l).attr("x2", W - m.r)
      .attr("y1", ys(uni)).attr("y2", ys(uni))
      .attr("stroke", MUTED).attr("stroke-dasharray", "4 3");
  }

  const idx = d3.range(nObs);
  plot.append("g").selectAll("circle")
    .data(idx.filter(i => !rec.outlier[i])).join("circle")
    .attr("cx", i => xs(res[i])).attr("cy", i => ys(pt[i]))
    .attr("r", 1.9).attr("fill", CALM).attr("fill-opacity", 0.82);
  plot.append("g").selectAll("circle")
    .data(idx.filter(i => rec.outlier[i])).join("circle")
    .attr("cx", i => xs(res[i])).attr("cy", i => ys(pt[i]))
    .attr("r", 3.2).attr("fill", CLAY).attr("fill-opacity", 1);

  const ring = hoverRing(plot);
  styleAxes(svg);
  return i => ring(i, i == null ? null : xs(res[i]), i == null ? null : ys(pt[i]));
}

/* Lorenz curve for the implied probabilities: observations sorted by weight,
   largest first, against cumulative share of the total. Uniform weighting is
   the diagonal, so the gap between curve and diagonal is the departure from
   equal weighting — the scalar in the readout is one point on this curve. */
function drawWeightLorenz(rec, pt) {
  const W = 640, H = 470, m = { l: 60, r: 14, t: 12, b: 42 };
  const svg = frame("#ptLorenz", W, H);

  // Sort indices by weight descending, then accumulate.
  const order = d3.range(nObs).sort((a, b) => pt[b] - pt[a]);
  const total = d3.sum(pt) || 1;
  const cum = new Float64Array(nObs);
  let run = 0;
  order.forEach((i, k) => { run += pt[i] / total; cum[k] = run; });

  // rank[i] is where observation i sits in the sorted order, so a hover on
  // any other panel can find its point on this curve.
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
  axisLabel(svg, (m.l + W - m.r) / 2, H - 8, "share of the sample, heaviest first");
  axisLabel(svg, 15, (m.t + H - m.b) / 2, "share of total weight", "middle", -90);

  const clip = "clip-ptlorenz";
  svg.append("clipPath").attr("id", clip).append("rect")
    .attr("x", m.l).attr("y", m.t)
    .attr("width", W - m.r - m.l).attr("height", H - m.b - m.t);
  const plot = svg.append("g").attr("clip-path", `url(#${clip})`);

  // Uniform reference: every observation holds 1/N, so the curve is y = x.
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

  // Mark where the contaminated observations run out, if any are present.
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

function weightsReadout(i, rec, gt, pt, res) {
  const el = d3.select("#weightsReadout");
  let outW = 0, nOut = 0;
  for (let k = 0; k < nObs; k++) if (rec.outlier[k]) { outW += pt[k]; nOut++; }
  const share = nOut ? outW / (nOut / nObs) : 0;   // 1 means "no downweighting"

  const summary =
    `<span class="dim">outliers hold</span> <b style="color:${CLAY}">${f2(100 * outW)}%</b>` +
    ` <span class="dim">of the weight, against</span> ${f1(100 * nOut / nObs)}%` +
    ` <span class="dim">of the sample &mdash; a factor of</span> <b>${f2(share)}</b>`;

  if (i == null) { el.html(summary); return; }
  el.html(
    `<span class="dim">obs</span> ${i + 1}` +
    `<span class="sep">|</span><span class="dim">gₜλ</span> ${f3(gt[i])}` +
    `<span class="sep">|</span><span class="dim">pₜ</span> ${d3.format(".3e")(pt[i])}` +
    ` <span class="dim">(${f2(pt[i] * nObs)}× uniform)</span>` +
    `<span class="sep">|</span><span class="dim">error</span> ${f3(res[i])}` +
    (rec.outlier[i] ? `<span class="sep">|</span><b style="color:${CLAY}">outlier</b>` : "")
  );
}

/* ===========================================================================
   view: likelihood surface
   =========================================================================== */

function renderLikelihood(rec) {
  // The export is row-major over th1; d3.contours wants the first axis varying
  // fastest, so transpose once here and th1 becomes the horizontal axis.
  const z = new Float64Array(nLL);
  for (let i1 = 0; i1 < n1; i1++)
    for (let i2 = 0; i2 < n2; i2++)
      z[i2 * n1 + i1] = rec.ll[i1 * n2 + i2];

  let iMax = 0;
  for (let k = 1; k < nLL; k++) if (z[k] > z[iMax]) iMax = k;
  const maxLL = z[iMax];
  const maxI1 = iMax % n1, maxI2 = (iMax - maxI1) / n1;

  const floor = maxLL - state.llWindow;
  const zc = Float64Array.from(z, v => Math.max(v, floor));

  drawLLSurface(zc, { maxLL, floor, maxI1, maxI2 });
  drawLLProfiles(z, { maxLL, floor, maxI1, maxI2 });

  const llTrue = bilinear(z, meta.beta_true.b0, meta.beta_true.b1);
  d3.select("#llReadout").html(
    `<span class="dim">grid max</span> <b>${f1(maxLL)}</b>` +
    ` <span class="dim">at</span> (${meta.th1[maxI1]}, ${meta.th2[maxI2]})` +
    `<span class="sep">|</span><span class="dim">at the truth</span> ` +
    `<b style="color:${CLAY}">${llTrue == null ? "—" : f1(llTrue)}</b>` +
    (llTrue == null ? "" :
      `<span class="sep">|</span><span class="dim">difference</span> ${f1(maxLL - llTrue)}`)
  );
}

/* Bilinear read of the surface at an arbitrary (th1, th2), for the truth
   marker and for the hover readout; the grid is evenly spaced in both. */
function bilinear(z, t1, t2) {
  const u = (t1 - meta.th1[0]) / d1, v = (t2 - meta.th2[0]) / d2;
  if (u < 0 || v < 0 || u > n1 - 1 || v > n2 - 1) return null;
  const i = Math.min(Math.floor(u), n1 - 2), j = Math.min(Math.floor(v), n2 - 2);
  const a = u - i, b = v - j;
  const at = (p, q) => z[q * n1 + p];
  const out = at(i, j) * (1 - a) * (1 - b) + at(i + 1, j) * a * (1 - b) +
              at(i, j + 1) * (1 - a) * b + at(i + 1, j + 1) * a * b;
  // The far corners of the grid are genuinely -Inf, and a corner touching the
  // cell makes the whole interpolation meaningless rather than merely small.
  return isFinite(out) ? out : null;
}

/* Bounding box of every contour ring, in grid-index units — the region
   where at least one threshold line actually falls. Used to zoom the axes
   to where there's something to see. */
function contourBBox(contours) {
  let x0 = Infinity, x1 = -Infinity, y0 = Infinity, y1 = -Infinity;
  contours.forEach(c => c.coordinates.forEach(poly => poly.forEach(ring =>
    ring.forEach(([x, y]) => {
      if (x < x0) x0 = x; if (x > x1) x1 = x;
      if (y < y0) y0 = y; if (y > y1) y1 = y;
    })
  )));
  return isFinite(x0) ? [x0, x1, y0, y1] : null;
}

function drawLLSurface(zc, info) {
  const W = 640, H = 470, m = { l: 54, r: 16, t: 12, b: 44 };
  const svg = frame("#llSurface", W, H);

  const levels = d3.ticks(info.floor, info.maxLL, 11).filter(v => v > info.floor);
  const contours = d3.contours().size([n1, n2]).thresholds(levels)(zc);

  // Zoom the axes to where the contours actually sit, so the corners of the
  // grid that are entirely at the floor value aren't shown as dead space.
  const bbox = contourBBox(contours);
  let dom1 = [meta.th1[0], meta.th1[n1 - 1]];
  let dom2 = [meta.th2[0], meta.th2[n2 - 1]];
  if (bbox) {
    const [gx0, gx1, gy0, gy1] = bbox;
    const padX = (gx1 - gx0) * 0.08 || 1;
    const padY = (gy1 - gy0) * 0.08 || 1;
    dom1 = [
      Math.max(meta.th1[0], meta.th1[0] + (gx0 - padX) * d1),
      Math.min(meta.th1[n1 - 1], meta.th1[0] + (gx1 + padX) * d1)
    ];
    dom2 = [
      Math.max(meta.th2[0], meta.th2[0] + (gy0 - padY) * d2),
      Math.min(meta.th2[n2 - 1], meta.th2[0] + (gy1 + padY) * d2)
    ];
  }
  const xs = d3.scaleLinear(dom1, [m.l, W - m.r]);
  const ys = d3.scaleLinear(dom2, [H - m.b, m.t]);

  // Contour coordinates arrive in grid-index units.
  const path = d3.geoPath(d3.geoTransform({
    point(x, y) { this.stream.point(xs(meta.th1[0] + x * d1), ys(meta.th2[0] + y * d2)); }
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

  plot.append("g").selectAll("path").data(contours).join("path")
    .attr("d", path)
    .attr("fill", c => fill(c.value))
    .attr("stroke", INK)
    .attr("stroke-opacity", 0.4)
    .attr("stroke-width", 0.95);

  const anchorX = xs(meta.th1[info.maxI1]), anchorY = ys(meta.th2[info.maxI2]);
  labelContours(plot, contours, xs, ys, anchorX, anchorY);

  const truth = [meta.beta_true.b0, meta.beta_true.b1];
  const mark = (cx, cy, colour, text, dy) => {
    plot.append("circle").attr("cx", cx).attr("cy", cy)
      .attr("r", 4.6).attr("fill", colour)
      .attr("stroke", PAPER).attr("stroke-width", 1.2);
    plot.append("text")
      .attr("x", Math.max(m.l + 26, Math.min(W - m.r - 26, cx)))
      .attr("y", cy + dy)
      .attr("text-anchor", "middle").attr("font-size", 11).attr("fill", colour)
      .attr("paint-order", "stroke").attr("stroke", PAPER).attr("stroke-width", 3.5)
      .text(text);
  };
  mark(anchorX, anchorY, INK, "max", -11);
  mark(xs(truth[0]), ys(truth[1]), CLAY, "true β", 18);

  svg.append("g").attr("transform", `translate(0,${H - m.b})`)
    .call(d3.axisBottom(xs).ticks(8).tickSize(4));
  svg.append("g").attr("transform", `translate(${m.l},0)`)
    .call(d3.axisLeft(ys).ticks(8).tickSize(4));
  axisLabel(svg, (m.l + W - m.r) / 2, H - 8, "$\\beta_0$");
  axisLabel(svg, 15, (m.t + H - m.b) / 2, "$\\beta_1$", "middle", -90);

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
      const v = bilinear(zc, t1, t2);
      cross.style("display", null);
      cx.attr("cx", px).attr("cy", py);
      d3.select("#llReadout").select(".probe").remove();
      d3.select("#llReadout").append("span").attr("class", "probe")
        .html(`<span class="sep">|</span><span class="dim">at</span> ` +
              `(${f2(t1)}, ${f2(t2)}) <span class="dim">ll</span> ` +
              `${v == null ? "—" : f1(v)}`);
    })
    .on("pointerleave", () => {
      cross.style("display", "none");
      d3.select("#llReadout").select(".probe").remove();
    });

  styleAxes(svg);
}

/* metR's labelled contours, rebuilt. Placing a label a fixed fraction of the
   way round each ring knots them all together on an elongated surface; taking
   the first crossing of the horizontal line through the maximum instead puts
   them in a readable row, because nested contours cross that line at
   increasing distance. Horizontal text with a halo: at such a crossing the
   contour is close to vertical, so rotating to the tangent is the worst
   possible angle. */
function labelContours(g, contours, xs, ys, x0, y0) {
  const layer = g.append("g");
  const px = ([gx, gy]) => [xs(meta.th1[0] + gx * d1), ys(meta.th2[0] + gy * d2)];
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
  const W = 420, H = 190, gap = 26, m = { l: 58, r: 12, t: 18, b: 30 };
  const svg = frame("#llProfiles", W, 2 * H + gap);

  const panels = [
    {
      label: `β₀ at β₁ = ${meta.th2[info.maxI2]}`,
      axis: meta.th1,
      v: d3.range(n1).map(i => Math.max(z[info.maxI2 * n1 + i], info.floor)),
      truth: meta.beta_true.b0,
      at: info.maxI1
    },
    {
      label: `β₁ at β₀ = ${meta.th1[info.maxI1]}`,
      axis: meta.th2,
      v: d3.range(n2).map(j => Math.max(z[j * n1 + info.maxI1], info.floor)),
      truth: meta.beta_true.b1,
      at: info.maxI2
    }
  ];

  panels.forEach((p, r) => {
    const g = svg.append("g").attr("transform", `translate(0,${r * (H + gap)})`);
    const xs = d3.scaleLinear(d3.extent(p.axis), [m.l, W - m.r]);
    // Clipped at the same window as the contours, so the two panels agree about
    // what counts as "near the peak".
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

    plot.append("line")
      .attr("x1", xs(p.truth)).attr("x2", xs(p.truth))
      .attr("y1", m.t).attr("y2", H - m.b)
      .attr("stroke", CLAY).attr("stroke-dasharray", "4 3");

    plot.append("path")
      .datum(p.v.map((y, i) => ({ x: p.axis[i], y })))
      .attr("fill", "none").attr("stroke", INK).attr("stroke-width", 1.5)
      .attr("d", d3.line().x(d => xs(d.x)).y(d => ys(d.y)));

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
    meta.pct.map((p, i) => ({ label: `$${+p}\\%$`, i })),
    d => d.i === state.iPct,
    d => { state.iPct = d.i; render(); });

  d3.select("#epsRange")
    .attr("max", nEps - 1).property("value", state.iEps)
    .on("input", function () { state.iEps = +this.value; render(); });
  d3.select("#simRange")
    .attr("max", nSim - 1).property("value", state.iSim)
    .on("input", function () { state.iSim = +this.value; render(); });

  // Per-view knobs, rebuilt on each render so their state always matches.
  const sk = d3.select("#sampleKnobs"); sk.selectAll("*").remove();
  sk.append("span").attr("class", "ctl").append("span").attr("class", "lbl").text("fits shown");
  sk.append("span").attr("id", "sampleLegend");
  methodLegend("#sampleLegend");

  const ek = d3.select("#estKnobs"); ek.selectAll("*").remove();
  ek.append("span").attr("class", "ctl").append("span").attr("class", "lbl").text("methods");
  ek.append("span").attr("id", "estLegend");
  methodLegend("#estLegend");
  checkbox(ek, "trim the outer 1% of estimates", () => state.trim, v => state.trim = v);

  const wk = d3.select("#weightsKnobs"); wk.selectAll("*").remove();
  const wc = wk.append("span").attr("class", "ctl");
  wc.append("span").attr("class", "lbl").text("fit");
  wc.append("span").attr("id", "gelTabs").attr("class", "tabs");
  tabGroup("#gelTabs", GEL_METHODS.map(m => ({ label: m.label, key: m.key })),
    d => d.key === state.gel,
    d => { state.gel = d.key; render(); });
  checkbox(wk, "log scale on the weights", () => state.logWeights, v => state.logWeights = v);

  const lk = d3.select("#llKnobs"); lk.selectAll("*").remove();
  const lc = lk.append("label").attr("class", "ctl");
  lc.append("span").attr("class", "lbl").text("contour window");
  lc.append("input").attr("type", "range").attr("min", 20).attr("max", 800)
    .attr("step", 20).property("value", state.llWindow)
    .on("input", function () { state.llWindow = +this.value; render(); });
  lc.append("span").attr("class", "val").attr("id", "llWindowVal").text(state.llWindow);
}

function syncKnobs() {
  d3.selectAll("#views .tab").classed("on", (d) => d.key === state.view);
  d3.selectAll("#pctTabs .tab").classed("on", d => d.i === state.iPct);
  d3.select("#epsRange").property("value", state.iEps);
  d3.select("#simRange").property("value", state.iSim);
  const e = meta.eps[state.iEps];
  d3.select("#epsVal").text(e > 0 ? `+${fg(e)}` : fg(e));
  d3.select("#simVal").text(`${meta.sim[state.iSim]} / ${nSim}`);
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
   replicate can be sent to someone. */
function writeHash() {
  const h = `#${state.view}/${meta.pct[state.iPct]}/${meta.eps_slug[state.iEps]}/${meta.sim[state.iSim]}`;
  if (location.hash !== h) history.replaceState(null, "", h);
}

function readHash() {
  const parts = location.hash.replace(/^#/, "").split("/");
  if (parts.length < 4) return;
  const [v, p, e, s] = parts;
  if (VIEWS.some(x => x.key === v)) state.view = v;
  const ip = meta.pct.indexOf(p);          if (ip >= 0) state.iPct = ip;
  const ie = meta.eps_slug.indexOf(e);     if (ie >= 0) state.iEps = ie;
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
    rec = await loadRep(state.iPct, state.iEps, state.iSim);
  } catch (err) {
    if (token !== renderToken) return;
    d3.select(`main .view[data-view="${state.view}"] .figure`)
      .html(`<p class="caption">Could not load ${repUrl(state.iPct, state.iEps, state.iSim)}.</p>`);
    return;
  }
  if (token !== renderToken) return;   // a knob moved while the fetch was in flight

  if (state.view === "sample") renderSample(rec);
  else if (state.view === "weights") renderWeights(rec);
  else if (state.view === "likelihood") renderLikelihood(rec);

  prefetch();
}

/* Warm the neighbours of the current cell, so dragging either slider one step
   is already resolved by the time the pointer stops. Failures are ignored:
   this is speculation, not a dependency. */
function prefetch() {
  const near = [
    [state.iPct, state.iEps, state.iSim + 1],
    [state.iPct, state.iEps, state.iSim - 1],
    [state.iPct, state.iEps + 1, state.iSim],
    [state.iPct, state.iEps - 1, state.iSim]
  ];
  for (const [p, e, s] of near) {
    if (e < 0 || e >= nEps || s < 0 || s >= nSim) continue;
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
  nEps = meta.eps.length;
  nSim = meta.sim.length;
  nObs = meta.n_obs;
  n1 = meta.th1.length;
  n2 = meta.th2.length;
  nLL = n1 * n2;
  d1 = (meta.th1[n1 - 1] - meta.th1[0]) / (n1 - 1);
  d2 = (meta.th2[n2 - 1] - meta.th2[0]) / (n2 - 1);

  // Open on a clean sample rather than an extreme one: eps = 0 if it is on the
  // grid, otherwise the middle of the sweep.
  const zero = meta.eps.indexOf(0);
  state.iEps = zero >= 0 ? zero : Math.floor(nEps / 2);

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
