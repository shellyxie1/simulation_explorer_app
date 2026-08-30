# Outlier Simulation Explorer

Four linked views over a contamination study of a simple linear model
`y = β₀ + β₁x + ε` with `β = (2, −3)`: a share of each sample is replaced by
outliers drawn around a mean error `ε̄_out`, and three estimators — OLS and two
exponentially tilted empirical likelihood fits — see the same data. The sweep is
3 contamination levels × 11 values of `ε̄_out` × 50 replicates × 1000
observations.

The site is at [fleverest.github.io/simulation_explorer_app](https://fleverest.github.io/simulation_explorer_app).

## What it shows

| view | question |
| --- | --- |
| one sample | Where did the outliers land in this replicate, and what did each estimator do with them? |
| sampling distributions | How do the 50 replicates spread out as contamination gets worse, and what does that cost in bias and RMSE? |
| weights | Which observations does the tilt `gₜλ` push down, and does the downweighting track distance from the true line? |
| likelihood surface | Where is the ETEL log-likelihood actually maximised, and how far is that from the truth? |

The three global knobs — contamination level, `ε̄_out`, replicate — are shared by
every view, and are written into the URL, so a particular sample can be linked
to directly: `#weights/05/-3/17`.

## Layout

```
docs/            the site — plain static files, no build step
  index.html
  site.css
  app.js         all four views
  d3.min.js      vendored (v7.9.0), so the site runs offline
  data/
    meta.json    factor levels, the theta grid, and the binary layout
    coefs.json   1650 x 3 coefficient pairs, loaded once
    rep/*.bin    one ~19 KB binary per replicate, fetched on demand
export_data.R    outputs/*.rds -> docs/data/   (base R only)
app.R            the original Shiny app, kept as the reference for what each
                 view plots; still runnable with shiny::runApp()
outputs/         the simulation output the export reads
```

## Running it

Any static file server will do:

```sh
cd docs && python3 -m http.server 8008
```

Opening `docs/index.html` straight off disk will not work: the replicate data is
fetched, and `file://` requests are blocked as cross-origin.

## Regenerating the data

After the simulation outputs change:

```sh
Rscript export_data.R
```

It reads `outputs/*.rds` and rewrites `docs/data/`. Nothing but base R is
needed, and it takes a few seconds.

Two things are worth knowing about the payload. Only the tilt `gₜλ` is
exported for each GEL fit, not the implied probabilities, because the second is
exactly `softmax` of the first and the site derives it. And the per-observation
columns go out as little-endian float32 binaries rather than JSON: as text they
would be roughly 75 MB, and the browser would parse a megabyte of digits every
time a slider moved. `meta.json` carries the byte offsets, so changing `n_obs`
or thinning the theta grid needs no change to `app.js`.

The whole of `docs/` is about 32 MB, against 244 MB for the Shinylive export it
replaces, and it renders without downloading a WebAssembly R.

## Deploying

GitHub Pages, *Settings → Pages → Deploy from a branch*, with the branch's
`/docs` folder as the source. There is no build step, so nothing else is needed.
