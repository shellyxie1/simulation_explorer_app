# Outlier Simulation Explorers

Two interactive simulation studies of how outliers move exponentially tilted
empirical likelihood (ETEL) estimates, published as one static site:

| study | model | what is contaminated | estimators |
| --- | --- | --- | --- |
| [linear](linear/) | `y = β₀ + β₁x + ε`, `β = (2, −3)` | a share of errors drawn around a shifted mean `ε̄_out` | OLS, ETEL from a biased start, ETEL from the truth |
| [lognormal](lognormal/) | `y ~ LN(μ, σ²)`, `(μ, σ²) = (0, 0.25)`, fitted through its first three raw moments | a share of observations with log-mean shifted by `δ` | ETEL from a biased start, ETEL from the truth |

Each study is 3 contamination levels × a sweep of the shift × 50 replicates of
1000 observations, and each explorer has the same four views: one sample, the
sampling distributions of the estimates, the ETEL weights, and the
log-likelihood surface. A fifth view, *Posterior means* (the second tab), shows kernel densities
of the BETEL posterior means over the 50 replications, clean against 10%
contaminated, from the thesis MCMC runs. The three global knobs are written into the URL, so a
particular sample can be linked to directly, for example
`lognormal/#weights/05/-3/17`.

The site is at [shellyxie1.github.io/simulation_explorer_app](https://shellyxie1.github.io/simulation_explorer_app).

## Layout

```
docs/                  the site — plain static files, no build step
  index.html           landing page linking the two explorers
  site.css             shared stylesheet
  d3.min.js            vendored (v7.9.0), shared, so the site runs offline
  linear/
    index.html, app.js
    data/              meta.json, coefs.json, postmean.json, rep/*.bin
  lognormal/
    index.html, app.js
    data/              meta.json, coefs.json, postmean.json, rep/*.bin
linear/                the linear-regression pipeline
  export_data.R        outputs/*.rds -> ../docs/linear/data/   (base R only)
  export_postmean.R    thesis MCMC draws -> ../docs/linear/data/postmean.json
  app.R                the original Shiny app, kept as the reference for each view
  prepare_lite_outputs.R
  outputs/             the simulation output the export reads
lognormal/             the log-normal pipeline
  00_*.R … 06_*.R      simulation, ETEL fitting, log-likelihood surface
  export_data.R        outputs/*.rds -> ../docs/lognormal/data/   (base R only)
  export_postmean.R    thesis MCMC means -> ../docs/lognormal/data/postmean.json
  app.R                the original Shiny app
  outputs/
```

## Running it

Any static file server will do, started from `docs/`:

```sh
cd docs && python3 -m http.server 8008
```

Then open <http://localhost:8008/>. Opening the HTML files straight off disk
will not work: the replicate data is fetched, and `file://` requests are
blocked as cross-origin.

## Regenerating the data

Each pipeline exports from inside its own folder:

```sh
cd linear    && Rscript export_data.R
cd lognormal && Rscript export_data.R
```

Both read `outputs/*.rds` and rewrite the matching `docs/*/data/`. Nothing but
base R is needed, and each takes a few seconds.

The payload design is the same in both. Only the tilt `λᵀgᵢ` is exported for
each ETEL fit, not the implied probabilities, because `log πᵢ = λᵀgᵢ −
logsumexp(λᵀg)` and the site derives it. The per-observation columns go out as
little-endian float32 binaries rather than JSON, and `meta.json` carries the
byte offsets, so a change in `n_obs` or the parameter grid needs no change on
the JavaScript side. The log-normal grid additionally has infeasible cells
(where the zero vector is outside the convex hull of the moment conditions);
they travel as `NaN` plus an explicit `hull` flag and are shaded on the surface.

The posterior means have one more export each, run the same way:

```sh
cd linear    && Rscript export_postmean.R
cd lognormal && Rscript export_postmean.R
```

These read the MCMC output from the thesis folder (`~/Desktop/thesis/...`,
set at the top of each script), not from `outputs/`, because the draws are too
large to keep here. Only the posterior means, one per replication, are written
to `docs/*/data/postmean.json`; the densities are drawn in the browser.

## Deploying

GitHub Pages, *Settings → Pages → Deploy from a branch*, with the branch's
`/docs` folder as the source. There is no build step, so nothing else is needed.
