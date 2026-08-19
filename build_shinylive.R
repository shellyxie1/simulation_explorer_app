# ============================================================================
# build_shinylive.R
#
# Run this ONCE (locally, in R/RStudio, NOT in the browser) to convert
# simulation_explorer_app/ into a static Shinylive site in ./docs/
#
# Prerequisites (install once):
#   install.packages("shinylive")
#   install.packages("httpuv")   # only needed to preview locally afterwards
# ============================================================================

# ---- 0. Make sure your data files are inside the app folder -----------------
# Shinylive bundles everything inside the app directory. Your app.R reads
# paths like "outputs/02_all_data_sims.rds" (relative to the app folder), so
# copy your outputs/ directory INTO simulation_explorer_app/ before exporting:
#
#   simulation_explorer_app/
#     app.R
#     outputs/
#       02_all_data_sims.rds
#       02_ols_results_wide.rds
#       04_gel_results_bias.rds
#       04_gel_results_bias2.rds
#       04_gel_results.rds
#       04_gel_results2.rds
#
# If you haven't done that yet, uncomment and run:
# dir.create("simulation_explorer_app/outputs", recursive = TRUE, showWarnings = FALSE)
# file.copy(
#   list.files("outputs", pattern = "\\.rds$", full.names = TRUE),
#   "simulation_explorer_app/outputs",
#   overwrite = TRUE
# )

# ---- 1. Export the app to a static site --------------------------------------
unlink("docs", recursive = TRUE)

shinylive::export(
  appdir  = ".",
  destdir = "docs"
)

# ---- 2. Preview it locally before deploying ----------------------------------
# This spins up a plain static file server (NOT a Shiny/R server) so you can
# check it works exactly as it would once hosted.
httpuv::runStaticServer("docs", port = 8008, browse = TRUE)

# ---- 3. Deploy -----------------------------------------------------------------
# GitHub Pages:
#   - Commit the docs/ folder to your repo (git add docs && git commit -m "Add shinylive site")
#   - Push, then in the repo's Settings > Pages, set source = "docs" folder on
#     your main branch. The site will be live at
#     https://<username>.github.io/<repo>/
#
# Netlify / other static hosts:
#   - Just upload/deploy the contents of docs/ as-is - no build step needed,
#     it's already fully static HTML/JS/WASM.
