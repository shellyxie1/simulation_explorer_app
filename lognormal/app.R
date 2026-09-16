# =============================================================================
# Log-normal Outlier Simulation Explorer
#
# Adapted from the linear-regression explorer. Differences:
#   * design is (pct, delta, sim); parameters are (mu, sigma2); no OLS
#   * fits come from the manual ETEL fitter (03_ETEL_truth.R / 03_ETEL_bias.R):
#       fit_ETEL = list(coef, loglike, conv, df = tibble(index, inner_pro, logpi))
#     so one file per method serves both the density and diagnostic tabs
#   * implied probabilities can be shown on the log scale (collapsed fits
#     underflow to 0 on the raw scale)
#   * LL surface comes from 06_ll_surface.rds with columns mu, sigma2, hull, ll;
#     infeasible cells are shaded and excluded from the contours
# =============================================================================

library(shiny)
library(bslib)
library(thematic)
library(tidyverse)
library(ggridges)
library(latex2exp)
library(metR)

thematic_shiny(font = "auto", bg = "white", fg = "black", accent = "#2C3E50")

# ============================== CONSTANTS ====================================

data_dir <- Sys.getenv("SIM_OUTPUT_DIR", unset = "outputs")

theta_true  <- c(mu = 0, sigma2 = 0.25)
gel_methods <- c("gel_bias", "gel")

outlier_cols <- c(`FALSE` = "steelblue", `TRUE` = "red")
outlier_labs <- c(`FALSE` = "Normal",    `TRUE` = "Outlier")

method_labels <- c(
  gel_bias = "ETEL (biased start)",
  gel      = "ETEL (start at truth)"
)

# ============================== LOAD DATA ====================================

read_output <- function(file) {
  path <- file.path(data_dir, file)
  if (!file.exists(path)) {
    stop(
      str_glue(
        "Missing data file: {path}\n",
        "Launch the app from the project root, or set SIM_OUTPUT_DIR."
      ),
      call. = FALSE
    )
  }
  read_rds(path)
}

all_data <- read_output("02_all_data_sims.rds")

gel_by_method <- list(
  gel_bias = read_output("04_gel_results_bias.rds"),
  gel      = read_output("04_gel_results.rds")
)

ll_surface <- read_output("06_ll_surface.rds")

# ============================== PREP =========================================

#' Extract (mu_hat, sigma2_hat) from a results tibble whose `coef` column is a
#' list of length-2 numerics (NULL where the fit failed).
extract_coefs <- function(results, method_label) {
  results |>
    mutate(
      mu_hat     = map_dbl(coef, \(cf) if (is.null(cf)) NA_real_ else cf[1]),
      sigma2_hat = map_dbl(coef, \(cf) if (is.null(cf)) NA_real_ else cf[2]),
      method     = method_label
    ) |>
    select(pct, delta, sim, mu_hat, sigma2_hat, loglike, conv, method)
}

all_coefs <- bind_rows(
  extract_coefs(gel_by_method$gel_bias, "gel_bias"),
  extract_coefs(gel_by_method$gel,      "gel")
) |>
  mutate(method = factor(method, levels = names(method_labels)))

as_choices <- function(x) {
  vals <- sort(unique(x))
  set_names(as.character(vals), format(vals, trim = TRUE))
}

pct_choices   <- as_choices(all_coefs$pct)      # character ("01", "05", "10")
delta_choices <- as_choices(all_coefs$delta)
sim_choices   <- as_choices(all_coefs$sim)
delta_numeric <- as.numeric(delta_choices)
delta_step    <- if (length(delta_numeric) > 1) min(diff(delta_numeric)) else 1

# ============================== PLOT HELPERS =================================

app_theme <- function(base_size = 14) {
  theme_bw(base_size = base_size) +
    theme(
      axis.text   = element_text(family = "serif", colour = "black"),
      axis.title  = element_text(colour = "black"),
      axis.ticks  = element_line(colour = "black"),
      plot.title  = element_text(face = "bold"),
      legend.position = "bottom"
    )
}

scale_outlier_colour <- function() {
  scale_colour_manual(values = outlier_cols, labels = outlier_labs, name = NULL, drop = FALSE)
}

get_gel_row <- function(method, pct_value, delta_value, sim_value) {
  gel_by_method[[method]] |>
    filter(pct == pct_value, delta == delta_value, sim == sim_value)
}

# ============================== UI ===========================================

ui <- page_navbar(
  id    = "nav",
  title = tagList(icon("chart-line"), "Log-normal Outlier Simulation Explorer"),
  theme = bs_theme(
    version      = 5,
    bootswatch   = "flatly",
    base_font    = font_google("Inter"),
    heading_font = font_google("Inter", wght = 600),
    primary      = "#2C3E50"
  ),
  fillable    = TRUE,
  collapsible = TRUE,

  sidebar = sidebar(
    width = 320,
    selectInput("pct", "Outlier %:", choices = pct_choices),

    conditionalPanel(
      "input.nav != 'coef'",
      selectInput("delta", "delta (log-scale shift):", choices = delta_choices),
      selectizeInput("sim", "Simulation replicate:", choices = NULL)
    ),

    conditionalPanel(
      "input.nav == 'coef'",
      radioButtons(
        "coef_param", "Parameter:",
        choices = c("\u03bc" = "mu", "\u03c3\u00b2" = "sigma2")
      ),
      checkboxGroupInput(
        "coef_methods", "Methods to show:",
        choices  = set_names(names(method_labels), method_labels),
        selected = names(method_labels)
      ),
      sliderInput(
        "coef_delta_range", "delta range:",
        min   = min(delta_numeric), max = max(delta_numeric),
        value = range(delta_numeric), step = delta_step
      )
    ),

    conditionalPanel(
      "input.nav == 'gtl' || input.nav == 'pt'",
      selectInput(
        "gel_method", "ETEL fit:",
        choices = set_names(gel_methods, method_labels[gel_methods])
      )
    ),

    conditionalPanel(
      "input.nav == 'pt'",
      checkboxInput("pt_log", "Log scale (log \u03c0\u1d62)", value = TRUE)
    ),

    conditionalPanel(
      "input.nav == 'll'",
      sliderInput(
        "ll_window", "Contour window (below max ll):",
        min = 50, max = 5000, value = 400, step = 50
      )
    ),

    hr(),
    uiOutput("tab_help")
  ),

  nav_panel(
    title = tagList(icon("table-cells"), "Original Data"), value = "data",
    card(full_screen = TRUE, plotOutput("data_scatter", height = "100%"))
  ),
  nav_panel(
    title = tagList(icon("chart-area"), "Coefficient densities"), value = "coef",
    card(full_screen = TRUE, plotOutput("coef_density", height = "100%"))
  ),
  nav_panel(
    title = tagList(icon("timeline"), "\u03bb\u1d40 g\u1d62"), value = "gtl",
    card(full_screen = TRUE, plotOutput("gtl_plot", height = "100%"))
  ),
  nav_panel(
    title = tagList(icon("scale-balanced"), "Implied Probabilities"), value = "pt",
    card(full_screen = TRUE, plotOutput("pt_plot", height = "100%"))
  ),
  nav_panel(
    title = tagList(icon("mountain-sun"), "LL Surface"), value = "ll",
    card(full_screen = TRUE, plotOutput("ll_surface_plot", height = "100%"))
  )
)

# ============================== SERVER =======================================

server <- function(input, output, session) {

  updateSelectizeInput(session, "sim", choices = sim_choices,
                       selected = sim_choices[1], server = TRUE)

  # pct stays character ("01"); delta and sim are numeric
  sel <- reactive({
    req(input$pct, input$delta, input$sim)
    list(
      pct   = input$pct,
      delta = as.numeric(input$delta),
      sim   = as.numeric(input$sim)
    )
  })

  caption <- reactive({
    s <- sel()
    str_glue("pct = {s$pct} \u00b7 delta = {s$delta} \u00b7 sim = {s$sim}")
  })

  output$tab_help <- renderUI({
    helpText(switch(
      input$nav %||% "data",
      data = "Simulated y for one replicate. Red points are the contaminated observations. Dashed line = true mean exp(\u03bc + \u03c3\u00b2/2).",
      coef = "Ridge densities of the ETEL estimates across replicates, one ridge per delta.",
      gtl  = "Inner products \u03bb\u1d40g\u1d62 per observation at the chosen ETEL estimate.",
      pt   = "ETEL implied probabilities. Dashed line = uniform weight 1/N. Collapsed fits are only visible on the log scale.",
      ll   = "Red = true (\u03bc, \u03c3\u00b2) (nearest grid node). Purple = grid maximum. Grey = infeasible (outside the convex hull)."
    ))
  })

  # ---- Tab: Original Data ---------------------------------------------------
  selected_replicate <- reactive({
    s <- sel()
    all_data |> filter(pct == s$pct, delta == s$delta, sim == s$sim)
  })

  output$data_scatter <- renderPlot({
    df <- selected_replicate()
    validate(need(nrow(df) > 0, "No data found for this combination."))

    ggplot(df, aes(x = index, y = y, colour = is_out)) +
      geom_point(alpha = 0.7, size = 2) +
      geom_hline(
        yintercept = exp(theta_true[["mu"]] + theta_true[["sigma2"]] / 2),
        linetype = "dashed", colour = "black"
      ) +
      scale_outlier_colour() +
      labs(
        title    = str_glue("Simulated data \u2014 {caption()}"),
        subtitle = "Dashed line = true mean",
        x = "Index", y = "y"
      ) +
      app_theme()
  }) |>
    bindCache(sel())

  # ---- Tab: Coefficient densities -------------------------------------------
  output$coef_density <- renderPlot({
    req(input$coef_param)
    validate(need(length(input$coef_methods) > 0, "Select at least one method."))

    param    <- input$coef_param
    x_var    <- if (param == "sigma2") "sigma2_hat" else "mu_hat"
    true_val <- theta_true[[param]]
    sym      <- if (param == "sigma2") "\\sigma^2" else "\\mu"

    plot_data <- all_coefs |>
      filter(
        pct    == input$pct,
        method %in% input$coef_methods,
        between(delta, input$coef_delta_range[1], input$coef_delta_range[2]),
        !is.na(.data[[x_var]])
      ) |>
      mutate(delta_f = fct_inseq(factor(delta)))

    validate(need(nrow(plot_data) > 0,
                  "No data for this selection \u2014 widen the range or pick a method."))

    ggplot(plot_data, aes(x = .data[[x_var]], y = delta_f, fill = method, colour = method)) +
      geom_density_ridges(alpha = 0.35, scale = 1.5, linewidth = 0.8) +
      geom_vline(xintercept = true_val, linetype = "dashed", colour = "black") +
      scale_fill_brewer(palette = "Set2", labels = method_labels, drop = FALSE) +
      scale_colour_brewer(palette = "Set2", labels = method_labels, drop = FALSE) +
      labs(
        title    = str_glue("Estimates of {if (param == 'sigma2') '\u03c3\u00b2' else '\u03bc'}, outliers = {input$pct}%"),
        subtitle = TeX(str_glue("Dashed line = true value (${sym} = {true_val}$)")),
        x        = TeX(str_glue("$\\hat{{{sym}}}$")),
        y        = TeX("$\\delta$"),
        fill = "Method", colour = "Method"
      ) +
      app_theme()
  })

  # ---- Shared per-observation tibble for the two diagnostic tabs -----------
  gel_obs <- reactive({
    req(input$gel_method)
    s   <- sel()
    row <- get_gel_row(input$gel_method, s$pct, s$delta, s$sim)
    validate(need(nrow(row) == 1, "No matching ETEL fit found for this combination."))
    validate(need(!is.null(row$fit_ETEL[[1]]), "This fit failed."))

    row$data[[1]] |>
      select(index, is_out, y) |>
      left_join(row$fit_ETEL[[1]]$df, by = "index")
  })

  fit_caption <- reactive({
    s   <- sel()
    row <- get_gel_row(input$gel_method, s$pct, s$delta, s$sim)
    cf  <- row$coef[[1]]
    str_glue("{method_labels[[input$gel_method]]}: ",
             "\u03bc\u0302 = {round(cf[1], 3)}, \u03c3\u0302\u00b2 = {round(cf[2], 3)}, ",
             "ll = {round(row$loglike, 1)}")
  })

  output$gtl_plot <- renderPlot({
    ggplot(gel_obs(), aes(x = index, y = inner_pro, colour = is_out)) +
      geom_point(alpha = 0.7, size = 2) +
      geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
      scale_outlier_colour() +
      labs(
        title    = str_glue("\u03bb\u1d40g\u1d62 \u2014 {caption()}"),
        subtitle = fit_caption(),
        x = "Observation index", y = TeX("$\\lambda^T g_i$")
      ) +
      app_theme()
  })

  output$pt_plot <- renderPlot({
    df <- gel_obs()
    n  <- nrow(df)

    if (isTRUE(input$pt_log)) {
      p <- ggplot(df, aes(x = index, y = logpi, colour = is_out)) +
        geom_hline(yintercept = log(1 / n), linetype = "dashed", colour = "grey40") +
        labs(y = TeX("$\\log \\pi_i$"), subtitle = str_glue("{fit_caption()} \u00b7 dashed = log(1/N)"))
    } else {
      p <- ggplot(df, aes(x = index, y = exp(logpi), colour = is_out)) +
        geom_hline(yintercept = 1 / n, linetype = "dashed", colour = "grey40") +
        labs(y = TeX("$\\pi_i$"), subtitle = str_glue("{fit_caption()} \u00b7 dashed = 1/N"))
    }

    p +
      geom_point(alpha = 0.7, size = 2) +
      scale_outlier_colour() +
      labs(title = str_glue("Implied probabilities \u2014 {caption()}"),
           x = "Observation index") +
      app_theme()
  })

  # ---- Tab: Log-likelihood surface -----------------------------------------
  output$ll_surface_plot <- renderPlot({
    s <- sel()

    ll_res <- ll_surface |>
      filter(pct == s$pct, delta == s$delta, sim == s$sim)

    validate(need(nrow(ll_res) > 0, "No log-likelihood surface found for this combination."))
    validate(need(any(!is.na(ll_res$ll)), "No feasible grid cells for this dataset."))

    pt_true <- ll_res |>
      filter(!is.na(ll)) |>
      slice_min((mu - theta_true[["mu"]])^2 + (sigma2 - theta_true[["sigma2"]])^2,
                n = 1, with_ties = FALSE)

    mx     <- ll_res |> slice_max(ll, n = 1, with_ties = FALSE)
    window <- input$ll_window
    marks  <- bind_rows(true = pt_true, max = mx, .id = "which")

    # Infeasible cells: shaded and set to the floor so contours close cleanly
    ll_res |>
      mutate(ll_clipped = pmax(coalesce(ll, mx$ll - window), mx$ll - window)) |>
      ggplot(aes(x = mu, y = sigma2)) +
      geom_tile(data = \(d) filter(d, !hull), fill = "grey85") +
      geom_contour(aes(z = ll_clipped), colour = "black", linewidth = 0.7,
                   breaks = MakeBreaks(binwidth = window / 10)) +
      geom_text_contour(aes(z = ll_clipped), size = 3, colour = "black", stroke = 0.2) +
      geom_point(data = marks, aes(colour = which), size = 4) +
      geom_text(data = marks, aes(colour = which, label = round(ll, 1)),
                size = 3.5, vjust = -1, show.legend = FALSE) +
      scale_colour_manual(
        values = c(true = "red", max = "purple"),
        labels = c(true = "True (\u03bc, \u03c3\u00b2) (nearest node)", max = "Grid maximum"),
        name   = NULL
      ) +
      labs(
        title = str_glue("ETEL log-likelihood \u2014 {caption()}"),
        x = TeX("$\\mu$"), y = TeX("$\\sigma^2$")
      ) +
      app_theme()
  }) |>
    bindCache(sel(), input$ll_window)
}

shinyApp(ui, server)