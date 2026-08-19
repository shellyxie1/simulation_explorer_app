library(shiny)
library(bslib)
library(thematic)
library(tidyverse)
library(ggridges)
library(latex2exp)
library(metR)

# Auto-style ggplot output (fonts/colors) to match the Bootstrap theme below.
# NOTE: if axis text/title colors don't come out black despite app_theme()
# setting them explicitly, thematic is likely overriding them - uncomment
# the fg = "black" argument below to pin foreground (text/line) color.
thematic_shiny(font = "auto")
# thematic_shiny(font = "auto", fg = "black")

# ============================== LOAD DATA ====================================
# All paths are relative to wherever this app is launched from — adjust if needed.

all_data          <- read_rds("outputs/02_all_data_sims.rds")
ols_results_wide  <- read_rds("outputs/02_ols_results_wide.rds")
# ols_results_long <- read_rds("outputs_lite/02_ols_results_long.rds")  # not currently used by any tab

# Coefficient-only GEL fits (used for the b0/b1 density tabs) — lightweight
# versions with beta0_hat/beta1_hat already extracted, fit_ETEL dropped.
gel_bias_coef <- read_rds("outputs/05_gel_coefs_bias.rds")
gel_coef      <- read_rds("outputs/05_gel_coefs.rds")

# Same fits' data, gt_lambda, and implied_prob (used for those two tabs) —
# lightweight versions with fit_ETEL dropped.
gel_bias_full <- read_rds("outputs/05_gel_full_bias.rds")
gel_full      <- read_rds("outputs/05_gel_full.rds")

# Log-likelihood surface over a (th1, th2) grid, per (pct, out_eps_mean, sim)
ll_surface <- read_rds("outputs/06_ll_surface.rds")

# NOTE: is_outlier is a column inside all_data itself (random per-replicate
# outlier assignment), so no separate outlier index lookup is needed.

beta_true <- c(2, -3)  # (Intercept), x

# Shared plot theme, applied to every plot in the app (matches the
# log-likelihood surface style: white background, black serif axis text).
app_theme <- function(base_size = 14) {
  theme_bw(base_size = base_size) +
    theme(
      axis.text   = element_text(family = "serif", colour = "black"),
      axis.title  = element_text(colour = "black"),
      axis.ticks  = element_line(colour = "black"),
      plot.title  = element_text(face = "bold")
    )
}

extract_gel_coefs <- function(results, method_label) {
  results %>%
    mutate(method = method_label) %>%
    select(pct, out_eps_mean, sim, beta0_hat, beta1_hat, method)
}

ols_coefs <- ols_results_wide %>%
  select(pct, out_eps_mean, sim, beta0_hat, beta1_hat) %>%
  mutate(method = "OLS")

all_coefs <- bind_rows(
  ols_coefs,
  extract_gel_coefs(gel_bias_coef, "gel_bias"),
  extract_gel_coefs(gel_coef, "gel")
) %>%
  mutate(method = factor(method, levels = c("OLS", "gel_bias", "gel")))

pct_choices        <- sort(unique(all_coefs$pct))
method_choices      <- levels(all_coefs$method)
eps_choices          <- sort(unique(all_coefs$out_eps_mean))
sim_choices          <- sort(unique(all_data$sim %||% ols_results_wide$sim))

gel_method_choices <- c("gel_bias", "gel")

get_gel_row <- function(method, pct_value, eps_value, sim_value) {
  source_df <- if (method == "gel_bias") gel_bias_full else gel_full
  
  source_df %>%
    filter(
      pct == pct_value,
      out_eps_mean == eps_value,
      sim == sim_value
    )
}

# ============================== UI ============================================

ui <- navbarPage(
  title = tagList(icon("chart-line"), "Outlier Simulation Explorer"),
  theme = bs_theme(
    version     = 5,
    bootswatch  = "flatly",
    primary     = "#2C3E50"
  ),
  collapsible = TRUE,
  header = tags$head(
    tags$style(HTML("
      .navbar-nav .nav-item { margin-right: 18px; }
      .navbar-nav .nav-link { padding-left: 4px; padding-right: 4px; }
    "))
  ),
  
  # ---- Tab 1: Original Data --------------------------------------------------
  tabPanel(
    tagList(icon("table-cells"), "Original Data"),
    sidebarLayout(
      sidebarPanel(
        selectInput("data_pct", "Outlier %:", choices = pct_choices, selected = pct_choices[1]),
        selectInput("data_eps", "out_eps_mean:", choices = eps_choices, selected = 0),
        selectInput("data_sim", "Simulation replicate:", choices = sim_choices, selected = sim_choices[1]),
        hr(),
        helpText("Shows the raw simulated (x, y) data for one replicate. Red points are the injected outliers.")
      ),
      mainPanel(
        card(
          full_screen = TRUE,
          card_body(plotOutput("data_scatter", height = "480px"))
        )
      )
    )
  ),
  
  # ---- Tab 2: Intercept (b0) -------------------------------------------------
  tabPanel(
    tagList(icon("chart-area"), "Intercept (\u03b2\u2080)"),
    sidebarLayout(
      sidebarPanel(
        selectInput("b0_pct", "Outlier %:", choices = pct_choices, selected = pct_choices[1]),
        checkboxGroupInput("b0_methods", "Methods to show:", choices = method_choices, selected = method_choices),
        sliderInput(
          "b0_eps_range", "out_eps_mean range:",
          min = min(eps_choices), max = max(eps_choices),
          value = c(min(eps_choices), max(eps_choices)), step = 0.5
        )
      ),
      mainPanel(
        card(
          full_screen = TRUE,
          card_body(plotOutput("b0_density", height = "680px"))
        )
      )
    )
  ),
  
  # ---- Tab 3: Slope (b1) ------------------------------------------------------
  tabPanel(
    tagList(icon("chart-area"), "Slope (\u03b2\u2081)"),
    sidebarLayout(
      sidebarPanel(
        selectInput("b1_pct", "Outlier %:", choices = pct_choices, selected = pct_choices[1]),
        checkboxGroupInput("b1_methods", "Methods to show:", choices = method_choices, selected = method_choices),
        sliderInput(
          "b1_eps_range", "out_eps_mean range:",
          min = min(eps_choices), max = max(eps_choices),
          value = c(min(eps_choices), max(eps_choices)), step = 0.5
        )
      ),
      mainPanel(
        card(
          full_screen = TRUE,
          card_body(plotOutput("b1_density", height = "680px"))
        )
      )
    )
  ),
  
  # ---- Tab 4: gt %*% lambda ----------------------------------------------------
  tabPanel(
    tagList(icon("timeline"), "gt \u00d7 lambda"),
    sidebarLayout(
      sidebarPanel(
        selectInput("gtl_method", "Method:", choices = gel_method_choices, selected = "gel_bias"),
        selectInput("gtl_pct", "Outlier %:", choices = pct_choices, selected = pct_choices[1]),
        selectInput("gtl_eps", "out_eps_mean:", choices = eps_choices, selected = 0),
        selectInput("gtl_sim", "Simulation replicate:", choices = sim_choices, selected = sim_choices[1]),
        hr(),
        helpText("gt %*% lambda per observation, for the chosen GEL fit. Red points are the injected outlier observations.")
      ),
      mainPanel(
        card(
          full_screen = TRUE,
          card_body(plotOutput("gtl_plot", height = "480px"))
        )
      )
    )
  ),
  
  # ---- Tab 5: Implied probabilities --------------------------------------------
  tabPanel(
    tagList(icon("scale-balanced"), "Implied Probabilities"),
    sidebarLayout(
      sidebarPanel(
        selectInput("pt_method", "Method:", choices = gel_method_choices, selected = "gel_bias"),
        selectInput("pt_pct", "Outlier %:", choices = pct_choices, selected = pct_choices[1]),
        selectInput("pt_eps", "out_eps_mean:", choices = eps_choices, selected = 0),
        selectInput("pt_sim", "Simulation replicate:", choices = sim_choices, selected = sim_choices[1]),
        hr(),
        helpText("GEL implied probabilities (pt) per observation. Dashed line = uniform weight 1/N. Red points are the injected outlier observations.")
      ),
      mainPanel(
        card(
          full_screen = TRUE,
          card_body(plotOutput("pt_plot", height = "480px"))
        )
      )
    )
  ),
  
  # ---- Tab 6: Log-likelihood surface -------------------------------------------
  tabPanel(
    tagList(icon("mountain-sun"), "LL Surface"),
    sidebarLayout(
      sidebarPanel(
        selectInput("ll_pct", "Outlier %:", choices = pct_choices, selected = pct_choices[1]),
        selectInput("ll_eps", "out_eps_mean:", choices = eps_choices, selected = 0),
        selectInput("ll_sim", "Simulation replicate:", choices = sim_choices, selected = sim_choices[1]),
        sliderInput("ll_window", "Contour window (below max ll):", min = 20, max = 800, value = 340, step = 20),
        hr(),
        helpText("Red point = true (\u03b2\u2080, \u03b2\u2081). Purple point = grid maximum log-likelihood. Contour values below (max \u2212 window) are clipped for readability.")
      ),
      mainPanel(
        card(
          full_screen = TRUE,
          card_body(plotOutput("ll_surface_plot", height = "550px"))
        )
      )
    )
  )
)

# ============================== SERVER =========================================

server <- function(input, output, session) {
  
  # ---- Tab 1: Original Data ---------------------------------------------------
  selected_replicate <- reactive({
    all_data %>%
      filter(
        pct == input$data_pct,
        out_eps_mean == as.numeric(input$data_eps),
        sim == as.numeric(input$data_sim)
      ) %>%
      mutate(row_id = row_number())
  })
  
  output$data_scatter <- renderPlot({
    df <- selected_replicate()
    validate(need(nrow(df) > 0, "No data found for this combination."))
    
    ggplot(df, aes(x = x, y = y, color = is_outlier)) +
      geom_point(alpha = 0.7, size = 3) +
      geom_abline(intercept = beta_true[1], slope = beta_true[2], linetype = "dashed", color = "black") +
      scale_color_manual(
        values = c(`FALSE` = "steelblue", `TRUE` = "red"),
        labels = c(`FALSE` = "Normal", `TRUE` = "Outlier"),
        name   = NULL
      ) +
      labs(
        title = paste0(
          "Simulated data — pct = ", input$data_pct,
          ", out_eps_mean = ", input$data_eps,
          ", sim = ", input$data_sim
        ),
        subtitle = "Dashed line = true regression line · red = injected outlier",
        x = "x", y = "y"
      ) +
      app_theme()
  })
  
  # ---- Shared plotting helper --------------------------------------------------
  method_colors <- c(OLS = "#2C3E50", gel_bias = "#E67E22", gel = "#8E44AD")
  
  make_density_plot <- function(pct_value, methods_selected, eps_range, param) {
    if (param == "b1") {
      x_var    <- "beta1_hat"
      true_val <- beta_true[2]
      x_lab    <- TeX("$\\hat{\\beta}_1$")
      subtitle <- TeX("Dashed line = true slope ($\\beta_1 = -3$)")
      title    <- paste0("Slope estimates, outliers = ", pct_value, "%")
    } else {
      x_var    <- "beta0_hat"
      true_val <- beta_true[1]
      x_lab    <- TeX("$\\hat{\\beta}_0$")
      subtitle <- TeX("Dashed line = true intercept ($\\beta_0 = 2$)")
      title    <- paste0("Intercept estimates, outliers = ", pct_value, "%")
    }
    
    plot_data <- all_coefs %>%
      filter(
        pct == pct_value,
        method %in% methods_selected,
        out_eps_mean >= eps_range[1],
        out_eps_mean <= eps_range[2]
      )
    
    validate(need(nrow(plot_data) > 0, "No data for this selection — widen the range or pick a method."))
    
    ggplot(plot_data, aes(x = .data[[x_var]], y = factor(out_eps_mean), fill = method, color = method)) +
      geom_density_ridges(alpha = 0.35, scale = 1.5, linewidth = 0.8) +
      geom_vline(xintercept = true_val, linetype = "dashed", color = "black") +
      scale_fill_manual(values = method_colors, drop = TRUE) +
      scale_color_manual(values = method_colors, drop = TRUE) +
      labs(
        title = title,
        subtitle = subtitle,
        x = x_lab,
        y = TeX("$\\bar{\\epsilon}_{out}$"),
        fill = "Method",
        color = "Method"
      ) +
      app_theme()
  }
  
  # ---- Tab 2: Intercept ---------------------------------------------------------
  output$b0_density <- renderPlot({
    validate(need(length(input$b0_methods) > 0, "Select at least one method."))
    make_density_plot(input$b0_pct, input$b0_methods, input$b0_eps_range, "b0")
  })
  
  # ---- Tab 3: Slope ---------------------------------------------------------------
  output$b1_density <- renderPlot({
    validate(need(length(input$b1_methods) > 0, "Select at least one method."))
    make_density_plot(input$b1_pct, input$b1_methods, input$b1_eps_range, "b1")
  })
  
  # ---- Shared helper: build a per-observation tibble for one GEL fit -----------
  build_gel_obs_df <- function(method, pct_value, eps_value, sim_value) {
    row <- get_gel_row(method, pct_value, as.numeric(eps_value), as.numeric(sim_value))
    validate(need(nrow(row) == 1, "No matching GEL fit found for this combination."))
    
    # Base x/y/is_outlier comes from all_data (same replicate, same row
    # order as when gt_lambda/implied_prob were computed — nest() preserves
    # within-group row order, so this lines up correctly).
    df <- all_data %>%
      filter(
        pct == pct_value,
        out_eps_mean == as.numeric(eps_value),
        sim == as.numeric(sim_value)
      )
    
    gt_lambda_vec    <- as.numeric(row$gt_lambda[[1]])
    implied_prob_vec <- as.numeric(row$implied_prob[[1]])
    
    df %>%
      mutate(
        row_id       = row_number(),
        gt_lambda    = gt_lambda_vec,
        implied_prob = implied_prob_vec
      )
  }
  
  # ---- Tab 4: gt x lambda -----------------------------------------------------
  output$gtl_plot <- renderPlot({
    df <- build_gel_obs_df(input$gtl_method, input$gtl_pct, input$gtl_eps, input$gtl_sim)
    
    ggplot(df, aes(x = row_id, y = gt_lambda, color = is_outlier)) +
      geom_point(alpha = 0.7, size = 3) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
      scale_color_manual(
        values = c(`FALSE` = "steelblue", `TRUE` = "red"),
        labels = c(`FALSE` = "Normal", `TRUE` = "Outlier"),
        name   = NULL
      ) +
      labs(
        title = paste0(
          "gt \u00d7 lambda — method = ", input$gtl_method,
          ", pct = ", input$gtl_pct,
          ", out_eps_mean = ", input$gtl_eps,
          ", sim = ", input$gtl_sim
        ),
        x = "Observation index", y = "gt \u00d7 lambda"
      ) +
      app_theme()
  })
  
  # ---- Tab 5: Implied probabilities --------------------------------------------
  output$pt_plot <- renderPlot({
    df <- build_gel_obs_df(input$pt_method, input$pt_pct, input$pt_eps, input$pt_sim)
    n_obs <- nrow(df)
    
    ggplot(df, aes(x = row_id, y = implied_prob, color = is_outlier)) +
      geom_point(alpha = 0.7, size = 3) +
      geom_hline(yintercept = 1 / n_obs, linetype = "dashed", color = "grey40") +
      scale_color_manual(
        values = c(`FALSE` = "steelblue", `TRUE` = "red"),
        labels = c(`FALSE` = "Normal", `TRUE` = "Outlier"),
        name   = NULL
      ) +
      labs(
        title = paste0(
          "Implied probabilities — method = ", input$pt_method,
          ", pct = ", input$pt_pct,
          ", out_eps_mean = ", input$pt_eps,
          ", sim = ", input$pt_sim
        ),
        subtitle = "Dashed line = uniform weight (1/N)",
        x = "Observation index", y = "Implied probability"
      ) +
      app_theme()
  })
  
  # ---- Tab 6: Log-likelihood surface -------------------------------------------
  output$ll_surface_plot <- renderPlot({
    
    ll_res <- ll_surface %>%
      filter(
        pct == input$ll_pct,
        out_eps_mean == as.numeric(input$ll_eps),
        sim == as.numeric(input$ll_sim)
      )
    
    validate(need(nrow(ll_res) > 0, "No log-likelihood surface found for this combination."))
    
    pt_true <- ll_res %>%
      filter(abs(th1 - beta_true[1]) < 1e-9, abs(th2 - beta_true[2]) < 1e-9)
    
    mx      <- ll_res %>% slice_max(ll, n = 1, with_ties = FALSE)
    max_ll  <- mx$ll
    window  <- input$ll_window
    
    ll_res %>%
      mutate(ll_clipped = pmax(ll, max_ll - window)) %>%
      ggplot(aes(x = th1, y = th2, z = ll_clipped)) +
      geom_contour(colour = "black", linewidth = 0.7) +
      geom_text_contour(size = 3, colour = "black", family = "serif") +
      geom_point(data = pt_true, aes(x = th1, y = th2),
                 colour = "red", size = 4, inherit.aes = FALSE) +
      geom_text(data = pt_true, aes(x = th1, y = th2,
                                    label = paste0("log-lik=", round(ll, 1))),
                colour = "red", size = 3.5, vjust = -1, family = "serif", inherit.aes = FALSE) +
      geom_point(data = mx, aes(x = th1, y = th2),
                 colour = "purple", size = 4, inherit.aes = FALSE) +
      geom_text(data = mx, aes(x = th1, y = th2,
                               label = paste0("log-lik=", round(ll, 1))),
                colour = "purple", size = 3.5, vjust = -1, family = "serif", inherit.aes = FALSE) +
      labs(
        title = paste0(
          "ETEL log-likelihood — pct = ", input$ll_pct,
          ", out_eps_mean = ", input$ll_eps,
          ", sim = ", input$ll_sim
        ),
        x = NULL, y = NULL
      ) +
      app_theme()
  })
}

shinyApp(ui, server)

