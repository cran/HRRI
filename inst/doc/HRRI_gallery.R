## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  collapse   = TRUE,
  comment    = "#>",
  fig.width  = 7,
  fig.height = 4.4,
  fig.align  = "center",
  dpi        = 150,
  out.width  = "100%",
  message    = FALSE,
  warning    = FALSE
)

if (requireNamespace("ggplot2", quietly = TRUE)) {
  old_theme <- ggplot2::theme_set(
    ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(
        panel.grid.minor = ggplot2::element_blank(),
        panel.grid.major = ggplot2::element_line(linewidth = 0.3,
                                                 colour = "#dde3e1"),
        axis.title       = ggplot2::element_text(colour = "#4a5451"),
        axis.text        = ggplot2::element_text(colour = "#4a5451"),
        strip.text       = ggplot2::element_text(face = "bold",
                                                 colour = "#1c2321"),
        plot.title       = ggplot2::element_text(face = "bold",
                                                 colour = "#1c2321"),
        legend.position  = "bottom"
      )
  )
}

## ----libraries----------------------------------------------------------------
library(HRRI)
library(ggplot2)
packageVersion("HRRI")

## ----simulate-----------------------------------------------------------------
PERTURB_START <- 12
PERTURB_END   <- 22

sim <- simulate_redox_holobiont(
  n_plot               = 2,
  n_depth              = 2,
  n_plant              = 3,
  n_time               = 40,
  p_micro              = 25,
  seed                 = 2026,
  scenario             = "flood_drain",
  n_cycles             = 1,
  disturbance_strength = 0.72,
  history_strength     = 0.55
)

nrow(sim$id)          # 2 x 2 x 3 x 40 = 480 observations
names(sim$latent_state)

## ----hidden-states, fig.height=5.2--------------------------------------------
ls_df <- as.data.frame(sim$latent_state)
one   <- with(sim$id, plot == "P1" & depth == "D1" & plant_id == "Plant1")

hidden <- data.frame(
  time  = sim$id$time[one],
  value = c(
    ls_df$Q_accept[one],
    ls_df$alpha_accept[one],
    ls_df$k_accept_h[one],
    ls_df$memory[one]
  ),
  state = factor(
    rep(c("Capacity  (Q)", "Connectivity  (α)",
          "Kinetics  (k)", "Memory  (M)"), each = sum(one)),
    levels = c("Capacity  (Q)", "Connectivity  (α)",
               "Kinetics  (k)", "Memory  (M)")
  )
)

ggplot(hidden, aes(time, value)) +
  annotate("rect", xmin = PERTURB_START, xmax = PERTURB_END,
           ymin = -Inf, ymax = Inf, fill = "#2f6b6b", alpha = 0.10) +
  geom_line(colour = "#2f6b6b", linewidth = 0.7) +
  facet_wrap(~ state, scales = "free_y", ncol = 2) +
  labs(x = "Time (days)", y = NULL,
       title = "The four hidden states during one flood-drain cycle",
       subtitle = "Shaded band = disturbance window") +
  theme(legend.position = "none")

## ----memory-components, fig.height=3.6----------------------------------------
mem <- data.frame(
  time  = rep(sim$id$time[one], 3),
  value = c(ls_df$memory[one],
            ls_df$plant_legacy[one],
            ls_df$micro_legacy[one]),
  component = factor(
    rep(c("Memory (total)", "Plant legacy (aerenchyma)",
          "Microbial legacy (community)"), each = sum(one)),
    levels = c("Memory (total)", "Plant legacy (aerenchyma)",
               "Microbial legacy (community)")
  )
)

ggplot(mem, aes(time, value, colour = component, linetype = component)) +
  annotate("rect", xmin = PERTURB_START, xmax = PERTURB_END,
           ymin = -Inf, ymax = Inf, fill = "#9a6a24", alpha = 0.10) +
  geom_line(linewidth = 0.7) +
  scale_colour_manual(values = c("#1c2321", "#8c4a2f", "#6b5b8a")) +
  scale_linetype_manual(values = c("solid", "dashed", "dotdash")) +
  labs(x = "Time (days)", y = "Legacy (0-1)", colour = NULL, linetype = NULL,
       title = "Memory decomposed across the holobiont")

## ----pipeline-----------------------------------------------------------------
res <- rri_pipeline_st(
  ROS_flux     = sim$ROS_flux,
  Eh_stability = sim$Eh_stability,
  micro_data   = sim$micro_data,
  id           = sim$id,
  time_col     = "time",
  group_cols   = c("plot", "depth", "plant_id"),
  mode         = "snapshot",
  reducer      = "per_domain",
  scaling      = "pnorm",
  direction_anchor_phys  = "FvFm",
  direction_anchor_soil  = "Eh",
  direction_anchor_micro = "ASV1"
)

scored <- attach_hrri_ids(res$row_scores, sim$id)
attr(scored, "id_alignment")
summary(scored$RRI)

## ----timeseries---------------------------------------------------------------
plot_rri_timeseries(
  sim, res,
  plot_id       = "P1",
  depth_id      = "D1",
  plant_id      = "Plant1",
  perturb_start = PERTURB_START,
  perturb_end   = PERTURB_END
)

## ----ternary, echo=TRUE, results="asis"---------------------------------------
## ggtern is a Suggests dependency. Loading it -- not drawing with it --
## patches ggplot2's element tree, and under ggplot2 >= 4.0.0 that patch makes
## every later ggplot in the session fail with
##   "The `tern.axis.ticks.length.major` theme element must be a <rel> object."
## Vignettes are built in one R session, so a requireNamespace() here would
## take the workflow vignette down with it. The ggplot2 version is therefore
## checked before ggtern is touched at all; try() alone is too late.
ggplot2_ok <- utils::packageVersion("ggplot2") < "4.0.0"
tern_ok <- ggplot2_ok &&
           requireNamespace("ggtern", quietly = TRUE) &&
           requireNamespace("viridis", quietly = TRUE)

if (tern_ok) {
  p_tern <- try(
    plot_RRI_ternary(res$row_scores_comp, point_size = 2.4,
                     show_centroid = TRUE),
    silent = TRUE
  )
  drawn <- !inherits(p_tern, "try-error") &&
           !inherits(try(print(p_tern), silent = TRUE), "try-error")
  if (!drawn) {
    cat("*The ternary plot could not be rendered: the installed **ggtern** is",
        "incompatible with this **ggplot2** version. The composition it would",
        "show is summarised numerically below.*\n\n")
  }
} else {
  cat("*The ternary panel is skipped here: **ggtern** is either not installed",
      "or not compatible with the installed **ggplot2**",
      sprintf("(%s).", utils::packageVersion("ggplot2")),
      "It is deliberately not loaded in that case, because loading it would",
      "break the remaining figures. The same composition is given numerically",
      "below.*\n\n")
}

## ----ternary-numeric----------------------------------------------------------
comp <- res$row_scores_comp[, c("Physio", "Soil", "Micro")]
round(colMeans(comp, na.rm = TRUE), 3)          # centroid
round(range(rowSums(comp, na.rm = TRUE)), 6)    # closure check: both 1

## ----state-space--------------------------------------------------------------
plot_rri_state_space(
  res,
  x_property = "Physio",
  y_property = "Soil",
  colour_by  = "RRI",
  group_cols = c("plot", "depth", "plant_id")
)

## ----recovery-----------------------------------------------------------------
rec <- rri_recovery_metrics(
  res           = res,
  id            = sim$id,
  time_col      = "time",
  group_cols    = c("plot", "depth", "plant_id"),
  perturb_start = PERTURB_START,
  perturb_end   = PERTURB_END,
  rri_col       = "RRI"
)

rec[1:4, c("plot", "depth", "plant_id", "baseline_rri", "depth_min_frac",
           "tau_lag", "overshoot_frac", "incomplete_return_frac",
           "displaced_plateau_flag", "fit_status")]

## ----recovery-map, fig.height=4.8---------------------------------------------
plot_rri_recovery_map(
  res           = res,
  id            = sim$id,
  rec           = rec,
  time_col      = "time",
  group_cols    = c("plot", "depth", "plant_id"),
  perturb_start = PERTURB_START,
  perturb_end   = PERTURB_END
)

## ----landscape, eval=requireNamespace("tidyr", quietly=TRUE) && requireNamespace("tidyselect", quietly=TRUE), fig.height=5----
## Name the metrics explicitly rather than relying on the function default.
## Older HRRI builds defaulted to A_norm / O_norm / tau_r, which
## rri_recovery_metrics() no longer produces; being explicit makes this chunk
## work against either version and documents which signatures are shown.
plot_rri_recovery_landscape(
  rec,
  metrics = intersect(
    c("depth_min_frac", "overshoot_frac", "I_norm", "k", "tau_lag", "t_half"),
    names(rec)
  ),
  order_by = "I_norm"
)

## ----properties, fig.height=5-------------------------------------------------
## soil_df is what makes Capacity available. Without it the Capacity axis is
## returned as NA and the radar shows a short spoke.
props <- rri_property_scores(res, rec = rec, soil_df = sim$soil_data)
props$property_table

plot_rri_properties(props, rri_value = mean(scored$RRI, na.rm = TRUE))

## ----validation, fig.height=4-------------------------------------------------
ok <- is.finite(scored$RRI) & is.finite(sim$latent_truth)
r  <- stats::cor(scored$RRI[ok], sim$latent_truth[ok])
cat("r(RRI, latent_truth) =", round(r, 3), "on", sum(ok), "observations\n")

ggplot(data.frame(truth = sim$latent_truth[ok], RRI = scored$RRI[ok]),
       aes(truth, RRI)) +
  geom_point(alpha = 0.25, size = 1.1, colour = "#2f6b6b") +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE,
              colour = "#8c4a2f", fill = "#8c4a2f", alpha = 0.12) +
  labs(x = "Simulator target  z(t)", y = "HRRI score",
       title = sprintf("Agreement with the prescribed target (r = %.3f)", r))

## ----own-data, eval=FALSE-----------------------------------------------------
# res <- rri_pipeline(
#   soil  = my_soil,      # Eh, pH, Fe pools, EAC/EDC ...
#   plant = my_plant,     # SPAD, Fv/Fm, ROL ...
#   micro = my_micro,     # ASV table or functional genes
#   id    = my_ids,       # plot, depth, plant_id, time
#   direction_anchor_soil = "Eh",
#   direction_anchor_phys = "FvFm"
# )

## ----restore-theme, include=FALSE---------------------------------------------
## theme_set() changes state that persists for the rest of the session.
## Vignettes build in their own process so nothing outside is affected, but
## restoring is the same courtesy CRAN asks for with par() and options().
if (exists("old_theme")) ggplot2::theme_set(old_theme)

## ----session------------------------------------------------------------------
sessionInfo()

