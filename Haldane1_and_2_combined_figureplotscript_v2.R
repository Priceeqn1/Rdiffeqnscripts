# ==============================================================================
# Numerical validation panel for:
# "Homeostatic Turnover, Repair, and the Demographic Form of Haldane's Principle"
#
# Panels:
#   A. Multiclass rare-mutation error: E_L(U) = O(U^2)
#   B. Binary threshold dynamics:
#        epsilon * U_c * tau * (1 - x(tau)) / zeta_rho -> 1
#   C. Binary repair regularization at U = U_c:
#        (1 - x^*) / sqrt(zeta_rho * nu / (epsilon * U_c)) -> 1
#
# Additional numerical check:
#   Multiclass O(U^2) convergence when repair acts at every damaged class,
#   for nu_i = 0.1 and nu_i = 0.5, i >= 1.
#
# Required package: deSolve
# Output:
#   numerical_validation_panel.pdf
#   numerical_validation_panel.png
#   multiclass_error_scaling.csv
#   all_class_repair_error_scaling.csv
#   threshold_dynamic_scaling.csv
#   threshold_repair_scaling.csv
#   numerical_validation_summary.txt
#
# August 16 2026
# Last updated: 9/19/2026
# Change Log
# Date            Change
# 8/16/2026       changed placement (h, w) of panel plot C's on-plot labels
# 8/30/2026       extended decimal places of slope labeled on the plot
# 9/19/2026       added five-point all-class repair convergence check
# 9/19/2026       removed upper bound of epsilon in accored with manuscript
#
# ==============================================================================

# Load required packages
library(deSolve)
 
# ------------------------------------------------------------------------------
# 1. User settings
# ------------------------------------------------------------------------------

output_dir <- "numerical_validation"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

pars <- list(
  n = 10L,
  epsilon = 0.5,
  zeta_rho = 0.8,
  zeta_gamma = 1.3
)

# Panel A: only first-step repair is active, so the numerical experiment
# isolates the leading rare-mutation mechanism in Corollary 1.
nu1_panel_A <- 0.10

# Small-U grid for the multiclass convergence test.
U_grid <- c(0.04, 0.02, 0.01, 0.005, 0.0025, 0.00125, 0.000625)

# Separate five-point grid used for the auxiliary all-class repair check
# reported in the manuscript. This is intentionally distinct from the
# seven-point Panel A grid above.
U_grid_all_class <- c(0.020, 0.010, 0.005, 0.002, 0.001)

# In the auxiliary check, every damaged class has the same repair rate.
nu_all_class_values <- c(0.10, 0.50)

# Equilibrium-integration settings for the multiclass system.
equilibrium_settings <- list(
  method = "lsoda",
  rtol = 1e-11,
  atol = 1e-13,
  maxsteps = 1e6,
  chunk = 100,
  max_chunks = 50L,
  derivative_tolerance = 1e-11
)

# Panel B: threshold trajectory grid.
threshold_times <- unique(c(
  0,
  seq(0.1, 10, length.out = 250),
  seq(10, 100, length.out = 300),
  seq(100, 1000, length.out = 350)
))

# Panel C: small positive repair values at U = U_c.
nu_grid <- c(
  1e-2, 5e-3, 3e-3, 2e-3, 1e-3,
  5e-4, 3e-4, 2e-4, 1e-4,
  5e-5, 3e-5, 2e-5, 1e-5,
  5e-6, 3e-6, 2e-6, 1e-6
)

# ------------------------------------------------------------------------------
# 2. Basic validation
# ------------------------------------------------------------------------------
# The primary biological regime discussed in the manuscript is
# 0 < epsilon < 1, but the analytical results do not require
# the upper bound epsilon < 1.

stopifnot(
  pars$n >= 1L,
  pars$epsilon > 0,
  #pars$epsilon < 1,
  pars$zeta_rho > 0,
  pars$zeta_rho < 1,
  pars$zeta_gamma > 1,
  pars$zeta_gamma > pars$zeta_rho,
  nu1_panel_A >= 0,
  all(U_grid > 0),
  all(U_grid < 1),
  all(U_grid_all_class > 0),
  all(U_grid_all_class < 1),
  all(nu_all_class_values >= 0),
  all(nu_grid > 0)
)

d <- pars$zeta_gamma - pars$zeta_rho
lambda <- pars$epsilon * d
U_c <- 1 - pars$zeta_rho / pars$zeta_gamma

# ------------------------------------------------------------------------------
# 3. Full nonlinear multiclass model with repair
# ------------------------------------------------------------------------------

feedback_ratio <- function(f, pars) {
  classes <- 0:pars$n
  w_rho <- sum(pars$zeta_rho^classes * f)
  w_gamma <- sum(pars$zeta_gamma^classes * f)
  
  if (!is.finite(w_rho) || w_rho <= 0) {
    stop("Mean relative fecundity became nonpositive or nonfinite.")
  }
  
  pars$epsilon * w_gamma / w_rho
}

make_multiclass_rhs <- function(U, pars, nu) {
  n <- pars$n
  
  if (length(nu) != n) {
    stop("nu must have length n, with nu[i] representing repair i -> i-1.")
  }
  
  function(time, state, parameters) {
    f <- state
    R <- feedback_ratio(f, pars)
    
    df <- numeric(n + 1L)
    
    # Class 0.
    df[1L] <-
      ((1 - U) * R - pars$epsilon) * f[1L] +
      nu[1L] * f[2L]
    
    # Classes 1,...,n-1.
    if (n > 1L) {
      for (i in 1:(n - 1L)) {
        df[i + 1L] <-
          U * pars$zeta_rho^(i - 1L) * R * f[i] +
          (
            (1 - U) * pars$zeta_rho^i * R -
              pars$epsilon * pars$zeta_gamma^i -
              nu[i]
          ) * f[i + 1L] +
          nu[i + 1L] * f[i + 2L]
      }
    }
    
    # Terminal class n: no mutation loss beyond the cap.
    df[n + 1L] <-
      U * pars$zeta_rho^(n - 1L) * R * f[n] +
      (
        pars$zeta_rho^n * R -
          pars$epsilon * pars$zeta_gamma^n -
          nu[n]
      ) * f[n + 1L]
    
    list(df)
  }
}

equilibrate_multiclass <- function(U, pars, nu, settings) {
  rhs <- make_multiclass_rhs(U, pars, nu)
  
  state <- c(1, rep(0, pars$n))
  names(state) <- paste0("f", 0:pars$n)
  
  converged <- FALSE
  derivative_norm <- Inf
  population_error <- Inf
  
  for (chunk_id in seq_len(settings$max_chunks)) {
    out <- deSolve::ode(
      y = state,
      times = c(0, settings$chunk),
      func = rhs,
      parms = NULL,
      method = settings$method,
      rtol = settings$rtol,
      atol = settings$atol,
      maxsteps = settings$maxsteps
    )
    
    state <- as.numeric(out[nrow(out), -1L])
    names(state) <- paste0("f", 0:pars$n)
    
    derivative <- rhs(0, state, NULL)[[1L]]
    derivative_norm <- max(abs(derivative))
    population_error <- abs(sum(state) - 1)
    
    if (
      derivative_norm < settings$derivative_tolerance &&
      population_error < 1e-9
    ) {
      converged <- TRUE
      break
    }
  }
  
  if (!converged) {
    warning(
      sprintf(
        "Equilibrium convergence criterion not met for U = %.8g; max |df/dtau| = %.3e",
        U, derivative_norm
      )
    )
  }
  
  list(
    frequency = state,
    converged = converged,
    derivative_norm = derivative_norm,
    population_error = population_error
  )
}

mutation_load <- function(f, pars) {
  classes <- 1:pars$n
  s <- pars$epsilon * (
    pars$zeta_gamma^classes - pars$zeta_rho^classes
  )
  
  sum(s * f[classes + 1L])
}

# ------------------------------------------------------------------------------
# 4. Panel A data: E_L(U) = O(U^2)
# ------------------------------------------------------------------------------

nu_panel_A <- rep(0, pars$n)
nu_panel_A[1L] <- nu1_panel_A

panel_A_rows <- lapply(
  U_grid,
  function(U) {
    eq <- equilibrate_multiclass(
      U = U,
      pars = pars,
      nu = nu_panel_A,
      settings = equilibrium_settings
    )
    
    L_numeric <- mutation_load(eq$frequency, pars)
    L_asymptotic <-
      pars$epsilon * U * lambda / (lambda + nu1_panel_A)
    
    data.frame(
      U = U,
      L_numeric = L_numeric,
      L_asymptotic = L_asymptotic,
      absolute_error = abs(L_numeric - L_asymptotic),
      error_over_U2 = abs(L_numeric - L_asymptotic) / U^2,
      derivative_norm = eq$derivative_norm,
      population_error = eq$population_error,
      converged = eq$converged
    )
  }
)

panel_A <- do.call(rbind, panel_A_rows)
panel_A <- panel_A[order(panel_A$U), ]

slope_fit <- stats::lm(
  log(panel_A$absolute_error) ~ log(panel_A$U)
)
panel_A_slope <- unname(stats::coef(slope_fit)[2L])

# Reference U^2 line, anchored at the smallest-U point.
U_reference <- sort(panel_A$U)
anchor_index <- which.min(panel_A$U)
reference_constant <-
  panel_A$absolute_error[anchor_index] / panel_A$U[anchor_index]^2
error_reference <- reference_constant * U_reference^2

# ------------------------------------------------------------------------------
# 4b. Auxiliary five-point check: repair at every damaged class
# ------------------------------------------------------------------------------

run_all_class_repair_check <- function(nu_value) {
  nu_all <- rep(nu_value, pars$n)
  
  rows <- lapply(
    U_grid_all_class,
    function(U) {
      eq <- equilibrate_multiclass(
        U = U,
        pars = pars,
        nu = nu_all,
        settings = equilibrium_settings
      )
      
      L_numeric <- mutation_load(eq$frequency, pars)
      
      # Corollary 1 depends at leading order only on nu_1. Here every
      # damaged class has the same repair rate, so nu_1 = nu_value.
      L_asymptotic <-
        pars$epsilon * U * lambda / (lambda + nu_value)
      
      data.frame(
        nu = nu_value,
        U = U,
        L_numeric = L_numeric,
        L_asymptotic = L_asymptotic,
        absolute_error = abs(L_numeric - L_asymptotic),
        error_over_U2 =
          abs(L_numeric - L_asymptotic) / U^2,
        derivative_norm = eq$derivative_norm,
        population_error = eq$population_error,
        converged = eq$converged
      )
    }
  )
  
  dat <- do.call(rbind, rows)
  dat <- dat[order(dat$U), ]
  
  fit <- stats::lm(
    log(dat$absolute_error) ~ log(dat$U)
  )
  
  list(
    data = dat,
    slope = unname(stats::coef(fit)[2L])
  )
}

all_class_checks <- lapply(
  nu_all_class_values,
  run_all_class_repair_check
)

all_class_repair <- do.call(
  rbind,
  lapply(all_class_checks, `[[`, "data")
)

all_class_slopes <- vapply(
  all_class_checks,
  `[[`,
  numeric(1),
  "slope"
)

names(all_class_slopes) <- sprintf(
  "nu_%.2f",
  nu_all_class_values
)

# ------------------------------------------------------------------------------
# 5. Panel B data: algebraic threshold relaxation
# ------------------------------------------------------------------------------

# At U = U_c and nu = 0, writing y = 1 - x gives
#   y' = -epsilon * U_c * y^2 / (zeta_rho + a y).
a <- 1 - pars$zeta_rho

threshold_y_rhs <- function(time, state, parameters) {
  y <- state[1L]
  dy <-
    -pars$epsilon * U_c * y^2 /
    (pars$zeta_rho + a * y)
  
  list(c(y = dy))
}

threshold_out <- deSolve::ode(
  y = c(y = 1),
  times = threshold_times,
  func = threshold_y_rhs,
  parms = NULL,
  method = "lsoda",
  rtol = 1e-11,
  atol = 1e-13,
  maxsteps = 1e6
)
threshold_out <- as.data.frame(threshold_out)

panel_B <- subset(threshold_out, time > 0)
panel_B$normalized_dynamic_scaling <-
  pars$epsilon * U_c * panel_B$time * panel_B$y / pars$zeta_rho

# ------------------------------------------------------------------------------
# 6. Exact binary equilibrium and Panel C repair scaling
# ------------------------------------------------------------------------------

binary_equilibrium <- function(U, nu, pars) {
  a <- 1 - pars$zeta_rho
  b <- pars$zeta_gamma - 1
  d <- pars$zeta_gamma - pars$zeta_rho
  
  cU <- d - U * b
  A <- pars$epsilon * cU + nu * a
  B <- pars$epsilon * (U + cU) + nu
  Delta <- B^2 - 4 * A * pars$epsilon * U
  
  if (Delta < -1e-12) {
    stop("Binary discriminant became negative beyond roundoff.")
  }
  Delta <- max(0, Delta)
  
  # Rationalized smaller root; numerically stable for small U.
  2 * pars$epsilon * U / (B + sqrt(Delta))
}

x_star_threshold <- vapply(
  nu_grid,
  function(nu) binary_equilibrium(U_c, nu, pars),
  numeric(1)
)

repair_asymptotic <- sqrt(
  pars$zeta_rho * nu_grid /
    (pars$epsilon * U_c)
)

panel_C <- data.frame(
  nu = nu_grid,
  x_star = x_star_threshold,
  one_minus_x_star = 1 - x_star_threshold,
  asymptotic = repair_asymptotic,
  normalized_repair_scaling =
    (1 - x_star_threshold) / repair_asymptotic
)
panel_C <- panel_C[order(panel_C$nu), ]

# ------------------------------------------------------------------------------
# 7. Write numerical tables
# ------------------------------------------------------------------------------

utils::write.csv(
  panel_A,
  file.path(output_dir, "multiclass_error_scaling.csv"),
  row.names = FALSE
)

utils::write.csv(
  panel_B,
  file.path(output_dir, "threshold_dynamic_scaling.csv"),
  row.names = FALSE
)

utils::write.csv(
  panel_C,
  file.path(output_dir, "threshold_repair_scaling.csv"),
  row.names = FALSE
)

utils::write.csv(
  all_class_repair,
  file.path(output_dir, "all_class_repair_error_scaling.csv"),
  row.names = FALSE
)

summary_lines <- c(
  "Numerical validation summary",
  "============================",
  "",
  sprintf("n = %d", pars$n),
  sprintf("epsilon = %.8g", pars$epsilon),
  sprintf("zeta_rho = %.8g", pars$zeta_rho),
  sprintf("zeta_gamma = %.8g", pars$zeta_gamma),
  sprintf("lambda = epsilon * (zeta_gamma - zeta_rho) = %.8g", lambda),
  sprintf("Panel-A nu_1 = %.8g", nu1_panel_A),
  sprintf("Binary U_c = %.12g", U_c),
  "",
  sprintf("Empirical slope of multiclass load error vs U = %.6f", panel_A_slope),
  sprintf(
    "Largest equilibrium derivative residual in Panel A = %.3e",
    max(panel_A$derivative_norm)
  ),
  sprintf(
    "Largest population-conservation residual in Panel A = %.3e",
    max(panel_A$population_error)
  ),
  sprintf(
    "Final normalized threshold dynamic quantity = %.8f",
    tail(panel_B$normalized_dynamic_scaling, 1L)
  ),
  sprintf(
    "Smallest-nu normalized repair quantity = %.8f",
    panel_C$normalized_repair_scaling[which.min(panel_C$nu)]
  ),
  sprintf(
    "Panel-A fit range: %.6g <= U <= %.6g (%d points)",
    min(panel_A$U),
    max(panel_A$U),
    nrow(panel_A)
  ),
  "Auxiliary all-class repair convergence check",
  sprintf(
    "All-class fit range: %.6g <= U <= %.6g (%d points)",
    min(U_grid_all_class),
    max(U_grid_all_class),
    length(U_grid_all_class)
  ),
  sprintf(
    "All-class repair slope, nu_i = 0.1 for all i >= 1 = %.6f",
    all_class_slopes[1L]
  ),
  sprintf(
    "All-class repair slope, nu_i = 0.5 for all i >= 1 = %.6f",
    all_class_slopes[2L]
  ),
  sprintf(
    "Largest equilibrium derivative residual in all-class check = %.3e",
    max(all_class_repair$derivative_norm)
  ),
  sprintf(
    "Largest population-conservation residual in all-class check = %.3e",
    max(all_class_repair$population_error)
  )
)

writeLines(
  summary_lines,
  con = file.path(output_dir, "numerical_validation_summary.txt")
)

# ------------------------------------------------------------------------------
# 8. Three-panel figure
# ------------------------------------------------------------------------------

draw_validation_panel <- function(device = c("pdf", "png")) {
  device <- match.arg(device)
  
  if (device == "pdf") {
    figure_file <- file.path(
      output_dir,
      "numerical_validation_panel.pdf"
    )
    grDevices::pdf(
      figure_file,
      width = 10.8,
      height = 4.3,
      useDingbats = FALSE
    )
  } else {
    figure_file <- file.path(
      output_dir,
      "numerical_validation_panel.png"
    )
    grDevices::png(
      figure_file,
      width = 2400,
      height = 1000,
      res = 220
    )
  }
  
  op <- graphics::par(
    mfrow = c(1, 3),
    mar = c(4.4, 4.7, 2.8, 1.0),
    oma = c(0.2, 0.2, 2.0, 0.2)
  )
  
  on.exit({
    graphics::par(op)
    grDevices::dev.off()
  }, add = TRUE)
  
  # --------------------------------------------------------------------------
  # Panel A: rare-U multiclass error scaling.
  # --------------------------------------------------------------------------
  
  graphics::plot(
    panel_A$U,
    panel_A$absolute_error,
    log = "xy",
    type = "b",
    pch = 16,
    lwd = 1.7,
    xlab = expression(U),
    ylab = expression(abs(L[num]^"*" - L[asym]^"*")),
    main = "A. Multiclass rare-mutation error"
  )
  
  graphics::lines(
    U_reference,
    error_reference,
    lty = 2,
    lwd = 1.5
  )
  
  graphics::legend(
    "topleft",
    legend = c(
      sprintf("Numerical error; slope = %.3f", panel_A_slope),
      expression(O(U^2))
    ),
    lty = c(1, 2),
    pch = c(16, NA),
    lwd = c(1.7, 1.5),
    bty = "n",
    cex = 0.82
  )
  
  # --------------------------------------------------------------------------
  # Panel B: normalized algebraic relaxation at the no-repair threshold.
  # Quantity should converge to 1.
  # --------------------------------------------------------------------------
  
  graphics::plot(
    panel_B$time,
    panel_B$normalized_dynamic_scaling,
    type = "l",
    lwd = 1.8,
    log = "x",
    xlab = expression(tau),
    ylab = expression(
      epsilon * U[c] * tau * (1 - x(tau)) / zeta[rho]
    ),
    main = "B. Threshold relaxation"
  )
  
  graphics::abline(
    h = 1,
    lty = 2,
    lwd = 1.5
  )
  
  graphics::legend(
    "bottomright",
    legend = c(
      "Numerical trajectory",
      "Asymptotic limit = 1"
    ),
    lty = c(1, 2),
    lwd = c(1.8, 1.5),
    bty = "n",
    cex = 0.82
  )
  
  # --------------------------------------------------------------------------
  # Panel C: normalized square-root repair displacement at U = U_c.
  # Quantity should converge to 1 as nu -> 0.
  # --------------------------------------------------------------------------
  
  graphics::plot(
    panel_C$nu,
    panel_C$normalized_repair_scaling,
    log = "x",
    type = "b",
    pch = 16,
    lwd = 1.7,
    xlab = expression(nu),
    ylab = expression(
      (1 - x^"*") /
        sqrt(zeta[rho] * nu / (epsilon * U[c]))
    ),
    main = "C. Repair regularization"
  )
  
  graphics::abline(
    h = 1,
    lty = 2,
    lwd = 1.5
  )
  
  graphics::legend(
    "topright",
    inset = c(0.02, 0.04),
    legend = c(
      "Exact binary equilibrium",
      "Asymptotic limit = 1"
    ),
    lty = c(1, 2),
    pch = c(16, NA),
    lwd = c(1.7, 1.5),
    bty = "n",
    cex = 0.82
  )
  
  graphics::mtext(
    bquote(
      epsilon == .(pars$epsilon) ~ "," ~
        zeta[rho] == .(pars$zeta_rho) ~ "," ~
        zeta[gamma] == .(pars$zeta_gamma)
    ),
    side = 3,
    outer = TRUE,
    line = -0.1,
    cex = 0.9
  )
  
  invisible(figure_file)
}

pdf_file <- draw_validation_panel("pdf")
png_file <- draw_validation_panel("png")

# ------------------------------------------------------------------------------
# 9. Console summary
# ------------------------------------------------------------------------------

cat("\nNumerical validation completed.\n")
cat("Output directory:", normalizePath(output_dir), "\n\n")
cat(paste(summary_lines, collapse = "\n"), "\n\n")
cat("Figure files:\n")
cat("  ", pdf_file, "\n", sep = "")
cat("  ", png_file, "\n", sep = "")

