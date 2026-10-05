# ============================================================================
# RepairMediatedPathFiltering_AllFigures_v4.R
#
# ~/Desktop/PathFiltering_v4_FigPlotter.R
#
# Numerical validation of repair-mediated path filtering for:
# "Repair-Mediated Path Filtering of Mutation Load
#  in Multi-Locus Homeostatic Populations"
#
# Compares the exact finite-U nonlinear suppression ratio
#
#     Q_k(U, nu) = L_k^*(U, nu) / L_k^*(U, 0)
#
# against the leading-order path-filtering prediction E_k,
#
#     Q_k(U, nu) = E_k + O(U).
#
# Because L_k^* = lambda_k f_k^* at fixed fitness landscape,
#
#     Q_k(U, nu) = f_k^*(U, nu) / f_k^*(U, 0).
#
# The no-repair denominator is evaluated from the exact analytical
# equilibrium rather than by a second numerical integration.
#
# Panels: k = 1, 2, 3
#
# Parameters:
#   n = 4, eps = 0.5, zeta_rho = 0.8, U = 0.01
#   zeta_gamma in {1.15, 1.30, 1.50}
#   nu in [0, 1]
#
# # June 28 2026
# v4
# Last Updated: October 4, 2026
# ============================================================================

library(deSolve)

# ---- Model parameters -------------------------------------------------------
n       <- 4L
eps     <- 0.5
zeta_r  <- 0.8
U       <- 0.01
T_end   <- 20000

zeta_g_vals <- c(1.15, 1.30, 1.50)
nu_vals     <- seq(0, 1, length.out = 41)

colors <- c("#2166ac", "#4dac26", "#d6604d")
labels <- expression(zeta[gamma] == 1.15,
                     zeta[gamma] == 1.30,
                     zeta[gamma] == 1.50)

# ---- Selection coefficient lambda_k ----------------------------------------
lam <- function(k, zg, zr = zeta_r, e = eps) {
  e * (zg^k - zr^k)
}

# ---- Leading mutation-source weight w_{jk} ---------------------------------
wjk <- function(j, k, zr = zeta_r, e = eps) {
  e * zr^j * choose(n - j, k - j)
}

# ---- Exact no-repair leading coefficients c_k^(0) --------------------------
compute_c0 <- function(zg) {
  d <- zg - zeta_r
  choose(n, 0:n) / d^(0:n)
}

# ---- Repair coefficients c_k(nu) -------------------------------------------
compute_cnu <- function(nu, zg) {
  cnu    <- numeric(n + 1)
  cnu[1] <- 1.0
  
  for (k in 1:n) {
    j_idx <- 0:(k - 1)
    numer <- sum(wjk(j_idx, k) * cnu[j_idx + 1])
    cnu[k + 1] <- numer / (lam(k, zg) + k * nu)
  }
  
  cnu
}

# ---- Analytical suppression E_k = c_k(nu) / c_k^(0) ------------------------
analytic_Ek <- function(k, nu, zg) {
  c0  <- compute_c0(zg)
  cnu <- compute_cnu(nu, zg)
  cnu[k + 1] / c0[k + 1]
}

# ---- Exact no-repair finite-U class equilibrium -----------------------------
fk_norepair_exact <- function(k, zg, Uval = U) {
  denom <- (1 - Uval) * zg - zeta_r
  
  if (denom <= 0) {
    stop(sprintf(
      paste0("No pristine-positive no-repair equilibrium: ",
             "(1-U) zeta_gamma - zeta_rho = %.6g <= 0."),
      denom
    ))
  }
  
  y <- Uval / denom
  choose(n, k) * y^k / (1 + y)^n
}

# ---- Exact birth rate into class k ------------------------------------------
birth_into_k <- function(k, f, R, zr, Uval) {
  total <- 0.0
  
  for (j in 0:k) {
    nmut    <- k - j
    n_avail <- n - j
    
    if (nmut > n_avail) next
    
    prob <- choose(n_avail, nmut) *
      Uval^nmut *
      (1 - Uval)^(n_avail - nmut)
    
    total <- total + R * zr^j * f[j + 1] * prob
  }
  
  total
}

# ---- ODE right-hand side ----------------------------------------------------
odes <- function(t, state, parms) {
  f    <- state
  zg   <- parms$zg
  zr   <- parms$zr
  nu   <- parms$nu
  Uval <- parms$U
  
  mean_g <- sum(zg^(0:n) * f)
  mean_r <- sum(zr^(0:n) * f)
  R      <- eps * mean_g / mean_r
  
  df <- numeric(n + 1)
  
  for (k in 0:n) {
    births     <- birth_into_k(k, f, R, zr, Uval)
    deaths     <- eps * zg^k * f[k + 1]
    repair_out <- k * nu * f[k + 1]
    repair_in  <- if (k < n) (k + 1) * nu * f[k + 2] else 0.0
    
    df[k + 1] <- births - deaths - repair_out + repair_in
  }
  
  list(df)
}

# ---- Integrate to steady state ----------------------------------------------
integrate_ss <- function(nu, zg, Uval = U,
                         ss_tol = 1e-10, mass_tol = 1e-10) {
  f_init <- c(1.0, rep(0.0, n))
  parms  <- list(zg = zg, zr = zeta_r, nu = nu, U = Uval)
  
  sol <- ode(
    y      = f_init,
    times  = c(0, T_end),
    func   = odes,
    parms  = parms,
    method = "lsoda",
    rtol   = 1e-10,
    atol   = 1e-12
  )
  
  f_end <- as.numeric(sol[2, 2:(n + 2)])
  
  rhs      <- odes(T_end, f_end, parms)[[1]]
  rhs_norm <- max(abs(rhs))
  
  if (rhs_norm > ss_tol) {
    warning(sprintf(
      paste0("Steady-state tolerance not met: ",
             "nu=%.4g, zeta_gamma=%.4g, ||df/dtau||_inf=%.3e"),
      nu, zg, rhs_norm
    ))
  }
  
  mass_error <- abs(sum(f_end) - 1.0)
  
  if (mass_error > mass_tol) {
    warning(sprintf(
      paste0("Simplex conservation tolerance not met: ",
             "nu=%.4g, zeta_gamma=%.4g, |sum(f)-1|=%.3e"),
      nu, zg, mass_error
    ))
  }
  
  f_end
}

# ---- Exact finite-U numerical suppression ratio Q_k -------------------------
numerical_Qk_from_state <- function(k, f_nu, zg) {
  f0_exact <- fk_norepair_exact(k, zg)
  
  if (f0_exact < 1e-30) return(NA_real_)
  
  f_nu[k + 1] / f0_exact
}

# ---- Internal consistency checks --------------------------------------------
for (zg in zeta_g_vals) {
  err <- max(abs(compute_cnu(0, zg) - compute_c0(zg)))
  
  if (err > 1e-10) {
    warning(sprintf(
      "Leading-coefficient consistency check failed for zeta_gamma=%.2f: %.3e",
      zg, err
    ))
  }
}

# ---- Compute analytical and numerical arrays --------------------------------
cat("Computing finite-U nonlinear validation...\n")

ana_list <- list()
num_list <- list()

for (zi in seq_along(zeta_g_vals)) {
  zg <- zeta_g_vals[zi]
  
  for (k in 1:3) {
    key <- paste0("k", k, "_z", zi)
    
    ana_list[[key]] <- sapply(
      nu_vals,
      function(nu) analytic_Ek(k, nu, zg)
    )
    
    num_list[[key]] <- numeric(length(nu_vals))
  }
  
  # One nonlinear integration per (nu, zeta_gamma) supplies k = 1,2,3.
  for (ni in seq_along(nu_vals)) {
    nu   <- nu_vals[ni]
    f_nu <- integrate_ss(nu, zg)
    
    for (k in 1:3) {
      key <- paste0("k", k, "_z", zi)
      num_list[[key]][ni] <- numerical_Qk_from_state(k, f_nu, zg)
    }
  }
  
  cat(sprintf("  zeta_gamma=%.2f: done\n", zg))
}

# ---- Three-panel main figure -------------------------------------------------
png(
  "cascade_figure_R.png",
  width = 9,
  height = 3.8,
  units = "in",
  res = 300
)

par(
  mfrow    = c(1, 3),
  mar      = c(4.2, 4.2, 2.2, 0.8),
  mgp      = c(2.5, 0.7, 0),
  cex.axis = 0.9,
  cex.lab  = 1.0
)

panel_labels <- c("(a)  k = 1", "(b)  k = 2", "(c)  k = 3")

for (ki in 1:3) {
  k <- ki
  
  plot(
    NA,
    xlim = c(0, 1),
    ylim = c(0, 1.05),
    xlab = expression("repair rate " * nu),
    ylab = if (ki == 1) expression(Q[k](U, nu)) else "",
    main = panel_labels[ki],
    font.main = 1,
    adj = 0
  )
  
  abline(h = 1.0, lty = 2, col = "black", lwd = 1.0)
  
  for (zi in seq_along(zeta_g_vals)) {
    key <- paste0("k", k, "_z", zi)
    col <- colors[zi]
    
    E_curve <- ana_list[[key]]
    Q_num   <- num_list[[key]]
    
    lines(nu_vals, E_curve, col = col, lwd = 2)
    points(nu_vals, Q_num, pch = 1, col = col, cex = 0.7, lwd = 1.0)
  }
  
  if (ki == 1) {
    legend(
      "topright",
      legend = labels,
      col    = colors,
      lty    = 1,
      lwd    = 2,
      bty    = "n",
      cex    = 0.78
    )
  }
}

dev.off()

cat("Figure saved to cascade_figure_R.png\n")

# =============================================================================
# Supplementary Figure S1: n-collapse check for fixed k
#
# The theorem predicts that for fixed k, the leading-order suppression factor
# E_k is independent of the total locus number n once n >= k. The finite-U
# nonlinear ratio
#
#   Q_k(U, nu; n) = L_k^*(U,nu;n) / L_k^*(U,0;n)
#
# need not be exactly n-independent at fixed U, but should lie close to the
# same analytical curve E_k(nu) when U is small.
# =============================================================================

# ---- Settings for the supplementary n-collapse check ------------------------
k_target        <- 3L
n_vals_sup      <- c(3L, 4L, 6L, 8L)
zg_sup          <- 1.30
U_sup           <- 0.01
nu_vals_sup     <- seq(0, 1, length.out = 41)

sup_colors <- c("#1b9e77", "#d95f02", "#7570b3", "#e7298a")
sup_pch    <- c(1, 2, 0, 5)

# ---- n-dependent helper functions -------------------------------------------
lam_sup <- function(k, zg, zr = zeta_r, e = eps) {
  e * (zg^k - zr^k)
}

wjk_sup <- function(j, k, n_total, zr = zeta_r, e = eps) {
  e * zr^j * choose(n_total - j, k - j)
}

compute_c0_sup <- function(n_total, zg) {
  d <- zg - zeta_r
  choose(n_total, 0:n_total) / d^(0:n_total)
}

compute_cnu_sup <- function(n_total, nu, zg) {
  cnu <- numeric(n_total + 1)
  cnu[1] <- 1.0
  
  for (k in 1:n_total) {
    j_idx <- 0:(k - 1)
    numer <- sum(wjk_sup(j_idx, k, n_total) * cnu[j_idx + 1])
    cnu[k + 1] <- numer / (lam_sup(k, zg) + k * nu)
  }
  
  cnu
}

analytic_Ek_sup <- function(k, nu, zg, n_total) {
  c0  <- compute_c0_sup(n_total, zg)
  cnu <- compute_cnu_sup(n_total, nu, zg)
  cnu[k + 1] / c0[k + 1]
}

fk_norepair_exact_sup <- function(k, n_total, zg, Uval) {
  denom <- (1 - Uval) * zg - zeta_r
  
  if (denom <= 0) {
    stop(sprintf(
      paste0("No pristine-positive no-repair equilibrium for n=%d: ",
             "(1-U) zeta_gamma - zeta_rho = %.6g <= 0."),
      n_total, denom
    ))
  }
  
  y <- Uval / denom
  choose(n_total, k) * y^k / (1 + y)^n_total
}

birth_into_k_sup <- function(k, f, R, zr, Uval, n_total) {
  total <- 0.0
  
  for (j in 0:k) {
    nmut    <- k - j
    n_avail <- n_total - j
    
    if (nmut > n_avail) next
    
    prob <- choose(n_avail, nmut) *
      Uval^nmut *
      (1 - Uval)^(n_avail - nmut)
    
    total <- total + R * zr^j * f[j + 1] * prob
  }
  
  total
}

odes_sup <- function(t, state, parms) {
  f       <- state
  n_total <- parms$n_total
  zg      <- parms$zg
  zr      <- parms$zr
  nu      <- parms$nu
  Uval    <- parms$U
  
  mean_g <- sum(zg^(0:n_total) * f)
  mean_r <- sum(zr^(0:n_total) * f)
  R      <- eps * mean_g / mean_r
  
  df <- numeric(n_total + 1)
  
  for (k in 0:n_total) {
    births     <- birth_into_k_sup(k, f, R, zr, Uval, n_total)
    deaths     <- eps * zg^k * f[k + 1]
    repair_out <- k * nu * f[k + 1]
    repair_in  <- if (k < n_total) (k + 1) * nu * f[k + 2] else 0.0
    
    df[k + 1] <- births - deaths - repair_out + repair_in
  }
  
  list(df)
}

integrate_ss_sup <- function(n_total, nu, zg, Uval,
                             ss_tol = 1e-10, mass_tol = 1e-10) {
  f_init <- c(1.0, rep(0.0, n_total))
  parms  <- list(
    n_total = n_total,
    zg      = zg,
    zr      = zeta_r,
    nu      = nu,
    U       = Uval
  )
  
  sol <- ode(
    y      = f_init,
    times  = c(0, T_end),
    func   = odes_sup,
    parms  = parms,
    method = "lsoda",
    rtol   = 1e-10,
    atol   = 1e-12
  )
  
  f_end <- as.numeric(sol[2, 2:(n_total + 2)])
  
  rhs      <- odes_sup(T_end, f_end, parms)[[1]]
  rhs_norm <- max(abs(rhs))
  if (rhs_norm > ss_tol) {
    warning(sprintf(
      paste0("Steady-state tolerance not met in n-collapse check: ",
             "n=%d, nu=%.4g, zeta_gamma=%.4g, ||df/dtau||_inf=%.3e"),
      n_total, nu, zg, rhs_norm
    ))
  }
  
  mass_error <- abs(sum(f_end) - 1.0)
  if (mass_error > mass_tol) {
    warning(sprintf(
      paste0("Simplex conservation tolerance not met in n-collapse check: ",
             "n=%d, nu=%.4g, zeta_gamma=%.4g, |sum(f)-1|=%.3e"),
      n_total, nu, zg, mass_error
    ))
  }
  
  f_end
}

numerical_Qk_sup_from_state <- function(k, f_nu, n_total, zg, Uval) {
  f0_exact <- fk_norepair_exact_sup(k, n_total, zg, Uval)
  
  if (f0_exact < 1e-30) return(NA_real_)
  
  f_nu[k + 1] / f0_exact
}

# ---- Compute supplementary data ---------------------------------------------
cat("Computing supplementary n-collapse check...\n")

sup_rows <- list()

for (n_total in n_vals_sup) {
  for (i in seq_along(nu_vals_sup)) {
    nu   <- nu_vals_sup[i]
    f_nu <- integrate_ss_sup(n_total, nu, zg_sup, U_sup)
    
    Qk_val <- numerical_Qk_sup_from_state(
      k = k_target,
      f_nu = f_nu,
      n_total = n_total,
      zg = zg_sup,
      Uval = U_sup
    )
    
    sup_rows[[length(sup_rows) + 1L]] <- data.frame(
      n = n_total,
      nu = nu,
      Qk = Qk_val
    )
  }
  
  cat(sprintf("  n=%d: done\n", n_total))
}

sup_data <- do.call(rbind, sup_rows)

# Analytical reference curve: theorem says E_k is independent of n for n >= k.
n_ref_sup <- max(n_vals_sup)
E_ref_sup <- sapply(
  nu_vals_sup,
  function(nu) analytic_Ek_sup(k_target, nu, zg_sup, n_ref_sup)
)

utils::write.csv(
  sup_data,
  file = "supplementary_ncollapse_k3.csv",
  row.names = FALSE
)

# ---- Supplementary figure ---------------------------------------------------
png(
  "supplementary_ncollapse_k3.png",
  width = 5.8,
  height = 4.4,
  units = "in",
  res = 300
)

par(
  mar = c(4.4, 4.6, 2.3, 0.8),
  mgp = c(2.6, 0.8, 0),
  cex.axis = 0.9,
  cex.lab  = 1.0
)

plot(
  NA,
  xlim = c(0, 1),
  ylim = c(0, 1.05),
  xlab = expression("repair rate " * nu),
  ylab = bquote(Q[.(k_target)](U, nu * ";" * n)),
  main = bquote("Supplementary: " ~ n * "-collapse check for " ~ k == .(k_target)),
  font.main = 1,
  cex.main = 0.75
)

abline(h = 1, lty = 2, lwd = 1.0, col = "black")

# Common analytic curve E_k
lines(nu_vals_sup, E_ref_sup, lwd = 2.6, col = "black")

for (ii in seq_along(n_vals_sup)) {
  n_total <- n_vals_sup[ii]
  df_n <- subset(sup_data, n == n_total)
  
  lines(df_n$nu, df_n$Qk, col = sup_colors[ii], lwd = 1.4)
  points(df_n$nu, df_n$Qk,
         pch = sup_pch[ii], col = sup_colors[ii], cex = 0.75, lwd = 1.0)
}

legend(
  "topright",
  legend = c("analytic  E_3", paste0("numerical  n=", n_vals_sup)),
  col    = c("black", sup_colors),
  lty    = c(1, rep(1, length(n_vals_sup))),
  lwd    = c(2.6, rep(1.4, length(n_vals_sup))),
  pch    = c(NA, sup_pch),
  bty    = "n",
  cex    = 0.82
)

mtext(
  bquote(
    epsilon == .(eps) ~ "," ~
      zeta[rho] == .(zeta_r) ~ "," ~
      zeta[gamma] == .(zg_sup) ~ "," ~
      U == .(U_sup)
  ),
  side = 3,
  line = 0.2,
  cex = 0.85
)

dev.off()

cat("Supplementary figure saved to supplementary_ncollapse_k3.png\n")
cat("Supplementary data saved to supplementary_ncollapse_k3.csv\n")

# =============================================================================
# Supplementary Figure S1: O(U) convergence of Q_k to E_k
#
# Tests the asymptotic prediction
#
#     Q_k(U, nu) = E_k + O(U)
#
# at the representative slice used in the manuscript:
#   epsilon = 0.5
#   zeta_rho = 0.8
#   zeta_gamma = 1.30
#   nu = 0.5
#   U in {0.02, 0.01, 0.005, 0.0025}
#   k = 1, 2, 3
#
# Outputs:
#   cascade_convergence_R.png
#   cascade_convergence_R.pdf
#   cascade_convergence_R.csv
# =============================================================================

cat("Computing supplementary O(U) convergence check...\n")

U_conv  <- c(0.02, 0.01, 0.005, 0.0025)
nu_conv <- 0.5
zg_conv <- 1.50

resid <- matrix(
  NA_real_,
  nrow = length(U_conv),
  ncol = 3
)

Q_conv <- matrix(
  NA_real_,
  nrow = length(U_conv),
  ncol = 3
)

E_conv <- matrix(
  NA_real_,
  nrow = length(U_conv),
  ncol = 3
)

for (ui in seq_along(U_conv)) {
  
  Uval <- U_conv[ui]
  
  f_nu <- integrate_ss(
    nu   = nu_conv,
    zg   = zg_conv,
    Uval = Uval
  )
  
  for (k in 1:3) {
    
    Qk <- f_nu[k + 1] /
      fk_norepair_exact(
        k,
        zg_conv,
        Uval = Uval
      )
    
    Ek <- analytic_Ek(
      k,
      nu_conv,
      zg_conv
    )
    
    Q_conv[ui, k] <- Qk
    E_conv[ui, k] <- Ek
    resid[ui, k]  <- abs(Qk - Ek)
  }
}

# ---- Log-log convergence slopes ---------------------------------------------

slopes <- sapply(
  1:3,
  function(k) {
    unname(
      coef(
        lm(
          log(resid[, k]) ~ log(U_conv)
        )
      )[2]
    )
  }
)

cat(
  sprintf(
    paste0(
      "Convergence slopes (k=1,2,3): ",
      "%.3f, %.3f, %.3f\n"
    ),
    slopes[1],
    slopes[2],
    slopes[3]
  )
)

# ---- Save convergence data --------------------------------------------------

convergence_data <- data.frame(
  U = U_conv,
  Q_k1 = Q_conv[, 1],
  E_k1 = E_conv[, 1],
  residual_k1 = resid[, 1],
  Q_k2 = Q_conv[, 2],
  E_k2 = E_conv[, 2],
  residual_k2 = resid[, 2],
  Q_k3 = Q_conv[, 3],
  E_k3 = E_conv[, 3],
  residual_k3 = resid[, 3]
)

utils::write.csv(
  convergence_data,
  file = "cascade_convergence_R.csv",
  row.names = FALSE
)

# ---- Draw Supplementary Figure S1 -------------------------------------------

draw_convergence_figure <- function(
    device = c("png", "pdf")
) {
  device <- match.arg(device)
  
  if (device == "png") {
    grDevices::png(
      "cascade_convergence_R.png",
      width = 5.8,
      height = 4.4,
      units = "in",
      res = 300
    )
  } else {
    grDevices::pdf(
      "cascade_convergence_R.pdf",
      width = 5.8,
      height = 4.4,
      useDingbats = FALSE
    )
  }
  
  op <- graphics::par(
    mar = c(4.4, 4.8, 1.2, 0.8),
    mgp = c(2.7, 0.8, 0),
    cex.axis = 0.9,
    cex.lab = 1.0
  )
  
  on.exit({
    graphics::par(op)
    grDevices::dev.off()
  }, add = TRUE)
  
  conv_cols <- c("#2166ac", "#4dac26", "#d6604d")
  conv_pch  <- c(16, 17, 15)
  conv_lty  <- c(1, 2, 3)
  
  graphics::plot(
    NA,
    xlim = range(U_conv),
    ylim = range(resid),
    log = "xy",
    xlab = expression(U),
    ylab = expression(abs(Q[k](U, nu) - E[k]))
  )
  
  for (k in 1:3) {
    
    ord <- order(U_conv)
    
    graphics::lines(
      U_conv[ord],
      resid[ord, k],
      col = conv_cols[k],
      lty = conv_lty[k],
      lwd = 1.8
    )
    
    graphics::points(
      U_conv[ord],
      resid[ord, k],
      col = conv_cols[k],
      pch = conv_pch[k],
      cex = 0.85
    )
  }
  
  # Slope-one reference line, anchored at the smallest-U k=1 residual.
  i_anchor <- which.min(U_conv)
  
  reference_constant <-
    resid[i_anchor, 1] /
    U_conv[i_anchor]
  
  U_reference <- sort(U_conv)
  
  graphics::lines(
    U_reference,
    reference_constant * U_reference,
    col = "grey40",
    lty = 2,
    lwd = 1.4
  )
  
  graphics::legend(
    "topleft",
    legend = c(
      sprintf("k = 1; slope = %.2f", slopes[1]),
      sprintf("k = 2; slope = %.2f", slopes[2]),
      sprintf("k = 3; slope = %.2f", slopes[3]),
      "slope 1 reference"
    ),
    col = c(conv_cols, "grey40"),
    lty = c(conv_lty, 2),
    pch = c(conv_pch, NA),
    lwd = c(1.8, 1.8, 1.8, 1.4),
    bty = "n",
    cex = 0.82
  )
  
  graphics::mtext(
    bquote(
      epsilon == .(eps) ~ "," ~
        zeta[rho] == .(zeta_r) ~ "," ~
        zeta[gamma] == .(zg_conv) ~ "," ~
        nu == .(nu_conv)
    ),
    side = 3,
    line = 0.1,
    cex = 0.85
  )
  
  invisible(NULL)
}

draw_convergence_figure("png")
draw_convergence_figure("pdf")

cat(
  paste0(
    "Supplementary convergence figure saved to ",
    "cascade_convergence_R.png and cascade_convergence_R.pdf\n"
  )
)

cat(
  "Supplementary convergence data saved to cascade_convergence_R.csv\n"
)

