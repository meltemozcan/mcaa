#' Distribution function (pdf) of a mixture of two normal distributions. 
#' \code{pnormmix} returns the cumulative probability of q or \eqn{1 - q} on the
#'  mixture normal distribution.
#' 
#' @param q A vector of quantiles.
#' @param mean1: Mean of the first normal distribution.
#' @param sd1: Standard deviation of the first normal distribution.
#' @param mean2: Mean of the second normal distribution.
#' @param sd2: Standard deviation of the second normal distribution.
#' @param pmix1: Mixing proportion for the first distribution. Should be a 
#'   number in the range (0, 1).
#' @param lower.tail: A logical scalar; if TRUE(default), probabilities are 
#' \eqn{P[X <= x]}; otherwise, \eqn{P[X > x]}. 
#' @return The output will be the cumulative probability of q or \eqn{1 - q} on 
#' the mixture normal distribution.
#' @examples
#' pnormmix(1, 0, 3.1, 1.7, 3.1, lower.tail = FALSE)
#' 

pnormmix <- function(q, mean1 = 0, sd1 = 1, mean2 = 0, sd2 = 1, pmix1 = 0.5, 
                     lower.tail = TRUE) {
  stopifnot(pmix1 > 0, pmix1 < 1)
  as.vector(c(pmix1, 1 - pmix1) %*% 
              sapply(q, pnorm, mean = c(mean1, mean2), sd = c(sd1, sd2), 
                     lower.tail = lower.tail))
}

#' Quantile function of a mixture of two normal distributions. 
#' \code{qnormmix} returns the quantile corresponding to p or \eqn{1 - q} on 
#' the mixture normal distribution.
#' 
#' @param p A vector of probabilities.
#' @param mean1: Mean of the first normal distribution.
#' @param sd1: Standard deviation of the first normal distribution.
#' @param mean2: Mean of the second normal distribution.
#' @param sd2: Standard deviation of the second normal distribution.
#' @param pmix1: Mixing proportion for the first distribution. Should be a 
#'   number in the range (0, 1).
#' @param lower.tail: A logical scalar; if TRUE(default), probabilities are 
#' \eqn{P[X <= x]}; otherwise, \eqn{P[X > x]}. 
#' @return The output will be the quantile corresponding to p or \eqn{1 - q} on 
#' the mixture normal distribution.
#' @examples
#' qnormmix(0.8, 0, 3.1, 1.7, 0.5, lower.tail = FALSE)
#' 

qnormmix <- function(p, mean1 = 0, sd1 = 1, mean2 = 0, sd2 = 1, pmix1 = 0.5, 
                     lower.tail = TRUE) {
  
  stopifnot(pmix1 > 0, pmix1 < 1, p >= 0, p <= 1)
  f <- function(x) (pnormmix(x, mean1, sd1, mean2, sd2, pmix1, 
                             lower.tail) - p)^2
  start <- as.vector(c(pmix1, 1 - pmix1) %*% 
                       sapply(p, qnorm, c(mean1, mean2), c(sd1, sd2), 
                              lower.tail = lower.tail))
  nlminb(start, f)$par
}

#' Helper function for computing the kernel for bivariate normal density. 
#' \code{.bvnorm_kernel} returns the kernel for bivariate normal density
#' 
#' @param x A normal distribution
#' @param mu_x: Mean of the normal distribution x.
#' @param y A normal distribution
#' @param sd_x: Standard deviation of the normal distribution x.
#' @param mu_y: Mean of the normal distribution y.
#' @param sd_y: Standard deviation of the normal distribution y.
#' @param cov_xy: covariance between x and y
#' 
#' @return The output will be the quantile corresponding to p or \eqn{1 - q} on 
#' the mixture normal distribution.
#' @examples
#' 

.bvnorm_kernel <- function(x, y, mu_x = 0, mu_y = 0, sd_x = 1, sd_y = 1, 
                           cov_xy = 0) {
  
  cor <- cov_xy / sd_x / sd_y
  numer <- (x - mu_x)^2 / sd_x^2 + (y - mu_y)^2 / sd_y^2 - 
    2 * cor * (x - mu_x) * (y - mu_y) / sd_x / sd_y
  numer / (1 - cor^2)
}

#' Plot contour for a bivariate normal distribution
#' \code{.bvnorm_kernel} returns the kernel for bivariate normal density
#' 
#' @param mean1: Mean of first normal distribution (on x-axis).
#' @param sd1: Standard deviation of first normal distribution.
#' @param mean2: Mean of second normal distribution (on y-axis).
#' @param sd2: Standard deviation of second normal distribution.
#' @param cor12: Correlation in the bivariate normal.
#' @param cov12: Covariance in the bivariate normal. If not input, compute the
#'          covariance using the correlation and the standard deviations.
#' @param density: Density level, i.e., probability enclosed by the ellipse.
#' @param length_out: Number of values on the x-axis and on the y-axis to be
#'               evaluated; default to 101.
#' @param bty: Argument passed to the `contour` function.
#' 
#' @return A plot showing the contour of the bivariate normal distribution on 
#'   a two-dimensional space.
#' @examples 
#' 
#'  

contour_bvnorm <- function(mean1 = 0, sd1 = 1, mean2 = 0, sd2 = 1, 
                           cor12 = 0, cov12 = NULL, 
                           density = .95, length_out = 101, 
                           bty = "L", 
                           ...) {
  # Error handling
  stopifnot(cor12 >= -1, cor12 <= 1)
  if (is.null(cov12)) cov12 <- cor12 * sd1 * sd2
  x_seq <- mean1 + seq(-3, 3, length.out = length_out) * sd1
  y_seq <- mean2 + seq(-3, 3, length.out = length_out) * sd2
  z <- outer(x_seq, y_seq, .bvnorm_kernel, mu_x = mean1, mu_y = mean2, 
             sd_x = sd1, sd_y = sd2, cov_xy = cov12)
  contour(x_seq, y_seq, z, levels = qchisq(density, 2), drawlabels = FALSE, 
          bty = bty, ...)
}

#' Computing summary statistics from a selection approach
#' \code{.partit_bvnorm} returns a table of selection accuracy indices
#' 
#' @param cut1: Cut score based on the latent score
#' @param cut2: Cut score based on the observed score
#' @param mean1: Mean of first normal distribution (on x-axis).
#' @param sd1: Standard deviation of first normal distribution.
#' @param mean2: Mean of second normal distribution (on y-axis).
#' @param sd2: Standard deviation of second normal distribution.
#' @param cor12: Correlation in the bivariate normal.
#' @param cov12: Covariance in the bivariate normal. If not input, compute the
#'          covariance using the correlation and the standard deviations.
#' 
#' @return A table of selection accuracy indices
#' @examples 
#' 
#'  

.partit_bvnorm <- function(cut1, cut2, mean1 = 0, sd1 = 1, mean2 = 0, sd2 = 1, 
                           cor12 = 0, cov12 = cor12 * sd1 * sd2) {
  Sigma <- matrix(c(sd1^2, cov12, cov12, sd2^2), nrow = 2)
  C <- pmnorm(c(cut1, cut2), c(mean1, mean2), Sigma)
  B <- pnorm(cut1, mean1, sd1) - C
  D <- pnorm(cut2, mean2, sd2) - C
  A <- 1 - B - C - D
  propsel <- A + B
  success_ratio <- A / propsel
  sensitivity <- A / (A + D)
  specificity <- C / (C + B)
  c(A, B, C, D, propsel, success_ratio, sensitivity, specificity)
}

#' Evaluate selection accuracy based on the MCAA Framework
#' \code{PartInvMulti_we} valuates partial measurement invariance using 
#'  an extension of Millsap & Kwok's (2004) approach
#' 
#' @param propsel: proportion of selection. If missing, computed using `cut_z`.
#' @param cut_z: pre-specified cutoff score on the observed composite. This 
#' argument is ignored when `propsel` has input.
#' @param weights_item: a vector of item weights
#' @param weights_latent: a  vector of latent factor weights
#' @param alpha_r: a vector of latent factor mean for the reference group.
#' @param alpha_f: (optional) a vector of latent factor mean for the focal group; 
#'            if no input, set equal to alpha_r.
#' @param psi_r: a matrix of latent factor variance for the reference group.
#' @param psi_f: (optional) a matrix of latent factor variance for the focal group; 
#'          if no input, set equal to psi_r.
#' @param lambda_r: a matrix of factor loadings for the reference group.
#' @param lambda_f: (optional) a matrix of factor loadings for the focal group; 
#'             if no input, set equal to lambda_r.
#' @param nu_r: a matrix of measurement intercepts for the reference group.
#' @param nu_f: (optional) a matrix of measurement intercepts for the focal group; 
#'          if no input, set equal to nu_r.
#' @param Theta_r: a matrix of the unique factor variances and covariances 
#'            for the reference group.
#' @param Theta_f: (optional) a matrix of the unique factor variances and 
#'            covariances for the focal group; if no input, set equal to Theta_r.
#' @param pmix_ref: Proportion of the reference group; 
#'            default to 0.5 (i.e., two populations have equal size).
#' @param plot_contour: logical; whether the contour of the two populations 
#'            should be plotted; default to TRUE.
#' @return The output will be  a list of four elements and a plot if plot_contour == TRUE:
#'     - propsel: echo the same argument as input.
#'     - cutpt_xi: cut point on the latent scale (xi).
#'     - cutpt_z: cut point on the observed scale (Z).
#'     - summary: A 8 x 2 table, with columns representing the reference
#'                and the focal groups, and the rows represent
#'                probabilities of true positive (A), false positive (B), 
#'                true negative (C), false negative (D); proportion selected, 
#'                success ratio, sensitivity, and specificity. 
#' @examples
#' # single dimension
#' PartInvMulti_we(propsel = .10,
#'                 weights_item = c(1, 0.9, 0.8, 1),
#'                 weights_latent = 0.9,
#'                 alpha_r = 0.5,
#'                 alpha_f = 0,
#'                 psi_r = 1,
#'                 lambda_r = c(.3, .5, .9, .7),
#'                 nu_r = c(.225, .025, .010, .240),
#'                 nu_f = c(.225, -.05, .240, -.025),
#'                 Theta_r = diag(.96, 4))
#' # multiple dimensions
#' lambda_matrix <- matrix(0,nrow = 5, ncol = 2)
#' lambda_matrix[1:2,1] <- c(.322, .655)
#' lambda_matrix[3:5,2] <- c(.398, .745, .543)
#' PartInvMulti_we(propsel = .05,
#'                 weights_item = c(1/4, 1/4, 1/6, 1/6, 1/6),
#'                 weights_latent = c(0.5, 0.5),
#'                 alpha_r = c(0, 0),
#'                 alpha_f = c(-0.3, 0.1),
#'                 psi_r = matrix(c(1, 0.5, 0.5, 1), nrow = 2),
#'                 lambda_r = lambda_matrix,
#'                 nu_r = c(.225, .025, .010, .240, .125),
#'                 nu_f = c(.225, -.05, .240, -.025, .125),
#'                 Theta_r = diag(1, 5),
#'                 Theta_f = c(1, .95, .80, .75, 1))
#' 

PartInvMulti_we <- function(propsel, cut_z = NULL, weights_item, weights_latent,
                            alpha_r, alpha_f = alpha_r, psi_r, psi_f = psi_r,
                            lambda_r, lambda_f = lambda_r,nu_r, nu_f = nu_r,
                            Theta_r, Theta_f = Theta_r, 
                            pmix_ref = 0.5, plot_contour = TRUE, ...) {
  # Error handling
  # stopifnot(length(alpha_r) != 1, length(alpha_f) != 1, length(psi_r) != 1, 
  #           length(psi_f) != 1)
  # check the dimensions of input parameters
  if (length(Theta_r) <= length(lambda_r)) Theta_r <- diag(Theta_r)
  if (length(Theta_f) <= length(lambda_f)) Theta_f <- diag(Theta_f)
  if (is.matrix(weights_item) == FALSE) weights_item <- t(weights_item)
  if (is.matrix(weights_latent) == FALSE) weights_latent <- t(weights_latent)
  if (is.matrix(alpha_r) == FALSE) alpha_r <- as.matrix(alpha_r)
  if (is.matrix(alpha_f) == FALSE) alpha_f <- as.matrix(alpha_f)
  library(mnormt)  # load `mnormt` package
  mean_zr <- c(weights_item %*% nu_r + weights_item %*% lambda_r %*% alpha_r)
  mean_zf <- c(weights_item %*% nu_f + weights_item %*% lambda_f %*% alpha_f)
  sd_zr <- c(sqrt(weights_item %*% lambda_r %*% psi_r %*% t(weights_item %*% lambda_r)+ weights_item %*% Theta_r %*% t(weights_item)))
  sd_zf <- c(sqrt(weights_item %*% lambda_f %*% psi_f %*% t(weights_item %*% lambda_f)+ weights_item %*% Theta_f %*% t(weights_item)))
  cov_z_xir <- c(weights_item %*% lambda_r %*% psi_r %*% t(weights_latent))
  cov_z_xif <- c(weights_item %*% lambda_f %*% psi_f %*% t(weights_latent))
  sd_xir <- c(sqrt(weights_latent %*% psi_r %*% t(weights_latent)))
  sd_xif <- c(sqrt(weights_latent %*% psi_f %*% t(weights_latent)))
  zeta_r <- c(weights_latent%*%alpha_r)
  zeta_f <- c(weights_latent%*%alpha_f)
  # if there is an input for selection proportion
  if (!missing(propsel)) {
    # and if there is an input for cut score
    if (!is.null(cut_z)) {
      warning("Input to `cut_z` is ignored.")
    }
    # compute the cut score using helper function qnormmix based on input selection
    # proportion
    cut_z <- qnormmix(propsel, mean_zr, sd_zr, mean_zf, sd_zf, 
                      pmix_ref, lower.tail = FALSE)
  } else if (!is.null(cut_z) & missing(propsel)) {
    # if missing selection proportion but has a cut score
    propsel <- pnormmix(cut_z, mean_zr, sd_zr, mean_zf, sd_zf, 
                        pmix_ref, lower.tail = FALSE)
  }
  cut_xi <- qnormmix(propsel, zeta_r, sd_xir,zeta_f, sd_xif,
                     pmix_ref, lower.tail = FALSE)
  # computing summary statistics using helper function .partit_bvnorm
  partit_1 <- .partit_bvnorm(cut_xi, cut_z, zeta_r, sd_xir, mean_zr, sd_zr, 
                             cov12 = cov_z_xir)
  partit_2 <- .partit_bvnorm(cut_xi, cut_z, zeta_f, sd_xif, mean_zf, sd_zf, 
                             cov12 = cov_z_xif)
  # selection indices for the focal group if latent dist matches the reference
  mean_zref <- as.numeric(weights_item %*% nu_f + weights_item %*% lambda_f %*% alpha_r)
  sd_zref <- as.numeric(sqrt(weights_item %*% lambda_f %*% psi_r %*% t(weights_item %*% lambda_r)+ weights_item %*% Theta_f %*% t(weights_item)))
  cov_z_xiref <- as.numeric(weights_item %*% lambda_f %*% psi_r %*% t(weights_latent))
  partit_1e2 <- .partit_bvnorm(cut_xi, cut_z, 
                               zeta_r, sd_xir, mean_zref, 
                               sd_zref, 
                               cov12 = cov_z_xiref)
  # result table
  dat <- data.frame("Reference" = partit_1, "Focal" = partit_2, 
                    "E_R(Focal)" = partit_1e2, 
                    row.names = c("A (true positive)", "B (false positive)", 
                                  "C (true negative)", "D (false negative)", 
                                  "Proportion selected", "Success ratio", 
                                  "Sensitivity", "Specificity"))
  # result plot
  p <- NULL
  if (plot_contour) {
    x_lim <- range(c(zeta_r + c(-3, 3) * sd_xir, 
                     zeta_f + c(-3, 3) * sd_xif))
    y_lim <- range(c(mean_zr + c(-3, 3) * sd_zr, 
                     mean_zf + c(-3, 3) *sd_zf))
    contour_bvnorm(zeta_r, sd_xir, mean_zr, sd_zr, cov12 = cov_z_xir, 
                   xlab = bquote("Latent Composite" ~ (zeta)), 
                   ylab = bquote("Observed Composite" ~ (italic(Z))), 
                   lwd = 2, col = "red", xlim = x_lim, ylim = y_lim, 
                   ...)
    contour_bvnorm(zeta_f, sd_xif, mean_zf, sd_zf, cov12 = cov_z_xif, 
                   add = TRUE, lty = "dashed", lwd = 2, col = "blue", 
                   ...)
    legend("topleft", c("Reference group", "Focal group"), 
           lty = c("solid", "dashed"), col = c("red", "blue"))
    abline(h = cut_z, v = cut_xi)
    x_cord <- rep(cut_xi + c(.25, -.25) * sd_xir, 2)
    y_cord <- rep(cut_z + c(.25, -.25) * sd_zr, each = 2)
    text(x_cord, y_cord, c("A", "B", "D", "C"))
    p <- recordPlot()
  }
  # return a list of results and the plot
  list(propsel = propsel, cutpt_xi = cut_xi, cutpt_z = cut_z, 
       summary = round(dat, 3), 
       ai_ratio = dat["Proportion selected", 3] / 
         dat["Proportion selected", 1], plot = p)
}

