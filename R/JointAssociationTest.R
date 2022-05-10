#' @title Jointly testing for main and interaction effects with splicing of the TE of one gene on the TE of another
#'
#' @description \code{jat} estimates the variance components and tests the splicing-TE interaction, as well as the fixed effect of TE jointly with a variance-component score test.
#' @details
#' This function \code{jat} jointly tested the fixed effect of TE and interaction effect of TE and splicing on other genes.
#' It was implemented based on the StructLMM method and used the GRAMMAR-GAMMA functions to estimate variance components (see reference).
#' Here we adopted a mixed-effects model in estimating the co-expression effects between pairs of genes while allowing sample-level variation as random effects,
#' with the covariance matrix of the random effects being proportional to sample-sample covariance of isoform composition.
#'
#' @param y gene expression level of the response gene Y. It is a standardized numeric vector.
#' @param x gene expression level of the regulator gene X. It is a standardized numeric vector.
#' @param Z covariates needed to be adjusted. It is a numeric vector or matrix.
#' @param IER intron excision ratio matrix of regulator gene X obtained from leafcutter or \link[ASTest]{JuncToIER}.
#' @param ... additional arguments.
#' @return A list containing:
#' \item{CALL}{function call.}
#' \item{data}{data used.}
#' \item{formula}{model under null hypothesis.}
#' \item{P}{the p-value of the joint effects of both fixed and random effects of the TE of one gene on the TE of another.}
#' \item{Q}{a score test statistic, which follows a mixed chi-square distribution.}
#' \item{rho}{relative contribution of the fixed effects and random effects in the joint association test.}
#' \item{residual}{residuals of the model under null hypothesis.}
#' \item{sigmau}{estimated variance component of the random intercept.}
#' \item{sigmae}{estimated variance component of the error term.}
#' \item{sigmaint}{estimated variance component of the random slope of TE.}
#' @references
#' Moore, R., Casale, F.P., Jan Bonder, M. et al. A linear mixed-model approach to study multivariate gene–environment interactions. Nat Genet 51, 180–186 (2019). \url{https://doi.org/10.1038/s41588-018-0271-0}
#'
#' Svishcheva, G., Axenovich, T., Belonogova, N. et al. Rapid variance components–based method for whole-genome association analysis. Nat Genet 44, 1166–1170 (2012). \url{https://doi.org/10.1038/ng.2410}
#' @export

jat <- function(y, x, Z, IER, ...) {
  ### check arguments ##
  ## 1.1 check data type ##
  if (!is.vector(y)) stop("y should be a vector.")
  y_m <- data.matrix(as.numeric(y))
  if (!is.vector(x)) stop("x should be a vector.")
  diagx <- diag(as.numeric(x))
  x_m <- data.matrix(as.numeric(x))
  if (!is.vector(Z)&!is.matrix(Z)) stop("Z should be a vector or a matrix.")
  ## 1.2 check length ##
  n <- nrow(y_m)
  if (length(x_m) != n) stop("x and y does not have the same number of samples.")
  if ((!is.null(Z)) & (nrow(Z) != n)) stop("Z and y does not have the same number of samples.")
  if ((!is.null(IER)) & (nrow(IER) != n)) stop("IER and y does not have the same number of samples.")
  if (ncol(IER) == 1) stop("need at least two IERs to perform variance-component score test.")
  ## 1.3 check IDs ##
  if (!all.equal(row.names(Z),rownames(IER))) stop("rownames of Z do not match rownames of IER")
  rownames(y_m) <- rownames(Z)
  colnames(y_m) <- "y"
  rownames(x_m) <- rownames(Z)
  colnames(x_m) <- "x"
  rownames(diagx) <- colnames(diagx) <- rownames(Z)
  colnames(IER) <- paste0("IER.",1:ncol(IER))
  ## 2.1 set up data and formula ##
  IER_s <- scale(IER)
  Sigma_IER <- cov(t(IER_s))
  Sigma_x <- diagx %*% Sigma_IER %*% diagx
  if(sd(Z[,1])!=0) {
    Z_i = cbind(intercept=1,Z)
  } else {
    Z_i = cbind(intercept=1,Z[,-1])
  }
  formula <- formula(paste0("y~",paste0(colnames(Z),collapse = "+")))

  data_input = data.frame(y = y_m, x_m, Z_i, IER_s)
  ## 2.2 estimate the parameters

  rho_obj <- RhoEstimation(y_m,x_m,Z_i,IER_s)
  rho <- rho_obj$rho
  var_int <- rho_obj$var_int
  # rho <- runif(1,0.8,1)
  jat_obj <- polygenic(formula, kin=Sigma_IER, data = data_input, quiet = T)
  InvH <- jat_obj$InvSigma
  res_asso <- jat_obj$residualY
  Py <- InvH %*% as.matrix(res_asso)
  P = InvH - InvH %*% Z_i %*% solve(t(Z_i) %*% InvH %*% Z_i) %*% t(Z_i) %*% InvH
  Sigma <- rho * Sigma_x + (1-rho) * x %*% t(x)

  sigmau <- 2 * tail(jat_obj$h2an$estimate,1)*tail(jat_obj$h2an$estimate,2)[1]
  sigmae <- 2 * tail(jat_obj$h2an$estimate,1)*(1-tail(jat_obj$h2an$estimate,2)[1])
  ## 2.3 calculate the score Q and p-value
  Q <- t(Py) %*% (Sigma) %*% Py
  lambda <- Re(eigen(P %*% Sigma)$values)
  lambda <- lambda[lambda>mean(lambda[lambda>=0])/100000]
  pval <- SKAT:::Get_PValue.Lambda(lambda,Q)$p.value

  result <- list(CALL = match.call(),
                 data = list(y=y,x=x,Z=Z,IER=IER),
                 formula = formula,
                 P = pval,
                 Q = Q,
                 rho = rho,
                 residual = res_asso,
                 sigmau = sigmau,
                 sigmae = sigmae,
                 sigmaint = var_int)
  return(result)
}


