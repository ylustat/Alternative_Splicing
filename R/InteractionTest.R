#' @title Modeling and testing for splicing-induced sample-level variation in TE-TE correlation
#'
#' @description \code{iat} estimates the variance components and tests the splicing-TE interaction with a variance-component score test.
#' @details
#' This function \code{iat} tests the interaction effect of TE and splicing on other genes.
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
#' \item{P}{the p-value of the TE-splicing interaction effect on the response gene expression.}
#' \item{Q}{a score test statistic, which follows a mixed chi-square distribution.}
#' \item{residual}{residuals of the model under null hypothesis.}
#' \item{beta}{fixed effect of gene X on gene Y under null hypothesis.}
#' \item{sigmau}{estimated variance component of the random intercept.}
#' \item{sigmae}{estimated variance component of the error term.}
#' @references
#' Moore, R., Casale, F.P., Jan Bonder, M. et al. A linear mixed-model approach to study multivariate gene–environment interactions. Nat Genet 51, 180–186 (2019). \url{https://doi.org/10.1038/s41588-018-0271-0}
#'
#' Svishcheva, G., Axenovich, T., Belonogova, N. et al. Rapid variance components–based method for whole-genome association analysis. Nat Genet 44, 1166–1170 (2012). \url{https://doi.org/10.1038/ng.2410}



iat <- function(y, x, Z, IER, ...) {
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
    Z_i = cbind(intercept=1,x_m,Z)
  } else {
    Z_i = cbind(intercept=1,x=x_m,Z[,-1])
  }
  formula <- formula(paste0("y~",paste0(colnames(Z_i), collapse = "+")))
  # formula1 <- formula(paste0("y~",paste0(colnames(Z)[-2],collapse = "+")))
  # formula2 <- formula(paste0("y~",paste0(colnames(Z),collapse = "+"),"+",paste0(colnames(IER),collapse = "+")))
  data_input = data.frame(y = y_m, Z_i, IER_s)

  ## 2.2 estimate the parameters
  iat_obj <- polygenic(formula, kin=Sigma_IER, data = data_input, quiet = T)
  InvH <- iat_obj$InvSigma
  res <- iat_obj$residualY
  Py <- InvH %*% as.matrix(res)
  P = InvH - InvH %*% Z_i %*% solve(t(Z_i) %*% InvH %*% Z_i) %*% t(Z_i) %*% InvH
  sigmau <- 2 * tail(iat_obj$h2an$estimate,1)*tail(iat_obj$h2an$estimate,2)[1]
  sigmae <- 2 * tail(iat_obj$h2an$estimate,1)*(1-tail(iat_obj$h2an$estimate,2)[1])
  beta <- iat_obj$h2an$estimate[2]

  ## 2.3 calculate the score Q and p-value
  Q <- t(Py) %*% Sigma_x  %*% Py
  lambda <- Re(eigen(P %*% Sigma_x)$values)
  lambda <- lambda[lambda>mean(lambda[lambda>=0])/100000]
  pval <- SKAT:::Get_PValue.Lambda(lambda,Q)$p.value
  result <- list(CALL = match.call(),
                 data = list(y=y,x=x,Z=Z,IER=IER),
                 formula = formula,
                 P = pval,
                 Q = Q,
                 residual = res,
                 beta = beta,
                 sigmau = sigmau,
                 sigmae = sigmae)
  return(result)
}
