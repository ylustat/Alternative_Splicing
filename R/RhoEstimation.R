#' @title Estimation of relative contribution of the fixed effects and random effects.
#'
#' @description \code{RhoEstimation} estimates the relative contribution of the fixed effects and random effects. It returns \eqn{\rho} which is informative in characterizing the relative contributions of main (fixed) effects and isoform-specific TE- TE effects (interaction/random effects) when studying regulation patterns of co-expressed genes.
#' @details
#' When \eqn{\rho = 0}, the joint association test is equivalent to the standard correlation test adjusting for IER and there is no splicing-induced sample-level variation in TE-TE effect. When ρ is large, the splicing-induced random slope heavily affects the variation of TE of gene Y.
#'
#' @param y gene expression level of the response gene Y. It is a standardized numeric vector.
#' @param x gene expression level of the regulator gene X. It is a standardized numeric vector.
#' @param Z covariates needed to be adjusted. It is a numeric vector or matrix.
#' @param IER intron excision ratio matrix of regulator gene X obtained from leafcutter or \link[ASTest]{JuncToIER}.
#' @param ... additional arguments.
#' @return A list containing:
#' \item{rho}{relative contribution of the fixed effects and random effects in the joint association test.}
#' \item{sigmaint}{estimated variance component of the random slope of TE.}
#' @export

RhoEstimation <- function(y,x,Z,IER) {
  # First get the random effects
  n <- nrow(Z)
  formula <- formula(paste0("y~","x+",paste0(colnames(Z),collapse = "+"),
                            "+",paste0(colnames(IER),collapse = "+")))
  data_input = data.frame(y = y, x = x, Z, IER)
  diagx <- diag(as.vector(x))
  colnames(diagx) <- rownames(diagx) <- rownames(x)
  Sigma_x <- diagx%*%cov(t(scale(IER)))%*%diagx
  rho_obj <- polygenic(formula, kin=Sigma_x, data = data_input, quiet = T)
  sigma_int <- 2 * tail(rho_obj$h2an$estimate,1)*
    tail(rho_obj$h2an$estimate,2)[1]
  P = diag(1,n) - matrix(1,n,n)/n
  tr <- function(M) {
    return(sum(diag(M)))
  }
  var_int = sigma_int * tr(P%*%Sigma_x)/(n-1)
  # Get the fixed effects
  formula <- formula(paste0("y~x+",paste0(colnames(Z), collapse = "+")))
  iat_obj <- polygenic(formula, kin=cov(t(scale(IER))), data = data_input, quiet = T)
  beta <- iat_obj$h2an$estimate[2]
  var_fix <- var(x) * beta^2
  rho = var_int/(var_int + var_fix)
  return(list(rho = as.numeric(rho),var_int = sigma_int))
}
