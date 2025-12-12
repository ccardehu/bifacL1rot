#' L1 rotation for Generalized Bi-factor Models
#'
#' @param A Input factor loading matrix (p x q)
#' @param Bstart Optional starting value for B matrix. Defaults to A.
#' @param Phi Optional starting correlation matrix. Defaults to identity.
#' @param rho Initial penalty parameter for augmented Lagrangian method. Default 1.
#' @param t Initial step size. Default 1/1000.
#' @param maxit.ou Maximum outer iterations. Default 5000.
#' @param maxit.in Maximum inner iterations. Default 300.
#' @param orthogonal Logical; if TRUE, constrains factors to be orthogonal. Default FALSE.
#' @param tol1 Convergence tolerance for parameter changes. Default 1e-6.
#' @param tol2 Tolerance for backtracking line search. Default 1e-4.
#' @param verbose Logical; print progress. Default TRUE.
#' @param v.every Print frequency (every v.every outer iterations). Default 10.
#' @param Lmax Clipping bound for Lagrange multipliers. Default 20.
#' @param c1 Multiplicative factor for rho increase (must be > 1). Default 1.05.
#' @param c2 Threshold for rho update (must be in (0,1)). Default 0.25.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{B} - Estimated factor loading matrix (follows a bi-factor structure)
#'   \item \code{Phi} - Estimated factor correlation matrix (follows a bi-factor structure)
#'   \item \code{rho.end} - Final value of rho
#'   \item \code{conv} - Logical; TRUE if converged before maxit.ou
#' }
#'
#' @examples
#' \dontrun{
#' set.seed(1234)
#' A <- matrix(0, 20, 4) # A 20 x 4 Factor loading matrix
#' A[,1] = runif(20, 1, 2) # Bi-factor structure, first factor
#' A[,2] = runif(20, 0.5, 1) * rbinom(20, 1, 0.25)
#' A[,3] = runif(20, 0.5, 1) * rbinom(20, 1, 0.25)
#' A[,4] = runif(20, 0.5, 1) * rbinom(20, 1, 0.25)
#'
#' # Create rotation matrix (via random matrix):
#' Ah =  array(rnorm(length(A), sd = .5), dim = dim(A))
#' Tr = eigen(t(Ah) %*% Ah)$vectors
#' Ah = (A)%*%(Tr)
#' res = bifacL1rot::bifactorL1(Ah)
#' Arot = res$B
#' Phi = res$Phi
#' }
#'
#' @export
#' @useDynLib bifacL1rot, .registration = TRUE
#' @importFrom Rcpp sourceCpp
bifactorL1 <- function(A, Bstart = NULL, Phi = NULL, rho = 1, t = 1/1000,
                maxit.ou = 5000, maxit.in = 300, orthogonal = FALSE,
                tol1 = 1e-6, tol2 = 1e-4, verbose = TRUE, v.every = 10,
                Lmax = 20, c1 = 1.05, c2 = 0.25) {

    # Input validation
    if (!is.matrix(A)) A <- as.matrix(A)
    if (!is.null(Bstart) && !is.matrix(Bstart)) Bstart <- as.matrix(Bstart)
    if (!is.null(Phi) && !is.matrix(Phi)) Phi <- as.matrix(Phi)

    result <- ALM_cpp(
        A = A,
        Bstart_ = Bstart,
        Phi_ = Phi,
        rho = rho,
        t = t,
        maxit_ou = as.integer(maxit.ou),
        maxit_in = as.integer(maxit.in),
        orthogonal = orthogonal,
        tol1 = tol1,
        tol2 = tol2,
        verbose = verbose,
        v_every = as.integer(v.every),
        Lmax = Lmax,
        c1 = c1,
        c2 = c2
    )

    return(result)
}
