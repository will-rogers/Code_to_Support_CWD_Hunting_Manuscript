M <- out$matrix #gathered from the altered cwd_det_model output
eigen.analysis(M, zero = T)
stable.stage(M)
Ma <- M[1:12,1:12]
eigen.analysis(Ma)
dim(M)
Ma <- M[13:24,13:24]

sensitivity <- function(A) {
  d <- eigen(A)$values   # eigen values
  w <- eigen(A)$vectors  # right eigen vectors
  v <- Conj(solve(w))    # complex conjugate of left eigen vectors
  # output of eigenvalues is decreasingly sorted.
  v[1,] %*% t(w[,1])
}

sensitivity(M)
