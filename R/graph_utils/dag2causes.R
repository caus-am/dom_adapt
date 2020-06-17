dag2causes <- function( G ) {
  # Calculates the causes-matrix for a DAG
  # tc[i,j] == 0 if i is not a causal ancestor of j
  # tc[i,j] == 1 if i is a causal ancestor of j
  # (this is basically the transitive closure, without the diagonal)

  library('RBGL')
  tc <- transitive.closure(as(G,'graphNEL'))
  tc <- t(wgtMatrix(tc))
  diag(tc) <- 0
  tc
}

