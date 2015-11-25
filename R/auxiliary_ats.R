# Auxiliary functions (adopted and modified from package "sensitivity"),
# available only locally:

ind.rep <- function (i, p) 
{
  (1:(p + 1)) + (i - 1) * (p + 1)
}

random.oat <- function(p, r, binf = rep(0, p), bsup = rep(1, p),
                       nl, design.step){
  B <- matrix(-1, nrow = p + 1, ncol = p)
  B[lower.tri(B)] <- 1
  delta <- design.step/(nl - 1)
  X <- matrix(nrow = r * (p + 1), ncol = p)
  for (j in 1:r) {
    D <- diag(sample(c(-1, 1), size = p, replace = TRUE))
    perm <- sample(p)
    P <- matrix(0, nrow = p, ncol = p)
    for (i in 1:p) {
      P[i, perm[i]] <- 1
    }
    x.base <- matrix(nrow = p + 1, ncol = p)
    for (i in 1:p) {
      x.base[, i] <- ((sample(nl[i] - design.step[i], size = 1) - 
                         1)/(nl[i] - 1))
    }
    X[ind.rep(j, p), ] <- 0.5 * (B %*% P %*% D + 1) %*% diag(delta) + 
      x.base
  }
  for (i in 1:p) {
    X[, i] <- X[, i] * (bsup[i] - binf[i]) + binf[i]
  }
  return(X)
}

################################################################################

haussdorf.distance <- function (x, set1, set2) 
{
  n1 <- length(set1)
  n2 <- length(set2)
  d <- matrix(nrow = n1, ncol = n2)
  for (i1 in 1:n1) {
    for (i2 in 1:n2) {
      d[i1, i2] <- sqrt(sum((x[set1[i1], ] - x[set2[i2], 
                                               ])^2))
    }
  }
  return(max(mean(apply(d, 1, min)), mean(apply(d, 2, min))))
}

################################################################################

kennard.stone <- function (dist.matrix, n) 
{
  out <- numeric(n)
  out[1] <- 1
  for (i in 2:n) {
    tmp <- dist.matrix[out, -out, drop = FALSE]
    out[i] <- (1:nrow(dist.matrix))[-out][which.max(apply(tmp, 
                                                          2, min))]
  }
  return(out)
}

################################################################################

morris.maximin <- function (x, r) 
{
  p <- ncol(x)
  R <- nrow(x)/(p + 1)
  d <- matrix(0, nrow = R, ncol = R)
  for (i in 1:(R - 1)) {
    for (j in (i + 1):R) {
      d[i, j] <- d[j, i] <- haussdorf.distance(x, ind.rep(i, 
                                                          p), ind.rep(j, p))
    }
  }
  kennard.stone(d, r)
}

#########################################################################

morris_ats <- function (model_matrix, pars, p, r, design, 
                        binf, bsup, scale, ...)
{
  r.max <- r
  if (!"type" %in% names(design)) {
    design$type <- "oat"
    warning("argument 'design$type' not found, set at 'oat'")
  }
  if (design$type == "oat") {
    if (!"levels" %in% names(design)) {
      stop("argument 'design$levels' not found")
    }
    nl <- design$levels
    if (length(nl) == 1) 
      nl <- rep(nl, p)
    if ("grid.jump" %in% names(design)) {
      jump <- design$grid.jump
      if (round(jump, 0) != jump) 
        stop("grid.jump must be integer")
      if (length(jump) == 1) 
        jump <- rep(jump, p)
    }
    else {
      jump <- rep(1, p)
      warning("argument 'design$grid.jump' not found, set at 1")
    }
  }
  else if (design$type == "simplex") {
    h <- design$scale.factor
  }
  else {
    stop("invalid argument design$type, waiting for \"oat\" or \"simplex\"")
  }
  if (length(binf) == 1) 
    binf <- rep(binf, p)
  if (length(bsup) == 1) 
    bsup <- rep(bsup, p)
  if (design$type == "oat") {
    X <- random.oat(p, r.max, binf, bsup, nl, jump)
  }
  else if (design$type == "simplex") {
    X <- random.simplexes(p, r.max, binf, bsup, h)
  }
  X.unique <- array(t(X), dim = c(p, p + 1, r.max))
  X.unique <- unique(X.unique, MARGIN = 3)
  X <- matrix(X.unique, ncol = p, byrow = TRUE)
  colnames(X) <- pars
  r.unique <- nrow(X)/(p + 1)
  if (r.unique < r.max) {
    warning(paste("keeping", r.unique, "repetitions out of", r.max))
  }
  r.max <- r.unique
  if (r < r.max) {
    ind <- morris.maximin(X, r)
    X <- X[sapply(ind, function(i) ind.rep(i, p)), ]
  }
  x <- list(model = model_matrix, factors = p, r = r, design = design, 
            binf = binf, bsup = bsup, scale = scale, X = X)
  response_ats(x, ...)
  tell_ats(x)
  return(x)
}

################################################################################

response_ats <- function (x, ...)
{
  id <- deparse(substitute(x))
  if (class(x$model) == "function") {
    Y <- x$model(x$X, ...)
  }
  else if (TRUE %in% (paste("predict.", class(x$model), sep = "") %in% 
                      methods(predict))) {
    Y <- predict(x$model, x$X, ...)
  }
  else {
    stop("The model isn't a function or does not have a predict method")
  }
  if (!class(Y) %in% c("numeric", "matrix")) {
    Y <- as.numeric(Y)
    warning("Conversion of the response to numeric")
  }
  x$Y <- Y
  assign(id, x, parent.frame())
}

################################################################################

ee.oat_ats <- function (X, Y) 
{
  p <- ncol(X)
  r <- nrow(X)/(p + 1)
  ee <- vector(mode = "list", length = r)
  for (i in 1:r) {
    j <- ind.rep(i, p)
    j1 <- j[1:p]
    j2 <- j[2:(p + 1)]
    ee[[i]] <- solve(X[j2, ] - X[j1, ], Y[j2,] - Y[j1,])
    t_num <- ncol(ee[[i]])
    colnames(ee[[i]]) <- paste0("traj", i, "_t", 1:t_num)
  }
  t_num <- ncol(ee[[1]])
  ee_matrix <- t(do.call(cbind, ee))
  nrow_ee_matrix <- nrow(ee_matrix)
  ee_by_time <- lapply(1:t_num, function(i){
    ee_matrix[seq(i, nrow_ee_matrix, r + 1), ]})
  names(ee_by_time) <- paste0("timepoint", 1:t_num)
  return(ee_by_time)
}

################################################################################

tell_ats <- function (x, Y = NULL, ...) 
{
  id <- deparse(substitute(x))
  if (!is.null(Y)) {
    x$Y <- Y
  }
  else if (is.null(x$Y)) {
    stop("Y not found")
  }
  X <- x$X
  Y <- x$Y
  if (x$scale) {
    Binf <- matrix(x$binf, nrow = nrow(X), ncol = length(x$binf), 
                   byrow = TRUE)
    Bsup <- matrix(x$bsup, nrow = nrow(X), ncol = length(x$bsup), 
                   byrow = TRUE)
    X <- (X - Binf)/(Bsup - Binf)
  }
  if (x$design$type == "oat") {
    x$ee <- ee.oat_ats(X, Y)
  }
  else if (x$design$type == "simplex") {
    x$ee <- ee.simplex(X, Y)
  }
  assign(id, x, parent.frame())
}

################################################################################

# out_ats <- function(x, times, pars){
#   mu <- lapply(x$ee, colMeans)
#   mu.star <- lapply(x$ee, abs)
#   mu.star <- lapply(mu.star, colMeans)
#   # mu.star <- lapply(x$ee, function(M){
#   #   apply(M, 2, function(x) mean(abs(x)))
#   # })
#   sigma <- lapply(x$ee, function(M){
#     apply(M, 2, sd)
#   })
#   out <- mapply(c, mu, mu.star, sigma, SIMPLIFY = TRUE)
#   out <- rbind(times, out)
#   rownames(out) <- c("time", paste0("mu_", pars), paste0("mu.star_", pars),
#                      paste0("sigma_", pars))
#   return(out)
# }
