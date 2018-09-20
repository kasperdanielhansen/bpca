prcomp_fsvd1 <- function(x, n = 3, i = 1L, method = c("qr", "svd", "exact")) {
    method <- match.arg(method)
    cx <- (x - rowMeans(x)) / rowSds(x)
    s <- fsvd(cx, k = n, method = method, i = i)
    out <- list(sdev = s$d / sqrt(ncol(x) - 1),
                rotation = s$u,
                x = sweep(s$v, MARGIN = 2, STATS = s$d, FUN = "*"))
    class(out) <- "prcomp"
    out
}
    

fsvd <- function(x, k, i = 1, p = 2, method = c("qr", "svd", "exact")) {
    method <- match.arg(method)
    l <- k + p
    n <- ncol(x)
    m <- nrow(x)
    tall <- TRUE
    if (l > ncol(x) && method != "exact") {
        stop("(k+p) is greater than the number of columns. Please decrease ",
             "the value of k.")
    }
    # NOTE: In the case a WIDE matrix is provided:
    if (m < n) {
        x <- t(x) # This is potentially expensive
        n <- ncol(x)
        m <- nrow(x)
        tall <- FALSE
    }

    # Construct G to be n x l and Gaussian
    G <- matrix(rnorm(n * l, 0, 1), nrow = n, ncol = l)

    if (method == "svd") {
        # Power method:
        H <- x %*% G # m x l matrix
        for (j in seq_len(i)) {
            H <- x %*% (crossprod(x, H))
        }
        # NOTE: We use a SVD to find an othogonal basis Q:
        # TODO: Remove this commented line if not needed
        # H = FF %*% Omega %*% t(S)
        svd <- svd(crossprod(H))
        FF <- svd$u # l x l
        omega <- diag(1/sqrt(svd$d)) # l x l
        S <- H %*% FF %*% omega # m x l
        # Define the orthogonal basis:
        Q <- S[, seq_len(k), drop = FALSE] # m x k
        # TODO: Remove these commented lines if not needed
        # TT <- t(x) %*% Q # n x k
        # TT <- t(TT)
        TT <- crossprod(Q, x)
    } else if (method == "qr") {
        # NOTE: Need to create a list of H matrices
        h.list <- vector("list", i + 1L)
        h.list[[1]] <- x %*% G
        for (j in seq(2, i + 1L)) {
            h.list[[j]] <- x %*% (crossprod(x, h.list[[j - 1L]]))
        }
        H <- do.call("cbind", h.list) # n x [(1+1)l] matrix
        # QR algorithm
        Q <- qr.Q(qr(H, 0))
        # TODO: Remove these commented lines if not needed
        # TT <- t(x)%*%Q # n x [(i+1)l]
        # TT <- t(TT)
        TT <- crossprod(Q, x)
    }
    if (method == "svd" || method == "qr") {
        svd <- svd(TT)
        u <- Q %*% svd$u[,seq_len(k), drop = FALSE]
        v <- svd$v[, seq_len(k), drop = FALSE]
        d <- svd$d[seq_len(k)]
    } else {
        # Exact SVD
        svd <- svd(x, nu = k, nv = k)
        u <- svd$u[, seq_len(k), drop = FALSE]
        v <- svd$v[, seq_len(k), drop = FALSE]
        d <- svd$d[seq_len(k)]
    }
    if (!tall) {
        uu <- v
        v <- u
        u <- uu
    }
    list(u = u, v = v, d = d)
}
