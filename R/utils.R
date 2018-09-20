describeMatrix <- function(mat) {
    cat(sprintf("class: %s\n", class(mat)[1]))
    cat(sprintf("dimension: %s x %s\n", nrow(mat), ncol(mat)))
    cat(sprintf("percent zero: %s\n", round(sum(mat == 0) / prod(dim(mat)) * 100, 1)))
    zero_rows <- sum(rowSums(mat != 0) == 0)
    cat(sprintf("number of zero rows: %s (percentage: %s)\n", zero_rows, round(zero_rows/nrow(mat) * 100, 1)))
    zero_cols <- sum(colSums(mat != 0) == 0)
    cat(sprintf("number of zero columns: %s (percentage: %s)\n", zero_cols, round(zero_cols/ncol(mat) * 100, 1)))
    invisible(NULL)
}

benchmark <- function(...) {
    exprs <- as.list(match.call(expand.dots = FALSE)$...)
    exprsnm <- sapply(exprs, function(e) paste(deparse(e), collapse = " "))
    timings <- sapply(exprs, function(ee) {
        time <- system.time( eval(ee) )
        time
    })
    cbind(data.frame(exprs = exprsnm), t(timings)[,1:3])
}
