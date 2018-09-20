# bpca


## Cross-product

Assume single-cell data with many more columns than rows. Input is
  x : m (rows, fixed) x n (columns, number of cells)

Want to compute
  tcrossprod(x) = x %*% t(x)
result has dimension
  m x m


```{r}
m <- 500     # Total number of rows
n <- 1000    # Total number of columns
X <- DelayedArray(matrix(rnorm(m*n), ncol=n))

# Regular crossprod
system.time({
  ref <- crossprod(as.matrix(X))
})

## Test 2
# Crossprod on DelayedArray of in memory matrix
system.time({
  out <- tcrossprod1
})


```
