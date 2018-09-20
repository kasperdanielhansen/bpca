library(bpca)
m <- 1000     # Total number of rows
n <- 100000    # Total number of columns
X <- DelayedArray(matrix(rnorm(m*n), ncol=n))
file.remove("bpca_test/X_disk.h5")
X_disk <- writeHDF5Array(X, chunkdim = c(nrow(X), 100),
                         filepath = "bpca_test/X_disk.h5")

# Regular crossprod
system.time({
  ref <- tcrossprod(as.matrix(X))
})

## Test 2
# Crossprod on DelayedArray of in memory matrix
system.time({
  out <- tcrossprod1(X)
})
system.time({
  out <- tcrossprod1(X_disk)
})



system.time({
  out <- tcrossprod1(X, columns_per_block = 1000)
})
system.time({
  out <- tcrossprod1(X_disk, columns_per_block = 1000)
})



system.time({
  out <- tcrossprod1(X, columns_per_block = ncol(X))
})

system.time({
  out <- tcrossprod1(X_disk, columns_per_block = ncol(X))
})


library(pryr)
library(TENxBrainData)
tenx <- TENxBrainData()
counts_in_memory <- as.matrix(counts(tenx[,5000 + 1:5000]))
counts_in_memory <- as.matrix(counts(tenx[,1:10000]))
counts_in_memory0 <- counts_in_memory[rowSums(counts_in_memory) > 0,]
object_size(counts_in_memory)
rm(counts_in_memory)
object_size(counts_in_memory0)
dim(counts_in_memory0)
sum(counts_in_memory0 == 0) / prod(dim(counts_in_memory0))

system.time({
    out <- tcrossprod(counts_in_memory0)
})
object_size(out)
system.time({
    out <- counts_in_memory0 %*% t(counts_in_memory0)
})
object_size(out)
library(help = Matrix)
system.time({
    sparse_counts_in_memory0 <- Matrix(counts_in_memory0, sparse = TRUE)
    out <- tcrossprod(sparse_counts_in_memory0)
})

rsums <- rowSums(counts_in_memory0)
rsds <- rowSds(counts_in_memory0)
scaled_counts_in_memory0 <- (counts_in_memory0 - rsums) / rsds
out_scaled <- tcrossprod(scaled_counts_in_memory0)

sparse_counts_in_memory0 <- Matrix(counts_in_memory0, sparse = TRUE)
