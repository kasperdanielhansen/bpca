setClass("DoubleChunk",
         representation(rowAccess = "HDF5Matrix",
                        colAccess = "HDF5Matrix",
                        rowSds = "numeric",
                        rowMeans = "numeric"))
dim.DoubleChunk = function(x) dim(x@rowAccess)
ncol.DoubleChunk = function(x) ncol(x@rowAccess)
nrow.DoubleChunk = function(x) nrow(x@rowAccess)

setMethod("%*%", signature(x = "DoubleChunk", y = "numeric"),
          function(x, y) {
    grid <- chunkGrid(x@rowAccess)
    blockResults <- lapply(seq_along(grid), function(ii) {
        block <- read_block(x@rowAccess, viewport = grid[[ii]])
        block %*% y
    }) # could be written as blockApply
    matrix((unlist(blockResults) - sum(y) * x@rowMeans) / x@rowSds, ncol = 1)
})
setMethod("%*%", signature(x = "numeric", y = "DoubleChunk"),
          function(x, y) {
    grid <- chunkGrid(y@colAccess)
    xs <- x / y@rowSds
    blockResults <- lapply(seq_along(grid), function(ii) {
        block <- read_block(y@colAccess, viewport = grid[[ii]])
        xs %*% block
    }) # could be written as blockApply 
    matrix(unlist(blockResults) - sum(xs * y@rowMeans), nrow = 1)
})


setClass("SingleChunk",
         representation(access = "HDF5Matrix",
                        rowSds = "numeric",
                        rowMeans = "numeric"))
dim.SingleChunk = function(x) dim(x@access)
ncol.SingleChunk = function(x) ncol(x@access)
nrow.SingleChunk = function(x) nrow(x@access)

setMethod("%*%", signature(x = "SingleChunk", y = "numeric"),
          function(x, y) {
    grid <- rowGrid(x@access)
    cat(sprintf("Ay (%s)\n", length(grid)))
    blockResults <- lapply(seq_along(grid), function(ii) {
        block <- read_block(x@access, viewport = grid[[ii]])
        as.vector(block %*% y)
    }) # could be written as blockApply
    matrix((unlist(blockResults) - sum(y) * x@rowMeans) / x@rowSds, ncol = 1)
})
setMethod("%*%", signature(x = "numeric", y = "SingleChunk"),
          function(x, y) {
    grid <- colGrid(y@access)
    cat(sprintf("yA (%s)\n", length(grid)))
    xs <- x / y@rowSds
    blockResults <- lapply(seq_along(grid), function(ii) {
        block <- read_block(y@access, viewport = grid[[ii]])
        as.vector(xs %*% block)
    }) # could be written as blockApply 
    matrix(unlist(blockResults) - sum(xs * y@rowMeans), nrow = 1)
})
