.tcrossprod1 <- function(x, columns_per_block = 100,
                            BPPARAM = SerialParam()) {
    ## Input (x) is a DelayedArray matrix of dimension
    ##    dim(x) = (m,n)
    ## This function computes
    ##   tcrossprod(x) = x %*% t(x)
    ## and assumes that we can have an (m,m) matrix in memory.
    ##
    ## It computes the crossproduct by
    ##    tcrossprod = sum_b tcrossprod(x[b])
    ## where x[b] is a block of columns.
    
    iterator <- function(x, grid = NULL) {
        ## grid <- DelayedArray:::.normarg_grid(grid, x)
        b <- 0L
        function() {
            if (b == length(grid))
                return(NULL)
            b <<- b + 1L
            viewport <- grid[[b]]
            read_block(x, viewport)
            ## if (!is.array(block))
            ##     block <- DelayedArray:::.as_array_or_matrix(block)
            ## attr(block, "from_grid") <- grid
            ## attr(block, "block_id") <- b
            ## block
        }
    }
    workhorse <- function(block, ...) {
        tcrossprod(as.matrix(block))
    }
    reducer <- function(x, y) {
        x + y
    }
    grid <- RegularArrayGrid(dim(x), c(nrow(x), columns_per_block))
    out <- bpiterate(ITER = iterator(x, grid),
                     FUN = workhorse,
                     REDUCE = reducer,
                     BPPARAM = BPPARAM)
}


## scale : sweep(x, 2L, scale, "/", check.margin = FALSE
## center: sweep(x, 2L, center, check.margin = FALSE
## But since this is tcrossprod and I assume row-based means, sd
## I can just do (x - center) / scale
