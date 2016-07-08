require(kernlab)
require(geigen)

supervised.kernel.pca <- function(X, Y, 
      kernel.x = c("linear", "poly", "gaussian", "laplacian", 
                   "sigmoid", "anova", "precomputed"),
      degree.x = NA, scale.x = NA, offset.x = NA, sigma.x = NA, order.x = NA, 
      kernel.y = c("linear", "poly", "gaussian", "laplacian", 
                   "sigmoid", "anova", "precomputed"),
      degree.y = NA, scale.y = NA, offset.y = NA, sigma.y = NA, order.y = NA)
{  
    if (nrow(X) != nrow(Y))
        stop("X and Y must be the same number of samples")
    
    if (kernel.x == "linear") {
        K <- kernelMatrix(kernel = vanilladot(), X)
    } else if (kernel.x == "poly") {
        K <- kernelMatrix(
            kernel = polydot(degree = degree.x, scale = scale.x, 
                             offset = offset.x), 
            X)
    } else if (kernel.x == "gaussian") {
        K <- kernelMatrix(kernel = rbfdot(sigma = sigma.x), X)
    } else if (kernel.x == "laplacian") {
        K <- kernelMatrix(kernel = laplacedot(sigma = sigma.x), X)
    } else if (kernel.x == "sigmoid") {
        K <- kernelMatrix(
            kernel = tanhdot(scale = scale.x, offset = offset.x), X)
    } else if (kernel.x == "anova") {
        K <- kernelMatrix(
            kernel = anovadot(sigma = sigma.x, degree = degree.x), X)
    } else if (kernel.x == "precomputed") {
        K <- X
    }
    if (kernel.y == "linear") {
        L <- kernelMatrix(kernel = vanilladot(), Y)
    } else if (kernel.y == "poly") {
        L <- kernelMatrix(
            kernel = polydot(degree = degree.y, scale = scale.y, 
                             offset = offset.y), 
            Y)
    } else if (kernel.y == "gaussian") {
        L <- kernelMatrix(kernel = rbfdot(sigma = sigma.y), Y)
    } else if (kernel.y == "laplacian") {
        L <- kernelMatrix(kernel = laplacedot(sigma = sigma.y), Y)
    } else if (kernel.y == "sigmoid") {
        L <- kernelMatrix(
            kernel = tanhdot(scale = scale.y, offset = offset.y), Y)
    } else if (kernel.y == "anova") {
        L <- kernelMatrix(
            kernel = anovadot(sigma = sigma.y, degree = degree.y), Y)
    } else if (kernel.y == "precomputed") {
        L <- Y
    }
    
    n.row <- nrow(K)
    H <- diag(1, n.row, n.row) - matrix(1, nrow = n.row, ncol = n.row) / n.row
    # decomposition of KHLHL
    eigen <- geigen(K %*% H %*% L %*% H %*% K, K)
    browser()
    values <- eigen$values[zapsmall(eigen$values) != 0]
    components <- matrix(eigen$vectors[, zapsmall(eigen$values) != 0], nrow = n.row)
    dimnames(components) <- list(colnames(K), 
                                 paste0("PC", seq_len(ncol(components))))
    pca <- list(components = components, eigenvalues = values, 
                xmatrix = X, 
                xkernel = list(kernel.x = kernel.x, degree.x = degree.x, 
                               scale.x = scale.x, offset.x = offset.x, 
                               sigma.x = sigma.x, order.x = order.x))
    class(pca) <- "skpca"
    return(pca)
}

transform.supervised.kernel.pca <- function(object, Xnew, n.components, 
                        kernel = c("inherited", "precomputed")) 
{
    if (n.components <= 1 || n.components > length(object$eigenvalues))
        stop(paste("n.components must be between 1 and the number of features", 
                   length(object$eigenvalues)))
    if (kernel == "inherited") {
        if (object$xkernel$kernel.x == "linear") {
            K <- kernelMatrix(kernel = vanilladot(), object$xmatrix, Xnew)
        } else if (object$xkernel$kernel.x == "poly") {
            K <- kernelMatrix(
                kernel = polydot(degree = object$xkernel$degree.x, 
                                 scale = object$xkernel$scale.x, 
                                 offset = object$xkernel$offset.x), 
                object$xmatrix, Xnew)
        } else if (object$xkernel$kernel.x == "gaussian") {
            K <- kernelMatrix(
                kernel = rbfdot(sigma = object$xkernel$sigma.x), 
                object$xmatrix, Xnew)
        } else if (object$xkernel$kernel.x == "laplacian") {
            K <- kernelMatrix(
                kernel = laplacedot(sigma = object$xkernel$sigma.x), 
                object$xmatrix, Xnew)
        } else if (object$xkernel$kernel.x == "sigmoid") {
            K <- kernelMatrix(
                kernel = tanhdot(scale = object$xkernel$scale.x, 
                                 offset = object$xkernel$offset.x), 
                object$xmatrix, Xnew)
        } else if (object$xkernel$kernel.x == "anova") {
            K <- kernelMatrix(
                kernel = anovadot(sigma = object$xkernel$sigma.x, 
                                  degree = object$xkernel$degree.x), 
                object$xmatrix, Xnew)
        } else if (object$xkernel$kernel.x == "precomputed") {
            K <- Xnew
        }
    }
    V <- object$components[, 1:n.components]
    newdata <- t(K) %*% V
    return(newdata)
}

print.supervised.kernel.pca <- function(x, print.x = FALSE, ...) {
    cat("Eigenvalues:\n")
    print(x$eigenvalues)
    invisible(x)
}