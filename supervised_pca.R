require(kernlab)
supervised.pca <- function(X, Y, 
                           kernel.y = c("linear", "poly", "gaussian", "laplacian", "sigmoid", "anova", "precomputed"),
                           degree.y = 3, scale.y = 1, offset.y = 1, sigma.y, order.y, center.x = TRUE, ...)
{  
    if (nrow(X) != nrow(Y))
        stop("X and Y must be the same number of samples")
    
    n.row <- nrow(X)
    n.col <- ncol(X)
    if (center.x == TRUE) {
        # Comparison of different ways to different ways to mean-center a matrix
        # http://gastonsanchez.com/how-to/2014/01/15/Center-data-in-R/
        center <- colMeans(X)
        X.adj <- X - rep(center, rep.int(n.row, n.col))
    } else {
        center <- rep(0, n.col)
        X.adj <- X
    }
    
    if (kernel.y == "linear") {
        L <- kernelMatrix(kernel = vanilladot(), Y, ...)
    } else if (kernel.y == "poly") {
        L <- kernelMatrix(kernel = polydot(degree = degree.y, scale = scale.y, offset = offset.y), Y, ...)
    } else if (kernel.y == "gaussian") {
        L <- kernelMatrix(kernel = rbfdot(sigma = sigma.y), Y, ...)
    } else if (kernel.y == "laplacian") {
        L <- kernelMatrix(kernel = laplacedot(sigma = sigma.y), Y, ...)
    } else if (kernel.y == "sigmoid") {
        L <- kernelMatrix(kernel = tanhdot(scale = scale.y, offset = offset.y), Y, ...)
    } else if (kernel.y == "anova") {
        L <- kernelMatrix(kernel = anovadot(sigma = sigma.y, degree = degree.y), Y, ...)
    } else if (kernel.y == "precomputed") {
        L <- Y
    }
    
    H <- diag(1, n.row, n.row) - matrix(1, nrow = n.row, ncol = n.row) / n.row
    eigen <- eigen(t(X.adj) %*% H %*% L %*% H %*% X.adj, symmetric = TRUE)
    sdev <- eigen$values/sqrt(max(1, n.row - 1))
    components <- eigen$vectors
    dimnames(components) <- list(colnames(X.adj), paste0("PC", seq_len(ncol(components))))
    pca <- list(components = components, sdev = sdev, center = center)
    class(pca) <- "spca"
    return(pca)
}

supervised.dual.pca <- function(X, Y, 
                           kernel.y = c("linear", "poly", "gaussian", "laplacian", "sigmoid", "anova", "precomputed"),
                           degree.y = 3, scale.y = 1, offset.y = 1, sigma.y, order.y, center.x = TRUE, ...)
{  
    if (nrow(X) != nrow(Y))
        stop("X and Y must be the same number of samples")
    
    n.row <- nrow(X)
    n.col <- ncol(X)
    if (center.x == TRUE) {
        # Comparison of different ways to different ways to mean-center a matrix
        # http://gastonsanchez.com/how-to/2014/01/15/Center-data-in-R/
        center <- colMeans(X)
        X.adj <- X - rep(center, rep.int(n.row, n.col))
    } else {
        center <- rep(0, n.col)
        X.adj <- X
    }
    
    if (kernel.y == "linear") {
        L <- kernelMatrix(kernel = vanilladot(), Y, ...)
    } else if (kernel.y == "poly") {
        L <- kernelMatrix(kernel = polydot(degree = degree.y, scale = scale.y, offset = offset.y), Y, ...)
    } else if (kernel.y == "gaussian") {
        L <- kernelMatrix(kernel = rbfdot(sigma = sigma.y), Y, ...)
    } else if (kernel.y == "laplacian") {
        L <- kernelMatrix(kernel = laplacedot(sigma = sigma.y), Y, ...)
    } else if (kernel.y == "sigmoid") {
        L <- kernelMatrix(kernel = tanhdot(scale = scale.y, offset = offset.y), Y, ...)
    } else if (kernel.y == "anova") {
        L <- kernelMatrix(kernel = anovadot(sigma = sigma.y, degree = degree.y), Y, ...)
    } else if (kernel.y == "precomputed") {
        L <- Y
    }
    
    H <- diag(1, n.row, n.row) - matrix(1, nrow = n.row, ncol = n.row) / n.row
    # decomposition of L
    svd.l <- svd(L)
    d.l <- diag(svd.l$d)
    l <- sqrt(d.l) %*% t(svd.l$v)
    # decomposition of psi
    svd <- svd(l %*% H %*% X.adj, nu = 0)
    components <- svd$v
    sdev <- svd$d / sqrt(max(1, n.row - 1))
    dimnames(components) <- list(colnames(X.adj), paste0("PC", seq_len(ncol(components))))
    pca <- list(components = components, sdev = sdev, center = center)
    class(pca) <- "spca"
    return(pca)
}

transform.supervised.pca <- function(object, Xnew, n.components) {
    if (ncol(Xnew)!=length(object$center)) {
        stop("data should have the same number of features as the training data")
    }
    if (n.components <= 1 || n.components > length(object$center))
        stop(paste("n.components must be between 1 and the number of features", length(object$center)))
    n.row <- nrow(Xnew)
    n.col <- ncol(Xnew)
    Xnew.adj <- Xnew - rep(object$center, rep.int(n.row, n.col))
    V <- object$components[, 1:n.components]
    newdata <- Xnew.adj %*% V
    return(newdata)
}

print.supervised.pca <- function(x, print.x = FALSE, ...) {
    cat("Standard deviations:\n")
    print(x$sdev)
    cat("\nPrincipal Components:\n")
    print(x$components)
    invisible(x)
}
