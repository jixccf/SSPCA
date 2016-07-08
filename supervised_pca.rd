\name{supervised.pca}
\title{Supervised PCA}

\description{
	Performs dual supervised principal component analysis on the given 
	data matrix and kernel matrix of target variable and returns the
	results as an object of class \code{pca}.
	The kernel matrix is implemented by \link[kernlab]{dots}.
}

\usage{
	supervised.pca <- function(X, Y, 
       kernel.y = c("linear", "poly", "gaussian", "laplacian", "sigmoid", 
       "anova", "precomputed"),
       degree.y = 3, scale.y = 1, offset.y = 1, sigma.y, order.y, center.x = TRUE, ...)
	transform.supervised.pca <- function(object, data, n.components)
}

\arguments{
	\item{X}{
	a numeric matrix. It could be a $n\times p$ training data matrix, 
	where $n$ is the number of samples and $p$ is the number of features. 
	Or it could be kernel matrix of training data.
	}
	\item{Y}{
	a $n\times l$ numeric matrix of target variable.
	}
	\item{Xnew}{
	a numeric matrix in which to apply the dimensionality reduction.
	}
	\item{n.components}{
	number of compents to keep.
	}
}
\value{
	The returned \code{spca} object of \code{supervised.pca} contains the following components.
	\item{components}{
	principal axes in feature space, representing the directions of maximum variance in the data
	}
	\item{sdev}{
	the standard deviations of the principal components
	}
	\item{center}{
	the means that were substracted.
	}
	The \code{transform.supervised.pca} returns a numeric matrix.
	The returned \code{skpca} object of \code{supervised.kernel.pca} contains the following components.
	\item{components}{
	principal axes in feature space, representing the directions of maximum variance in the data
	}
	\item{eigenvalues}{
	the eigenvalues of the principal components
	}
	\item{xmatrix}{
	the input \bold{X} matrix
	}
	\item{xkernel}{
	the prameters of kernel function of X
	}
	The \code{transform.supervised.kernel.pca} returns a numeric matrix.

}
\reference{
	E. Barshan, A. Ghodsi, Z. Azimifar, M. Jahromi
	Supervised principal component analysis: Visualization, classification and regression on subspaces and submanifolds
	Pattern Recognition 44:1357--1371, 2011.
}
\examples{
	set.seed(1)
	data("iris")
	y <- with(iris["Species"], model.matrix(~ Species + 0))
	index <- c(1:nrow(iris))
	train.index <- sample(index, 105)
	train.X <- as.matrix(iris[train.index, 1:4])
	train.Y <- y[train.index, ]
	test.X <- as.matrix(iris[-train.index, 1:4])
	test.Y <- y[-train.index, ]

	# supervised pca
	iris.pca <- supervised.pca(train.X, train.Y, kernel.y = "linear")
	transform.train <- transform.supervised.pca(iris.pca, data = train.X, n.components = 2)
	transform.test <- transform.supervised.pca(iris.pca, data = test.X, n.components = 2)
	plot(transform.train[,1], transform.train[,2], pch = 20,col = iris[train.index, 5])
	points(transform.test[,1], transform.test[,2], pch = 24, col = iris[-train.index, 5])

	# compare non-supervised pca to prcomp
	supervised.pca(train.X, diag(105), kernel.y = "precomputed")
	prcomp(train.X)
	
	# kernel supervised PCA
	iris.kernel.pca <- supervised.kernel.pca(train.X, train.Y, kernel.x = "gaussian", sigma.x = 10, kernel.y = "linear")
	transform.train.kernel <- transform.supervised.kernel.pca(iris.kernel.pca, Xnew = train.X, n.components = 2, "inherited")
	transform.test.kernel <- transform.supervised.kernel.pca(iris.kernel.pca, Xnew = test.X, n.components = 2, "inherited")
	plot(transform.train.kernel[,1], transform.train.kernel[,2], pch = 20,col = iris[train.index, 5])
	points(transform.test.kernel[,1], transform.test.kernel[,2], pch = 24, col = iris[-train.index, 5])

	set.seed(1)
	data("Sonar")
	y <- with(Sonar["Class"], model.matrix(~ Class + 0))
	index <- c(1:nrow(Sonar))
	train.index <- sample(index, 140)
	train.X <- as.matrix(Sonar[train.index, 1:4])
	train.Y <- y[train.index, ]
	test.X <- as.matrix(Sonar[-train.index, 1:4])
	test.Y <- y[-train.index, ]

	# supervised pca
	Sonar.pca <- supervised.pca(train.X, train.Y, kernel.y = "linear")
	transform.train <- transform.supervised.pca(Sonar.pca, data = train.X, n.components = 2)
	transform.test <- transform.supervised.pca(Sonar.pca, data = test.X, n.components = 2)
	plot(transform.train[,1], transform.train[,2], pch = 20,col = Sonar[train.index, 61])
	points(transform.test[,1], transform.test[,2], pch = 24, col = Sonar[-train.index, 61])

	# compare non-supervised pca to prcomp
	supervised.pca(train.X, diag(140), kernel.y = "precomputed")
	prcomp(train.X)
	
	# kernel supervised PCA
	Sonar.kernel.pca <- supervised.kernel.pca(train.X, train.Y, kernel.x = "gaussian", sigma.x = 1, kernel.y = "linear")
	transform.train.kernel <- transform.supervised.kernel.pca(Sonar.kernel.pca, Xnew = train.X, n.components = 2, "inherited")
	transform.test.kernel <- transform.supervised.kernel.pca(Sonar.kernel.pca, Xnew = test.X, n.components = 2, "inherited")
	plot(transform.train.kernel[,1], transform.train.kernel[,2], pch = 20,col = Sonar[train.index, 61])
	points(transform.test.kernel[,1], transform.test.kernel[,2], pch = 24, col = Sonar[-train.index, 61])
}
