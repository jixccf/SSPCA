set.seed(1)
require(mlbench)

# iris --------------------------------------------------------------------------------
data("iris")
y <- with(iris["Species"], model.matrix(~ Species + 0))
index <- c(1:nrow(iris))
train.index <- sample(index, 105)
train.X <- as.matrix(iris[train.index, 1:4])
train.Y <- y[train.index, ]
test.X <- as.matrix(iris[-train.index, 1:4])
test.Y <- y[-train.index, ]

# supervised pca (eigen)
iris.pca <- supervised.pca(train.X, train.Y, kernel.y = "linear")
transform.train <- transform.supervised.pca(iris.pca, Xnew = train.X, n.components = 2)
transform.test <- transform.supervised.pca(iris.pca, Xnew = test.X, n.components = 2)
plot(transform.train[,1], transform.train[,2], pch = 20,col = iris[train.index, 5])
points(transform.test[,1], transform.test[,2], pch = 24, col = iris[-train.index, 5])

# supervised dual pca (svd)
iris.dual.pca <- supervised.dual.pca(train.X, train.Y, kernel.y = "linear")
transform.train <- transform.supervised.pca(iris.pca, Xnew = train.X, n.components = 2)
transform.test <- transform.supervised.pca(iris.pca, Xnew = test.X, n.components = 2)
plot(transform.train[,1], transform.train[,2], pch = 20,col = iris[train.index, 5])
points(transform.test[,1], transform.test[,2], pch = 24, col = iris[-train.index, 5])

# compare to prcomp
supervised.pca(train.X, diag(105), kernel.y = "precomputed")
supervised.dual.pca(train.X, diag(105), kernel.y = "precomputed")
prcomp(train.X)

# kernel supervised PCA
iris.kernel.pca <- supervised.kernel.pca(train.X, train.Y, kernel.x = "gaussian", sigma.x = 20, kernel.y = "linear")
transform.train.kernel <- transform.supervised.kernel.pca(iris.kernel.pca, Xnew = train.X, n.components = 2, "inherited")
transform.test.kernel <- transform.supervised.kernel.pca(iris.kernel.pca, Xnew = test.X, n.components = 2, "inherited")
plot(transform.train.kernel[,1], transform.train.kernel[,2], pch = 20,col = iris[train.index, 5])
points(transform.test.kernel[,1], transform.test.kernel[,2], pch = 24, col = iris[-train.index, 5])


# Sonar -----------------------------------------------------------------------------------------
data("Sonar")
y <- with(Sonar["Class"], model.matrix(~ Class + 0))
index <- c(1:nrow(Sonar))
train.index <- sample(index, 142)
train.X <- as.matrix(Sonar[train.index, 1:60])
train.Y <- y[train.index, ]
test.X <- as.matrix(Sonar[-train.index, 1:60])
test.Y <- y[-train.index, ]

# supervised pca (eigen)
Sonar.pca <- supervised.pca(train.X, train.Y, kernel.y = "linear")
transform.train <- transform.supervised.pca(Sonar.pca, Xnew = train.X, n.components = 2)
transform.test <- transform.supervised.pca(Sonar.pca, Xnew = test.X, n.components = 2)
plot(transform.train[,1], transform.train[,2], pch = 20,col = Sonar[train.index, 61])
points(transform.test[,1], transform.test[,2], pch = 24, col = Sonar[-train.index, 61])

# supervised dual pca (svd)
Sonar.dual.pca <- supervised.dual.pca(train.X, train.Y, kernel.y = "linear")
transform.train <- transform.supervised.pca(Sonar.dual.pca, Xnew = train.X, n.components = 2)
transform.test <- transform.supervised.pca(Sonar.dual.pca, Xnew = test.X, n.components = 2)
plot(transform.train[,1], transform.train[,2], pch = 20,col = Sonar[train.index, 61])
points(transform.test[,1], transform.test[,2], pch = 24, col = Sonar[-train.index, 61])

# kernel supervised PCA
Sonar.kernel.pca <- supervised.kernel.pca(train.X, train.Y, kernel.x = "gaussian", sigma.x = 5, kernel.y = "linear")
transform.train.kernel <- transform.supervised.kernel.pca(Sonar.kernel.pca, Xnew = train.X, n.components = 2, "inherited")
transform.test.kernel <- transform.supervised.kernel.pca(Sonar.kernel.pca, Xnew = test.X, n.components = 2, "inherited")
plot(transform.train.kernel[,1], transform.train.kernel[,2], pch = 20,col = Sonar[train.index, 61])
points(transform.test.kernel[,1], transform.test.kernel[,2], pch = 24, col = Sonar[-train.index, 61])
