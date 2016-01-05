data 

design <- matrix(0, nrow=36, ncol=36)
design

X <- read.csv("~/Desktop/designMatrix.csv", header=FALSE)


dig <- diag(1, nrow=36, ncol=36)

nm <- c("u", paste("g", 1:9, sep=""), paste("s", 1:36, sep=""))
X <- cbind(X, dig)
X <- as.matrix(X)
fit <- lm(sample$Yield ~ X)
names(X) <- nm
