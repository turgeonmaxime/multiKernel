computeLeastSq <- function(Y, Kmat, alpha, Z, B) {
  0.5 * norm(Y - Kmat %*% alpha - Z %*% B, type = "F")
}

# Stripped down version of caret::createFolds
createFold <- function(Y, K = 5) {
  cuts <- floor(nrow(Y)/K)
  if (cuts < 2) 
    cuts <- 2
  if (cuts > 5) 
    cuts <- 5
  breaks <- unique(quantile(Y[,1], probs = seq(0, 1, length = cuts)))
  y <- cut(Y[,1], breaks, include.lowest = TRUE)
  if (K < length(y)) {
    y <- factor(as.character(y))
    numInClass <- table(y)
    foldVector <- vector(mode = "integer", length(y))
    for (i in 1:length(numInClass)) {
      min_reps <- numInClass[i]%/%K
      if (min_reps > 0) {
        spares <- numInClass[i]%%K
        seqVector <- rep(1:K, min_reps)
        if (spares > 0) 
          seqVector <- c(seqVector, sample(1:K, spares))
        foldVector[which(y == names(numInClass)[i])] <- sample(seqVector)
      }
      else {
        foldVector[which(y == names(numInClass)[i])] <- sample(1:K, 
                                                               size = numInClass[i])
      }
    }
  }
  out <- split(seq(along = Y[,1]), foldVector)
  names(out) <- paste("Fold", gsub(" ", "0", format(seq(along = out))),
                      sep = "")
  return(out)
}