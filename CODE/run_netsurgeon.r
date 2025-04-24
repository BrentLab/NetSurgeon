args <- commandArgs(trailingOnly = TRUE)
regulatoryNetworkMatrixFile <- toString(args[1])
startGoalStateDEVectorFile <- toString(args[2])
regulatorGeneNamesFile <- toString(args[3])
targetGeneNamesFile <- toString(args[4])
outputDirectory <- toString(args[5])
overexpressionFlag  <- as.integer(args[6])
outputFileName <- "intervention_scores.txt"

regulatoryNetworkMatrix <- as.matrix(read.table(regulatoryNetworkMatrixFile))
startGoalStateDEVector <- scan(startGoalStateDEVectorFile, what=numeric())
regulatorGeneNames <- scan(regulatorGeneNamesFile, what=character())
targetGeneNames <- scan(targetGeneNamesFile, what=character())

regulatoryNetworkMatrix <- regulatoryNetworkMatrix / max(abs(regulatoryNetworkMatrix))

if(overexpressionFlag == 1) {
	regulatoryNetworkMatrix <- regulatoryNetworkMatrix * -1
}

startGoalStateDEVector <- startGoalStateDEVector / max(abs(startGoalStateDEVector))

## By default this function adds autoregulation, but autoregulation can be removed by setting score to 0
modifyAutoRegulation <- function(network, tfs, orfs, score=max(abs(network))) {
  selected <- cbind(match(intersect(tfs, orfs), tfs), match(intersect(tfs, orfs), orfs))
  network[selected] <- score
  network
}

measureDirectEnrichment <- function(M, tfs, orfs, deVect, binStart, binStop, binStep) {
  M.sign <- sign(M)
  M <- abs(M)
  bins <- seq(binStart, binStop, by=binStep)
  cutoffs <- quantile(M, probs=1-c(bins/length(M)))
  tf.inds <- match(tfs, orfs)
  right_direction <- matrix(nrow=length(tfs), ncol=length(bins), data=1)
  wrong_direction <- matrix(nrow=length(tfs), ncol=length(bins), data=1)

  white <- sum(abs(sign(deVect))) ## Num DE Genes
  scaling_factor <- (white / sum(abs(deVect)))
  black <- length(deVect) - white
  de_pred_signs <- t(t(M.sign) * sign(deVect))

  for (i in 1:length(cutoffs)) {
    cutoff <- cutoffs[i]
    M.cut <- M>=cutoff
    target.count <- apply(M.cut,1,sum)

    for(j in 1:length(tfs)) {

      target_inds <- which(M.cut[j,])
      
			draws <- (sum(abs(deVect[target_inds])) * scaling_factor) + target.count[j] - sum(abs(sign(deVect[target_inds])))

      ## Find the amount of changing regulated genes that are regulated in the right direction
      target_pos_inds <- target_inds[which(de_pred_signs[j,target_inds] == 1)]
      successes <- sum(abs(deVect[target_pos_inds])) * scaling_factor
			right_direction[j,i] <- phyper(successes - 1, white, black, draws, lower.tail=FALSE)

      ## Find the amount of changing regulated genes that are regulated in the wrong direction
      target_neg_inds <- target_inds[which(de_pred_signs[j,target_inds] == -1)]
      successes <- sum(abs(deVect[target_neg_inds])) * scaling_factor
      wrong_direction[j,i] <- phyper(successes - 1, white, black, draws, lower.tail=FALSE)
    }
  }

  right_direction_enrichment_max <- -log10(apply(right_direction, 1, min) + 1e-100)

  results <- cbind(right_direction_enrichment_max)
  row.names(results) <- tfs
  colnames(results) <- c("intervention_score")
  results
}

binStart <- 4000
binStop <- 40000
binStep <- 4000
interventionScores <- measureDirectEnrichment(modifyAutoRegulation(regulatoryNetworkMatrix, regulatorGeneNames, targetGeneNames, 0), regulatorGeneNames, targetGeneNames, startGoalStateDEVector, binStart, binStop, binStep)

interventionScores <- interventionScores[order(interventionScores[,1], decreasing=T),]

write.table(as.table(interventionScores),file=file.path(outputDirectory,outputFileName),row.names=FALSE,quote=FALSE,sep='\t')
