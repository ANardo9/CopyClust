library(devtools)
library(ggplot2)

load_all()

load("~/Desktop/Caldas Lab/Caldas Lab/Datatest_CopyClust.RData")

features = CC_format(Y, reference_genome = "hg19", probes = 26050)

results = CopyClust(features)

results_numeric = data.frame(label = as.integer(results$IntClust_Label))

ggplot(data = results_numeric, aes(x=factor(label))) +
  geom_histogram(stat = "count") +
  xlab("IntClust")

#####

X <- data.frame(ID=Y$ID, chrom=Y$chr, loc.start=Y$loc.start, loc.end=Y$loc.end,
                num.mark=Y$num.mark, seg.mean=Y$seg.mean)
X$chrom[which(X$chrom=="X")] <- 23
features <- CC_format(X, reference_genome="hg19", probes=26050)
results <- CopyClust(features)

#####

