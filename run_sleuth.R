library(sleuth)
library(dplyr)

### MAYBE SEE IF THIS WORKS ON SAMPLE INPUT DATA FIRST ###

stab <- read.table("sleuth_input.txt", header=TRUE) # read in table describing kallisto output

so <- sleuth_prep(stab) # initialize sleuth object

so <- sleuth_fit(so, ~condition, "full") # fit a model comparing the two conditions

so <- sleuth_fit(so, ~1, "reduced") # fit reduced model to compare in likelihood ratio test

so <- sleuth_lrt(so, "reduced", "full") # perform likelihood ratio test for differential expression

sleuth_table <- sleuth_results(so, "reduced:full", "lrt", show_all=FALSE) # extract test results

sleuth_significant <- dplyr::filter(sleuth_table, qval <= 0.05) |> dplyr::arrange(pval) # filter most significant results, sort by pval

head(sleuth_significant, n=10) # print top 10 transcripts

write.table(sleuth_significant, file="fdr05_results.txt", quote=FALSE, row.names=FALSE) # write most significant transcripts to file

head(dplyr::select(sleuth_significant, target_id, pval, qval), n=10) # just show transcript, pval, and qval