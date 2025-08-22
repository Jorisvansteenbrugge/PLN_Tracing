library(vcfR)
library(ggplot2)
library(ape)
library(magrittr)


vcf_file <- "PLN_international_200kb_rename_nonMissing.filtered.phased.vcf.gz"
vcf <- read.vcfR(vcf_file)

nl_rename <- read.table("NL_rename.csv", sep=';') |>
nl_rename <- setNames(as.list(nl_rename[[1]]), nl_rename[[2]])

# import genotype field ("GT") as strings ( "0|1", "1|1", etc.)
gt <- extract.gt(vcf, element = "GT", as.numeric = FALSE)
# 'gt' is a matrix where rows = variants, columns = samples.
# e.g.,: gt[1, "SampleA"] = "0|1"



# Clean up sample names
colnames(gt) %<>%
  sapply(function(x){
    snames = strsplit(x,split="_") |>unlist()
    return(snames[1])
  }) %>%
  as.character()

colnames(gt) %<>%
  sapply(function(x){
    if(x %in% names(nl_rename)){
      return(nl_rename[[x]])
    }
    else {
      return(x)
    }
  }) %>%
  as.character

samples <- colnames(gt)
n_samples <- length(samples)
n_variants <- nrow(gt)

# ============================================================
# 3) DIPLOID ANALYSIS
# ============================================================

# Convert phased haplotypes ("0|0", "0|1", "1|1") to a single score (e.g, 0/1/2)
# (0 = homo ref, 1 = heterozygoot, 2 = homo alt)
diploid_geno_to_numeric <- function(g) {
  # bv. "0|1" -> c("0","1") -> sum(0,1) = 1
  # "1|1" -> c("1","1") -> sum(1,1) = 2
  # "0|0" -> sum(0,0) = 0
  alleles <- strsplit(g, "\\|")[[1]]
  return(sum(as.numeric(alleles)))
}

# Create a matrix (n_samples x n_variants) in which each cell is 0,1 or 2
diploid_mat <- matrix(NA, nrow = n_samples, ncol = n_variants)
rownames(diploid_mat) <- samples  # row-names is sample names
colnames(diploid_mat) <- rownames(gt)  # variant IDs

for (i in seq_len(n_samples)) {
  sample_id <- samples[i]
  # gt[, sample_id] is a vector of length n_variants
  diploid_mat[i, ] <- sapply(gt[, sample_id], diploid_geno_to_numeric)
}

rownames(diploid_mat) %<>%
  sapply(function(x){
    snames = strsplit(x,split="_") |>unlist()
    return(snames[1])
  }) |>
  as.character()

# Hierarchical clustering
dist_mat_dip <- dist(diploid_mat, method = "manhattan")
hc_dip <- hclust(dist_mat_dip, method = "complete")


pdf(file="haplotypes_merged.pdf", width = 12)
plot(hc_dip, main = "Dendrogram Merged Haplotypes)",
     xlab = "", sub = "")
dev.off()



# ============================================================
# 4) HAPLOTYPE-ANALYSE (twee rijen per sample)
# ============================================================
# Idea: "0|1" -> haplotypeA=0, haplotypeB=1
# build a matrix (2*n_samples x n_variants):
#   - row1 = SampleA_hap1
#   - row2 = SampleA_hap2
#   - row3 = SampleB_hap1
#   - row4 = SampleB_hap2
#   ...
# only 0 or 1 per haplotype.

haplo_mat <- matrix(NA, nrow = 2 * n_samples, ncol = n_variants)

hap_names <- character(2 * n_samples)

# Convert "0|1" -> (0,1), "1|0" -> (1,0) etc.
split_haplotypes <- function(g) {
  alleles <- strsplit(g, "\\|")[[1]]
  return(as.numeric(alleles))  # e.g., c(0,1)
}

row_index <- 1

for (s in seq_len(n_samples)) {
  sample_id <- samples[s]

  gt_vec <- gt[, sample_id]  # length = n_variants, e.g., "0|0", "0|1", ...

  # Split each genotype into two haplotype-alleles
  # Create a list of c(alleleA, alleleB), repeated for n_variants
  hap_list <- lapply(gt_vec, split_haplotypes)

  # hapA = vector with first allele of each SNP
  # hapB = vector with second allele of each SNP
  hapA <- sapply(hap_list, `[`, 1)
  hapB <- sapply(hap_list, `[`, 2)


  sample_id <- unlist(strsplit(sample_id,split="_"))[1]

  haplo_mat[row_index, ] <- hapA
  hap_names[row_index] <- paste0(sample_id, "_hapA")


  haplo_mat[row_index + 1, ] <- hapB
  hap_names[row_index + 1] <- paste0(sample_id, "_hapB")

  row_index <- row_index + 2
}

rownames(haplo_mat) <- hap_names
colnames(haplo_mat) <- rownames(gt)
# Now we have a matrix with (2 * n_samples) rows, and (haplotypes x n_variants) columns.

# Clustering (Manhattan + average, b.v.)
dist_mat_hap <- dist(haplo_mat, method = "manhattan")
hc_hap <- hclust(dist_mat_hap, method = "complete")


pdf(file="haplotypes_split.pdf", width = 24)
plot(hc_hap, main = "Dendrogram Split haplotypes",
     xlab = "", sub = "")
dev.off()


hap_phylo <- as.phylo(hc_hap)
write.tree(hap_phylo, file = "haplotype_clusters.tree")
