library(ggplot2)
library(dplyr)
library(patchwork)
library(plotly)
library(htmlwidgets)




Add_Pop_Meta <- function(data) {
  return(
    case_when(
      startsWith(data$FID, "Car") ~ "Iceland",
      startsWith(data$FID, "U175754CFGIAB") ~ "GIAB",
      startsWith(data$FID, "30028") ~ "Netherlands (No mutation)",
      TRUE ~ "Netherlands"
    )
  )
}

plot_mds <- function(data, label = "") {
  p12 <- ggplot(data, aes(x = C1, y = C2)) +
  geom_point(aes(color = Population), size = 2) +
  theme_minimal() +
  labs(title = paste("MDS-Plot", label), x = "Dim 1", y = "Dim 2")

p13 <- ggplot(data, aes(x = C1, y = C3)) +
  geom_point(aes(color = Population), size = 2) +
  theme_minimal() +
  labs(title = paste("MDS-Plot", label), x = "Dim 1", y = "Dim 3")

p23 <- ggplot(data, aes(x = C2, y = C3)) +
  geom_point(aes(color = Population), size = 2) +
  theme_minimal() +
  labs(title = paste("MDS-Plot", label), x = "Dim 2", y = "Dim 3")

  p12 + p13 + p23

}



mds_data <- read.table("cluster_3.mds", header=T)
mds_subset <- read.table("subset_cluster_3.mds", header=T)
mds_LD <- read.table("clusterLD.mds", header=T)

mds_data$Population <- Add_Pop_Meta(mds_data)
mds_subset$Population <- Add_Pop_Meta(mds_subset)
mds_LD$Population <- Add_Pop_Meta(mds_LD)


full <- plot_mds(mds_data, "200kb region")
subset <- plot_mds(mds_subset, "20kb region")
LD <- plot_mds(mds_LD, "LD correction")


png(filename='MDS.png', width=1000, height=400)
full / su<bset
dev.off()


png(file="MDS_LD.png", width=1000,height=400)
LD
dev.off()

fst_data <- read.table("fst.fst", header = T)
fst_subset <- read.table("subset_fst.fst", header=T)
fst_ld <- read.table("LD.fst", header=T)

fst_full <- fst_data |>
  ggplot() +
  geom_point(aes(x=POS, y=FST))+
  ggtitle("Fixation Index 200kb region") +
  theme_bw()

fst_subset <- fst_subset |>
  ggplot() +
  geom_point(aes(x=POS, y=FST)) +
  ggtitle("Fixation Index 20kb region") +
  theme_bw()

fst_LD <- fst_ld |>
  ggplot() +
  geom_point(aes(x=POS, y=FST)) +
  ggtitle("Fixation Index 200kb LD") +
  theme_bw()

png(filename = "FST.png", width=1000, height=500)
fst_full / fst_subset
dev.off()




######

plt <-ggplot(mds_LD, aes(x = C2, y = C3, text=IID)) +
  geom_point(aes(color = Population), size = 2) +
  theme_minimal() +
  labs(title = paste("MDS-Plot", "C2&C3"), x = "Dimensie 2", y = "Dimensie 3")
