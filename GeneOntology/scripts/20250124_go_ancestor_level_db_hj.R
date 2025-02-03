######################
# database construct #
# hojinlee ###########
# 2025.01.24 #########
######################

library(dplyr)
library(fgsea)
library(ggplot2)
library(GO.db)

setwd("/Users/hojin/Dropbox/project/enrichmate/20250124/")

# Setting environment for analysis ----

# 1. human - GOBP ----
gogmt <- gmtPathways("/Users/hojin/Dropbox/project/enrichmate/20250124/GeneOntology/v2024.1.Hs/GOBP/c5.go.bp.v2024.1.Hs.symbols.gmt")
goid_terms_df <- read.table("/Users/hojin/Dropbox/project/enrichmate/20250124/GeneOntology/v2024.1.Hs/GOBP/20250124_c5.go.bp.v2024.1.Hs.term.id_mapping_hj.txt",header = T, sep = "\t") # this file is made from json file

go_ancestor <- as.list(GOBPANCESTOR) # from GO.db package
names(gogmt) %>% length() # 7647
which(goid_terms_df$GOID %in% names(go_ancestor)) %>% length() # 7636

## do not run this step ----
# this step is time consuming, so the results object was already saved as object.

df <- data.frame()
for (go_tmp in names(go_ancestor)) {
  ancestor_tmp <- go_ancestor[[go_tmp]]
  df_tmp <- data.frame(TARGET_GOID = rep(go_tmp, length(ancestor_tmp)+1),
                       ANCESTOR_GOID = c(ancestor_tmp, go_tmp)) # ANCESTOR + TARGET GOID
  df <- rbind(df, df_tmp)
}

df2 <- data.frame()
for (go_tmp in unique(df$ANCESTOR_GOID)) {
  gotest <- GOTERM[[go_tmp]]
  pathway <- gotest@Term
  df_tmp <- data.frame(ANCESTOR_GOID = go_tmp, ANCESTOR_GOTERM = pathway)
  df2 <- rbind(df2, df_tmp)
}

df <- merge(df, df2, by = "ANCESTOR_GOID", all.x = T)

saveRDS(df, file = "/Users/hojin/Dropbox/project/enrichmate/20250124/GeneOntology/v2024.1.Hs/GOBP/20250124_c5.go.bp.v2024.1.Hs_ANCESTOR_df_hj.rds")

df <- data.frame()
for (level_tmp in seq(1,18)) {
  goid_tmp <- GOxploreR::Level2GOTermBP(level = level_tmp) # Level2GOTermBP
  df_tmp <- data.frame(LEVEL = rep(level_tmp, length(goid_tmp)),
                       GOID = goid_tmp)
  df <- rbind(df, df_tmp)
}

saveRDS(df, file = "/Users/hojin/Dropbox/project/enrichmate/20250124/GeneOntology/v2024.1.Hs/GOBP/20250124_c5.go.bp.v2024.1.Hs_LEVEL_df_hj.rds")

# 2. mouse - GOBP ----
gogmt <- gmtPathways("/Users/hojin/Dropbox/project/enrichmate/20250124/GeneOntology/v2024.1.Mm/GOBP/m5.go.bp.v2024.1.Mm.symbols.gmt")
goid_terms_df <- read.table("/Users/hojin/Dropbox/project/enrichmate/20250124/GeneOntology/v2024.1.Mm/GOBP/20250124_m5.go.bp.v2024.1.Mm.term.id_mapping_hj.txt",header = T, sep = "\t") # this file is made from json file

go_ancestor <- as.list(GOBPANCESTOR) # from GO.db package
names(gogmt) %>% length() # 7713
which(goid_terms_df$GOID %in% names(go_ancestor)) %>% length() # 7704

## do not run this step ----
# this step is time consuming, so the results object was already saved as object.

df <- data.frame()
for (go_tmp in names(go_ancestor)) {
  ancestor_tmp <- go_ancestor[[go_tmp]]
  df_tmp <- data.frame(TARGET_GOID = rep(go_tmp, length(ancestor_tmp)+1),
                       ANCESTOR_GOID = c(ancestor_tmp, go_tmp)) # ANCESTOR + TARGET GOID
  df <- rbind(df, df_tmp)
}

df2 <- data.frame()
for (go_tmp in unique(df$ANCESTOR_GOID)) {
  gotest <- GOTERM[[go_tmp]]
  pathway <- gotest@Term
  df_tmp <- data.frame(ANCESTOR_GOID = go_tmp, ANCESTOR_GOTERM = pathway)
  df2 <- rbind(df2, df_tmp)
}

df <- merge(df, df2, by = "ANCESTOR_GOID", all.x = T)

saveRDS(df, file = "/Users/hojin/Dropbox/project/enrichmate/20250124/GeneOntology/v2024.1.Mm/GOBP/20250124_m5.go.bp.v2024.1.Mm_ANCESTOR_df_hj.rds")

df <- data.frame()
for (level_tmp in seq(1,18)) {
  goid_tmp <- GOxploreR::Level2GOTermBP(level = level_tmp, organism = "Mouse") # Level2GOTermBP
  df_tmp <- data.frame(LEVEL = rep(level_tmp, length(goid_tmp)),
                       GOID = goid_tmp)
  df <- rbind(df, df_tmp)
}

saveRDS(df, file = "/Users/hojin/Dropbox/project/enrichmate/20250124/GeneOntology/v2024.1.Mm/GOBP/20250124_m5.go.bp.v2024.1.Mm_LEVEL_df_hj.rds")
