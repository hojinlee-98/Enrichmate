###############
# Enrichmate
# hojinlee
# 2025.02.03
###############

### function ### -----

enrichmate_sim_mat <- function(goids = NULL,
                               sim_method = "Rel",
                               hsGO = NULL) {
  
  sim_method <- "Rel"
  set.seed(1234); w <- GOSemSim::mgoSim(goids, goids, semData = hsGO, measure = sim_method, combine = NULL) # w is similarity matrix. 
  na_goid <- names(which(is.na(w[1,]))) # select NA goid.
  goids <- goids[which(!goids %in% na_goid)]
  set.seed(1234); w <- GOSemSim::mgoSim(goids, goids, semData = hsGO, measure = sim_method, combine = NULL) # w is similarity
  
  return(w)
}

enrichmate_outlier_plot <- function(w = NULL,
                                    min_cluster = 3) {
  
  cutoff <- seq(0.6, 0.9, by = 0.01)
  cutoff = cutoff[cutoff >= 0.5 & cutoff <= 1]
  s1 = s2 = s3 = NULL
  
  for (i in seq_along(cutoff)) {
    set.seed(1234)
    cl = binary_cut(w, cutoff = cutoff[i])
    s1[i] = difference_score(w, cl)
    tb = table(cl)
    s2[i] = length(tb)
    s3[i] = sum(tb < min_cluster)
    
  }
  
  df1 <- data.frame(cutoff, s2)
  colnames(df1) <- c("method", "value")
  df2 <- data.frame(cutoff, s3)
  colnames(df2) <- c("method", "value")
  df1$type <- "All sizes"
  df2$type <- paste("Size", "<", min_cluster)
  df <- rbind(df1, df2)
  
  p <- df %>% ggplot(aes(x = method, y = value, color = type)) +
    geom_point() +
    theme_classic() +
    xlab("Cutoff") +
    ylab("# of clusters") +
    ggtitle("Outlier detection") +
    scale_color_discrete("Type of clusters") +
    theme(aspect.ratio = 0.5,
          plot.title = element_text(hjust = 0.5, size = 15)) + 
    scale_color_manual(values = c("darkgrey", "red3"))
  
  return(p) # select outlier cutoff using this plot.
  
}

enrichmate_clustering <- function(cutoff = 0.84,
                                  w = NULL,
                                  k_vals = NULL,
                                  outlier_detect = T,
                                  min_cluster = 3,
                                  goid_terms_df = NULL,
                                  goid_level_df = NULL) {
  
  ## assign pseudo-name for outliers ----
  # if you want to not assign outliers, change the outlier_detect to FALSE.
  
  if (outlier_detect == T) {
    
    ## binary clustering -----
    set.seed(1234); cl <- binary_cut(mat = w, cutoff = cutoff, partial = F)
    min_term_numb <- 1
    #set.seed(1234); cluster_plot <- ht_clusters(mat = w, cl = cl, draw_word_cloud = F, min_term = min_term_numb) # in this heatmap, small clusters will be considered as outliers.
    bc_res <- data.frame(GOID = colnames(w), cluster = cl) 
    outlier_clusters <- bc_res %>% dplyr::group_by(cluster) %>% dplyr::summarise(n = n()) %>% dplyr::filter(n < min_cluster) %>% dplyr::pull(cluster)
    
    set.seed(1234); pseudoname_outliers <- sample(10000:20000, size = length(outlier_clusters), replace = F)
    pseudoname_outliers_df <- data.frame(pseudoname = pseudoname_outliers, cluster = outlier_clusters) 
    pseudoname_outliers_df <- merge(pseudoname_outliers_df, bc_res, by = "cluster", all.x = T)
  } else {
    pseudoname_outliers_df <- NULL
  }
  
  ## hclust ward.D ----
  
  set.seed(1234); dmat <- stats::as.dist(1 - w) # distance matrix.
  set.seed(1234); hc <- stats::hclust(dmat, method = "ward.D") # h.clust using ward.D.
  
  hc_res_list <- list()
  for (k_val in k_vals) {
    
    set.seed(1234); hclust_wardD <- stats::cutree(hc, k = k_val)
    hc_res_tmp <- data.frame(k = rep(k_val, length(hclust_wardD)),
                             GOID = names(hclust_wardD),
                             cluster = hclust_wardD,
                             row.names = NULL)
    
    hc_res_tmp <- hc_res_tmp[hc$order, ] # ordering
    
    hc_res_list[[paste0("hclust_", k_val)]] <- hc_res_tmp
  }
  
  hc_res_df <- do.call(rbind, hc_res_list)
  rownames(hc_res_df) <- NULL
  
  final_res_df <- merge(hc_res_df, goid_terms_df, by = "GOID", all.x = T)
  final_res_df <- merge(final_res_df, goid_level_df, by = "GOID", all.x = T) # get LEVEL information
  final_res_df <- final_res_df %>% dplyr::arrange(cluster) # arrange cluster
  
  my_list <- list("pseudoname_outliers_df" = pseudoname_outliers_df,
                  "hc_res_df" = hc_res_df,
                  "final_res_df" = final_res_df)
  
  return(my_list)
  
}

enrichmate_representative_terms <- function(final_res_df = NULL,
                                            k_val = NULL,
                                            outlier_detect = T,
                                            goid_ancestor_level_df = NULL,
                                            pseudoname_outliers_df = NULL,
                                            representative_term_level_cutoff = 1,
                                            GO_explain = 2
) {
  
  goid_ancestor_level_df <- goid_ancestor_level_df %>% dplyr::filter(ANCESTOR_LEVEL >= representative_term_level_cutoff)
  
  if (is.null(pseudoname_outliers_df)) {
    
    outlier_detect = F
    
  } else {
    
    if (nrow(pseudoname_outliers_df) == 0) {
      
      outlier_detect = F
      
    }
  }
  
  if (outlier_detect == T) {
    final_res_df <- final_res_df %>% dplyr::filter(!(GOID %in% pseudoname_outliers_df$GOID))
  }
  
  df <- data.frame()
  ans_df <- data.frame()
  first_rep_c_goid_df <- data.frame()
  for (k_tmp in k_val) {
    final_res_df_tmp <- final_res_df %>% dplyr::filter(k == k_tmp)
    
    
    for (c_tmp in sort(unique(final_res_df_tmp$cluster))) {
      final_res_df_c_tmp <- final_res_df_tmp %>% dplyr::filter(cluster == c_tmp)
      ancestor_df_tmp <- goid_ancestor_level_df %>% dplyr::filter(TARGET_GOID %in% final_res_df_c_tmp$GOID)
      ancestor_df_tmp$k <- k_tmp
      ancestor_df_tmp$cluster <- c_tmp
      ans_df <- rbind(ans_df, ancestor_df_tmp)
      
      seed_c_tmp <- final_res_df_c_tmp
      remained_goid <- final_res_df_c_tmp$GOID
      
      
      while (TRUE) {
        
        explained_go <- floor(nrow(seed_c_tmp)/GO_explain)
        ancestor_sub_df_tmp <- goid_ancestor_level_df %>% dplyr::filter(TARGET_GOID %in% seed_c_tmp$GOID)
        
        ancestor_sub_df_filtered <- ancestor_sub_df_tmp %>%
          dplyr::group_by(ANCESTOR_GOID, ANCESTOR_LEVEL) %>%
          dplyr::summarise(n=n()) %>%
          dplyr::arrange(desc(n)) %>%
          dplyr::filter(n > explained_go)
        
        if (nrow(ancestor_sub_df_filtered) != 0) {
          
          max_ancestor_tmp <- ancestor_sub_df_filtered %>%
            dplyr::ungroup(ANCESTOR_GOID, ANCESTOR_LEVEL) %>%
            dplyr::filter(ANCESTOR_LEVEL == max(ANCESTOR_LEVEL)) %>%
            dplyr::filter(n == max(n))
          
          first_rep_c_goid <- ancestor_sub_df_tmp %>%
            dplyr::filter(ANCESTOR_GOID %in% max_ancestor_tmp$ANCESTOR_GOID) %>%
            dplyr::pull(TARGET_GOID) %>% unique()
          
          # get the common ancestor terms and child term 
          first_rep_c_goid_df_tmp <- ancestor_sub_df_tmp %>% dplyr::filter(TARGET_GOID %in% first_rep_c_goid &
                                                                             ANCESTOR_GOID %in% max_ancestor_tmp$ANCESTOR_GOID)
          first_rep_c_goid_df_tmp$cluster <- c_tmp
          first_rep_c_goid_df <- rbind(first_rep_c_goid_df, first_rep_c_goid_df_tmp) # ancestor and child dataframe
          
          remained_goid <- seed_c_tmp %>%
            dplyr::filter(!GOID %in% first_rep_c_goid) %>%
            dplyr::pull(GOID)
          
          seed_c_tmp <- seed_c_tmp %>% dplyr::filter(GOID %in% remained_goid)
          
          df_tmp <- data.frame(k = rep(k_tmp, nrow(max_ancestor_tmp)),
                               cluster = rep(c_tmp, nrow(max_ancestor_tmp)),
                               max_level = rep(unique(max_ancestor_tmp$ANCESTOR_LEVEL), nrow(max_ancestor_tmp)),
                               n = rep(unique(max_ancestor_tmp$n), nrow(max_ancestor_tmp)),
                               total_gobp = rep(nrow(final_res_df_c_tmp), nrow(max_ancestor_tmp)),
                               n_total_gobp = rep(paste(c(unique(max_ancestor_tmp$n), nrow(final_res_df_c_tmp)), collapse = " out of "), nrow(max_ancestor_tmp)),
                               ANCESTOR_GOID = max_ancestor_tmp$ANCESTOR_GOID)
          
          df <- rbind(df, df_tmp)
  
          if (length(remained_goid) == 1) {
            df_tmp <- data.frame(k = rep(k_tmp, nrow(seed_c_tmp)),
                                 cluster = rep(c_tmp, nrow(seed_c_tmp)),
                                 max_level = rep(".", nrow(seed_c_tmp)),
                                 n = rep(1, nrow(seed_c_tmp)),
                                 total_gobp = rep(nrow(final_res_df_c_tmp), nrow(seed_c_tmp)),
                                 n_total_gobp = rep(paste(c(1, nrow(final_res_df_c_tmp)), collapse = " out of "), nrow(seed_c_tmp)),
                                 ANCESTOR_GOID = seed_c_tmp$GOID)
            df <- rbind(df, df_tmp)
            
            break
          } else if (length(remained_goid) == 0) {
            break
          }
          
        } else {
          
          df_tmp <- data.frame(k = rep(k_tmp, nrow(seed_c_tmp)),
                               cluster = rep(c_tmp, nrow(seed_c_tmp)),
                               max_level = rep(".", nrow(seed_c_tmp)),
                               n = rep(1, nrow(seed_c_tmp)),
                               total_gobp = rep(nrow(final_res_df_c_tmp), nrow(seed_c_tmp)),
                               n_total_gobp = rep(paste(c(1, nrow(final_res_df_c_tmp)), collapse = " out of "), nrow(seed_c_tmp)),
                               ANCESTOR_GOID = seed_c_tmp$GOID)
          
          df <- rbind(df, df_tmp)
          
          break
        }
        
      }
    }
  }
  
  rep_df <- df
  
  return(rep_df)
  
}

enrichmate_go_heatmap_summary <- function(w = NULL,
                                          k_val = NULL,
                                          outlier_detect = T,
                                          heatmap_filename = "test.pdf",
                                          heatmap_width = 20,
                                          heatmap_height = 20,
                                          filename1 = "total_GOBP.xlsx",
                                          filename2 = "representative_term.xlsx",
                                          pseudoname_outliers_df = NULL,
                                          ancestor_annotation = T,
                                          representative_term_df = NULL,
                                          font_size = 3,
                                          plot_pdf = T,
                                          goid_terms_df = NULL,
                                          rep_df = NULL) {
  
  col <- my_palette <- colorRampPalette(c("white","red"))(n = 300)
  col_fun <- colorRamp2(seq(0, quantile(w[w > 0], 0.975), length = length(c("white","red"))), c("white","red"))
  set.seed(1234); dmat <- stats::as.dist(1 - w) # distance matrix.
  set.seed(1234); hc <- stats::hclust(dmat, method = "ward.D") # h.clust using ward.D.
  hclust_wardD <- stats::cutree(hc, k = k_val)
  
  hc_cluster_df <- data.frame(GOID = names(hclust_wardD), cluster = hclust_wardD, row.names = NULL)
  hc_cluster_df <- hc_cluster_df[hc$order, ] # ordering 
  rownames(hc_cluster_df) <- NULL
  
  if (is.null(pseudoname_outliers_df)) {
    
    outlier_detect = F
    
  } else {
    
    if (nrow(pseudoname_outliers_df) == 0) {
      
      outlier_detect = F
      
      }
  }
  
  # if you do not want to assign outliers, change outlier_detect
  if (outlier_detect == T) {
    
    pseudoname_outliers_df <- pseudoname_outliers_df %>% dplyr::select(-c("cluster")) %>% dplyr::rename("cluster" = pseudoname) 
    
    hc_cluster_df <- hc_cluster_df %>% dplyr::filter(!(GOID %in% pseudoname_outliers_df$GOID))
    hc_cluster_df$order <- as.numeric(rownames(hc_cluster_df))
    bc_cluster_df <- pseudoname_outliers_df %>% dplyr::arrange(cluster) %>% dplyr::select(c("GOID"))
    bc_cluster_df <- pseudoname_outliers_df %>% dplyr::arrange(cluster)
    
    bc_cluster_tb <- table(pseudoname_outliers_df$cluster)
    bc_cluster_numb <- length(bc_cluster_tb)
    bc_cluster_orig_name <- as.numeric(names(bc_cluster_tb))
    
    max_cl <- max(hc_cluster_df$cluster)
    cl_mold <- data.frame(new_cluster = seq(max_cl+1, max_cl+bc_cluster_numb),
                          cluster = bc_cluster_orig_name)
    
    bc_cluster_df <- merge(pseudoname_outliers_df, cl_mold, by = "cluster", all.x = T) %>% dplyr::select(c("GOID", "new_cluster")) %>% dplyr::rename("cluster" = new_cluster)
    bc_cluster_df$order <- max(as.numeric(rownames(hc_cluster_df))) + as.numeric(rownames(bc_cluster_df))
    
    final_order <- rbind(hc_cluster_df, bc_cluster_df)
    #final_order <- hc_cluster_df
    rownames(final_order) <- final_order$GOID
    
    outliers_rm_goids <- rownames(w)[!rownames(w) %in% pseudoname_outliers_df$GOID]
    w <- w[outliers_rm_goids, outliers_rm_goids]
    final_order <- final_order %>% dplyr::filter(cluster <= k_val)
    
    total_clusters <- final_order %>% dplyr::group_by(cluster) %>% dplyr::summarise(n=n()) %>% dplyr::select(cluster) %>% table() %>% names()
    orig_order <- data.frame(GOID = colnames(w), orig_order = seq(1, nrow(w)), row.names = colnames(w))
    orig_order <- orig_order %>% dplyr::filter(!GOID %in% pseudoname_outliers_df$GOID)
    orig_order <- orig_order[final_order$GOID, ]
    
  } else if (outlier_detect == F) {
    hc_cluster_df$order <- as.numeric(rownames(hc_cluster_df))
    
    max_cl <- max(hc_cluster_df$cluster)
    final_order <- hc_cluster_df
    rownames(final_order) <- final_order$GOID
    
    final_order <- final_order %>% dplyr::arrange(order)
    total_clusters <- final_order %>% dplyr::group_by(cluster) %>% dplyr::summarise(n=n()) %>% dplyr::select(cluster) %>% table() %>% names()
    orig_order <- data.frame(GOID = colnames(w), orig_order = seq(1, nrow(w)), row.names = colnames(w))
    orig_order <- orig_order[final_order$GOID, ]
  }
  
  # Summary results 
  final_order_summary <- merge(final_order, goid_terms_df, by = "GOID", all.x = T) %>% dplyr::arrange(order)
  clusters_orig <- final_order_summary$cluster[(!duplicated(final_order_summary$cluster))]
  c_df <- data.frame(clusters_orig = clusters_orig, clustsers_anno = seq(1, length(clusters_orig)))
  final_order_summary <- merge(final_order_summary, c_df, by.x = "cluster", by.y = "clusters_orig", all.x = T)
  final_order_summary$GOTERM <- gsub(x = final_order_summary$GOTERM, pattern = "GOBP_", replacement = "") %>%
    gsub(x = ., pattern = "_", replacement = " ") %>% tolower()
  
  # total GO terms with the information
  WriteXLS(final_order_summary, ExcelFileName = filename1) # this dataframe includes GO terms and cluster information.
  
  #return(list(final_order_summary, rep_df))
  final_order_summary <- merge(final_order_summary, rep_df, by.x = "cluster", by.y = "cluster", all.x = T)
  rep_results_anno <- final_order_summary %>% dplyr::select(c("clustsers_anno", "max_level", "ANCESTOR_GOID", "n_total_gobp", "n", "total_gobp"))
  rep_results_anno <- rep_results_anno[!duplicated(rep_results_anno[c("clustsers_anno", "max_level", "ANCESTOR_GOID", "n_total_gobp")]),]
  rep_results_anno <- rep_results_anno %>% dplyr::filter(!is.na(ANCESTOR_GOID)) # filtering out outliers
  
  rep_results_anno_final <- data.frame()
  for (i in seq(1:nrow(rep_results_anno))) {
    rep_results_anno_tmp <- rep_results_anno[i,]
    
    GOTERM_tmp <- goid_ancestor_level_df %>%
      dplyr::filter(ANCESTOR_GOID == rep_results_anno_tmp$ANCESTOR_GOID) %>%
      dplyr::slice(1) %>%
      dplyr::pull(ANCESTOR_GOTERM)
    rep_results_anno_tmp$GOBP_ANCESTOR_GOTERM <- GOTERM_tmp
    rep_results_anno_final <- rbind(rep_results_anno_final, rep_results_anno_tmp)
  }
  
  rep_results_anno_final <- rep_results_anno_final %>% dplyr::arrange(clustsers_anno, desc(n))
  
  # final results
  WriteXLS(rep_results_anno_final, ExcelFileName = filename2) # this dataframe is the final results. put the table with final heatmap.
  
  # heatmap ----
  # frequency for clusters
  cl_vec <- c()
  freq_cl_vec <- c()
  for (anode in unique(final_order$cluster)) {
    cl_tmp <- as.numeric(table(final_order$cluster)[as.character(anode)])
    cl_vec <- append(cl_vec, cl_tmp)
    term_sum <- length(final_order$cluster)
    freq_cl_tmp <- cl_tmp/term_sum
    freq_cl_vec <- append(freq_cl_vec, freq_cl_tmp)
  }
  
  col_fun <- colorRamp2(seq(0, quantile(w[w > 0], 0.975), length = length(c("white","red"))), c("white","red"))
  
  if (plot_pdf == T) {
    pdf(heatmap_filename, width = heatmap_width, height = heatmap_height) 
  }
  
  if (ancestor_annotation == T) {
    
    p1 <- Heatmap(mat = w, col = col_fun, name = "Similarity",
                  row_order = orig_order$orig_order,
                  column_order = orig_order$orig_order,
                  border = "#404040", row_title = NULL,
                  show_column_dend = F,
                  show_row_dend = F,
                  show_row_names = F,
                  show_column_names = F,
                  row_dend_width = unit(2, "cm"),
                  width = unit(20, "cm"),
                  height = unit(20, "cm"),
                  left_annotation = rowAnnotation(ggplot1 = anno_empty(height = unit(20, "cm"), width = unit(1, "cm")),
                                                  show_annotation_name = FALSE, gp = gpar(col = "white")),
                  right_annotation = rowAnnotation(ggplot2 = anno_empty(border = FALSE, height = unit(20, "cm"), width = unit(15, "cm")),
                                                   show_annotation_name = FALSE, gp = gpar(col = "white")
                  ))
    
  } else if (ancestor_annotation == F) {
    
    p1 <- Heatmap(mat = w, col = col_fun, name = "Similarity",
                  row_order = orig_order$orig_order,
                  column_order = orig_order$orig_order,
                  border = "#404040", row_title = NULL,
                  show_column_dend = F,
                  show_row_dend = F,
                  show_row_names = F,
                  show_column_names = F,
                  row_dend_width = unit(2, "cm"),
                  width = unit(20, "cm"),
                  height = unit(20, "cm"),
                  left_annotation = rowAnnotation(ggplot1 = anno_empty(height = unit(20, "cm"), width = unit(1, "cm")),
                                                  show_annotation_name = FALSE, gp = gpar(col = "white")))
    
  }
  
  print(p1)
  
  freq_cl_cumsum <- cumsum(freq_cl_vec)
  freq_cl_cumsum <- freq_cl_cumsum[-length(freq_cl_cumsum)]
  
  gap_filler <- 0.0001
  decorate_heatmap_body("Similarity", {
    for (freq_cl in freq_cl_cumsum) {
      grid.lines(c(freq_cl + gap_filler, freq_cl + gap_filler), c(0, 1), gp = gpar(lty = 1, lwd = 1.5))
    }
  })
  
  freq_cl_cumsum_rev <- cumsum(rev(freq_cl_vec))
  freq_cl_cumsum_rev <- freq_cl_cumsum_rev[-length(freq_cl_cumsum_rev)]
  
  decorate_heatmap_body("Similarity", {
    for (freq_cl in freq_cl_cumsum_rev) {
      grid.lines(c(1, 0), c(freq_cl + gap_filler, freq_cl + gap_filler), gp = gpar(lty = 1, lwd = 1.5))
    }
  })
  coor_df <- data.frame(x1 = 0, x2 = 5,
                        y1 = c(0,freq_cl_cumsum_rev), y2 = c(freq_cl_cumsum_rev, 1))
  
  if (outlier_detect == T) {
    p <- ggplot() +
      geom_rect(data=coor_df,
                mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2, fill=as.factor(y1)),
                color= "black",
                alpha=1) +
      scale_x_continuous(expand = c(0,0)) +
      scale_y_continuous(expand = c(0,0)) +
      theme_void() +
      theme(panel.background = element_blank(),
            panel.grid.major= element_blank(),
            panel.grid.minor = element_blank(),
            axis.line = element_blank(),
            axis.ticks = element_blank(),
            axis.title = element_blank(),
            axis.text = element_blank(),
            legend.position = "none",
            plot.background = element_blank()) +
      labs(x=NULL, y=NULL) +
      geom_text(data = coor_df,
                aes(x=(x1+x2)/2, y=y1+((y2-y1)/2), label = rev(c(paste0("c", total_clusters))), size = 3))
    
  } else if (outlier_detect == F) {
    
    p <- ggplot() +
      geom_rect(data=coor_df,
                mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2, fill=as.factor(y1)),
                color= "black",
                alpha=1) +
      scale_x_continuous(expand = c(0,0)) +
      scale_y_continuous(expand = c(0,0)) +
      theme_void() +
      theme(panel.background = element_blank(),
            panel.grid.major= element_blank(),
            panel.grid.minor = element_blank(),
            axis.line = element_blank(),
            axis.ticks = element_blank(),
            axis.title = element_blank(),
            axis.text = element_blank(),
            legend.position = "none",
            plot.background = element_blank()) +
      labs(x=NULL, y=NULL) +
      geom_text(data = coor_df,
                aes(x=(x1+x2)/2, y=y1+((y2-y1)/2), label = rev(c(paste0("c", total_clusters[1:k_val]))), size = 3))
    
  }
  
  if (ancestor_annotation == T) {
    
    #return(rep_results_anno_final)
    
    representative_term_df_plot <- rep_results_anno_final %>%
      dplyr::arrange(clustsers_anno, desc(n)) %>%
      dplyr::group_by(clustsers_anno) %>%
      dplyr::filter(n == max(n))
    
    #return(representative_term_df_plot)
    
    representative_term_df_plot_final <- data.frame()
    
    for (c_tmp in sort(unique(representative_term_df_plot$clustsers_anno))) {
      
      representative_term_df_plot_tmp <- representative_term_df_plot %>%
        dplyr::filter(clustsers_anno == c_tmp)
      
      if (representative_term_df_plot_tmp[1,]$n != 1) {
        
        plot_anno_tmp <- paste(representative_term_df_plot_tmp$GOBP_ANCESTOR_GOTERM, collapse = "\n")
        representative_term_df_plot_tmp$plot_anno <- plot_anno_tmp 
        
        representative_term_df_plot_final <- rbind(representative_term_df_plot_final, representative_term_df_plot_tmp)
        
      } else if (representative_term_df_plot_tmp[1,]$n == 1) {
        
        if (nrow(representative_term_df_plot_tmp) < 4) {
          
          plot_anno_tmp <- paste(representative_term_df_plot_tmp$GOBP_ANCESTOR_GOTERM, collapse = "\n")
          representative_term_df_plot_tmp$plot_anno <- plot_anno_tmp 
          
          representative_term_df_plot_final <- rbind(representative_term_df_plot_final, representative_term_df_plot_tmp)
          
        } else {
          
          #plot_anno_tmp <- paste(representative_term_df_plot_tmp$GOBP_ANCESTOR_GOTERM, collapse = "\n")
          representative_term_df_plot_tmp$plot_anno <- "" 
          
          representative_term_df_plot_final <- rbind(representative_term_df_plot_final, representative_term_df_plot_tmp)
          
        }
        
      }
      
    }
    
    representative_term_df_plot_final <- representative_term_df_plot_final[!duplicated(representative_term_df_plot_final[c("clustsers_anno", "n")]),]
    
    #return(representative_term_df_plot_final)
    
    # representative_term_df_plot <- rep_results_anno_final %>%
    #   dplyr::arrange(clustsers_anno, desc(n)) %>%
    #   dplyr::group_by(clustsers_anno) %>%
    #   dplyr::filter(n == max(n)) %>%
    #   dplyr::group_by(clustsers_anno) %>% 
    #   dplyr::mutate(plot_anno = paste(GOBP_ANCESTOR_GOTERM, collapse = "\n")) %>%
    #   dplyr::mutate(plot_anno = case_when(n == 1 ~ "",
    #                                       T ~ plot_anno))
    
    #representative_term_df_plot <- representative_term_df_plot[!duplicated(representative_term_df_plot[c("clustsers_anno", "n")]),]
    
    p2 <- ggplot() +
      geom_rect(data=coor_df,
                mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2), fill=NA, color=NA,
                alpha=1) +
      scale_x_continuous(expand = c(0,0)) +
      scale_y_continuous(expand = c(0,0)) +
      theme_void() +
      theme(panel.background = element_blank(),
            panel.border = element_blank(),
            panel.grid.major= element_blank(),
            panel.grid.minor = element_blank(),
            axis.line = element_blank(),
            axis.ticks = element_blank(),
            axis.title = element_blank(),
            axis.text = element_blank(),
            legend.position = "none",
            plot.background = element_blank()) +
      labs(x=NULL, y=NULL) +
      geom_text(data = coor_df, hjust = 0, lineheight = 0.7, size = font_size,
                aes(x=0, y=y1+((y2-y1)/2), label = rev(c(representative_term_df_plot_final$plot_anno)), color = as.factor(y1)))
    
  } else if (ancestor_annotation == F) {
    
    p2 <- NULL
    
  }
  
  if (ancestor_annotation == T) {
    
    decorate_annotation("ggplot1", {
      vp = current.viewport()$name
      print(p, vp = vp)
    })
    
    decorate_annotation("ggplot2", {
      vp = current.viewport()$name
      print(p2, vp = vp)
    })
    
  } else if (ancestor_annotation == F) {
    
    decorate_annotation("ggplot1", {
      vp = current.viewport()$name
      print(p, vp = vp)
    }) 
    
  }
  
  if (plot_pdf == T) {
    dev.off()
  }
  
  my_list <- list("total_GOBP_order" = final_order_summary, "rep_results_anno_final" = rep_results_anno_final)
  
  return(my_list)
  
}
