## 

    # rm(list=ls())
    source("analysis/ms__van__spreadsheeting__shortcuts.R")
    set.seed(2205)
    
    library(rtk)
    library(parallel)

    head(mgdat)
        
    
## rarefaction   ========================================================

      #   summary(rowSums(mgfeat))
      #   # ggplot(data=mgdat, aes(x=reorder(sample, SeqDepth), y=SeqDepth, color=substrate)) +
      #   #   geom_point()
      #   
      #   # need ot pick a cutting point for the data - 7100 drops 4 samples
      #   rtk_usable_depth <- "7100"
      #   rtk_its <- 100
      # 
      # 
      #   # locally manageable
      #   mgrar <- rtk(mgfeat,
      #                margin = 1,
      #                repeats = rtk_its,
      #                threads = 7,
      #                depth = seq(100, 10000, by = 1000),
      #                ReturnMatrix = rtk_its)
      #   # saveRDS(mg_rar, "output/ms__som__rar_100-10000_100.RDS")
      # 
      #   # bad cuve
      #   plot(mgrar, div = "shannon", groups = mgdat$substrate)
      #   
      # ## set div at 7K - lines fairly stable therafter
      #   d7k <- t(sapply(mgrar[[rtk_usable_depth]][["divvs"]], function(aa){   # aa <- mgrar[["9100"]][["divvs"]][[15]]
      #     bb_dmed <- lapply( aa[2:7], function(aaa){ median(aaa)})
      #     unlist(c( aa[[1]], bb_dmed))
      #   } ))
      #   rownames(d7k) <- d7k[,1] ; d7k <- t(apply(d7k[,-1], 1, as.numeric))
      #   colnames(d7k) <- c("richness", "shannon", "simpson", "invsimpson", "chao1", "evenness")
      #   head(d7k)
      #   ## fix names / order!


      ## faith's PD on 5k 100 replicates  ---------------------------------------
      mgtree <- readRDS("output/ms__van_phangoGTR_NJ_tree.RDS")
      # 
      ## if using RTK
          # rtk_mats <- mgrar[[rtk_usable_depth]][["raremat"]]
          # 
          # faith_df <- mclapply( rtk_mats, function(aa){   # aa <- rtk_mats[[3]]
          #   bb <- apply(aa, 1, function(aaa){   # aaa <- aa[10,]
          #     abdiv::faith_pd( as.vector(aaa), mgtree)
          #   })
          # },  mc.cores = 7)
          # length(faith_vec <- colMeans(do.call("rbind", faith_df)))
          # 
          # # get missing samples - use the full sample directly, as this is all that would be possible at this depth
          # miss_faith <- sapply( mgdat[ !(mgdat$sample %in% names(faith_vec)), "sample"], function(aa){
          #   abdiv::faith_pd( mgfeat[ aa, ], mgtree)
          # })
          # 
          # # make a matrix to stick back together, more reliable than matchin names
          # all_faiths <- c(miss_faith, faith_vec)
          # af_mat <- matrix(all_faiths, nrow = length(all_faiths), dimnames = list(names(all_faiths),"faith_pd"))
          # head(af_mat)
          # 
          # ## bring together in correct order
          # mgdd <- data.frame( mgdat[ rownames(mgdat) , ],
          #                     d7k[ rownames(mgdat) , ],
          #                     "faith_pd" = af_mat[ rownames(mgdat) , ] )
          # head(mgdd)
          # dim(mgdd)
          # saveRDS(mgdd, "output/ms__van__metadata_63-19_A-Div.RDS")
          # mgdd <- readRDS("output/ms__van__metadata_63-36_A-Div.RDS")

    
    ##  vegan diversity ---------------------------------

      # mgdd <- mgdat
      # 
      # # # we suspect
      # # mgdd[ "D15" , "timepoint"] <- "T0"
      # 
      # mgdd$shan <- vegan::diversity( mgfeat[ rownames(mgdd) , ], index = "shannon")
      # mgdd$invs <- vegan::diversity( mgfeat[ rownames(mgdd) , ], index = "invsimpson")
      # mgdd$rich <- apply(mgfeat[ rownames(mgdd) , ], 1, function(aa) sum( aa > 0 ) )
      # 
      # ## if not using rtk - go basic
      # mgdd$faith_pd <- unlist( mclapply( 1:nrow(mgfeat), function(aa){   # aa <- 12
      #   abdiv::faith_pd( as.vector( mgfeat[ rownames(mgdd),][aa,] ), mgtree)
      # }, mc.cores = 7, mc.preschedule = TRUE))
      # 
      # saveRDS(mgdd, "output/ms__van__metadata_58-12_A-Div.RDS")
      # 
      mgdd <- readRDS("output/ms__van__metadata_58-12_A-Div.RDS")
      
      
## alpha testing   -------------------------------------
    
    # can use models for this if residuals are norm, n'es pas?
    head(mgdd)
    
    divs <- c("rich", "shan", "invs", "faith_pd","SeqDepth") # , "chao1", 
    mdivs <- mgdd[ , c("van_type", "sample", divs) ] 
    head(mdivs)
    
    # not particularly parametric
    # md_m <- reshape2::melt(mdivs)
    # ggplot(md_m, aes(x = value)) +
    #   facet_wrap(variable~, scale = "free") +
    #   geom_density(aes(fill = van_type), alpha = 0.5)
    
    
## alpha div plots ================================================================
    
    # bees? not violins, sadly..
    # col_gend <-nr_col[c(13,14)]   # show_col(nr_col)
    
    ad1 <-  ggplot(mgdd, aes(x = timepoint, fill = diet, y  = rich)) +
      # facet_grid(.~diet, space = "free") +
      geom_boxplot( size = 0.9, outlier.shape = 21) + # aes(fill = diet), 
      geom_point(aes(fill = diet), position = position_jitterdodge(jitter.width = 0.25), shape = 21, size = 1.5, alpha = 0.6) +
      # ggbeeswarm::geom_beeswarm(aes(fill = diet)) + 
      scale_fill_manual(values = nr_col) +
      labs(y = "species richness", x = "") +
      theme( axis.text.x = element_text(angle = 0))
    
    ad2 <-  ggplot(mgdd, aes(x = timepoint, fill = diet, y  = shan)) +
      # facet_grid(.~diet, space = "free") +
      scale_fill_manual(values = nr_col) +
      geom_boxplot( size = 0.9,  outlier.shape = 21) +
      geom_point(aes(fill = diet), position = position_jitterdodge(jitter.width = 0.25), shape = 21, size = 1.5, alpha = 0.6) +
      # ggbeeswarm::geom_beeswarm(aes(fill = diet)) + 
      labs(y = "Shannon Diversity", x = "") +
      theme( axis.text.x = element_text(angle = 0))
    
    ad3 <-  ggplot(mgdd, aes(x = timepoint, fill = diet, y  = invs)) +
      # facet_grid(.~diet, space = "free") +
      scale_fill_manual(values = nr_col) +
      geom_boxplot( size = 0.9,  outlier.shape = 21) +
      geom_point(aes(fill = diet), position = position_jitterdodge(jitter.width = 0.25), shape = 21, size = 1.5, alpha = 0.6) +
     # ggbeeswarm::geom_beeswarm(aes(fill = diet)) + 
      labs(y = "Inverse Simpsons", x = "") +
      theme( axis.text.x = element_text(angle = 0))
    
    ad4 <-  ggplot(mgdd, aes(x = timepoint, fill = diet, y  = faith_pd)) +
      # facet_grid(.~diet, space = "free") +
      scale_fill_manual(values = nr_col) +
      # coord_cartesian(ylim = c(4, 5.25)) +
      geom_boxplot( size = 0.9, outlier.shape = 21) +
      geom_point(aes(fill = diet), position = position_jitterdodge(jitter.width = 0.25), shape = 21, size = 1.5, alpha = 0.6) +
      # ggbeeswarm::geom_beeswarm(aes(fill = diet)) + 
      labs(y = "Faith's Phylogenetic Diversity", x = "") +
      theme( axis.text.x = element_text(angle = 0))
    
    ad5 <-  ggplot(mgdd, aes(x = timepoint, fill = diet, y  = log10(SeqDepth))) +
      # facet_grid(.~diet, space = "free") +
      scale_fill_manual(values = nr_col) +
      coord_cartesian(ylim = c(4, 5.25)) +
      geom_boxplot( size = 0.9, outlier.shape = 21) +
      geom_point(aes(fill = diet), position = position_jitterdodge(jitter.width = 0.25), shape = 21, size = 1.5, alpha = 0.6) +
       # ggbeeswarm::geom_beeswarm(aes(fill = diet)) + 
      labs(y = "Sequencing Depth (log10 reads / sample)", x = "") +
      theme( axis.text.x = element_text(angle = 0))

    
    # ggpubr::ggarrange(ad1, ad2, ad3, ad4, ad5, common.legend = TRUE, legend = "right", ncol = 5, nrow = 1 )
    ad_plot <- ggpubr::ggarrange(ad1, ad2, ad3, ad4, ad5, common.legend = TRUE, legend = "right", ncol = 5, nrow = 1 )
    ad_plot
    
    ggsave(ad_plot, filename = "vis/ms__van_alpha-divx5.png", device = "png", width = 10, height = 6, bg = "white")
    
    
    # # differential abundance
    # # 
    # # DMM
    # # 
    # # network
    
