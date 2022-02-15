
## ms__som__DA__LMM-EMM-plots_feats_boxplot

    rm(list=ls())
    
    library(ggplot2)
    library(vegan)
    source("analysis/background_code/R__fns_jfg/fn_definitions.R")
    source("analysis/ms__van__spreadsheeting__shortcuts.R")

    
    
## ============================================================

    test_df <- readRDS("output/ms__van__DA_EMM-12__kA-0.005-0.1.RDS")

    test_df$feat <- as.vector(test_df$feat)
    test_df$id <- as.vector(test_df$id)
    test_df$facet_by <- mgtax[ test_df$id , "c" ]
    
    head(test_df)

    test_df$stable <-  gsub("^ ", "", gsub(" $", "", gsub("\\.", "", paste( test_df$timepoint, test_df$diet), perl = TRUE ), perl = TRUE ), perl = TRUE )
    test_out <- test_df[ , c( "feat", "id", "facet_by", "contrast", "stable", "diet", "timepoint", "estimate", "t.ratio", "p.value", "FDR")]

        
    ## LOOK AT THE OUTPUT WITH YOUR EYES
    View(test_out)
    View(dplyr::filter(test_out, FDR < 0.05))
    
    # out of interest
    sum(test_out$FDR < 0.01)  # 228
    sum(test_out$FDR < 0.05)  # 305

    
    # get abundance values
    ms_ids <- unique(test_out$id)
    head(ms_idm <- reshape2::melt(as.matrix( mg_clr[ , ms_ids ])))
    colnames(ms_idm) <- c("sample", "id_melt", "clr")


## =============================================================================================

    ## complicate - three groups
      
      head(test_out)
      # dim(test_out_d0 <- dplyr::filter(test_out, grepl("D0", test_out$contrast), FDR < 0.05))
      dim(test_out_d0 <- dplyr::filter(test_out, FDR < 0.05))
    
        ##   
        ##   
        ##   h e e r e   -    i s h 
        ##   
        ##   
          #     dim(test_out_d13inter <- dplyr::filter(test_out, grepl("D0", test_out$contrast), FDR < 0.05))
          #     dim(test_out_d16inter <- dplyr::filter(test_out, grepl("D0", test_out$contrast), FDR < 0.05))
          #     dim(test_out_d13_d16 <- dplyr::filter(test_out, grepl("D0", test_out$contrast), FDR < 0.05))
          #     test_out_gc <- dplyr::filter(test_out, !grepl("D", test_out$contrast))
          # # check all correct  -  ... in what sense?.. Sensible?
          # # View( unique(rbind( test_out_fm, test_out_gc, test_out_ms) )) #[,c(4,9,10)]) ) 
    
    
    ## ind box plots  
  
      # p filter
      # test_out_fm_1 <- dplyr::filter(test_out_fm, FDR < 0.1)     # not 0.05!
      # test_out_gc_05 <- dplyr::filter(test_out_gc, FDR < 0.05)
      
      

##   control   ## ----------------------------------------------------------------------------
    
      head( test_out_d0)
      
                # # comparison( id( samples( condition() ) ) ) 
                #   test_fm_box <- do.call("rbind", lapply( unique(test_out_fm_1$stable), function(aa){   # aa <- " GDX"
                #       bb_df <- dplyr::filter(test_out_fm_1, stable == aa)
                #       bb_df <- rbind( bb_df, bb_df)
                #       bb_df$stable <- c("GDX", "Control")
                #       cc_id <- unique(bb_df$id)
                #       dd_ab <- do.call("rbind",lapply( cc_id, function(aaa){ ms_idm[ grep(aaa, ms_idm$id_melt) , ]  }))
                #       #
                #       # only show samples of interest in comparison, otherwise will plot eg F-GDX in M-Ctrl plot
                #       test_cond <- unlist(strsplit( unique(bb_df$stable), split = " "))
                #       test_samp <- dplyr::filter( mgdat, exp2 %in% test_cond)$sample     # substrate %in% test_cond, 
                #       dd_ab <- dplyr::filter(dd_ab, sample %in% test_samp )
                #       #
                #       dd_ab$stable <- aa
                #       dd_ab$FDR <- bb_df[ sapply( dd_ab$id_melt, function(aaaa) match( aaaa, bb_df$id)), "FDR"]   #grep(aaaa, bb_df$id)) , "FDR"]
                #       dd_ab$t.ratio <- bb_df[ sapply( dd_ab$id_melt, function(aaaa) match( aaaa, bb_df$id)), "t.ratio"]   #grep(aaaa, bb_df$id)) , "FDR"]
                #       #
                #       dd_ab$exp2 <- dplyr::recode("Control (diestrus)" = "Control", mgdat[ dd_ab$sample, "exp.group"])
                #       dd_ab$sex <- mgdat[ dd_ab$sample, "sex"]
                #       dd_ab$substrate <- mgdat[ dd_ab$sample, "substrate"]
                #       #
                #       dd_ab$plot_by <- mgtax[ as.character(dd_ab$id_melt), "g"]
                #       # dd_ab$facet_by <- mgtax[ as.character(dd_ab$id_melt), "c"]
                #       dd_ab$facet_by <- ifelse( dd_ab$t.ratio > 0, "higher in Female", "higher in Male")
                #       #
                # 
                #       dd_ab
                #     }))
                #   head(test_fm_box)
                #   # # View(test_fm_box)
    

    # do plot  -  the range of FDR here is 0.01 - 0.001 - not hugely worth the complication
      # boxplot_fm_st <- 
      ggplot( test_out_d0, aes(fill = diet, x = reorder(feat, t.ratio), y= t.ratio)) +    # , alpha = FDR < 0.05
        # facet_grid(facet_by~stable, scale = "free_y", space = "free_y", switch = "y") +
        facet_grid(facet_by~timepoint, scale = "free_y", space = "free_y") +
        coord_flip() +
        geom_hline(yintercept = 0, lty = 2, colour = "grey35", size = 0.25) +
        geom_boxplot(outlier.shape = NA, width = 0.75, size = 0.9) +
        # geom_point(shape = 21, size = 1.5, alpha = 0.4, position = position_jitterdodge(jitter.width =0.2) ) +
        scale_fill_manual(values = var_cols, limits = unique(test_out_d0$diet) ) + 
        # scale_alpha_manual(values = c(0.1, 0.5)) +
        theme(
          # axis.line.y = element_line(size = 1),
          axis.text.x = element_text(angle = 0, size = 9),
          legend.position = "right",
          # plot.margin = unit(c(0,2,0,0), "cm"),   # 2.2
          plot.subtitle = element_text(size = 8),
          panel.grid.major.y = element_line(colour = "grey70", size = 0.10),
          strip.background = element_blank(),
          strip.text.y = element_text(angle = 0,size = 14), #element_blank() #,
          strip.text.x = element_text(size = 14),
          text=element_text(size = 13)
        ) +
        labs(y = "log-ratio to mean sample abundance" ,  #"ratio of modelled abundances",
             x = "",
             title = "pairwise comparison of Gender: Female v. Male ",
             subtitle = "substrate samples only; FDR threshold of 0.1, prevalence threshold of 0.5% in 10% of samples (95 ASVs)",
             NULL)
    
  
  ## 
  ##  well, that's wrong! you'll need a pplot for each sub comparison. 
  ## 
      
  
      # boxplot_fm_st

      # ggsave(boxplot_exp_st_fm.f.m, filename = "vis/ms__som_boxplDA-stool-genus_exp_fm.f.m.png", device = "png", width = 8.5, height = 6)
      # # ggsave(boxplot_exp_st_spec, filename = "vis/ms__som_boxplDA-stool-species_exp.png", device = "png", width = 12, height = 7)
      # ggsave(boxplot_fm_st, filename = "vis/ms__som_boxplDA-stool_gender__FDR-0.1.png", device = "png", width = 9, height = 4.5)
          
        
        