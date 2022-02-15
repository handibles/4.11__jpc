## ms__som__chunk__abundance_plots.R

    ## can catch if tax colnames dont match
    ## can catch if some samples empty
    ## can catch if no hclust

    # rm(list = ls()) 
  
    source("analysis/ms__van__spreadsheeting__shortcuts.R")

  
    plot_by <- "g"
    plot_n <- 15
    facet_by <- "p"
    facet_n <- 3       # currently doing nothing
    colour_by <- "o"
    feats <- mgfeat_ra
    hier <- mgtax
    moniker <- "other"
    na_term <- "unkn."
    
    ab_col <- sort(nr_col)

        
## crush to plottable ------------------------------------------

    
    top_Nplot <- names( sort(sapply( unique(hier[ , plot_by]), function(aa){    # aa <- hier[ 10 , plot_by]
      sum( feats[ , rownames(hier)[grep(aa, hier[,plot_by]) ] ], na.rm = TRUE)
    }), decreasing = TRUE)[ 1:plot_n])
    
    top_Nfacet <- names( sort(sapply( unique(hier[ , facet_by]), function(aa){
      sum( feats[ , rownames(hier)[grep(aa, hier[,facet_by]) ] ], na.rm = TRUE)
    }), decreasing = TRUE)[ 1:facet_n])

    glomf <- t(do.call("rbind", lapply( top_Nplot, function(aa){  # aa <- "UCG-008"     aa <- hier[ 10 , plot_by]
      if( sum( hier[ , plot_by] == aa) == 1 ){
      feats[ , hier[,plot_by]==aa ]
      }else{
      rowSums(feats[ , grep(aa, hier[,plot_by]) ]) }
    }))  )
    glomfo <- cbind( glomf, rowSums(feats) - rowSums(glomf))
    colnames(glomfo) <- c( top_Nplot, moniker)
    class(glomfo)
    dim(glomfo) ; glomfo[1:5, 1:5]


  ## hclust on the subset of data presented
    mgdat[ rownames(glomfo) , "new_hclust"] <- hclust( as.dist( vegdist(glomfo, method = "bray")), method = "ward.D2" )$order
    plot( mgdat$hcluster, mgdat$new_hclust, main = "clustering of total v. top 15 are much less different!")
    
  ## plot  
  
    mg_m <- reshape2::melt( as.matrix(glomfo[,]))
    mg_m$van_type <- mgdat[ mg_m$Var1 , "van_type"]
    mg_m$diet <- mgdat[ mg_m$Var1 , "diet"]
    mg_m$timepoint <- mgdat[ mg_m$Var1 , "timepoint"]
    mg_m$hcluster <- factor(mgdat[ mg_m$Var1 , "hcluster"], levels = mgdat$hclust)
    # mg_m$hcluster <- factor(mgdat[ mg_m$Var1 , "new_hclust"], levels = mgdat$hclust)
    mg_m$plot_var <- hier[ match( mg_m$Var2, hier[ , plot_by]) , plot_by]
    mg_m$facet_var <- hier[ match( mg_m$Var2, hier[ , plot_by]) , facet_by]
    mg_m[ is.na(mg_m[ , "plot_var"]) , "plot_var" ] <- moniker
    mg_m[ is.na(mg_m[ , "facet_var"]) , "facet_var" ] <- moniker
    head(mg_m)

    abund <- ggplot( mg_m, aes(fill = plot_var, x = hcluster, y = value)) +
      facet_grid(facet_var ~ timepoint + diet, space = "free", scale = "free") +
      geom_col(color = "black") +
      labs( xlab = "samples") +
      scale_fill_manual(values = sample(ab_col, length(nr_col)), "genus") + 
      theme(
        axis.text.x = element_text(angle = 90),
        strip.text.y = element_text(angle = 0)
      )
    abund
    
    
    ggsave(
      abund +
        labs( title = "16S populations by dietary condition (VanDat)",
              y = "relative abundance (%)",
              x = "Ward D2 clustered samples, Bray-Curtis"),
      filename = "vis/ms__van_basicAbundance__genus.png", device = "png", width = 12, height = 8)


## shade ranks
    # 
    # hier_plot <- mgtax[ mgtax[ , plot_by] %in% top_Nplot , ]
    # plot_cols <- shade_ranks2(TAXA = hier_plot, COLOUR_BY = colour_by, SHADE_BY = plot_by,
    #              MONIKER = moniker, NA_TERM = na_term)
    # scales::show_col(plot_cols)
    # 
    # plot_cols_sh <- plot_cols[ names(plot_cols) %in% c(top_Nplot, moniker, na_term) ]
    # 
    # mg_m$plot_var <- factor(mg_m$plot_var, levels = names(plot_cols))
    # ggplot( mg_m, aes(fill = plot_var, x = hcluster, y = value)) +
    #   facet_grid(facet_var ~ gender+exp.group, space = "free", scale = "free") +
    #   scale_fill_manual(values = plot_cols_sh) +
    #   geom_col(color = "black") +
    #   labs( xlab = "samples") +
    #   theme(
    #     axis.text.x = element_text(angle = 90),
    #     strip.text.y = element_text(angle = 0)
    #   )
    
    
    