      
##    collection of metadata

      
      # rm(list=ls())
      
      # source("../../SeqBiome/sb_4.11_master/analysis/background_code/R__fns_jfg/fn_definitions.R")
      library("ggplot2")
      library("vegan")
      
      print("  + + +   loading   ===   V A N E S S A ' S   ===   data   + + +")
      
      
      dim( mgdat <- readRDS("output/ms__van__metadata_58_7.RDS")[ , ] )   #  -1
      dim( mgfeat <- ( readRDS("output/ms__van__mgfeat.RDS")[ , ] ))   # -1
      dim( mgtax <- readRDS("output/ms__van__mgtax.RDS") )
      
      dim(mgfeat_ra <- t(apply( mgfeat, 1, function(aa) aa/sum(aa))))


      bc_dist <- vegdist( mgfeat_ra, method = "bray")
      ja_dist <- vegdist( mgfeat_ra, method = "jaccard")
      
      mgdat$hcluster <- hclust( as.dist(bc_dist), method = "ward.D2")$order

      ## CLR
      mg_czm <- zCompositions::cmultRepl( mgfeat,  label=0,  method="CZM")    # samples must be ROWS ; returns samples as rows
      mg_clr <- t( apply(mg_czm, 1, function(x){log(x) - mean(log(x))}) )    # samples as ROWS ; returns samplas as rows - but should be gmean?
      dim(mg_clr)
      saveRDS(mg_clr, "output/ms__van__mgfeat-CLR_nonGMP.R")
      mg_clr <- readRDS("output/ms__van__mgfeat-CLR_nonGMP.R")

      
  ## plotting info  =========
      
      # colours, factors, whathavveyou  
      nr_col <- c("#1f78b4", "#33a02c", "#e31a1c", "#ff7f00", "#6a3d9a", "#FF1493",
                  "#b15928", "#737f33", "#8B008B", "#32fbd8", "#fdbf6f",
                  # RColorBrewer::brewer.pal(5, "Spectral")[1],
                  "#b2df8a", "#fb9a99", "#d9e627", "#EE82EE", "#DEB887",
                  "#a6cee3" )
      # scales::show_col(nr_col)  
      # # NA,
      # "white",
      # RColorBrewer::brewer.pal(5, "Spectral")[2:5])
      
      gg_color_hue <- function(n) {
        hues = seq(15, 375, length = n + 1)
        hcl(h = hues, l = 65, c = 100)[1:n]
      }
      var_cols <- gg_color_hue(6)
      names(var_cols) <- c("LF", "HF", "AHF", "D0", "D13", "D16")
      
      
## theme variables    -------------------------------------------------------------------
      
      print("applying theme_update in van_spreadsheeting_shortcuts")
      theme_set( theme_minimal())
      theme_update(
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.title.y = element_text(size = 14, colour = "grey20"),
        axis.title.x = element_text(size = 14, colour = "grey20"),
        axis.line = element_line(colour = "grey80", size = 0.2),
        #
        legend.box = "horizontal",
        legend.position = "right",
        legend.text = element_text(size = 9),
        legend.key.size = unit(0.75, "cm"),
        #
        panel.grid.major.x = element_line(colour = "lightskyblue3", size = 0.1),
        panel.grid.major.y = element_line(colour = "lightskyblue3", size = 0.1),
        panel.grid.minor.x = element_line(colour = "lightskyblue3", size = 0.1),
        panel.grid.minor.y = element_line(colour = "lightskyblue3", size = 0.1),
        panel.spacing.x = unit(1, "lines"),
        #
        plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"),  # if plotting class using super
        plot.tag = element_text(size = 25, face = "bold", margin = margin()),  #  margin = margin(l = 10)),  # if plotting class using super
        plot.title = element_text(size = 18, face = "bold"),  #  margin = margin(l = 10)),  # if plotting class using super
        plot.subtitle = element_text(size=14),
        plot.background = element_rect(fill = "white", colour = "white"),
        #
        strip.background = element_blank(),
        strip.text.x = element_text(angle = 0, size = 16),   # Top
        strip.text.y = element_text(angle = 0, size =14),
        NULL
      )
      
      
##  F U N C T I O N S   ===============================
      
      # one FN  - note use of scores
      ord_plot <- function(vegout, ax_lab = "axis", pdatf = mgdat, ptax = mgtax, pfeat = mgfeat, 
                           DO_SPEC = FALSE, DO_VARS = FALSE, alt_title = NULL, SHADE = FALSE, OMIT = TRUE, samples_are_rows = TRUE, 
                           tax_cent_thresh = 4, filt_asv = TRUE, taxon_level = "c"){
        
        # vegout <- ja_nmds_ft
        # pdatf <- mgdat
        # ptax <- mgtax
        # pfeat <- mgfeat
        # 
        
        ## testers
            # alt_title <- NULL
            # ax_lab <- "axis"
            # vegout =  ja_nmds_ra # feat_cah
            # pdatf = mgdat
            # alt_title <- "paed_cf mgfeat_ra BC"
            # pfeat = mgfeat
            # ptax = mgtax #[ -c(1,2), ]
            # OMIT = TRUE
            # samples_are_rows = FALSE
            # DO_SPEC = TRUE
            # DO_VARS = FALSE
            # SHADE <- FALSE
            # tax_cent_thresh = 3
            # taxon_level = "c"
            # filt_asv = TRUE
      
        
        require(ggrepel)
        require(ggvegan)
      
        ## caveats mf        
        if(DO_VARS){
          print("won't do VAR _and_ SPEC - setting to SPEC = FALSE")
          DO_SPEC <- FALSE
          }       
        if(DO_SPEC & ( !("species" %in% names(scores(vegout))) && !("species" %in% names(vegout)) ) ){ 
          stop("\n   + + +    asked for SPEC but no species in chosen ordivanion\n    + + +    see scores(vegout) and ?life.choices")
          }       
        if(DO_VARS & !("cca" %in% class(vegout) ) ){ 
          stop("\n   + + +    asked for VARS but no variables in chosen ordivanion (it's not constrained)\n    + + +    see class(vegout) and ?life.choices")
          }       
        
        if(!samples_are_rows){ pfeat <- t(pfeat) }    
        
        ptax <- as.matrix( ptax )
        
        
        ## handle vegan PCoA ELSE vegan RDA/CCA
        if( any ("matrix" %in% class(vegout) | grepl("MDS", class(vegout))) ){
          RD1 <- paste0( ax_lab, " 1")
          RD2 <- paste0( ax_lab, " 2")
          ord_df = data.frame(data.frame(pdatf[rownames(scores(vegout)),], stringsAsFactors = FALSE), 
                              "axis1" = scores(vegout)[,1],
                              "axis2" = scores(vegout)[,2])
        }else{      
          if( is.null( vegout[["CCA"]])){
            RD1 <- paste0( ax_lab, ' 1: ', round((vegout[["CA"]][["eig"]][1] / vegout[[ "tot.chi" ]])*100, 0),'%')
            RD2 <- paste0( ax_lab, ' 2: ', round((vegout[["CA"]][["eig"]][2] / vegout[["tot.chi" ]])*100, 0),'%')
          }else{
            RD1 <- paste0( ax_lab, ' 1: ', round((vegout[["CCA"]][["eig"]][1] / vegout[[ "tot.chi" ]])*100, 0),'%')
            RD2 <- paste0( ax_lab, ' 2: ', round((vegout[["CCA"]][["eig"]][2] / vegout[[ "tot.chi" ]])*100, 0),'%')
          }
          ord_df = data.frame(data.frame(pdatf, stringsAsFactors = FALSE), 
                              # mg_richness,
                              "axis1" = scores(vegout)$sites[,1],
                              "axis2" = scores(vegout)$sites[,2])
        }
        
        ##  !!!!!!!!!!!!  
        # ord_shape <- c(1,2,  16,17) ; names(ord_shape) <- c("F mucosa", "M mucosa",  "F stool", "M stool")
        # ord_df$plot_shape <- paste(ord_df$researcher, ord_df$Genotype)
        
        if(DO_SPEC){
          if( DO_SPEC & !(any(grepl("MDS", class(vegout)))) ){
            asv_cent <- data.frame(scores(vegout)$species)
          }else if( DO_SPEC & (any(grepl("MDS", class(vegout)))) ){
            asv_cent <- data.frame(vegout$species)
          }
          colnames(asv_cent) <- c("asvAxis1", "asvAxis2")
          # asv_cent$fill_var <- ptax[ rownames(asv_cent), "c"]    # broken
          asv_cent$fill_var <- sapply( rownames(asv_cent), function(aa){ ifelse( aa %in% rownames(ptax), ptax[ aa, taxon_level ], "non-defined taxon") })   # aa <- rownames(asv_cent)[30]
          asv_cent$size_var <- colMeans(pfeat[ , rownames(asv_cent)])
          ## better still to use conf ints, these v. small ad limited
          taxon_cent <- data.frame(
            aggregate( asvAxis1 ~ fill_var, FUN = mean, data = asv_cent)[ , ],
            "asvAxis2" = aggregate( asvAxis2 ~ fill_var, FUN = mean, data = asv_cent)[ , 2],
            "label_var" = apply(aggregate( . ~ fill_var, FUN = length, data = asv_cent)[ , 1:2], 1, function(aa){paste0(aa[[1]], " (n=", aa[[2]], ")" )}),
            "size_var" = aggregate( size_var ~ fill_var, FUN = sum, data = asv_cent)[ , 2]
          )
      
          tax_count <- sapply(unique(asv_cent$fill_var), function(aa){sum(asv_cent$fill_var == aa)}) 
          taxon_cent <- taxon_cent[ taxon_cent$fill_var %in% names( tax_count[ tax_count > tax_cent_thresh]) , ]
          tax_removed <- names( tax_count[ tax_count <= tax_cent_thresh])
          tax_removed_sub <- paste("omitted", length(tax_removed), "taxa from labels as incidence <", 4,":\n", paste(abbreviate(tax_removed, minlength = 7,strict = TRUE, dot = TRUE), collapse = ", "))
          print(tax_removed_sub)
          
          if( filt_asv ){
            # asv_cent$fill_var <- ifelse( !(asv_cent$fill_var %in% tax_removed), "other", asv_cent$fill_var)   # rename - doesn't work as banjaxes fill/colour overlap, but why?...
            asv_cent <- dplyr::filter( asv_cent, !(fill_var %in% tax_removed))   # remove altogether
            print( paste0("   + + +    hiding all features with incidence below ", tax_cent_thresh))
          }else{
            print("   + + +    showing all features, even if not in centroids - will clutter up legend\n   + + +    consider \\'filt_asv = TRUE\\'")
          }
          
        }
          
        
        if(is.null(alt_title)){plot_title <-  deparse(substitute(vegout)) }else{plot_title <- alt_title}
        
        
        ## plot limits
        plot_lim <- c(
          min(unlist(scores(vegout)), na.rm = TRUE),
          max(unlist(scores(vegout)), na.rm = TRUE)
        )
        
        
      ## ===    P L O T S   ========================================================================================== ##
        
        ## ---- A  
        pre_o_plot1 <- ggplot(ord_df, aes(x = axis1, y = axis2)) +  
          geom_hline(yintercept = 0, colour = "grey80")  +   # weight is the wrong arg
          geom_vline(xintercept = 0, colour = "grey80")  +
          geom_path( aes(group = sample), colour = "grey50", alpha = 0.1, size = 0.5) +
          coord_fixed(ratio = 1, xlim = plot_lim, ylim = plot_lim )
        
        ## ---- B  
        if(DO_SPEC){
          # pre_o_plot2 <- pre_o_plot1 + geom_point(data = asv_cent, aes(x = asvAxis1, y = asvAxis2, colour=fill_var, fill=fill_var, size = size_var), alpha = 0.8, shape = 21)
          pre_o_plot2 <- pre_o_plot1 + 
            geom_point(aes(x = axis1, y = axis2, shape = van_type), colour = "grey30", size = 1.5, alpha = 0.15) +
            geom_point(data = asv_cent, aes(x = asvAxis1, y = asvAxis2, colour = fill_var, fill = fill_var, size = size_var), shape = 20, alpha = 0.3) +    #, size = 1.5
            stat_ellipse(data = asv_cent, aes(x = asvAxis1, y = asvAxis2, colour=fill_var, fill=fill_var, size = size_var), geom = "polygon", size=0.25, level=0.5, alpha = 0.15) +
            geom_point(data = taxon_cent, aes(x = asvAxis1, y = asvAxis2, fill = fill_var, size = size_var), alpha = 0.8, shape = 21, colour = "black") +
            ggrepel::geom_label_repel(data = taxon_cent,
                                      aes(x = asvAxis1, y = asvAxis2, fill = NULL), 
                                      label = taxon_cent$label_var,
                                      size = 3,
                                      force =2,
                                      box.padding = 0.75,
                                      direction = "both",
                                      max.overlaps = Inf,
                                      min.segment.length = 0,
                                      segment.size = 0.5,
                                      hjust = 0.7,
                                      alpha = 0.85,
                                      label.r = unit(0.2, "cm"))
        }else{
          pre_o_plot2 <- pre_o_plot1 + 
            geom_point(aes(fill=van_type, shape=researcher, colour=van_type), size = 4, stroke=0.8) +    
            stat_ellipse(aes(x = axis1, y=axis2, fill=van_type, colour=van_type, lty=researcher), geom="polygon", size=0.25, level=0.8 , alpha=0.08, show.legend = TRUE) +   #, colour=event, colour=van_type
            geom_point(aes(fill=van_type, shape=researcher, colour=van_type), size = 4, stroke=0.8)
        }
        
        ## ---- C  is for centroids
        if("cca" %in% class(vegout) && !is.null(scores(vegout)$centroids) && !all(is.na(scores(vegout)$centroids))){    
          env_cent <- data.frame(scores(vegout)$centroids)
          colnames(env_cent) <- c("centAxis1", "centAxis2")
          
          if(DO_SPEC){
            
            pre_o_plot3 <- pre_o_plot2 +
              geom_point(data = env_cent,
                         aes(x = centAxis1, y = centAxis2, fill = ), 
                         size = 4,
                         alpha = 0.8,
                         shape = 13) #+
            # guides( fill = "none",
            #         colour  ="none")
            # ggrepel::geom_text_repel(data = env_cent,
            #                           aes(x = centAxis1, y = centAxis2, fill = NULL),
            #                           label = rownames(env_cent),
            #                           force_pull   = 0, # do not pull toward data points
            #                           nudge_y      = 7.5,
            #                           direction    = "x",
            #                           angle        = 0,
            #                           hjust        = 0,
            #                           segment.size = 0.05,
            #                           alpha = 0.7,
            #                           label.r = unit(0.2, "cm"))
            
            if(!OMIT){ pre_o_plot3 <- pre_o_plot3 + labs(subtitle = tax_removed_sub) }
            
          }else{
            
            pre_o_plot3 <- pre_o_plot2 +
              ggrepel::geom_label_repel(data = env_cent,
                                        aes(x = centAxis1, y = centAxis2, fill = NULL), 
                                        label = rownames(env_cent),
                                        max.overlaps = Inf,
                                        size = 3,
                                        alpha = 0.7,
                                        label.r = unit(0.2, "cm"))
          }
          
        }else{ 
          pre_o_plot3 <- pre_o_plot2
        }
        
        ## ----------------- D
        if(DO_VARS){
          if( "species" %in% names(scores(vegout)) ){
          arrow_mult <- ggvegan:::arrowMul( summary(vegout)$biplot,
                                            rbind(summary(vegout)$sites, summary(vegout)$species))
          }else if(!("species" %in% names(scores(vegout))) ){
          arrow_mult <- ggvegan:::arrowMul( summary(vegout)$biplot,
                                            summary(vegout)$sites )
          }
          #
          dd_cents <- dplyr::filter( fortify(vegout), Score == "centroids")
          ee_vars <- dplyr::filter( fortify(vegout), Score == "biplot", !(Label %in% dd_cents$Label))
          colnames(ee_vars)[3:8] <- paste0("Axis", 1:6)
          #
          ee_vars[,3:8] <- ee_vars[,3:8]*arrow_mult
          
          # cover with annot_rect to dampen samples  
          if(SHADE){   pre_o_plot3 <- pre_o_plot3 + annotate(geom = "rect", xmin=-10,xmax=10, ymin=-10, ymax=10, alpha = 0.5, fill = "white") }
          
          pre_o_plot4 <- pre_o_plot3 +
            geom_segment(data = ee_vars, aes(xend = Axis1, x = 0 , yend = Axis2, y = 0), colour = ifelse(grepl("fa", ee_vars$Label), "red", "blue"), arrow = arrow(length = unit(0.01, "npc")), lineend = "butt", linejoin= "bevel" ) +
            ggrepel::geom_text_repel(data = ee_vars, aes(x = Axis1, y = Axis2),
                                     label = ee_vars$Label,
                                     # colour = ifelse(grepl("fa", ee_vars$Label), "red", "blue"),
                                     max.overlaps = Inf,
                                     min.segment.length = 0,   # always draw
                                     direction = "both",
                                     box.padding = 0.5,
                                     force = 2,
                                     force_pull = 0,
                                     segment.size = 0.1,
                                     size = 2.5
            )
          
        }else{
          pre_o_plot4 <- pre_o_plot3
        }
        
        
        o_plot <- pre_o_plot4 +
          # scale_shape_manual("researcher x Genotype", values=ord_shape) +    # c(22,21)) sq/circ
          guides(lty = "none") +
          labs(
            title = plot_title )+
          # subtitle = "Ordivanion of community composition (Bray-Curtis PCoA, unconstrained)") +
          # alpha = guide_legend(override.aes = list(size = 6), title = "p < 0.05", nrow = 2),
          labs( x = RD1, y = RD2) 
        
        
        
        # (o_plot  )   #+ coord_fixed(ratio = 1)
        return( o_plot )
      }
      
      
      
      
      ## ===========================================================
      
