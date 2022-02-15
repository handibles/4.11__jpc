      
      # rm(list=ls())
      
      source("analysis/ms__van__spreadsheeting__shortcuts.R")
      set.seed(2205)
      
      library(ggpubr)
      
      ## gotta be timepoints!

      
##  BC NMDS
      
      mg_nmds <- metaMDS( mgfeat_ra, distance = "bray", try = 200)
    
      mgdat_ed <- mgdat
      mgdat_ed$van_type <- mgdat_ed$diet
      
      bc_nmds <- ord_plot( mg_nmds, pdatf = mgdat_ed, alt_title = "BC-NMDS of Vanessa's rat microbiome samples")  
      
      bc_nmds_spec <- ord_plot( mg_nmds, pdatf = mgdat_ed, DO_SPEC = TRUE, alt_title = "BC-NMDS of Vanessa's rat microbiota")  
      
      # ggsave( bc_nmds, width = 12, height = 8, filename = "vis/ms__van_bc_ndms_samp.png")
      # ggsave( bc_nmds_spec, width = 9, height = 6, filename = "vis/ms__van_bc_ndms_spec.png")      

##  BC PCoA
      
      mg_pcoa <- cmdscale(  bc_dist)
    
      bc_pcoa <- ord_plot( mg_pcoa, pdatf = mgdat_ed, alt_title = "BC-PCoA of Vanessa's rat microbiome")  
  
      ## nope!    
      # ord_plot( mg_pcoa, DO_SPEC = TRUE)  
      
      # ggsave( bc_pcoa, width = 12, height = 8, filename = "vis/ms__van_bc_pcoa_samp.png")      

            
##  CA
      
      mg_cca <- cca( mgfeat_ra ~ 1, data = mgdat)
      
      ord_cca <- ord_plot( mg_cca, pdatf = mgdat_ed, alt_title = "CCA of Vanessa's rat microbiome samples")  
      
      ord_cca_spec <- ord_plot( mg_cca, pdatf = mgdat, DO_SPEC = TRUE, alt_title = "CCA of Vanessa's rat microbiota")
      
      # ggsave( ord_cca, width = 12, height = 8, filename = "vis/ms__van_cca_samp.png")      
      # 
      # ggsave( ord_cca_spec, width = 12, height = 8, filename = "vis/ms__van_cca_spec.png")      
      
      