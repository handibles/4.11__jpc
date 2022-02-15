
## de-chimera, de-contaminate   =================================================================

    rm(list=ls())
  
    library('ggplot2')
    library('dada2')
    library('decontam')
            # library('scales')        # hue_pal
            # library('RColorBrewer') # brewer.pal
                        # library('reshape2')     # melt 
                        # library('stringr')
                                          # library('vegan')
                                          # library('dada2')        # sequencetable 
                                          # library('knitr')        # kable
                                          # library('ShortRead')    # sread
                                        # library('ggpubr')       # ggarrange
                                          # library('zCompositions')       # cMultRepl
                                          # library('robCompositions')     # gmean_sum
                                          # 

## gets, get fna   ---------------------------------------------------------------
    
    # exfil <- "/data/Food/primary/R0602_microsupport/ms__natval/5__outputs/"
    exfil <- "output/"
    d2tax <- readRDS( paste0( exfil, "ms__xxx__taxa_D2-assSpec_fig-0-0.RDS"))
    d2feat <- readRDS( paste0(exfil, "ms__xxx__feat_D2-seqtab.nochim.fin_fig-0-0.RDS"))
    
    d2tax_seqs <- as.vector( d2tax[ , 8]) 
    names(d2tax_seqs) <- rownames(d2tax)

    # !!!        
    rownames(d2tax) <- gsub("d2aS", "asv_", rownames(d2tax))
    
    d2feat[1:10, 1:7]
    d2tax[1:10, 1:7]
    dim(d2tax) ; dim(d2feat)
          

### How long are my reads?  ------------------------------------------------------------------------

    ## Know your primers, and your adapters.
    ## serious issue with adapter contamination outside of 400-405 and 410-420 peaks.
    ## Huge improvement post-cutadapt (as much as feasible, passing QC), but still contaminated. 
    ## Cut remainder by length, taxon
    
    # replace NA in taxonomy
    d2tax[ is.na(d2tax)] <- "unkn."
    
    ## split to STUDIES    
    d2feats_N <- d2feat[ grepl( "^A|^B|^C|^D|^E", rownames(d2feat), perl = TRUE) , ]
    d2feats_V <- d2feat[ grepl( "Mic|16S", rownames(d2feat), perl = TRUE) , ]
    # remove singletons :: 3829 / 1265
    dim(d2tax_N <- d2tax[ colSums(d2feats_N) > 1 , ])
    dim(d2tax_V <- d2tax[ colSums(d2feats_V) > 1 , ])
    dim(d2feats_N <- d2feats_N[ , colSums(d2feats_N) > 1 ])
    dim(d2feats_V <- d2feats_V[ , colSums(d2feats_V) > 1 ])

        
    # png( paste0( exfil, "ms__xxx_ASV_lengthXsubstrate_n100.png"), width = 1000, height = 800)
    par(mfrow = c(3,2))
    hist(sapply( d2tax[,8], function(x) nchar(x) ), breaks = 500,
         main = "vast majority of ASVs within expected range\n(400-430bp)")
    
    hist( sapply( d2tax[ colSums(d2feat)  > 100,8], function(x) nchar(x) ), breaks = 500,
          main = "filtering (abu. >100) retains 99.1% of abundances \nwhile removing low-copy ASVs")
    round(sum(d2feat[ , colSums(d2feat)  > 100]) / sum(d2feat), digits = 3)
    
    
    ## plot nseq  /  sum(feat) -  like hist but better
    phy_col <- RColorBrewer::brewer.pal(12, "Paired")
    names(phy_col) <- unique(d2tax[,2])
    # par(mfrow = c(1,2))
    
    
    plot( sapply( d2tax_N[,8], nchar), log10(colSums(d2feats_N)), col = phy_col[ d2tax_N[,2]], pch = 18, main = "total abundances of Natalia's\n ASVs (n = 4656) by length", xlab = "16S amplicon length", ylab = "log10 total abundance", xlim = c(250, 430), ylim = c(0,6) )
    legend("topleft", legend = names(phy_col), col = phy_col, pt.bg = phy_col, pch = c(22), pt.cex = 2, cex =0.8, bty = "n" )
    
    plot( sapply( d2tax_V[,8], nchar), log10(colSums(d2feats_V)), col = phy_col[ d2tax_V[,2]], pch = 18, main = "total abundances of Vanessa's\n ASVs (n = 1969) by length", xlab = "16S amplicon length", ylab = "log10 total abundance", xlim = c(250, 430), ylim = c(0,6) )
    legend("topleft", legend = names(phy_col), col = phy_col, pt.bg = phy_col, pch = c(22), pt.cex = 2, cex =0.8, bty = "n" )
    # dev.off()    
    
    
    # par(mfrow = c(1,2))
    
    fil_feat_V <- d2feats_V[ , colSums(d2feats_V)  > 100]
    fil_tax_V <- d2tax[ colnames(fil_feat_V) , ]
    fil_feat_N <- d2feats_N[ , colSums(d2feats_N)  > 100]
    fil_tax_N <- d2tax[ colnames(fil_feat_N) , ]
    
    plot( sapply( fil_tax_N[,8], nchar), log10(colSums(fil_feat_N)), col = phy_col[ fil_tax_N[,2]], pch = 18, main = "abundances of Natalia's\n ASVs (ab>100) by length", xlab = "16S amplicon length", ylab = "log10 total abundance", xlim = c(389, 430), ylim = c(1.8,6) )
    legend("topleft", legend = names(phy_col), col = phy_col, pt.bg = phy_col, pch = c(22), pt.cex = 2, cex =0.8, bty = "n" )
    
    plot( sapply( fil_tax_V[,8], nchar), log10(colSums(fil_feat_V)), col = phy_col[ fil_tax_V[,2]], pch = 18, main = "abundances of Vanessa's\n ASVs (ab>100) by length", xlab = "16S amplicon length", ylab = "log10 total abundance", xlim = c(389, 430), ylim = c(1.8,6) )
    legend("topleft", legend = names(phy_col), col = phy_col, pt.bg = phy_col, pch = c(22), pt.cex = 2, cex =0.8, bty = "n" )

    # dev.off()

        
    # check sequence complexity - is there bimodalism, reflecting real and artefactual sources of variation?  
    # --->  no.
    # png(paste0( exfil, "ms__xxx_ASV_complexityXlength.png"), width = 800, height = 600)
    hist(seqComplexity(d2tax[,8]),breaks=500, xlim = c(10,16), ylim = c(0,60), 
         main = "total sequence information complexity")
    hist(seqComplexity(d2tax[ nchar(d2tax[,8]) < 398, 8]), breaks = 500, xlim = c(10,16), ylim = c(0,60), 
         main = "complexity of spurious seqs")
    hist(seqComplexity(d2tax[ nchar(d2tax[,8]) > 398, 8]), breaks = 500, xlim = c(10,16), ylim = c(0,60),
         main = "complexity of real seqs")
    plot( seqComplexity(d2tax[ ,8]), nchar(d2tax[,8]), col = phy_col[ d2tax[,2]], 
        main = "note ~uniform dist of real/spur. and complexity" )
    # dev.off()           
    
    
## remove by dada2 taxonomy ---------------------------------------------------------------
    
    par(mfrow=c(1,1)) # dev.off()

  # crucial reference, check out the bad stuf: mainly unkn. at d/p, or Rickettsia-Mitochondria/Chloroplasts
    bad_range <- c(0:399) # , 405:418, 422:424)
    dim(badt <- d2tax[ nchar(d2tax[,8]) %in% bad_range , ])   # 2139
    badt[,8] <- nchar(badt[,8])
    # some real stuff, way out of range, exclude
    apply(badt[,1:7], 2, unique)
    # View(badt)
    
  # from above, keep 400;426 -  cherry-picking length too biased ; c(400:404, 419,420,421, 425,426)
    good_range <- c(400:428)  #
  # then remove taxonomic nonsense
    dim(subt <- d2tax[ nchar(d2tax[,8]) %in% good_range , ])   # 1498
    dim(subt <- dplyr::filter(subt, f!="Mitochondria", o!="Chloroplast", p!="unkn.", d!="unkn."))

    
    subt[,9] <- nchar(subt[,8])
    apply(subt[,1:7], 2, unique)
    dim(subf <- d2feat[ , rownames(subt)])
    
    
      ## same same?
          # subf_N <- subf[ grepl( "^A|^B|^C|^D|^E", rownames(subf), perl = TRUE) , ]
          # subf_V <- subf[ grepl( "Mic|16S", rownames(subf), perl = TRUE) , ]
          # # remove singletons :: 3829 / 1265
          # dim(subt_M <- subt[ colSums(subf_N) > 1 , ])
          # dim(subt_S <- subt[ colSums(subf_V) > 1 , ])
          # dim(subf_N <- subf_N[ , colSums(subf_N) > 1 ])
          # dim(subf_V <- subf_V[ , colSums(subf_V) > 1 ])
          # 
          # # png( paste0( exfil, "ms__xxx_ASV_deadapted_detax_split.png"), width = 800, height = 600)
          # par(mfrow=c(1,2))
          # plot( subt_M[,9], log10(colSums(subf_N)), col = phy_col[ subt_M[,2]], pch = 18, main = "total abundances of Natalia's\n ASVs by length", xlab = "16S amplicon length", ylab = "log10 total abundance" , ylim = c(0,6))
          # legend("topleft", legend = names(phy_col), col = phy_col, pt.bg = phy_col, pch = c(22), pt.cex = 2, cex =0.8, bty = "n" )
          # plot( subt_S[,9], log10(colSums(subf_V)), col = phy_col[ subt_S[,2]], pch = 18, main = "total abundances of Vanessa's\n ASVs by length", xlab = "16S amplicon length", ylab = "log10 total abundance", ylim = c(0,6) )
          # legend("topleft", legend = names(phy_col), col = phy_col, pt.bg = phy_col, pch = c(22), pt.cex = 2, cex =0.8, bty = "n" )
          # # dev.off()
    
    # View(subt)
    
    
## doublecheck - remove by decipher taxonomy  ---------------------------------------------
    
    d2tax_d <- readRDS( paste0( exfil, "ms__xxx__taxa_DEC-IdTax_fig-0-0.RDS") )
    rownames(d2tax_d) <- gsub("deci", "asv_", rownames(d2tax_d))
    
  # crucial comparison
    baddec <- d2tax_d[ !(rownames(d2tax_d) %in% rownames(subt)) , ]
    baddec[,9] <- nchar( baddec[,9])
    # View(baddec)
    
  ## decipher picks up all <399 as "unk._Root", consistent across df
    subdec <- d2tax_d[ rownames(subt) , ]
    dim(subdec <- subdec[ subdec$domain !="unk._Root" , ])
    dim(subdec <- subdec[ subdec$phylum !="unk._Bacteria" , ])
    
    length(deci_keep_asvs <- rownames(subdec))
    dim(d2feat <- d2feat[ , deci_keep_asvs])
    

## plotting data  ----------------------------------    
    
    # mgdat <- readRDS("output/ms__xxx__metadata_65-11.RDS")
    mgdat <- read.table("input/ms__xxx_metadat.csv", sep = ",", header = TRUE) ; 

    
  ## fix names  ~~~~~~~   
    rownames(mgdat) <- mgdat$sample
    rownames(d2feat) <- gsub("_S.*", "", gsub("\\-", "_", rownames(d2feat)), perl = TRUE)
    all( rownames(d2feat) %in% rownames(mgdat) )
    all(rownames(mgdat) %in%  rownames(d2feat) )
    head(mgdat)
    # env colours...    
    # taxa colours...
    

## decontam -----------------------------------------    

    mgdat$IsNeg <- ifelse( grepl("Neg", mgdat$sample), "NegCtrl", as.character(mgdat$researcher) )
    mgdat$SeqDepth <- rowSums( d2feat[ rownames(mgdat) , ])
    ggplot(data=mgdat, aes(x=reorder(sample, SeqDepth), y=SeqDepth, color=researcher)) + 
      geom_point()
      
    # library(decontam)
  # # Method A - Frequency :: DNA Concentrations from Teagasc? 
    # dat_0 <- mgdat[ mgdat$IsNeg != "NegCtrl" , ]
    # ft_0 <- d2feat[ dat_0[,"sample"] , ]
    # contamdf.freq <- isContaminant(ft_0, method="frequency", conc=dat_0$dna_ngul, threshold = 0.1)
    # head(contamdf.freq)
    # sum(contamdf.freq$contaminant)

    # # nice! only plots contaminants
    # plot_frequency(ft_0, (colnames(ft_0)[contamdf.freq$contaminant]), conc=dat_0$dna_ngul ) + 
    #   xlab("DNA Concentration (raw extract qubit)")
    
  # Method B - Prevalence :: only one neg ctrl so not really going to help much here.
    contamdf.prev <- isContaminant(d2feat, method="prevalence", neg=(mgdat$IsNeg=="NegCtrl"))
    table(contamdf.prev$contaminant)
    # 56 true, 3534 false 

  # decontam
    decontam_keep_asv <- rownames(contamdf.prev)[!contamdf.prev$contaminant]
      
          
## WRITE OUT DATA ===============================================================

    dim( mgdat <- dplyr::filter(mgdat, IsNeg!="NegCtrl")[, -5] )
    saveRDS(mgdat, "output/ms__xxx__metadata_122_5.RDS")
    d2feat <- d2feat[ mgdat$sample, ]
    
    mgfeat <- d2feat[ , intersect( decontam_keep_asv, deci_keep_asvs)]
    mgtax_s <- d2tax[ intersect( decontam_keep_asv, deci_keep_asvs) , ]
    mgtax_d <- d2tax_d[ intersect( decontam_keep_asv, deci_keep_asvs) , c(1:7,9)]
    colnames(mgtax_d) <- c(colnames(mgtax_d)[1:7], "seq")
    
    saveRDS(mgfeat, "output/ms__xxx__mgfeat.RDS")
    saveRDS(mgtax_s, "output/ms__xxx__mgtax.RDS")
    saveRDS(mgtax_d, "output/ms__xxx__mgtax_decipher.RDS")
        
    
    samp_thresh <- 1000
    spec_thresh <- 1
    print( paste0("  + + +   remove samples/tax   ===   B E L O W   ", samp_thresh, " / ", spec_thresh, "   ===    + + +"))

    
  ## Natalia
    
    # # some 0 / near-0 samples, but comparison with van shows that this is a researcher effect. Sucks for them, but saves time. ⎄
    # plot( rowSums(mgfeat), col = c(2:3)[ as.factor(mgdat$researcher)])
    
    mgdat_N <- dplyr::filter(mgdat, researcher == "Natalia")
    mgdat_N$timepoint <- gsub(".*(T.*)", "\\1", mgdat_N$sample)
    mgfeat_N <- mgfeat[ rownames(mgdat_N), ]
    
    mgfeat_N <- mgfeat_N[ rowSums(mgfeat_N) > samp_thresh , ]
    mgfeat_N <- mgfeat_N[ , colSums(mgfeat_N) > spec_thresh ]
    mgdat_N <- mgdat_N[ rownames(mgfeat_N) , ]
    d2tax_Ns <- mgtax_s[ colnames(mgfeat_N) , ]
    d2tax_Nd <- mgtax_d[ colnames(mgfeat_N) , ]
    
    dim(mgfeat_N)
    dim(mgdat_N)
    
    # View(mgdat_N) ; dim(mgdat_N)
    saveRDS(mgdat_N,  "output/ms__nat__metadata_57_6.RDS")
    saveRDS(mgfeat_N, "output/ms__nat__mgfeat.RDS")
    saveRDS(d2tax_Ns, "output/ms__nat__mgtax.RDS")
    saveRDS(d2tax_Nd, "output/ms__nat__mgtax_decipher.RDS")
    
    
  ## Vanessa
    
    mgdat_V <- dplyr::filter(mgdat, researcher == "Vanessa")
    mgdat_V$timepoint <- gsub(".*(D.*)_.*", "\\1", mgdat_V$sample)
    mgdat_V$diet <- gsub("Mic.*", "C", gsub(".*_(.*F)_?.*", "\\1", mgdat_V$sample) )
    mgfeat_V <- mgfeat[ rownames(mgdat_V), ]
    
    mgfeat_V <- mgfeat_V[ rowSums(mgfeat_V) > samp_thresh , ]
    mgfeat_V <- mgfeat_V[ , colSums(mgfeat_V) > spec_thresh ]
    mgdat_V <- mgdat_V[ rownames(mgfeat_V) , ]
    d2tax_Vs <- mgtax_s[ colnames(mgfeat_V) , ]
    d2tax_Vd <- mgtax_d[ colnames(mgfeat_V) , ]
    
    dim(mgfeat_V)
    dim(mgdat_V)
    
    # View(mgdat_V) ; dim(mgdat_V)
    saveRDS(mgdat_V, "output/ms__van__metadata_58_7.RDS")
    saveRDS(mgfeat_V, "output/ms__van__mgfeat.RDS")
    saveRDS(d2tax_Vs, "output/ms__van__mgtax.RDS")
    saveRDS(d2tax_Vd, "output/ms__van__mgtax_decipher.RDS")
    
    
##  nreads summary  -  plot track object  -------------------------------------------------------
    
  ## update for decontamination
    
    # track <- as.data.frame(readRDS('output/ms__xxx__track_d2-FIGARO.RDS'))
    # track$substrate <- gsub("(^.).*", "\\1", rownames(track))
    # track$substrate <- dplyr::recode( track$substrate, "M"="mucosa", "S"="stool")
    # track$substrate <- ifelse(grepl("[MS]0_.*", rownames(track), perl = TRUE), "neg ctrl", track$substrate) 
    # head(track)
    # 
    # z.rm <- reshape2::melt(track)   # saveRDS(z.rm , '../6_PHYLO_Seqtab/tracK_z.rm.RDS') # z.rm <- readRDS()
    # highest_count <- round( max(track[,1:6]), digits = -4)
    # 
    # # outlier labels: https://stackoverflow.com/a/43659981
    # 
    # ggplot(z.rm, aes(y=value, x=variable,  fill=as.factor(substrate) ,  colour=as.factor(substrate) )) +
    #   geom_boxplot(outlier.shape=19 , colour='black') + 
    # #  geom_point(width=0.2, shape=21,  size=2.5,  aes(colour='black' , fill=as.factor(matrix), alpha=0.6 )) + #note factoring-by-source requires use of the aes arg
    #   scale_fill_manual( "substrate",values=c("#1e8bc3",  "#C30000",  "#ffb90f", 'darkslategray') ) +         
    #   scale_color_manual("substrate",values=c("#1e8bc3",  "#C30000",  "#ffb90f", 'darkslategray') ) +
    #   # scale_y_continuous(limits = c(0, highest_count*1.05)) +
    #   scale_y_continuous(limits = c(0, 500000)) +
    #   labs(title='Read-counts through pipeline' ) + 
    #   xlab('pipeline Step') + 
    #   ylab('number of reads') + 
    #   theme_classic() +
    #   theme(text = element_text(size = 14),
    #         axis.text.x = element_text(angle = 90, vjust = -0.005))
    # 
    # # knitr::kable(track[,1:6])
