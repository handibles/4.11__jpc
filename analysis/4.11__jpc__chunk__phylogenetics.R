## ms__van__phylogenetics

  # rm( list = ls() )

  library('phangorn')
  library('phytools')   # source of persistent message?
  library('DECIPHER')
  
  source("analysis/ms__van__spreadsheeting__shortcuts.R")


# Make a New Tree from ASV sequences   ------------------------------

  # theres got to be a better way to be doing this, or a tleast a more informed way.
    seqs <- mgtax[,8] ; names(seqs) <- rownames(mgtax) ; head(seqs)    
  
  ## DECIPHER  -  get seqs from tax_table instead of from SeqTab
    alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA, processors = 9)

    gc()
    
  ## Phangorn  -  better done in DECIPHER?
    # print('Phangorn: - Using Max-Likelihood, JC69 and then GTR substitution models, Neighbour-Joining, GTR .. read up Jamie!')
    phang.align <- phyDat(as(alignment, "matrix"), type="DNA")
    dm <- dist.ml(phang.align)
    treeNJ <- NJ(dm) # Note, tip order != sequence order
    fit = pml(treeNJ, data=phang.align)
    fitGTR <- update(fit, k=4, inv=0.2)
    fitGTR <- optim.pml(fitGTR,                             ##     < ! >  >40min with 2k ASV!
                        model="GTR",                        ##            >10min with 1k ASV!
                        optInv=TRUE,                        ##            >30min with 940 ASV (adj)
                        optGamma=TRUE,
                        rearrangement = "stochastic",
                        control = pml.control(trace = 0)
    )
    
    saveRDS(fitGTR, 'output/ms__van_phangoGTR_NJ.RDS')
    saveRDS(fitGTR$tree, 'output/ms__van_phangoGTR_NJ_tree.RDS')

