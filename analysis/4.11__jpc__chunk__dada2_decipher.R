

  ## ref - Nov 2021
    # wget https://zenodo.org/record/4587955/files/silva_nr99_v138.1_train_set.fa.gz & 
    # wget https://zenodo.org/record/4587955/files/silva_species_assignment_v138.1.fa.gz & 
    # wget http://www2.decipher.codes/Classification/TrainingSets/SILVA_SSU_r138_2019.RData ; mv ./SILVA_SSU_r138_2019.RData ./decipher_SILVA_SSU_r138_2019.RData


  ##  from illumina protocol
    # F :: 5' 341F :: GCC TAC GGG NGG CWG CAG 
    # R :: 3' 805R :: GGA CTA CHV GGG TAT CTA ATC C


  ## and possibly 
    # conda activate jfg_r   # teaghpc


    library(dada2)  #; packageVersion("dada2")
    library(DECIPHER) #; packageVersion("DECIPHER")
  
    path <- "/home/jfg/data/4.11__jpc"               # project base point  ~~figFilt output, not cutadapt
    outpath <- "/home/jfg/data/4.11__jpc/4__dada2"
    list.files(path)
    list.files(outpath)

      
## params -------------------------------------------------------------

  # learnErrors - multithread = TRUE  
    le_nbases <- 3e8   # 1e9 = 42 samples, memfault
  
  # dada2 - multithread = TRUE 
    d_pool = FALSE    # "pseudo"  # TRUE
  
  # mergePairs
    mp_minOverlap = 10
    mp_maxMismatch = 0
  
  
  # removeBimerasDenovo  -  multithread = TRUE
    rbd_method <- "consensus"
  
  # assignTaxonomy
    at_trainset <- "~/data/ref/silva_nr99_v138_train_set.fa.gz"
  
  # assignSpecies - 100% identity
    as_trainset <- "~/data/ref/silva_species_assignment_v138.fa.gz"
    as_allowMultiple <- TRUE
  
  # DECIPHER  -  multithread = TRUE
    decipher_traindat <- "~/data/ref/SILVA_SSU_r138_2019.RData"
  
  # GTDB? - it is there, on the DECIPHER resources page.
  
  
## ====================================================================================== ##  
  
    ## the figaro step (& filterandTrim) are done from the CLI, and save in the 2__filt_d2 dir 
    
  ## figaro output
    fnFs <- sort(list.files( paste0(path, "/2__filt_d2"), pattern="R1.fastq.gz", full.names = TRUE))
    fnRs <- sort(list.files( paste0(path, "/2__filt_d2"), pattern="R2.fastq.gz", full.names = TRUE))
    head(fnFs) ; tail(fnRs)
  
    ##
    ##
    ##   intermediate step missing - plot the parameters chosen (figaro_diagnostic)
    ##
    ##
    

  # Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
    head( sample.names <- gsub(".*2__filt_d2\\/(.*)_d2fat.*", "\\1", fnFs, perl = TRUE) )

  # filt seq outputs
    filtFs <- gsub("2__filt_d2", "4__dada2", fnFs)
    filtRs <- gsub("2__filt_d2", "4__dada2", fnRs)
    
        
## filter sequences  ----------------------------------------------------
  
  ## manual filter : perform the FAT step natively.
        #   out <- filterAndTrim(fnFs, filtFs,
        #                        fnRs, filtRs,
        #                        maxN = 0,
        #                        maxEE = c(5,5),
        #                        truncQ = 2,
        #                        trimLeft = c(18,22),    # *
        #                        # truncLen = nfat_truncLen,    # *
        #                        rm.phix = TRUE,
        #                        compress = TRUE,
        #                        multithread=10 )
        # ## track pieces missing
    
  ## replaced in assembly with fAT-figaro step
    # - probably better to call figaro in here though... fack. 
    # - ...explain why? 
  
  ## write :    # figaro_to_r() ?

  
## get errors  ----------------------------------------------------------
    
    errF <- learnErrors(filtFs, multithread=8, randomize = TRUE, nbases = le_nbases )
    errR <- learnErrors(filtRs, multithread=8, randomize = TRUE, nbases = le_nbases )
    # plotErrors(errF, nominalQ=TRUE)
    # plotErrors(errR, nominalQ=TRUE)
    
    save.image() ; gc()
  
  
## dodo dada   ----------------------------------------------------------
  
  ## missing samples ?
    # sample.names <- unique(gsub("(.*)_L001.*", "\\1", list.files("~/data/4.11__jpc_raw"), perl = TRUE))
    # filt.names <- gsub("(.*)_F_.*", "\\1", list.files("~/data/4.11__jpc/2__filt_d2/", pattern = "_F_"), perl = TRUE)
    # 1
    # # length(filt.names) ; length(sample.names)
    # 
    # filt.names[ !sample.names  %in% filt.names]  

    dadaFs <- dada(filtFs, err=errF, multithread=TRUE, pool = d_pool) 
    dadaRs <- dada(filtRs, err=errR, multithread=TRUE, pool = d_pool)
    
    dadaFs[[10]]
    
    mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE,
                          minOverlap = mp_minOverlap,
                          maxMismatch = mp_maxMismatch)
    mergers[[1]]
    saveRDS(mergers, paste0(path, "/4__dada2/4.11__jpc__fAT-5-5__noFig__D2-mergers.RDS"))
    
    seqtab <- makeSequenceTable(mergers)
    dim(seqtab)  # jpc : 83 X 8925                xxx: 124 * 92827
    
    table(nchar(getSequences(seqtab)))
    hist(nchar(getSequences(seqtab)), breaks = 100)

    
  ## sensible length-based step
  
    
## == = ==   == = ==   == = ==   == = ==   == = ==   == = ==   == = ==   == = ==   == = ==   == = ==   == = ==   == = ==  
##   T A X O N O M I E S   =============================================================================================  
## == = ==   == = ==   == = ==   == = ==   == = ==   == = ==   == = ==   == = ==   == = ==   == = ==   == = ==   == = ==  
  
  
## slay 'meras   ---------------------------------------------------------
  
    seqtab.nochim <- removeBimeraDenovo(seqtab, method = rbd_method, multithread=TRUE, verbose=TRUE)
    dim(seqtab.nochim)  # xxx : 124 * 5559
    
    sum(seqtab.nochim)/sum(seqtab)   # 84% retained
  # par(mfrow=c(1,2))
  # hist(nchar(getSequences(seqtab.nochim)), breaks = 100)
  
  
## de-UCHIME ------------------------------------------------------------
      
  ## broken    
    
    # make fna 
    seqtab_seqs <- colnames(seqtab.nochim)
    fna_outfile <- paste0( path, "/5__outputs/4.11__jpc_seqtab_fig-0-0.fna")
    FIRST <- TRUE
    for (i in seq_along(seqtab_seqs)) {
      if (i == 2)  FIRST <- FALSE
          fna_outdata <- paste0(">seqtab_uchim__", i,"\n", seqtab_seqs[[i]], "\n")  # stringr::str_pad(i, 4, pad="0")
          if (FIRST == TRUE) {
            write(fna_outdata,file = fna_outfile, sep="", append = FALSE)
          } else {
            write(fna_outdata,file = fna_outfile, sep="", append = TRUE)
          }
    }

    
    # ## hugely inefficient, retune, as 80% of asvs unidientified (possibly reflects adapter seqs etc)  -   from 8_UCHIME.stripChim.R v0.1, and from PHYLO_16S_pipeline.sh
    # system( "usearch_64 -uchime_ref ~/Dropbox/MicroSupp/4.11__jpc/output/4.11__jpc_seqtab_fig-0-0.fna -db ~/data/ref/chimera_slayer_gold_db/rRNA16S_gold_v20110519.fasta -uchimeout ~/Dropbox/MicroSupp/4.11__jpc/output/4.11__jpc__uc_out_fig-0-0.uchime -strand plus -nonchimeras ~/Dropbox/MicroSupp/4.11__jpc/output/4.11__jpc__uc_good_fig-0-0.fasta &> ~/Dropbox/MicroSupp/4.11__jpc/output/4.11__jpc__uchimera_removal_fig-0-0.log", intern = FALSE)
    # ls()
    # uc <- read.table("~/Dropbox/MicroSupp/4.11__jpc/output/4.11__jpc__uc_out_fig-0-0.uchime")
    # uchimeras <- as.character(uc[ which( uc[,17] == "Y") , 2])
    # 
    # seqtab.nochim.uc <- seqtab.nochim[ , -as.numeric(gsub("seqtab_uchim__", "", uchimeras)) ]
    
    seqtab.nochim.fin <- seqtab.nochim
    
        
## track  - need to re-figure this for the additional CA / Fig steps   -------------------
    
        # getN <- function(x) sum(getUniques(x))
        # track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
        
        # # figfat_out <- do.call("rbind", lapply( list.files("~/data/4.11__jpc/1__qc/", pattern = ".RDS", full.names = TRUE), function(aa){readRDS(aa)}) )
        # track <- cbind( 
        #                 out, 
        #                 # figfat_out, 
        #                 sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
        # colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
        # rownames(track) <- sample.names
        # track <- as.data.frame(track)
        # track[,"substrate"] <- gsub("(.)\\d*_.*", "\\1", rownames(track))
        # track[,"substrate"] <- dplyr::recode(track$substrate, "M"="mucosa", "S"="stool")
        # head(track)
    
    save.image()
    gc()
  
    
##  make d2feat obj  ----------------------------------------------------
    
    # parse out, grab seqs, set names
    d2feat <- seqtab.nochim.fin
    mf_seqs <- colnames(d2feat)
    names(mf_seqs) <- paste0("asv__", stringr::str_pad(1:length(mf_seqs), 4, pad="0") )
    colnames(d2feat) <- paste0("asv__", stringr::str_pad(1:ncol(d2feat), 4, pad="0") )
    rownames(d2feat) <- gsub("(.*)_d2fat.*", "\\1", rownames(d2feat), perl = TRUE)

    dim(d2feat) ; length(mf_seqs)
      
  
## 100 %  taxonomic assignment --------------------------------------------
  
    d2_ranks <- c("d", "p", "c", "o", "f", "g", "s", "seq")
    
    taxa <- assignTaxonomy(seqtab.nochim.fin,
                           at_trainset,
                           multithread=6,
                           verbose = TRUE)
    taxa_sp <- addSpecies(taxa, 
                          as_trainset, 
                          allowMultiple = as_allowMultiple) ## 100% matching
    
    taxadf <- data.frame( taxa, "seq" = as.vector(getSequences(seqtab.nochim.fin)), stringsAsFactors = FALSE)
    colnames(taxadf) <- d2_ranks[1:6]
    rownames(taxadf) <- paste0("d2aT_", stringr::str_pad(1:nrow(taxadf), 4, pad="0"))
    
    taxa_spdf <- data.frame( taxa_sp, "seq" = as.vector(getSequences(seqtab.nochim.fin)), stringsAsFactors = FALSE)
    colnames(taxa_spdf) <- d2_ranks
    rownames(taxa_spdf) <- paste0("d2aS_", stringr::str_pad(1:nrow(taxa_spdf), 4, pad="0"))
    
    save.image()
    gc()
  
  
## or, DECIPHER taxonomic assignment ---------------------------------------
    
    # specific to converting to matrix
    dec_ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species", "seq")
    
    # worked out to be near exact same format
    # dna <- DNAStringSet(getSequences(seqtab.nochim.fin)) # Create a DNAStringSet from the ASVs
    dna <- DNAStringSet(as.vector(taxa_spdf[ , 8])) 
    load( decipher_traindat )
    
    gc() 
    
    ## if we can subset DB here to be the same as the amplified region, would be amaze and  == Q2
    
    ids <- IdTaxa(dna,
                  trainingSet,
                  minDescend = 1,  # default is 0.98, changed to 1 with respect to specification of 100% identity above
                  strand="top",
                  processors=NULL,
                  verbose=FALSE) # use all processors
    
    taxid <- data.frame(t(sapply(ids, function(x) {
      m <- match(dec_ranks, x$rank)
      taxa <- x$taxon[m]
      taxa <- gsub("unclassified_", "unk._", taxa)   # check out the startsWith fn!
      taxa
    })),
    as.vector(getSequences(seqtab.nochim.fin)),
    stringsAsFactors = FALSE)
    colnames(taxid) <- dec_ranks
    rownames(taxid) <- paste0("deci_", stringr::str_pad(1:nrow(taxid), 4, pad="0"))
    head(taxid[,1:6])
    
    dev.off()
    
    save.image()
    
  
##   output   ==========================================================
  
    exfil <- "~/Dropbox/4.11/4.11__jpc/output/4.11__jpc__overlap20"
    
    # feats
    saveRDS( d2feat, paste0(exfil, "4.11__jpc__feat_D2-seqtab.nochim.fin_fig-0-0.RDS") )
    
    # D2 tax and tax-spec method
    saveRDS( taxadf, paste0(exfil, "4.11__jpc__taxa_D2-assTax_fig-0-0.RDS"))
    saveRDS( taxa_spdf, paste0(exfil, "4.11__jpc__taxa_D2-assSpec_fig-0-0.RDS"))
    # data.frame( taxadf[20:40,4] , taxa_spdf[20:40,4])
  
    # DECIPHER tax method
    saveRDS( taxid, paste0(exfil, "4.11__jpc__taxa_DEC-IdTax_fig-0-0.RDS"))
    
    # # save the track object
    # saveRDS(track, "output/4.11__jpc__track_d2-FIGARO_fig-0-0.RDS")
    save.image( paste0(exfil, "4.11__jpc__dada2_deciph__fig-0-0.RData"))
    # load("output/4.11__jpc__overlap20/4.11__jpc__dada2_deciph__fig-0-0.RData")
    
    ## hand off to decontam   ...
  
