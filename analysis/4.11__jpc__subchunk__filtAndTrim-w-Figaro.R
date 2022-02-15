#! /usr/bin/Rscript --vanilla --default-package=rjson,dada2

    library("rjson")
    library("dada2")
    
    ## parse figaro output via JSON to RDS, combine with script parameters, then fAT each sample separately

    ## first pre-UCHIME implementation gave twin peaks of ASVs, circa 250 ~ 400-420. 250 peak is of low abundance,
    ## set penalty to decrease maxEE
	
	
##   get args   ===========================================

    
    ## testers
      # bash_sampn <- "AMX1B"     # $TEST
      # bash_indir <- "/home/jfg/data/4.11__jpc/2__cutad"   # $CUTAD
      # bash_outdir <- "/home/jfg/data/4.11__jpc/"           # $PROJ
      # bash_Fprim <- 18                                                            # miseq standard
      # bash_Rprim <- 21                                                            # miseq standard
      # bash_maxEE_penalty <- c(
      #   0 ,
      #   0
      # )
    
    
    bash_args <- commandArgs(trailingOnly = TRUE)  
    # bash_args <- c("S23_S56", "~/data/raw/ms__som_raw", "~/data/ms__som", 18, 21)
    bash_sampn <- bash_args[1]
    bash_indir <- bash_args[2]
    bash_outdir <- bash_args[3]
    bash_Fprim <- as.numeric( bash_args[4] )
    bash_Rprim <- as.numeric( bash_args[5] )
  # to harsh out "twin peaks" in length(nseq)
    bash_maxEE_penalty <- c( as.numeric( bash_args[6] ) ,  as.numeric( bash_args[7] ) )
    
    
##  figaro out  ===========================================

  ## test case
    eg_json <- "/home/jfg/data/4.11__jpc/1__qc/figout/AMX1B/trimParameters.json"
    figout <- do.call("rbind", lapply( rjson::fromJSON(file = eg_json ), function(aa){  unlist(aa) }) )
    best_fig <- unlist(figout[ which(figout[,6] == max( figout[,6])) , ] )
    
    figout <- do.call("rbind", lapply( rjson::fromJSON(file = paste0( bash_outdir, "/1__qc/figout/", bash_sampn, "/trimParameters.json") ), function(aa){  unlist(aa) }) ) 
    best_fig <- unlist(figout[ which(figout[,6] == max( figout[,6]))[1] , ] )  # careful you only grab one "best" in case of ties!


##   dada2::filtAndTrim  =================================
    
    fnFs <- paste0( bash_indir, "/", bash_sampn, "_cutad_R1.fastq.gz")   # had ot strip the "001" from the "...001.fastq.gz" for figaro - revise at some point
    fnRs <- paste0( bash_indir, "/", bash_sampn, "_cutad_R2.fastq.gz")
    filtFs <- paste0( bash_outdir, "/2__filt_d2/", bash_sampn, "_d2fat_R1.fastq.gz")
    filtRs <- paste0( bash_outdir, "/2__filt_d2/", bash_sampn, "_d2fat_R2.fastq.gz")
    
    fat_truncQ <- 2
    fat_truncLen <- c(best_fig[1], best_fig[2])
    fat_maxEE <- c(best_fig[3], best_fig[4]) - bash_maxEE_penalty
    fat_trimLeft <- c(bash_Fprim, bash_Rprim)
    
    out <- filterAndTrim(fnFs, filtFs,
                         fnRs, filtRs,
                         maxN = 0,
                         maxEE = fat_maxEE,    # *
                         truncQ = fat_truncQ,
                         trimLeft = fat_trimLeft,    # *
                         truncLen = fat_truncLen,    # *
                         rm.phix = TRUE,
                         compress = TRUE,
                         # multithread=11)   # already parallel
                         multithread=1 )
    # for track  
    saveRDS( c( bash_sampn, out), paste0(bash_outdir, "/1__qc/dada2_filtAndTrim_out_for_track_", bash_sampn, ".RDS"))
    print( c( bash_sampn, out)  )  
    
    