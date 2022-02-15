#!/usr/bin/Rscript --vanilla --default-package=rjson

  ##  fallow, as function wrapped into filterAndTrim call, but left for autonomy

      # ## test case
      #   # eg_json <- "/home/handibles/Dropbox/MicroSupp/ms__som/output/som__figaro/M6_S6/trimParameters.json"
      #   # figout <- do.call("rbind", lapply( rjson::fromJSON(file = eg_json ), function(aa){  unlist(aa) }) ) 
      #   # unlist(figout[ which(figout[,6] == max( figout[,6])) , ], use.names = FALSE )
      #   
      # ## get BEST figaro sumamry from dir @ arg[1]:
      #   library("rjson")
      #   bash_args <- commandArgs(trailingOnly = TRUE)  
      #   bash_path <- bash_args[1]
      #   bash_sample <- bash_args[2]
      #   
      #   figout <- do.call("rbind", lapply( rjson::fromJSON(file = paste0( bash_path, "/", bash_sample, "/trimParameters.json") ), function(aa){  unlist(aa) }) ) 
      #   saveRDS(( unlist(figout[ which(figout[,6] == max( figout[,6])) , ] ) ), paste0("output/som_figaro/", bash_sample, "/figParams.RDS") )

    