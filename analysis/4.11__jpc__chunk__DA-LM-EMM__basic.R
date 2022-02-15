## ms__van__chunk__diff_ab
      

      rm(list=ls())
      source("analysis/background_code/R__fns_jfg/fn_definitions.R")
      source("analysis/ms__van__spreadsheeting__shortcuts.R")
      set.seed(2205)
      
      library(lme4)
      library(car)
      library(emmeans)
      library(e1071)
      
      

##  outline comparisons   =================================================================


## JKEANE  ----------------------------------------------------
        # 
        # AHF = Akkermansia High-Fat protective diet, HF = high far, LF = low fat
        # expectation is that Ak will be protective against metabolic syndrome outcomes
        # 
        # Consideration: D0:D13, D0:D16 comparisons, and then D13:D16 comparisons
        # 
      # -------------------------------------------------------
    

## prep  ---------------------------------------------------------------------

      # filter at 0.5% in 10%
      dim(mg_ka <- t(k_A( t(mgfeat_ra), k = 0.005, A = 0.1)))
      # filter at all?...
      dim(clr_df <- mg_clr[ rownames(mgdat), colnames(mg_ka)]) 
      
      ## not sure, but possibly need to set the factor levels such that D0:C is first as the control for comparisons
      mgdat$timepoint <- factor( mgdat$timepoint, levels = c( "D0", "D13", "D16"))
      mgdat$diet <- factor( mgdat$diet, levels = c( "C", "LF", "HF", "AHF"))

## explore ----------------------------------------------------------------
      
    ##  . . . 
      
    ## check model assumptions    -------------------------------------------
          
          # # 1 - independence of samples  - judge for self
          # 
          # # 2 - linearity
          # 
          # # 4 - homoskedasticity
          # 
          # # 3 - non co-linear FEs
          # 
          # # 5  -normal resids  -  ehh..
          # 
          # # 6 - no extreme values / outliers -  none are extreme
      
      
## test  ----------------------------------------------------------------------
  
    ## need to compare the different timepoints:  D0:D13, D0:D16 comparisons, and then D13:D16 comparisons
      
      
      
      test_feats <- colnames(clr_df)  # [1:10]
      
      # test_df <- do.call("rbind", parallel::mclapply( test_feats, function(aa){     # aa <- test_feats[30]
      test_df <- do.call("rbind", lapply( test_feats, function(aa){     # aa <- test_feats[30]
      
          feat_dat <- mgdat
          t_feat <- clr_df[ , aa]
                                  # feat_dat$binAbs <- factor(ifelse(t_feat > 0, "pres", "abs"))
          
        ## attempt to model sparseness using RanEf (see subchunk) hugely ineffective, deprecated
        ## basic but effective
          feat_test <- lm( t_feat ~ diet*timepoint, data=feat_dat )
          
        ## model-fit diagnostics
          # hist(residuals(feat_test))
          shap.test <- shapiro.test(residuals(feat_test))
          lev.test  <- leveneTest(residuals(feat_test), group=feat_dat$diet)   #  change?
          
          
        # see subchunk_EMM_config. original EMM_feat version used pairs in the exact same manner   
          (emm.int <- emmeans(feat_test, specs = ~ diet + timepoint))
          (contr.int <- as.data.frame( contrast(emm.int, 
                                                combine = TRUE,        # cobine text outp to a df outp
                                                simple = list("diet", "timepoint"),
                                                type = "response",
                                                interaction =  FALSE,  # if true will estimate just the EFFECT OF THE INTERACTION, == s:e:s, a single term
                                                method = "pairwise",
                                                by = NULL,             # "blocks" to work in, as removed/sep from testing - set to NULL
                                                #
                                                offset = NULL,
                                                name = NULL,
                                                adjust = "none",
                                                ratios = NULL,  
                                                reverse = TRUE   ) ) )
      
          bb <- cbind( "feat" = mgtax[ aa, 6],
                       "id" = aa,
                       as.data.frame(contr.int),
                       "shapiro_w" = shap.test$statistic,
                       "shapiro_p" = shap.test$p.value,
                       "levene_F" = lev.test[1,2],
                       "levene_p" = lev.test[1,3])
          
          cc <- bb[ !apply(bb, 1, function(aaa){ any( is.na(aaa)) }) , ]
          colnames(cc)[5] <- "contrast"

            
        ## D0-D13/D16 comparisons, we have a single Control case - use the appropriate function
          dd_ctrl <- summary(emmeans(feat_test, trt.vs.ctrl~diet+timepoint, adjust="none"), infer=c(TRUE,TRUE))$contrasts
          dd_ctrl <- dd_ctrl[ !apply(dd_ctrl, 1, function(aaa){ any( is.na(aaa)) }) &
                                 !apply(dd_ctrl, 1, function(aaa){ grepl( "D16", aaa[ "contrast"] ) }) , ]

                    
        ## combine the first and second step EMM objects with the correct names, text yr Dad, tell him he's a good fella.           
          rbind(
            cc,
            cbind( "feat" = mgtax[ aa, 6],
                         "id" = aa,
                         "timepoint" = unlist( gsub(".*(D..).*", "\\1", dd_ctrl$contrast, perl = TRUE)),
                         "diet" = unlist( gsub("^(.*F)\\s.*", "\\1", dd_ctrl$contrast, perl = TRUE)),
                         as.data.frame(dd_ctrl)[ , -c(5,6)],
                         "shapiro_w" = shap.test$statistic,
                         "shapiro_p" = shap.test$p.value,
                         "levene_F" = lev.test[1,2],
                         "levene_p" = lev.test[1,3])
          )


          
        } ))
        # }, mc.cores  = 7 ))
      
        head(test_df)
        # View(test_df)

      
      ## unitque instances of comparison
        uniq_comp <- expand.grid( unique(test_df$contrast), unique(test_df$timepoint))
        test_df_fdr <-  do.call( "rbind", 
                                  apply( uniq_comp, 1, function(aa){  # aa <- uniq_comp[ 10 ,]
                  bb_df <- dplyr::filter( test_df, timepoint == unlist(aa["Var2"]), contrast == unlist(aa["Var1"]) )
                  bb_df$FDR <- p.adjust( bb_df$p.value, method = "fdr")
                  bb_df
                }) )
        
        
        saveRDS(test_df_fdr, "output/ms__van__DA_EMM-12__kA-0.005-0.1.RDS")
        # test_df <- readRDS("output/ms__van__DA_EMM-12__kA-0.005-0.1.RDS")
        # write.table(test_df, file = "output/ms__mouseVan_DA_LMM-EMM-12_filtered_DiffAb", sep = "\t")
        
  
        # save.image("output/ms__van__DA-EMM-12_test_output.RData")
        # # load("output/ms__van__DA-EMM-12_test_output.RData")
        

