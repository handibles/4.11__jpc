


# ## F U N C T I O N A B L E  ======
# 
#     # melt FEAT+TAXA ; take VARIABLE: repeat each value ntaxa times to append to df
#       #
#     # note 2x fixing sort-by-order  w. jdh$order
#       jjj <- melt(cbind(j_crush, j[, jdh$order]), variable_name = "Sample")    
#       jjj <- cbind(jjj,  
#                  unlist( lapply( metadata[jdh$order, "sample"], function(m) 
#                      rep(m, nrow(taxo) )
#                          )
#                        )
#                      )
#       # head(jjj)
#       colnames(jjj)[ncol(jjj)] <- "Sample"

# ## ================================    


## ================================================================================================================================ ##
 # -------------------------------------------------------------------------------------------------------------------------------- #



## parallel compositional information UNIFRAC  ;  see abdiv package but not parallel

infofrac_mc <- function(abfeat, abtree, mc=6){ 
  abmat<- as.dist(as.matrix(do.call("rbind", 
                   mclapply( 1:ncol(abfeat), function(aa){  # aa <- 10
                     print("those samples had best be columns this time!")  # no print :(
                     aa_v <- abfeat[,aa]
                     apply(abfeat, 2, function(aaa){  # aaa <- bb_m[,3]
                       abdiv::information_unifrac(aa_v, aaa, tree = abtree)
                     }) }, mc.cores = mc) ), nrow = ncol(abfeat)))
return(abmat)
}





## fns of others :: tpq's CLR transform via sweep()

tpq.clr <- function(x) {
  #  samples are ROWS
  if( any( x == 0) ){ print(" use zCompositions::cmultRepl() to remove 0's") }
  sweep(log(x), 1, rowMeans(log(x)), "-")
        }


# ------------------------------------------------------------------------


## bin numericla vector


bin_vec <- function(
  vec = data_vector,
  bwidth = 10){
  
  print("assuming all values are posiitve - regut if not")
  vec <- as.vector(vec)
  bwidth <- bwidth
  bins <- seq(floor(min(vec)), ceiling(max(vec)), bwidth)
  sapply(vec, function(aa) {   # aa <- vec[300]
    bins[which( aa < bins)[1]]
  } )
}



# simple CLR! isn't gmean, which seems to be the norm but disagrees with the TEACHINGS

clr <- function(aa){    print(paste("samples must be ROWS ; dim =", dim(aa)[1], "x", dim(aa)[2])) 
                           t(apply(aa, 1, function(bb) log(bb) - log(mean(bb)) )) 
                        }


setlist <- function(set="all", excl_fns=TRUE){

  ## Needs: 
    # - way to take random text input, i.e. not assigned to set
    # - better way to decide whether to use knitr or base..
    # could add switch to choose how sorted
  
    
  # stuff to check
  if(all(set=="all")){
    stuff <- ls(envir = .GlobalEnv)
  }else{
    stuff <- unlist(lapply(set, function(aa){   # aa <- set[2]
      ls(pattern=aa, envir = .GlobalEnv)
    }) )
  }
  stuff <- unlist(stuff)
  
  # check the stuff
  setl <-  do.call("rbind",   
                   lapply(stuff, function(aa){          # aa <-stuff[2]
                     bb <- mget(aa, envir = .GlobalEnv )
                     if(is.null( dim(bb[[1]])) ){
                       cc <- c(aa, class(bb[[1]]), length(bb[[1]]), "0" )
                     }else{
                       cc <- c(aa, class(bb[[1]]), dim(bb[[1]])[1], dim(bb[[1]])[2])
                       cc
                     } })
  )
  
  # format, print
    # sort by added dims
    setl <- setl[ order(as.numeric(setl[,3]) + as.numeric(setl[,4]), decreasing=TRUE)   , ]
    colnames(setl) <- c("obj", "class", "nrows", "ncols")

  # kick out funcitons if you dont want em
    if(excl_fns){ setl <- setl[ !(grepl("function", setl[,2])) , ] }
        
    if("knitr" %in% installed.packages()){
      knitr::kable(setl)
    }else{
      print(setl)
    }
    # ## somehow, but who knows      
    # tryCatch(, finally = print(setl))
}


# --------------------------------------------------------------------------------------------------------------------------------

ill_gains <- function(bp_output= 15000000000,  # 15Gb, miseq default
                      paired=TRUE,             # PE def
                      readlength=250,          # 150bp def
                      nsamples=36,
                      ncontrols=0,
                      phiX=0.15){               # amplicon
  
  ## calculator for the illumina chemistries  
  
  options(scipen=99999999)
  # nextseq <- 120000000000
  # miseq <- 15000000000
  # nanoseq <- 1000000000
  # 
  # phiX <- 0.15
  # nsamples <- 36
  # ncontrols <- 4
  # readlength <- 300
  
  if(paired){read_n <- 2}else{read_n <- 1}
  SEQ <- bp_output

  outbp <- floor(SEQ - (SEQ*phiX))
  bp_persamp <- floor(outbp / (nsamples))
  outread <- floor((outbp/(readlength*read_n)))    # ie N paired if paired; N single if not
  read_persamp <- floor(outread / (nsamples))
  persamp_ctrl <- floor(outread / (nsamples+ncontrols))
  cat(paste("For a platform of", format(bp_output, big.mark=',', trim=TRUE), "basepairs max output:\ntotal bases output (less PhiX of",
            phiX,
            "%):\n",
            format(outbp, big.mark=',', trim=TRUE),
            "\nbase pairs per sample (n =",
            format(nsamples, big.mark=',', trim=TRUE),
            "):\n",
            format(bp_persamp, big.mark=',', trim=TRUE),
            "\ntotal reads out:\n ",
            format(outread, big.mark=',', trim=TRUE),
            "\nreads per sample (n =",
            format(nsamples, big.mark=',', trim=TRUE),
            "):\n ",
            format(read_persamp, big.mark=',', trim=TRUE),
            "\nreads per sample (n =",
            format(nsamples, big.mark=',', trim=TRUE),
            "; nctrl =",
            format(ncontrols, big.mark=',', trim=TRUE),
            "):\n",
            format(persamp_ctrl, big.mark=',', trim=TRUE)) )
  
}


# ------------------------------------------------------------------------------------------


## k over A

k_A <- function(
  counts,     # count table, taxa are ROWS
  k=0,        # raw counts
  A=0         # prevalence
){
  print("k_A: assuming samples are COLUMNS, and na.rm = TRUE")
  # parachute
  data <- counts     
  # to % as necessary
  if(k < 1 & any(colSums(counts)>1)){ counts <- apply(counts, 2, function(a) a/sum(a))}   
  
  if(A < 1){
    filt <- apply(counts, 1, function(x) sum(x>k, na.rm = TRUE)/ncol(counts) > A )    # x <- counts[20,]
  }else{
    filt <- apply(counts, 1, function(x) sum(x, na.rm = TRUE) > k )
    # filt <- apply(j, 1, function(x) sum(x>k) > 1 )    # wtf was this?? j?? >1 presumably a clumsy any() ???
  }
  return(data[ filt, ])
}


## crush
# if RANK ABUND < k_A :  crush RANK TAXON

crush_ranks2 <- function(FEATS=feats,
                         TAX=tax,
                         RANK=colnames(TAX)[3],
                         k=0.001,
                         A=0.1,
                         MONIKER="etc.",
                         NA_TERM="unknown",
                         P=NULL){

  # print("note: crush ranks won't handle 1-column DFs")
  
  print("note: assuming samples are COLS")
  dim(FEATS)
  
  if(!is.null(P) & class(P)=="phyloseq"){FEATS <- otu_table(P) ; TAX <- tax_table(P)}       
  
  # cull NAs - unnecessary, hackish, but brain is cold
  TAX <- as.matrix(TAX)
  TAX[is.na(TAX)] <- NA_TERM
  
  FEATS <- data.frame(FEATS, stringsAsFactors = FALSE)
  TAX <- data.frame(TAX, stringsAsFactors = FALSE)
  
  
  # SUBSET o at RANK (keep feat names)
  unq <- unique(TAX[, RANK])
  
  # collect subset
  T_crushed <- (lapply(unq, function(u){     # u <- unq[8]
    FEATS_u <- FEATS[ TAX[,RANK]==u , ]
    TAX_u <- TAX[ rownames(FEATS_u) , ]
    if( (sum(colSums(FEATS_u) > k) / ncol(FEATS_u)) < A){ TAX_u[ rownames(FEATS_u), RANK ] <- c(MONIKER) }
    TAX_u
  } )) 
  T_crushed <- do.call(rbind.data.frame, T_crushed)
  T_crushed <- T_crushed[rownames(TAX) , ]   # return same order
  print( paste0("only ", length(unique(T_crushed[ , RANK])), " ranks in crushed ",RANK ))
  return(T_crushed)
}


# -----------------

## Bray Curtis of one sample against a vector of samples
oneWayBC <- function(
  i = vector,
  df = datafrm,
  exclude=TRUE){
  
  ## phew... define BC for one sample versus a vector
  # as per vegan expresses as (A+B-2*J)/(A+B) but == .
  # set as RA prior to here
  
  ## returns vector of length == ncol(df)
  
  owBC <-    apply(df, 2, function(j){
                ij <- rbind(i,matrix(j,nrow = 1))
                Cij <- sum(apply(ij, 2, function(pairs){ if(all(pairs > 0)){min(pairs)}else{ 0 }   }) )    # min abund if feat in both samples, else 0
                Si <- sum(i)
                Sj <- sum(j)
                1-((2*Cij)/(Si+Sj))     # bray-curtis dissim
              })
  
  # in case compared with self (i.e. BC=0), exclude
  if(exclude & any(owBC == 0)){ 
    owBC <- owBC[ owBC > 0]
    print(paste("excluded identical sample at column", as.character((1:ncol(df))[owBC==0]), "as presumably is same as input i"))
    }
  return(owBC)
}


### =========================================================================================== ###
## ===================== ##   R A N K    O P E R A T I O N S   ## ===============================##
## ---------------------                                          -------------------------------##

rank_ranks <- function(
  p = NULL,
  o = abund,
  t = categ, 
  rank_above = colnames(t)[2],
  mean_thresh = NULL,
  tot_thresh = NULL,
  # trim_last = TRUE,
  digits = 3,
  mc.cores = 6,
  summarise = c('basic', 'percent', 'raw')
  ){
  
  #### ----------------------------
      #   o <- abu[1:500, 1:20]
      #   t <- abu_tax[1:500, ]
      #   rank_above <- colnames(t)[2]
      #   mean_thresh <- 0
      #   tot_thresh = NULL
      #   p <- NULL
      #   digits <- 3
      #   summarise <-  c('basic', 'percent', 'raw')
      #   options(scipen = 999)
      #   mc.cores <- 6
  #### ----------------------------
  
  if(!is.null(p)){
    require('phyloseq')
    o <- data.frame(otu_table(p), stringsAsFactors = FALSE)
    t <- data.frame(tax_table(p), stringsAsFactors = FALSE)
    
  }else{
    
    o <- data.frame(o, stringsAsFactors = FALSE)      ## still making assumptions about orientation, TAXA_ARE_ROWS
    t <- data.frame(t, stringsAsFactors = FALSE)
  }
  
  t[is.na(t)]  <- 'unknown'   # remove NAs, which are problematic at lower ranks
  
  # remove the final column as huge job and uninformative?
  
  # don't feed relative abundances (percentages) into this: distorts summary
  if (any(colSums(o) == 1)){ 
    
    return(print('input is in relative abundances (%) - need raw (read) abundances. (or table not orientated properly and sum(ASV)=1 - try transposing table with "t( )" )'))
    
  } else{
    

    unq <- as.vector( unique(t[ , rank_above]) )
    
    # original but slow...
    ## make & populate basic table - thank you sapply!
    # m <- t(sapply( unq, function(u){      #  u <- unq[10]
    #   z <- o[ (t[,rank_above] == u) , ]
    #   m1 <- sum(z)                               # u total as reads
    #   m2 <- sum(z)/ncol(o)                       # avg. reads / sample
    #   m3 <- sum(z)/sum(o)                        # u total as RA
    #   m4 <- mean(apply(z,2,sum)/colSums(o))      # avg. rank RA / sample
    #   m5 <- mean(as.matrix(z/colSums(o)))        # mean(RA)
    #   m6 <- sd(as.matrix(z/colSums(o)))          # sdev(RA)
    #   m7 <- mean(apply(z,1,sum)/sum(o))          # mean(SumASV/SumTot%)
    #   m8 <- sd(apply(z,1,sum)/sum(o))            # sdev(SumASV/SumTot%)
    #   c( m1, m2, m3, m4, m5, m6, m7, m8 )
    # })
    # )
    
    ## gone par:   -  seems fine lol
    m_par <- (parallel::mclapply( unq, function(u){      #  u <- unq[10]
      z <- o[ (t[,rank_above] == u) , ]
      m1 <- sum(z)                               # u total as reads
      m2 <- sum(z)/ncol(o)                       # avg. reads / sample
      m3 <- sum(z)/sum(o)                        # u total as RA
      m4 <- mean(apply(z,2,sum)/colSums(o))      # avg. rank RA / sample
      m5 <- mean(as.matrix(z/colSums(o)))        # mean(RA)
      m6 <- sd(as.matrix(z/colSums(o)))          # sdev(RA)
      m7 <- mean(apply(z,1,sum)/sum(o))          # mean(SumASV/SumTot%)
      m8 <- sd(apply(z,1,sum)/sum(o))            # sdev(SumASV/SumTot%)
      c( m1, m2, m3, m4, m5, m6, m7, m8 )
    }, mc.cores = mc.cores)
    )
    mm_par <- do.call("rbind", m_par)
    rownames(mm_par) <- unq   # head(mm)
    
    m <- mm_par  
    ##
    
    
    colnames(m) <- c('rank_sum' ,'rank_mean' ,'rank_TotRA' ,'rank_meanSampRA' ,'ASV_meanSampRA' ,'ASV_sdSampRA' ,'ASV_meanTotRA' ,'ASV_sdTotRA' )
    
    
    ## per-sample columns :: summarise raw, relab 
    n <- NULL
    n_ra <- NULL
    
    if('raw' %in% summarise){  
      n <- t(sapply(unq, function(u){        # use of SAPPLY avoids setting matrix, cleaner
        colSums(o[ t[,rank_above]==u , ])    # works so well as ONLY RETURNING ONE OPERAITON
      }))
      # n <- cbind( 'RAW' = rep('|', nrow(n)) , n )     #non-numeric argument to mathematical function
    }
    
    if('percent' %in% summarise){  
      n_ra <- t(sapply(unq, function(u){
        z <- o[ t[,rank_above]==u , ]
        (colSums(z)/colSums(o))*100
      })) 
      # n <- cbind( 'PERCENT' = rep('|', nrow(n)) , n )    # non-numeric argument to mathematical function
    }
    
    # bring together, sort
    m <- cbind(m, n, n_ra)  
    m <- m[order(m[ , "rank_sum"], m[ , "rank_TotRA"], decreasing = TRUE) , ] 
    
    ## threhsold
    if( !is.null(mean_thresh) ){
      print(paste0('discarding the levels: ' , paste( rownames(m)[ m[,"rank_meanSampRA"]<mean_thresh ] , collapse = ',') ) )
      # if(mean_thresh < 1){mean_thresh <- mean_thresh*100}
      m <- m[ m[,"rank_meanSampRA"] > mean_thresh , ] }
    
    if( !is.null(tot_thresh) ){
      print(paste0('discarding the levels: ' , paste( rownames(m)[ m[,"rank_TotRA"]<tot_thresh ] , collapse = ',') ) )
      # if(tot_thresh < 1){tot_thresh <- tot_thresh*100}
      m <- m[ m[,"rank_TotRA"] > tot_thresh , ] }
    
  print("__still__ not using aggregate() :( ")
  
    m <- round(m, digits = digits)
    return(m)
  } 
}


# --------------------------------------------------------------------------------------------------------------------------------

## glom taxa - should match the dims of the below

# march 2020
glom_tax <- function(
  TAX=tax,
  RANK=colnames(tax)[ncol(tax) - 1],
  NA_TERM="unknown"){
  
  # simplification: don't use NAs as rownames you freak
  TAX[ is.na(TAX) ] <- NA_TERM
  
  uniq_glom <- unique(TAX[ , RANK])
  tax_glom <- (unique(lapply(uniq_glom, function(aa) unlist(TAX[ match(aa,TAX[,RANK]) , 1:match(RANK,colnames(TAX))]  )) ))   # d <- "Streptococcus"
  tax_glom <- do.call("rbind", tax_glom)
  rownames(tax_glom) <- uniq_glom
  
  return(tax_glom)

}


# --------------------------------------------------------------------------------

# glom_counts  -  just flatten at a given rank
# glom_counts <- function(
#                         O=feats,
#                         T=tax,
#                         RANK=colnames(T)[5]){      
#     unq <- unique(T[, RANK])
#     O_glommed  <- lapply(unq, function(u){  # u <- unq[21]
#         O_u <- O[ T[,RANK]==u , ]
#         T_u <- T[ rownames(O_u) , ]
#         unq <- unique(T[, RANK])
#         if(!(is.null(dim(O_u)))){
#           O_g <- t(matrix(colSums(O_u)))
#         }else{
#           O_g <- as.vector(sum(O_u))
#           }
#         colnames(O_g) <- colnames(O_u)
#         if(!(is.null(dim(O_g)))){
#           rownames(O_g) <- u
#           # O_g <- t(matrix(colSums(O_u)))
#         }else{
#           names(O_g) <- u
#         }
#         
#         O_g
#         })
#         O_glommed <- do.call(rbind.data.frame, O_glommed)
#     return(O_glommed)
#     } 

## how many times have we written this function?...

    ## march 2020
glom_counts <- function(
    abundances = asv,
    identities = tax,
    glomrank = "genus",
    NA_TERM = "unknown"){
    
    # abundances = asv
    # identities = tax
    # glomrank = "genus"
    # NA_TERM = "unknown"
    
    if(!( glomrank %in% colnames(identities))){
      print("glomrank not in the hierarchy supplied")
      stop()
    }  
    
    identities[is.na(identities)] <- NA_TERM    # kill NA
    
    uniq_glom <- unique(identities[ , glomrank])       # uniques
    print( "psst: assuming samples are ROWS" )
    
    asv_glom <- sapply(uniq_glom, function(a){                   # a <- NA            a <- uniq_glom[10]
      if(is.na(a)){
        rowSums(abundances[ , is.na(identities[,glomrank])])    # catch pesky NAs - hopefully redundant now, but keep
      }else{
        b <- abundances[ , which(identities[,glomrank] == a) ]  # catch ASVs
        if(is.null(dim(b))){
          as.vector(b)                                          # catch n=1 taxa
        }else{
          c <- as.vector(rowSums( b ))
        }
      }
    })
    rownames(asv_glom) <- rownames(abundances)
    
    return(asv_glom)
}


# --------------------------------------------------------------------------------------------------------------------------------

## colour/shade

shade_ranks2 <- function(
  TAXA = NULL,    
  COLOUR_BY =  colnames(taxo)[4],
  SHADE_BY =  colnames(taxo)[6],
  MONIKER = "etc.",
  NA_TERM = "n/a",
  MON_COL = "grey50",
  NA_COL = "grey75",
  PADDING = 2,
  # NAMED_COLS = NULL,                #pre-made named colour vector
  SHOW=TRUE,
  colour_and_shade = c("both", "colour", "shade"),     # this bit still whack
  P=NULL,  #phylo
  #
  lum = 65,         # 45
  chrom = 100,      # 175
  h.start = 0,
  h = c(0, 360) + 15 ){

  # ====
  if(is.null(TAXA) & is.null(P)){cat("No input. feed me!") ; return(NULL)}
  
  if(!is.null(P) & class(P)=="phyloseq"){TAXA <- tax_table(P)}
  require(scales)
  
  print("looks like shade_ranks2 doesnt take phylo objects  :( ")
  
        # if( !is.null(NAMED_COLS)){    ## pass a pre-made named colour vector
        #   tax_cols <- NAMED_COLS
        #   unq <- names(tax_cols)
        # }else{
  
  unq <- unique(TAXA[ , COLOUR_BY])
  unq <- unq[ !(unq %in% c(MONIKER, NA_TERM, NA)) ] 
  unq <- c(unq,  MONIKER)    # add back in last to catch taxa orphaned by CRUSHing
  tax_cols <- hue_pal( l = lum, c = chrom, h=h, h.start = h.start)(length(unq))
  names(tax_cols) <- unq

  shades <- lapply(unq, function(u){
    unq.shades <- unique(TAXA[ (TAXA[, COLOUR_BY] == u) , SHADE_BY])
    unq.shades <- unq.shades[ !(unq.shades %in% c(MONIKER, NA_TERM, NA)) ] 
    shade.ramp <- colorRampPalette(c(tax_cols[u == names(tax_cols)], "white"))(length(unq.shades) + PADDING)[1:length(unq.shades)]
    names(shade.ramp) <- unq.shades
    shade.ramp <- unlist(shade.ramp)
  })
  shades <- unlist(shades)

  # only add misc terms once
  tax_cols <- tax_cols[ !(names(tax_cols) %in% c(MONIKER, NA_TERM, NA)) ] 
  shades <- shades[ !(shades %in% c(MONIKER, NA_TERM, NA)) ] 
  MONIKER_NA <- c(MON_COL, NA_COL)
  names(MONIKER_NA) <- c(MONIKER, NA_TERM)
  
  print("colour and shade is still broken")
  if(colour_and_shade == "colour"){ 
    tax_cols <- c(tax_cols, MONIKER_NA) 
  }else if(colour_and_shade=="shade"){shades <- c(shades, MONIKER_NA)
  }else {tax_cols <- c(tax_cols, MONIKER_NA, shades) }
  
  return(tax_cols)    # show_col(tax_cols)
}


## glom ranks. ...


###################################################################################
# ____________________   G  G   P L O T   F N s   _________________________________

# guess what? you can functionalise ggplot calls, how unexpected. Key is aes_string() 


###   E O E   H A R D C O D E D 

  # assumes ranks already crushed and shaded, for melted phylos
  jfg_plot <- function(
    YRDATA = z.rm5,   #  a melted object
    RANK = "genus",
    TAX_COL = g_col,
    TAG= "A", 
    ...){             # working?
    
    # define for meltyface
    YRDATA[ , RANK] <- factor(YRDATA[ , RANK], levels= rev(c( names(TAX_COL[ !(names(TAX_COL) %in% c("various", "unknown"))]), "various", "unknown") ))
    
    ggplot( YRDATA, aes(x=ID_Sample,y=Abundance), aes_string(fill = RANK)) +
      facet_grid(Description~., scales='free', space='free', switch="y") +                                  # scale/size control
      geom_bar(aes_string(fill = RANK, color = NULL ), stat ="identity", position="stack", width=1, color='black', size=0.05) +      # , color='black'
      theme(axis.text.x = element_text(angle = 0, )) +  # rotate
      scale_y_continuous(breaks=seq(0,1, by=0.10)) +
      scale_fill_manual(values = (TAX_COL), na.value='white')  +
      theme_transparent() +   # put theme call first so doesn't overwrite other theme-calls
      theme(plot.margin=unit(c(0,0,0,0), "null"),
            # panel.grid.major.x = element_line(colour='grey45', size = 0.3),
            # panel.grid.major.y = element_line(colour='white', size = 0),
            # panel.border = element_rect(colour = "grey45", fill=NA, size=1),
            axis.line=element_blank(),
            axis.text.y =element_blank(),
            axis.ticks.y =element_blank(),
            strip.text.y = element_text(size = 10, angle = 180),
            plot.title = element_text(size=16),
            plot.tag = element_text(size=20),
            legend.position = 'right',
            legend.box = "vertical") +
      guides(fill = guide_legend(ncol = 1, reverse=T)) +
      labs(title=paste0(RANK," level"), tag=TAG) +
      coord_flip(ylim=c(0,1))
    
  }


  ## E O E   H A R D C O D E D    # =================================================

# assumes ranks already crushed and shaded
ejl_plot <- function(
  a = otu_table(eoe_crush),
  t = tax_table(eoe_crush),
  m = sample_data(eoe_crush),
  RANK = "genus",
  TAX_COL = e_col,
  TAG= "A"){
  
  # RANK <- "family"
  # TAX_COL <- e_col_f
  
  # get each Description's total
  a_sub <- sapply(unique(m$Description), function(d){    # d <- "PRE"
    desc_sub <- a[m$Description == d , ]
    sapply(unique(t[ , RANK]), function(g){  # g <- "Alloprevotella"
      (sum(desc_sub[ , t[ , RANK] == g])/nrow(desc_sub))*100   # conditional on matched row/colsnames :: sum(rownames(t) == colnames(a))
    })
  })
  colnames(a_sub) <- unique(m$Description)
  
  
  # melt m & a
  a_sub_m <- melt(a_sub)
  colnames(a_sub_m) <- c(RANK, "Description","Abundance")
  a_sub_m[ , RANK] <- factor(a_sub_m[ , RANK], levels= rev(c( names(TAX_COL[ !(names(TAX_COL) %in% c("various", "unknown"))]), "various", "unknown") ))
  
  # Control
  
  # recode variables for plotting
  a_sub_m[ , "Description"] <- plyr::mapvalues(a_sub_m[ , "Description"], from=c( "POST_FED", "POST_PPIs", "POST_STEROIDS", "PRE", "CONTROL"), to=c("EoE post-FED (n=9)", "EoE post-PPI (n=9)", "EoE post-STC (n=9)", "EoE baseline (n=28)","Control (n=10)"))
  a_sub_m[ , "Description"] <- factor(a_sub_m[ , "Description"], levels = (c("EoE post-STC (n=9)", "EoE post-PPI (n=9)", "EoE post-FED (n=9)", "EoE baseline (n=28)","Control (n=10)")))
  
  
  # plot each condition (i.e. 5 rows) by the colours defined by EJL
  # annotate with the appropriate stuff
  ggplot(a_sub_m, aes(x=Description, y=Abundance), aes_string(fill=RANK)) +
    geom_bar(aes_string(fill = RANK, color = NULL ), stat ="identity", position="stack", width=0.75 , color='black', size=0.15) +      #
    # theme_transparent() +
    theme(axis.text.x = element_text(angle = 0, )) +  # rotate   
    scale_y_continuous(breaks=seq(0,100, by=10)) +
    scale_fill_manual(values = (TAX_COL), na.value='pink')  +
    theme(panel.grid.major.x = element_line(colour='grey45', size = 0.3),
          panel.grid.major.y = element_line(colour='white', size = 0), 
          panel.border = element_rect(colour = "grey45", fill=NA, size=1),
          plot.margin=unit(c(0,0,0,0), "null"),
          axis.line=element_blank(),
          strip.text.y = element_text(size = 10, angle = 0),
          plot.title = element_text(size=16),
          plot.tag = element_text(size=20),
          legend.position = 'right', 
          legend.box = "vertical") +
    guides(fill = guide_legend(ncol = 1, reverse=T)) +              # change to 5
    labs(title=paste0(RANK," level"), tag=TAG) +
    xlab("           POST                PRE   normal") +   # could but dont
    ylab("Relative Abundance (%)") +
    coord_flip(ylim=c(1,99))    
}

