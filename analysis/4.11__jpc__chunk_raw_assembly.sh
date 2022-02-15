
## assembly for 16S , more specific

## where is the place I put all my stuff?

											# # # place & summon basemount
												# # basemount $BASE
											
											# # # check 

		


## welcome to teagasc! working from asp    ====================================================
    
    # setup / auth basespace bs & bs-cp :: https://developer.basespace.illumina.com/docs/content/documentation/cli/cli-overview
    # sudo bash -c "$(curl -L https://basemount.basespace.illumina.com/install)"

    ## from illumina protocol
          
      # F :: TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG[CCTACGGGNGGCWGCAG]       [V3-V4]
      # R :: GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG[GACTACHVGGGTATCTAATCC]  [V3-V4]

      # Rev-comp F, on R2 :: [CTGCWGCCNCCCGTAGG]CTGTCTCTTATACACATCTGACGCTGCCGACGA      [V3-V4]
      # Rev-comp R, on R1 :: [GGATTAGATACCCBDGTAGTC]CTGTCTCTTATACACATCTCCGAGCCCACGAGAC  [V3-V4]

      # F :: 5' 341F :: CC TAC GGG NGG CWG CAG   -   S-D-Bact-0341-b-S-17
      # R :: 3' 805R :: GA CTA CHV GGG TAT CTA ATC C  -  	S-D-Bact-0785-a-A-21

    ## expected product:
      # 464 (805 - 341)
    
    # primers removed? nope.
      # removed ILLUMINACLIP:/home/jfg/bin/trimmomatic/adapters/NexteraPE-PE.fa:2:30:10 
      # removed the follwoing as per reviewer's recommendations :: 


      # conda create -n mqc_env multiqc
    ## cutadapt to deal with NAstorm
      # conda create -n cutad_env cutadapt
      # conda activate cutad_env
    ## figaro
      # git clone https://github.com/Zymo-Research/figaro.git
      # cd figaro
      # python3 setup.py bdist_wheel
      # pip3 install --force-reinstall dist/*.whl


## 0 ::   init   ## ===============================================================================

		# variables
		DB=~/data/ref/
		# JPDAT=/data/Food/primary/R0602_microsupport
		# MSOUT=/data/Food/analysis/R0602_microsupport/jamie.fitzgerald

    ## first tranche only
		JPDAT=/home/jfg/data/raw/4.11__jpc_1
		RAW=$JPDAT/4.11__jpennycook__raw	

	## PROJ is our generic term
		PROJ=~/data/4.11__jpc

		QC=$PROJ/1__qc
		FIGOUT=$QC/figout
		CUTAD=$PROJ/2__cutad
		FILT=$PROJ/2__filt_d2
		D2=$PROJ/4__dada2
		OUT=$PROJ/5__outputs
		MAT=$PROJ/Materials
		
		TEST=AMX1B
		
		## pipe home? 
		# ..........
		
		## there are many instances below where modules are loaded, either using module, or conda. 
		#  these have beel left in line, but could be moved here
		## similarly, many folders are made for QC - these could also be moved to the below mkdir step, but are currently in-line

		
		if [ -d $PROJ/Materials ] ; then  mv $PROJ/Materials $PROJ/Materials_$(date "+%b_%d_%Y") ;  fi
		# in order
		mkdir $RAW $PROJ $QC $FIGOUT $CUTAD $FILT $D2 $OUT $MAT 
    

##  1 ::   get data   ## ===================================================================    

    ## not baseclear! 
    
    ## gdrive
      # https://drive.google.com/drive/folders/1ol1PK139iKQ_ftC0yAfGv2tDf_7dR1fW?usp=sharing  #  tranche 1
      # https://drive.google.com/drive/folders/17RxjICDGut6GBx8I3pMNUMdRYNSD2Exu?usp=sharing  #  tranche 2


##  2 ::   QC raw data   ## ====================================================================

		module load fastqc multiqc
		mkdir $QC/qc_raw $QC/qc_raw_multi    
		parallel -j 11 fastqc {} -o $QC/qc_raw ::: $RAW/*fastq.gz

		# conda activate mqc_env    # matplotlibprobs
		multiqc $QC/qc_raw -o $QC/qc_raw_multi
		# cp -r $QC/qc_multi_raw ~/Dropbox/MicroSupp/ms__som/output/
		# scp -i ~/.ssh/id_rsa_ucc -r $QC/qc_multi_raw jamie@dunsh.ucc.ie:/home/jamie/   # scp -P 55555 -r ~/* jfg@143.239.154.14:/home/jfg/Dropbox/teag_sync/

		# conda deactivate
  
    
##  3 ::   CUTADAPT   ## ======================================================================

	## SP actually preferes TrimGalore - encodes illumina sequences
					
		## Cutadapt removes instances where the adapter (not primer) is in our sequence, as well as other tasks
		#  most common case is read-through at the 3' end - we spcify this (-a) and it's RC (-A) to cut it out.
		#  that is why there is no -g/-G below by default - we do not expect any surplus at the 5' end.
		#  note: sequences removed can be much longer than the adapter on its own
		## Primer removal specified and carried out in DADA2
		
		## note also that for FIGARO to work, we need uniform read legntsh per F and R reads (not necess _b/w_ F and R) 
		## hence -l and -m, based on length distrib from raw multiQC
		
		pigz -d $RAW/*fastq.gz
    rename 's/_001//' $RAW/*fastq  # -m

	# cutadapt
		# module load cutadapt
		conda activate cutad_env
		ls $RAW/*R1.fastq | sed -E "s/.*raw\/\/(.*)_R.*/\1/" | parallel -j 22 cutadapt -m 245 -l 245 -e 0.2 --discard-trimmed -a CTGTCTCTTATACACATCTCCGAGCCCACGAGAC -A CTGTCTCTTATACACATCTGACGCTGCCGACGA -o $CUTAD/{}_cutad_R1.fastq -p $CUTAD/{}_cutad_R2.fastq $RAW/{}_R1.fastq $RAW/{}_R2.fastq  >> $QC/qc_trim_param_testing.txt
		# overhead
		# conda deactivate
		pigz $CUTAD/*fastq
		# rm $RAW/*fastq    #  maybe don't

	# check
		mkdir $QC/qc_cutad
		parallel -j 22 fastqc {} -o $QC/qc_cutad ::: $CUTAD/*gz
		mkdir $QC/qc_cutad_multiqc
		# conda activate mqc_env
		cd $QC/qc_cutad_multiqc
		multiqc $QC/qc_cutad $QC/qc_cutad_multiqc
		# conda deactivate
		# cp -r $QC/qc_cutad_multiqc ~/Dropbox/MicroSupp/ms__som/output/
		scp -i ~/.ssh/id_rsa_ucc -r $QC/qc_cutad_multiqc jamie@dunsh.ucc.ie:/home/jamie/   # scp -P 55555 -r ~/* jfg@143.239.154.14:/home/jfg/Dropbox/teag_sync/

    
    ## gawjuss

    
##  4 ::   FIGARO   ## ======================================================================

	## infrastructure
    conda create -n figar_env
    conda activate figar_env
    conda install -c conda-forge r-base=3.6.1
    which R
    R            # use r cmd ? 
    install.packages("devtools")
    library("devtools")
    install_github("benjjneb/dada2")
    conda install -c anaconda python
    pip3 install numpy
    pip3 install scipy
    pip3 install matplotlib
    # get figaro reading on data, then gather for filterAndTrim via script
    cd ~/bin
    git clone https://github.com/Zymo-Research/figaro.git
    cd figaro
    pip3 install -r requirements.txt

		mkdir $QC/figout      
		# figaro needs paired reads in their own folder
		ls $CUTAD/*R1.fastq.gz | sed -E "s/.*cutad\/(.*)_cutad.*/\1/" | parallel -j 11  "mkdir $CUTAD/4.11__jpc_{} ; mv $CUTAD/{}*fastq.gz  $CUTAD/4.11__jpc_{}/ ; python3 ~/bin/figaro/figaro/figaro.py -i $CUTAD/4.11__jpc_{} -o $QC/figout/{} -a 426 -f 18 -r 22 "
		mv $CUTAD/4.11__jpc*/*gz $CUTAD
		rmdir $CUTAD/4.11__jpc*
		
	## copy 4.11__jpc__subchunk__filtAndTrim-w-Figaro to $MAT/chunk_figaro.R
		# basic and stupid
		

##  5 ::   DADA2   ## ======================================================================

	## now, implement filtAndTrim in parallel w. figaro's output. Note penalty parameter (6&7, set to 0,0 at first), but really should be calling FIGARO from R..
		# on hpc, try which or env below

    # hope you;ve put your figaro chunk in $MAT....
		conda activate figar_env
		ls $CUTAD/*R1.fastq.gz | sed -E "s/.*cutad\/(.*)_cutad.*/\1/" | parallel -j 8 "/usr/bin/Rscript ~/Dropbox/4.11/4.11__jpc/analysis/4.11__jpc__subchunk__filtAndTrim-w-Figaro.R {} $CUTAD $PROJ 18 21 0 0"


    # module load fastqc multiqc
		mkdir $QC/qc_fig $QC/qc_fig_multi    
		parallel -j 11 fastqc {} -o $QC/qc_fig ::: $FILT/*fastq.gz
		# conda activate mqc_env    # matplotlibprobs
		multiqc $QC/qc_fig -o $QC/qc_fig_multi


	## then implement DADA2 flow in chunk__dada2_decipher.R   >> >  >   > 
		# ...
		# ...
		# ...
		## gotta be a better step off method than this...


##  6 ::   check DADA data   ## ===================================================================    

	## check the dada2 output

		mkdir $QC/qc_trimmed $QC/qc_trimmed_multi    
		parallel -j 11 fastqc {} -o $QC/qc_trimmed ::: $FILT/*fastq.gz
		conda activate mqc_env
		multiqc $QC/qc_trimmed -o $QC/qc_trimmed_multi
		conda deactivate
		cp -r $QC/qc_trimmed_multi ~/Dropbox/MicroSupp/4.11__jpc/output/qc_trimmed_multi


##  7 ::   DADA stuffs   ## ==================================================================

    ## see the r script


##  X ::   trim data   ## ==================================================================

      # primer sequences are at the 5'&3' ends of the 464bp product. 
      # when split to paired reads, the F primer is at the __start/5'__ of the F read,
      # and R primer is at the __start/5'__ of the R read
      # so you must trim the primers from the start in either case
      
      ## could also ask  - do we need to do QC here instead of within D2? 
            ### YES ALWAYS
        # need to run to exclude adapter sequences etc - after figaro?
            ### that functionality is also present in Trimmomatic under ILLUMINACLIP
        # and if doing so, use cutadapt?
            ### SP prefers CA, that's about it.
				### actually no he prefers TrimGalore
      
      ## excluded from current run
            # 
            #   cd $RAW
            #   for i in $(ls *1.fastq.gz | sed -e "s/\(.*\)_R.*/\1/g") ;
            #     do java -jar /home/jfg/bin/trimmomatic/trimmomatic.jar PE -threads 11 \
            #     ${i}_R1_001.fastq.gz ${i}_R2_001.fastq.gz \
            #     $TRIM/trimmed_$i.R1.fastq.gz $TRIM/2.1_lost/lost_$i.R1.fastq.gz \
            #     $TRIM/trimmed_$i.R2.fastq.gz $TRIM/2.1_lost/lost_$i.R2.fastq.gz \
            #     HEADCROP:20 \
            #     MINLEN:150 ;
            #   done > $MAT/4.11__jpc__2_trim-raw.stout
            # 
            #   ## remove if not necess, can be better scrutinised in DADA2. Loss from MINLEN is ~25%, SLIDWIN = ~5%; LEAD/TRAIL=<5%
            #     # ILLUMINACLIP:/home/jfg/bin/trimmomatic/adapters/TruSeq3-PE.fa:2:30:10 \
            #     # SLIDINGWINDOW:4:15 \
            #     # CROP:280 \
            #     # LEADING:24 \
            #     # TRAILING:24 \
            # 
            # ## evaluate trimm
            # 
            #   mkdir $QC/trimmed_qc
            #   parallel -j 11 fastqc {} -o $QC/trimmed_qc ::: $TRIM/*fastq.gz
            #   multiqc $QC/trimmed_qc -o $QC/trimmed_qc_multi
            #   cp -r $QC/trimmed_qc_multi ~/Dropbox/MicroSupp/4.11__jpc/output/trimmed_qc_multi

    ## alt :: use cutadapt
            # python3 -m pip install --user --upgrade cutadapt
            # cutadapt -a 

    
