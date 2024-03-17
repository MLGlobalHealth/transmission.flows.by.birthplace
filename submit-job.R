require(data.table)

hmc_chains_n <- 4

## important note:
# the combination of stanModelFile and job_tag should be unique for each analysis. job_tag should be the same for MSM and HSX model
# all outputs with stanModelFile-job_tag are assumed to be several HMC chains run in parallel

#	function to make PBS header
make.PBS.header <- function(hpc.walltime=23, hpc.select=1, hpc.nproc=8, hpc.mem= "80gb", hpc.load= "module load anaconda3/personal\nsource activate src", hpc.q=NA, hpc.array=1 )
{
  pbshead <- "#!/bin/sh"
  tmp <- paste("#PBS -l walltime=", hpc.walltime, ":59:00", sep = "")
  pbshead <- paste(pbshead, tmp, sep = "\n")
  tmp <- paste("#PBS -l select=", hpc.select, ":ncpus=", hpc.nproc,":ompthreads=", hpc.nproc,":mem=", hpc.mem, sep = "")
  pbshead <- paste(pbshead, tmp, sep = "\n")
  pbshead <- paste(pbshead, "#PBS -j oe", sep = "\n")
  if(hpc.array>1)
  {
    pbshead	<- paste(pbshead, "\n#PBS -J 1-", hpc.array, sep='')
  }
  if(!is.na(hpc.q))
  {
    pbshead <- paste(pbshead, paste("#PBS -q", hpc.q), sep = "\n")
  }
  pbshead	<- paste(pbshead, hpc.load, sep = "\n")
  pbshead
}

# input args
if(1)
{
  hpc.nproc.cmdstan <- 5
  args <- data.table(
    cmdstan_dir = '/apps/cmdstan/2.33.0',
    source_dir = '~/git/transmission.flows.by.birthplace',
    in_dir='/rds/general/project/ratmann_roadmap_data_analysis/live',
    out_dir= '/rds/general/project/ratmann_roadmap_data_analysis/live/transmission_sources',
    script_make_pairs= 'scripts/formulate-Amsterdam-pairs.R',
    script_file= 'scripts/run-stan.R',
    script_converting_file = "scripts/stan-convert-csv-to-rda.r",
    #stanModelFile = 'mm_sigHierG_bgUnif_piVanilla_220408b', # vanilla model
    #stanModelFile = 'mm_sigHierG_bgUnif_piReg_230111b', # covariate model
    stanModelFile = 'mm_bgUnif_piGP_221027b_hpc', # 2D HSGP model
    #stanModelFile = 'mm_bgUnif_pi1DGP_Ams_230224', # 1D HSGP model
    #stanModelFile = 'mm_bgUnif_pi1DGP_Ams_230224b', # 2 * 1D HSGP model
    analysis= 'analysis_220713',
    pairs_dir = 'bs_tsi_MSM-2010_2022',
    job_tag = 'bs_tsi_MSM-2010_2022',
    UndiagStanModelFile='undiagnosed_211102',
    job_tag_undiagnosed = 'cohort_2010_2015',
    trsm = 'MSM',
    time_period='2010-2022',
    clock_model = '/rds/general/project/ratmann_roadmap_data_analysis/live/transmission_sources/molecular_clock/hierarchical',
    hmc_stepsize= 0.25,
    hmc_num_samples= 2000,
    hmc_num_warmup= 500,
    cmdstan = 1L,
    seed = 42,
    chain = 1,
    local = 0,
    m1 = 24,
    m2 = 24,
    B = 1.2,
    bs_tsi = 0, # for bootstrap resampling TSIs
    bs_phylo = 0 # for bootstrap resampled alignments
  )
}

if(1)
{
  tmp <- data.table(chain=1:hmc_chains_n)
  tmp[, seed:= c(256081,483137,136323,756305)]
  set(args, NULL, colnames(tmp), NULL)
  tmp[, dummy:= 1L]
  args[, dummy:= 1L]
  args <- merge(args, tmp, by='dummy')
  set(args, NULL, 'dummy', NULL)
}

#if(args$bs_tsi[1]==1){
  tsi_seed <- abs(round(rnorm(1) * 1e6)) # set same seed across chains for resampling pairs
#}

# make commands
cmds <- vector('list', nrow(args))
for(i in seq_len(nrow(args)))
{
  #	general housekeeping
  # first generate pairs
  cmd    <- ''
  #	general housekeeping
  cmd    <- paste0(cmd,"CWD=$(pwd)\n")
  cmd    <- paste0(cmd,"echo $CWD\n")
  tmpdir.prefix	<- paste0('src_',format(Sys.time(),"%y-%m-%d-%H-%M-%S"))
  tmpdir    <- paste0("$CWD/",tmpdir.prefix)
  cmd    	<- paste0(cmd,"mkdir -p ",tmpdir,'\n')
  #	generate data set and run if not using cmdstan
  cmd     <- paste0( cmd, 'echo "----------- Generating input data: ------------"\n')
  tmp     <- paste0('Rscript ', file.path(args$source_dir[i],args$script_make_pairs[i]),
                    ' -source_dir ', args$source_dir[i],'',
                    ' -indir ', args$in_dir[i],'',
                    ' -outdir ', tmpdir,'',
                    ' -analysis ', args$analysis[i],'',
                    ' -clock_model ', args$clock_model[i],'',
                    ' -jobtag "', args$job_tag[i],'"',
                    ' -trsm "', args$trsm[i],'"',
                    ' -seed ', tsi_seed,'',
                    ' -bs_tsi ', args$bs_tsi[i],' ',
                    ' -bs_phylo ', args$bs_phylo[i],' '
  )
  cmd    	<- paste0(cmd, tmp, '\n')

  # then run stan for each HMC chain
  #cmd    	<- paste0(cmd,"CWD=$(pwd)\n")
  #cmd    	<- paste0(cmd,"echo $CWD\n")
  #tmpdir.prefix	<- paste0('src_',format(Sys.time(),"%y-%m-%d-%H-%M-%S"))
  #tmpdir    <- paste0("$CWD/",tmpdir.prefix)
  #cmd    	<- paste0(cmd,"mkdir -p ",tmpdir,'\n')
  #	generate data set and run if not using cmdstan
  cmd     <- paste0( cmd, 'echo "----------- Generating input data: ------------"\n')
  tmp     <- paste0('Rscript ', file.path(args$source_dir[i],args$script_file[i]),
                   ' -source_dir "', args$source_dir[i],'"',
                   ' -stanModelFile "', args$stanModelFile[i],'"',
                   #' -seed ', args$seed[i],
                   ' -chain ', args$chain[i],
                   ' -indir ', args$in_dir[i],'',
                   ' -outdir ', tmpdir,'',
                   ' -jobtag "', args$job_tag[i],'"',
                   ' -trsm "', args$trsm[i],'"',
                   ' -cmdstan ', args$cmdstan[i],
                   ' -pairs_dir ', args$pairs_dir[i],
                   ' -clock_model ', args$clock_model[i],
                   ' -time_period ', args$time_period[i],
                   ' -m1 ', args$m1[i],
                   ' -m2 ', args$m2[i],
                   ' -B ', args$B[i],
                   ' -local ', args$local[i]
  )
  cmd    	<- paste0(cmd, tmp, '\n')
  #	if using cmdstan
  if(args$cmdstan[i]==1)
  {
    cmd <- paste0(cmd, 'echo "----------- Building Stan model file: ------------"\n')
    #	clean up any existing model code
    cmd <- paste0(cmd, 'rm ', file.path('$CWD',paste0(args$stanModelFile[i],'.*')), ' \n')
    #	copy stan model file
    cmd	<- paste0(cmd, 'cp -R ',file.path(args$source_dir[i], 'stan_model_files',paste0(args$stanModelFile[i],'.stan')),' .\n')
    #	build model
    cmd <- paste0(cmd, 'cd ', args$cmdstan_dir[i], '\n')
    cmd <- paste0(cmd, 'make STAN_THREADS=TRUE ', file.path('$CWD',args$stanModelFile[i]), ' \n')
    cmd <- paste0(cmd, 'cd $CWD\n')
    #	set up env variables
    #cmd <- paste0( cmd, 'JOB_DIR=$(ls -d "',tmpdir,'"/*/)\n') # lists all directories in tmpdir
    cmd <- paste0( cmd, 'JOB_DIR=$(ls -d "',tmpdir,'/mm"*/)\n') # only find directories starting with mm (prefix to stan model name) i.e. not pairs directory
    cmd <- paste0( cmd, 'JOB_DIR=${JOB_DIR%?}\n')
    cmd <- paste0( cmd, 'JOB_DIR_NAME=${JOB_DIR##*/}\n')
    cmd <- paste0( cmd, 'SCRIPT_DIR=',args$source_dir[i],'\n')
    cmd <- paste0( cmd, 'IN_DIR=',args$in_dir[i],'\n')
    cmd <- paste0( cmd, 'PAIRS_DIR=',args$pairs_dir[i],'\n')
    cmd <- paste0( cmd, 'CLOCK_DIR=',args$clock_model[i],'\n')
    cmd <- paste0( cmd, 'ANALYSIS=',args$analysis[i],'\n')
    cmd <- paste0( cmd, 'TRSM=',args$trsm[i],'\n')
    cmd <- paste0( cmd, 'TIME_PERIOD=',args$time_period[i],'\n')
    cmd <- paste0( cmd, 'STAN_MODEL_FILE=',args$stanModelFile[i],'\n')
    cmd <- paste0( cmd, 'STAN_DATA_FILE=$(find ', tmpdir, ' -name "*cmdstanin.R")\n')
    cmd <- paste0( cmd, 'STAN_INIT_FILE=$(find ', tmpdir, ' -name "*cmdstaninit.R")\n')
    cmd <- paste0( cmd, 'STAN_OUT_FILE=', file.path('$JOB_DIR','${JOB_DIR##*/}_stanout.csv'),' \n')
    #	run model
    cmd <- paste0( cmd, 'echo "----------- env variables are: ------------"\n')
    cmd <- paste0( cmd, 'echo $JOB_DIR\n')
    cmd <- paste0( cmd, 'echo $JOB_DIR_NAME\n')
    cmd <- paste0( cmd, 'echo $STAN_DATA_FILE\n')
    cmd <- paste0( cmd, 'echo $STAN_OUT_FILE\n')
    cmd <- paste0( cmd, 'echo $IN_DIR\n')
    cmd <- paste0( cmd, 'echo $ANALYSIS\n')
    cmd <- paste0( cmd, 'echo "----------- Starting Stan sampling: ------------"\n')
    #
    tmp <- paste0( './',args$stanModelFile[i],' ',
                   'sample num_samples=',args$hmc_num_samples[i],' num_warmup=',args$hmc_num_warmup[i],' save_warmup=0 thin=1 ',
                   'adapt delta=0.95 ',
                   'algorithm=hmc engine=nuts max_depth=15 stepsize=',args$hmc_stepsize[i],' ',
                   'data file=$STAN_DATA_FILE ',
                   'init=$STAN_INIT_FILE ',
                   'random seed=',args$seed[i],' ',
                   'output file=$STAN_OUT_FILE' )
    cmd <- paste0(cmd, tmp, '\n')
    # convert csv to rdata
    cmd		<- paste0( cmd, 'echo "----------- Converting Stan output to RDA file: ------------"\n')
    tmp		<- paste0('Rscript ', file.path(args$source_dir[i],args$script_converting_file[i]),
                   ' -csv_file "', "$STAN_OUT_FILE",'"',
                   ' -rda_file "', file.path('$JOB_DIR','${JOB_DIR##*/}_stanout.RData'),'"'
    )
    cmd		<- paste0(cmd, tmp, '\n')
  }

  #	general housekeeping
  cmd 	<- paste0( cmd, 'echo "----------- Copy files to out directory: ------------"\n')
  tmpdir2	<- file.path(args$out_dir[i], paste0(args$stanModelFile[i],'-',args$job_tag[i]))
  if(i==1)
  {
    dir.create(tmpdir2)
  }
  cmd		<- paste0(cmd,"mkdir -p ",tmpdir2,'\n')
  cmd		<- paste0(cmd, 'cp -R --no-preserve=mode,ownership "', tmpdir,'"/* ', tmpdir2,'\n')
  cmd		<- paste0(cmd, 'chmod -R g+rw ', tmpdir2,'\n')

  # create post-processing shell script for central analyses

  if(i==1)
  {
    cmd2 <- make.PBS.header( hpc.walltime=23,
                             hpc.select=1,
                             hpc.nproc=48,
                             hpc.mem= "124gb",
                             hpc.load= "module load anaconda3/personal\nsource activate src",
                             hpc.q=NA,
                             hpc.array= 1)
    cmd2 <- paste0(cmd2,'\n')
    # set up env variables
    cmd2 <- paste0(cmd2,'SCRIPT_DIR=',args$source_dir[i],'\n',
                   'IN_DIR=',args$in_dir[i],'\n',
                   'OUT_DIR=',tmpdir2,'\n',
                   'JOB_TAG=',args$job_tag[i],'\n',
                   'STAN_MODEL_FILE=',args$stanModelFile[i],'\n',
                   'CHAINS=', hmc_chains_n,'\n',
                   'UNDIAGNOSED=', paste0(args$out_dir[i],args$UndiagStanModelFile[i],args$job_tag_undiagnosed[i]),'\n',
                   'JOB_TAG_UNDIAGNOSED=', args$job_tag_undiagnosed[i],'\n',
                   'TRSM=', args$trsm[i], '\n',
                   'ANALYSIS=', args$analysis[i], '\n',
                   'LOCAL=', args$local[i], '\n',
                   'OVERWRITE=0\n'
    )
    # save posterior samples
    tmp <- paste0('Rscript ', file.path('$SCRIPT_DIR','scripts','post-processing-save-posterior-samples.R'),
                  ' -source_dir $SCRIPT_DIR -stanModelFile $STAN_MODEL_FILE -outdir $OUT_DIR -job_tag $JOB_TAG -numb_chains $CHAINS -trsm $TRSM')
    cmd2 <- paste0(cmd2,tmp,'\n')
    tmp <- paste0('Rscript ', file.path('$SCRIPT_DIR','scripts','post-processing-assess-mixing.R'),
                  ' -source_dir $SCRIPT_DIR -stanModelFile $STAN_MODEL_FILE -outdir $OUT_DIR -job_tag $JOB_TAG -local $LOCAL')
    cmd2 <- paste0(cmd2,tmp,'\n')
    tmp <- paste0('Rscript ', file.path('$SCRIPT_DIR','scripts','post-processing-sampling-prob-incident-cases.R'),
                  ' -source_dir $SCRIPT_DIR -stanModelFile $STAN_MODEL_FILE -indir $IN_DIR -outdir $OUT_DIR -job_tag $JOB_TAG -undiagnosed $UNDIAGNOSED -job_tag_undiag $JOB_TAG_UNDIAGNOSED -analysis $ANALYSIS')
    cmd2 <- paste0(cmd2,tmp,'\n')
    tmp <- paste0('Rscript ', file.path('$SCRIPT_DIR','scripts','post-processing-posterior-sources.R'),
                  ' -source_dir $SCRIPT_DIR -indir $IN_DIR -outdir $OUT_DIR -stanModelFile $STAN_MODEL_FILE -job_tag $JOB_TAG -analysis $ANALYSIS')
    cmd2 <- paste0(cmd2,tmp,'\n')
    # write submission file
    post.processing.file <- file.path(tmpdir2, 'post_processing.sh')
    cat(cmd2, file=post.processing.file)
    # set permissions
    Sys.chmod(post.processing.file, mode='644')
  }
  #	schedule post-processing
  cmd		<- paste0( cmd, 'echo "----------- Post-processing: ------------"\n')
  tmp		<- paste("if [ $(find ",tmpdir2," -name '*_stanout.RData' | wc -l) -ge ",hmc_chains_n," ]; then\n",sep='')
  cmd		<- paste(cmd,tmp,sep='')
  post.processing.file <- file.path(tmpdir2, 'post_processing.sh')
  cmd 	<- paste0(cmd, '\tcd ', dirname(post.processing.file),'\n')
  cmd 	<- paste0(cmd,'\tqsub ', basename(post.processing.file),'\n')
  cmd		<- paste0(cmd,"fi\n")
  cmd		<- paste(cmd, "rm -rf $CWD/", basename(args$source_dir[i]),'\n',sep='')
  cat(cmd)
  cmds[[i]]	<- cmd
}
if(args$cmdstan[1]==0)
{
  pbshead <- make.PBS.header(	hpc.walltime=23,
                              hpc.select=1,
                              hpc.nproc=1,
                              hpc.mem= "30gb",
                              hpc.load= "module load anaconda3/personal\nsource activate src",
                              hpc.q=NA,
                              hpc.array= length(cmds) )
}
if(args$cmdstan[1]==1)
{
  pbshead <- make.PBS.header(	hpc.walltime=23,
                              hpc.select=1,
                              hpc.nproc=hpc.nproc.cmdstan,
                              hpc.mem= paste0(hpc.nproc.cmdstan*9,'gb'),
                              hpc.load= paste0("module load cmdstan/2.33.0 anaconda3/personal\nsource activate src\nexport STAN_NUM_THREADS=",hpc.nproc.cmdstan,"\nexport TBB_CXX_TYPE=gcc\nexport CXXFLAGS+=-fPIE"),
                              hpc.q=NA,
                              hpc.array= length(cmds) )
}

#	make array job
for(i in seq_len(nrow(args)))
{
  cmds[[i]] <- paste0(i,')\n',cmds[[i]],';;\n')
}
cmd		<- paste0('case $PBS_ARRAY_INDEX in\n',paste0(cmds, collapse=''),'esac')
cmd		<- paste(pbshead,cmd,sep='\n')

#	submit job
outfile		<- gsub(':','',paste("src",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),'sh',sep='.'))
outfile		<- file.path(args$out_dir[1], outfile)
cat(cmd, file=outfile)
cmd 		<- paste("qsub", outfile)
cat(cmd)
cat(system(cmd, intern= TRUE))
