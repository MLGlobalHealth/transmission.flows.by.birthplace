
## preamble ----
require(data.table)

if (1)
{
  args_dir <- list(
    source_dir = '~/Documents/GitHub/transmission.flows.by.birthplace',
    #indir = '~/Box\ Sync/Roadmap/source_attribution',
    indir = '~/Box\ Sync/Roadmap',
    analysis = 'analysis_220713',
    pairs.dir = '/Users/alexb/Documents/GitHub/transmission.flows.by.birthplace/out_Amsterdam/agegps_updated_criteria_MSM-2010_2022',
    #outdir = '/Users/alexb/Documents/GitHub/source.attr.with.infection.time.fork/out_Amsterdam/mm_bgUnif_piGP_221027b-agegps_TE16_MSM-2010_2022-1665619',
    outdir = '/Users/alexb/Documents/GitHub/transmission.flows.by.birthplace/out_Amsterdam/mm_bgUnif_piGP_221027b-agegps_sensanalysis_210216_MSM-618873',
    clock_model = '/Users/alexb/Box Sync/Roadmap/source_attribution/molecular_clock/hierarchical',
    stanModelFile = 'mm_bgUnif_piGP_221027b',
    scenario = 15,
    reps = 1,
    rep = 1,
    simulate_data = T,
    job_tag = 'agegps_sensanalysis_210216_MSM-2010_2022'
  )
}

#	function to make PBS header
make.PBS.header <- function(hpc.walltime=23, hpc.select=1, hpc.nproc=8, hpc.mem= "80gb", hpc.load= "module load anaconda3/personal\nsource activate bpm", hpc.q=NA, hpc.array=1 )
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


cmd2 <- make.PBS.header( hpc.walltime=23,
                         hpc.select=1,
                         hpc.nproc=48,
                         hpc.mem= "124gb",
                         hpc.load= "module load anaconda3/personal\nsource activate bpm",
                         hpc.q=NA,
                         hpc.array= 1)
cmd2 <- paste0(cmd2,'\n')
# set up env variables
cmd2 <- paste0(cmd2,'SCRIPT_DIR=',args$source_dir[i],'\n',
               'IN_DIR=',args$in_dir[i],'\n',
               'OUT_DIR=',tmpdir2,'\n',
               'JOB_TAG=',args$job_tag[i],'\n',
               'STAN_MODEL_FILE=',args$stanModelFile[i],'\n',
               'ANALYSIS=',args$analysis[i],'\n',
               'NUMB_CHAINS=', max(args$chain),'\n',
               'OVERWRITE=0\n',
               'TRSM=', args$trsm[i],'\n',
               'START_Y=', args$start_d[i],'\n',
               'END_Y=', args$end_d[i],'\n',
               'UNDIAG_JOB=', args$undiag_job[i],'\n'
)
# save posterior samples
tmp <- paste0('Rscript ', file.path('$SCRIPT_DIR','scripts','check-convergence.R'),
              ' -stanModelFile $STAN_MODEL_FILE -out_dir $OUT_DIR -job_tag $JOB_TAG -numb_chains $NUMB_CHAINS -source_dir $SCRIPT_DIR -trsm $TRSM')
cmd2 <- paste0(cmd2,tmp,'\n')
# postprocessing assess HMC mixing
tmp <- paste0('Rscript ', file.path('$SCRIPT_DIR','scripts','make-panel-plot-nonB-cases.R'),
              ' -stanModelFile $STAN_MODEL_FILE -out_dir $OUT_DIR -job_tag $JOB_TAG')
cmd2 <- paste0(cmd2,tmp,'\n')
# postprocessing assess HMC mixing
tmp <- paste0('Rscript ', file.path('$SCRIPT_DIR','scripts','make-plot-prevalence-unsuppressed-flows.R'),
              ' -stanModelFile $STAN_MODEL_FILE -out_dir $OUT_DIR -job_tag $JOB_TAG')
cmd2 <- paste0(cmd2,tmp,'\n')
# postprocessing make plot  of transmission pars
tmp <- paste0('Rscript ', file.path('$SCRIPT_DIR','scripts','make-plot-stratified-flows-heatmap-temporal-trends.R'),
              ' -stanModelFile $STAN_MODEL_FILE -out_dir $OUT_DIR -job_tag $JOB_TAG -source_dir $SCRIPT_DIR')
cmd2 <- paste0(cmd2,tmp,'\n')
# postprocessing importations
tmp <- paste0('Rscript ', file.path('$SCRIPT_DIR','scripts','make-table-characteristics.R'),
              ' -stanModelFile $STAN_MODEL_FILE -analysis $ANALYSIS -in_dir $IN_DIR -out_dir $OUT_DIR -job_tag $JOB_TAG -trsm $TRSM -source_dir $SCRIPT_DIR')
cmd2 <- paste0(cmd2,tmp,'\n')
# postprocessing posterior predictive check
tmp <- paste0('Rscript ', file.path('$SCRIPT_DIR','scripts','make-table-cases-sequenced-sources.R'),
              ' -stanModelFile $STAN_MODEL_FILE -analysis $ANALYSIS -in_dir $IN_DIR -out_dir $OUT_DIR -job_tag $JOB_TAG -trsm $TRSM -source_dir $SCRIPT_DIR')
cmd2 <- paste0(cmd2,tmp,'\n')

# write submission file
post.processing.file <- file.path(tmpdir2, 'post_processing.sh')
cat(cmd2, file=post.processing.file)
# set permissions
Sys.chmod(post.processing.file, mode='644')


#	schedule post-processing
cmd		<- paste0( cmd, 'echo "----------- Post-processing: ------------"\n')
tmp		<- paste("if [ $(find ",tmpdir2," -name '*_stanout.RData' | wc -l) -ge ",max( args$chain )," ]; then\n",sep='')
cmd		<- paste(cmd,tmp,sep='')
post.processing.file <- file.path(tmpdir2, 'post_processing.sh')
cmd 	<- paste0(cmd, '\tcd ', dirname(post.processing.file),'\n')
cmd 	<- paste0(cmd,'\tqsub ', basename(post.processing.file),'\n')
cmd		<- paste0(cmd,"fi\n")
cmd		<- paste(cmd, "rm -rf $CWD/", basename(args$source_dir[i]),'\n',sep='')
cat(cmd)
cmds[[i]]	<- cmd

pbshead <- make.PBS.header(	hpc.walltime=23,
                            hpc.select=1,
                            hpc.nproc=hpc.nproc.cmdstan,
                            hpc.mem= paste0(hpc.nproc.cmdstan*9,'gb'),
                            hpc.load= paste0("module load cmdstan/2.33.0 anaconda3/personal\nsource activate bpm\nexport STAN_NUM_THREADS=",hpc.nproc.cmdstan,"\nexport TBB_CXX_TYPE=gcc\nexport CXXFLAGS+=-fPIE"),
                            hpc.q=NA,
                            hpc.array= length(cmds) )

#	submit job
outfile		<- gsub(':','',paste("bpm",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),'sh',sep='.'))
outfile		<- file.path(args$out_dir[1], outfile)
cat(cmd, file=outfile)
cmd 		<- paste("qsub", outfile)
cat(cmd)
cat(system(cmd, intern= TRUE))

