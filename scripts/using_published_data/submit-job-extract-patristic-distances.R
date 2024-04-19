require(data.table)
require(stringr)

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
    source_dir = '~/git/transmission.flows.by.birthplace',
    in.dir = '/rds/general/project/ratmann_roadmap_data_analysis/live/analysis_220713',
    out.dir = '/rds/general/project/ratmann_roadmap_data_analysis/live/transmission_sources/patristic_distances',
    #out.dir = 'data_Ams/analysis_220713/trees'
    trsm = 'AmsMSM'#,
    #REP = '000'
  )
  out.dir <- args$out.dir
}

if(1)
{
  tmp <- data.table(REP=0:2)
  tmp[, REP:= as.character(str_pad(REP, 3, pad = "0"))]
  set(args, NULL, colnames(tmp), NULL)
  tmp[, dummy:= 1L]
  args[, dummy:= 1L]
  args <- merge(args, tmp, by='dummy')
  set(args, NULL, 'dummy', NULL)
}

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
                    ' -indir ', args$in.dir[i],'',
                    ' -out.dir ', tmpdir,'',
                    ' -trsm ', args$trsm[i],'',
                    ' -REP ', args$REP[i],''
  )
  cmd    	<- paste0(cmd, tmp, '\n')

  #	general housekeeping
  cmd 	<- paste0( cmd, 'echo "----------- Copy files to out directory: ------------"\n')
  cmd		<- paste0(cmd, 'cp -R --no-preserve=mode,ownership "', tmpdir,'"/* ', out.dir,'\n')
  cmd		<- paste0(cmd, 'chmod -R g+rw ', out.dir,'\n')

  cmds[[i]]	<- cmd
}

pbshead <- make.PBS.header(hpc.walltime=23,
                            hpc.select=1,
                            hpc.nproc=1,
                            hpc.mem= "30gb",
                            hpc.load= "module load anaconda3/personal\nsource activate phylo",
                            hpc.q=NA,
                            hpc.array= length(cmds) )

#	make array job
for(i in seq_len(nrow(args)))
{
  cmds[[i]] <- paste0(i,')\n',cmds[[i]],';;\n')
}
cmd		<- paste0('case $PBS_ARRAY_INDEX in\n',paste0(cmds, collapse=''),'esac')
cmd		<- paste(pbshead,cmd,sep='\n')

#	submit job
outfile		<- gsub(':','',paste("pds",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),'sh',sep='.'))
outfile		<- file.path(args$out.dir[1], outfile)
cat(cmd, file=outfile)
cmd 		<- paste("qsub", outfile)
cat(cmd)
cat(system(cmd, intern= TRUE))
