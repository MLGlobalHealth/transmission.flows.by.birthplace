cat(" \n -------------------------------- \n \n Running save-posterior-samples.r \n \n -------------------------------- \n")

suppressMessages(library(data.table, quietly = TRUE))
suppressMessages(library(dplyr, quietly = TRUE))
suppressMessages(library(abind, quietly = TRUE))


`%notin%` <- Negate(`%in%`)

# for testing
if(0){
	args_dir <- list()
	args_dir[['stanModelFile']] <- 'mm_sigHierG_bgUnif_piReg_221123'
	args_dir[['outdir']] <- '/rds/general/project/ratmann_roadmap_data_analysis/live/source_attribution/mm_sigHierG_bgUnif_piReg_221123-simulations_network_500truepairs_prop_subsample100pct-scenario_1-464913'
	args_dir[['job_tag']] <- 'simulations_network_500truepairs_prop_subsample100pct'
	args_dir[['numb_chains']] <- 4
	args_dir[['source_dir']] <- '~/git/source.attr.with.infection.time'
	args_dir[['trsm']] <- 'MSM'
}

# save args for report before loading those from running session
args_line <-  as.list(commandArgs(trailingOnly=TRUE))
if(length(args_line) > 0)
{
  stopifnot(args_line[[1]]=='-source_dir')
  stopifnot(args_line[[3]]=='-stanModelFile')
	stopifnot(args_line[[5]]=='-outdir')
	stopifnot(args_line[[7]]=='-job_tag')
	stopifnot(args_line[[9]]=='-numb_chains')
	stopifnot(args_line[[11]]=='-trsm')

	args_dir <- list()
	args_dir[['source_dir']] <- args_line[[2]]
	args_dir[['stanModelFile']] <- args_line[[4]]
	args_dir[['outdir']] <- args_line[[6]]
	args_dir[['job_tag']] <- args_line[[8]]
	args_dir[['numb_chains']] <- as.integer(args_line[[10]])
	args_dir[['trsm']] <- args_line[[12]]
}

# start script
args_dir[['work_dir']] <- getwd()
source(file.path(args_dir$source_dir, 'R', 'functions.R'))

cat(" \n --------------------------------  with arguments -------------------------------- \n")
str(args_dir)

dout <- data.table(F=list.files(args_dir$outdir, pattern='*_stanout.RData$', recursive=TRUE, full.name=TRUE))
cat(paste("\n", nrow(dout),"/",args_dir$numb_chains, "chains are finished \n"))

#	load all input variables for this analysis run
z <- load( gsub('pbs_stanout.RData','pbs_stanin.RData',dout[1,F]) )
str(args)

outfile.base <- unique( dout[, file.path(dirname(dirname(F)), paste0(args_dir$stanModelFile,'-',args_dir$job_tag))] )
stopifnot(length(outfile.base)==1 )

cat(" \n -------------------------------- load: fit -------------------------------- \n")
#	reading job output, merge separate stanfits into one consolidated stanfit object
rf <- vector('list', nrow(dout))
median_lp_ <- vector('numeric', nrow(dout))
for(i in seq_len(nrow(dout)))
{
	cat('Loading output in ',dout[i,F],'\n')
	z <- load(dout[i,F])
	stopifnot('fit' %in% z)
	median_lp_[i] = median(rstan:::extract(fit)$lp__)
	#if(all(rstan::summary(fit)$summary[,1] == 0) & all(is.na(rstan::summary(fit)$summary[,2]))) next
	rf[[i]] <- fit
}
fit <- rstan:::sflist2stanfit(rf[lapply(rf,length)>0])
re <- rstan::extract(fit)

cat(" \n -------------------------------- save: fit -------------------------------- \n")
file = paste0(outfile.base,'-stanout-fit.RDS')
cat("\n save file:", file)
io_saveRDS(fit, args_dir[['work_dir']], dirname(file), basename(file), check_if_saved_n=10)
fit <- NULL
rf <- NULL
gc()

cat(" \n -------------------------------- save: pars -------------------------------- \n")
file =paste0(outfile.base,'-stanout-pars.RDS')
cat("\n save file:", file)
cat("\n saving objects:", names(re))
io_saveRDS(re, args_dir[['work_dir']], dirname(file), basename(file), check_if_saved_n=10)
gc()

cat(" \n -------------------------------- save: tpair_prob gq -------------------------------- \n")
tpair_prob_w_gqs = c("tpair_prob_w")
if(all(tpair_prob_w_gqs %in% names(re))){
	file = paste0(outfile.base,'-stanout-tpairprobw-gqs.RDS')
	io_saveRDS(re[tpair_prob_w_gqs], args_dir[['work_dir']], dirname(file), basename(file), check_if_saved_n=10)
	re[tpair_prob_w_gqs] <- NULL
}

cat(" \n -------------------------------- \n \n End post-processing.r \n \n -------------------------------- \n")

