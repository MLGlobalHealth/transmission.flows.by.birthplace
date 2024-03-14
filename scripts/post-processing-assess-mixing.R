# assess-mixing.R
# 
###############################################################################

cat(" \n -------------------------------- \n \n Running assess-mixing.R \n \n -------------------------------- \n")

suppressMessages(library(rstan, quietly = TRUE))
suppressMessages(library(data.table, quietly = TRUE))
suppressMessages(library(dplyr, quietly = TRUE))
suppressMessages(library(bayesplot, quietly = TRUE))
suppressMessages(library(hexbin, quietly = TRUE))

# for dev purposes
if(0){
  args_dir <- list(
	  source_dir = '~/Documents/GitHub/transmission.flows.by.birthplace',
	  indir = '~/Box\ Sync/Roadmap/source_attribution',
	  outdir = '/Users/alexb/Documents/GitHub/transmission.flows.by.birthplace/out_Amsterdam/mm_bgUnif_piGP_221027b-agegps_sensanalysis_210216_MSM-618873',
	  stanModelFile = 'mm_bgUnif_piGP_221027b',
	  job_tag = 'agegps_sensanalysis_210216_MSM',
	  local = 1
  )
}

args_line <-  as.list(commandArgs(trailingOnly=TRUE))
if(length(args_line) > 0) 
{
	stopifnot(args_line[[1]]=='-source_dir')
	stopifnot(args_line[[3]]=='-stanModelFile')
	stopifnot(args_line[[5]]=='-outdir')
	stopifnot(args_line[[7]]=='-job_tag')
	stopifnot(args_line[[9]]=='-local')
	
	args_dir <- list()
	args_dir[['source_dir']] <- args_line[[2]]
	args_dir[['stanModelFile']] <- args_line[[4]]
	args_dir[['outdir']] <- args_line[[6]]
	args_dir[['job_tag']] <- args_line[[8]]
	args_dir[['local']] <- args_line[[10]]
} 

## start script
cat(" \n -------------------------------- with post-processing arguments -------------------------------- \n")
str(args_dir)

outfile.base <- paste0(args_dir$out_dir, "/",
                       args_dir$stanModelFile , "-", args_dir$job_tag)

file <- paste0(outfile.base,'-stanout-fit.RDS')
cat("\n read RDS:", file)
fit <- readRDS( file )

#
# extract neff 
#

cat("\n ----------- calculate pars with small neff: start ----------- \n")
summary.par = summary(fit)$summary
neff <- as.numeric(summary.par[, which(colnames(summary.par) == "n_eff")])
Rhat <- summary.par[, which(colnames(summary.par) == "Rhat")]

bound <- 500
if(sum(neff < bound,na.rm = TRUE) >0){
  pars.with.small.neff <- summary.par[which(neff < bound),,drop=FALSE]
}else{
  pars.with.small.neff <- summary.par[order(neff)[1:10],]
}

cat("\n Write to file",paste0(outfile.base,'-pars-with-small-neff.rds'))
saveRDS(pars.with.small.neff, file=paste0(outfile.base,'-pars-with-small-neff.rds'),version = 2)

# find worst parameter
wstpar <- rownames(summary.par)[which.min(summary.par[-length(rownames(summary.par)),'n_eff'])]

pars.small.neff <- summary.par[order(neff)[1:10],]
lp <- grepl("lp__",rownames(pars.small.neff))
pars.small.neff <- pars.small.neff[-lp,]
write.csv(pars.small.neff,file=paste0(outfile.base,'-pars-with-smallest-neff.csv'))
saveRDS(pars.small.neff,file=paste0(outfile.base,'-pars-with-smallest-neff.rds'),version = 2)
cat("\n ----------- calculate pars with small neff: end ----------- \n")

#
# Sampler diagnostics
#
cat("\n ----------- report sampler diagnostics: start ----------- \n")
sampler_params <- rstan::get_sampler_params(fit, inc_warmup = FALSE)
sampler_diagnostics <- data.table()
for (i in colnames(sampler_params[[1]])) {
  tmp <- data.table(t(sapply(sampler_params, function(x) quantile(x[, i],probs = c(0.025,0.5,0.975)))))
  tmp[, diagnostics:=i ]
  tmp[, chain:= seq_len(length(sampler_params))]
  sampler_diagnostics <- rbind(sampler_diagnostics, tmp)
}
saveRDS(sampler_diagnostics,file=paste0(outfile.base,'-sampler_diagnostics.rds'),version = 2)
cat("\n ----------- report sampler diagnostics: end ----------- \n")


#
# Monitor posterior distribution 
#
fit.target.pars <- c('logit_y_mix_0','gpscale','lscale[1]','lscale[2]',"lp__")
target.pars <- names(fit)[ grepl(paste(paste0('^',fit.target.pars),collapse = '|'),names(fit)) ]
summary <- rstan::monitor(rstan::extract(fit, pars=target.pars, permuted = FALSE, inc_warmup = FALSE))
saveRDS(summary, file=paste0(outfile.base,'-summary-pars.rds'),version = 2)


#
# Trace plots
#
cat("\n ----------- make trace plots: start ----------- \n")
color_scheme_set("mix-blue-red")
  p <- rstan::traceplot(fit, pars=target.pars, inc_warmup=TRUE, ncol = 1)
  ggsave(p, file = paste0(outfile.base, "-HMC-trace.png"), w=10, h=length(target.pars)*2)
  p <- rstan::traceplot(fit, pars='lp__', inc_warmup=TRUE, ncol = 1)
  ggsave(p, file = paste0(outfile.base, "-HMC-log_likelihood.png"), w=10, h=8)
  p <- rstan::traceplot(fit, pars=wstpar, inc_warmup=TRUE, ncol = 1)
  ggsave(p, file = paste0(outfile.base, "-HMC-trace-worst_par.png"), w=8, h=3)
  cat("\n ----------- make trace plots: end ----------- \n")

#
# Pairs plots
#
cat("\n ----------- make pairs plots: start ----------- \n")
  p <-	mcmc_pairs(rstan::extract(fit, pars=target.pars, permuted=FALSE, inc_warmup=FALSE), diag_fun = "dens", off_diag_fun = "hex")
  ggsave(p, file = paste0(outfile.base, "-HMC-pairs.png"), w=length(target.pars)*2, h=length(target.pars)*2)
cat("\n ----------- make pairs plots: end ----------- \n")

cat(" \n -------------------------------- \n \n End assess-mixing.R \n \n -------------------------------- \n")
