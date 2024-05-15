
## preamble ----
require(data.table)
require(bayesplot)
require(hexbin)
require(knitr)
require(ggplot2)
require(rstan)  # run Stan from R
require(cmdstanr)
require(ggsci)
require(scales)
require(grid)
require(viridis)
require(loo)

if (1)
{
  args <- list(
    source_dir = '~/Documents/GitHub/transmission.flows.by.birthplace',
    indir = '~/Box\ Sync/Roadmap/source_attribution',
    outdir = '/Users/alexb/Documents/GitHub/transmission.flows.by.birthplace/out_Amsterdam/mm_bgUnif_piGP_221027b-published_data_MSM-2010_2021-860418',
    stanModelFile = 'mm_bgUnif_piGP_221027b',
    job_tag = 'published_data_MSM-2010_2021',
    local = 1
  )
}

## command line parsing if any
args_line <-  as.list(commandArgs(trailingOnly=TRUE))
if(length(args_line) > 0)
{
  stopifnot(args_line[[1]]=='-source_dir')
  stopifnot(args_line[[3]]=='-stanModelFile')
  stopifnot(args_line[[5]]=='-outdir')
  stopifnot(args_line[[7]]=='-job_tag')
  stopifnot(args_line[[9]]=='-local')

  args <- list()
  args[['source_dir']] <- args_line[[2]]
  args[['stanModelFile']] <- args_line[[4]]
  args[['outdir']] <- args_line[[6]]
  args[['job_tag']] <- args_line[[8]]
  args[['local']] <- args_line[[10]]
}
args

cat(" \n --------------------------------  with arguments -------------------------------- \n")

## read stanin
cat('\nReading Stan input data...')
infile.stanin <- list.files(args$outdir, pattern=paste0('_stanin.RData$'), recursive=TRUE)[1]
stopifnot(length(infile.stanin)>0)
stopifnot(length(infile.stanin)<=1)
tmp <- load(file.path(args$outdir, infile.stanin))
stopifnot(c('args','stan_data')%in%tmp)

outfile.base <- paste0(args$outdir, "/",
                       args$stanModelFile , "-", args$job_tag)

if(args$local==1){
  file <- paste0(outfile.base,'-fitted_stan_model.rds')
}else{
  file <- paste0(outfile.base,'-stanout-fit.RDS')
}
cat("\n read RDS:", file)
model_fit <- readRDS(file = file)

# Check convergence and diagnostics ----

cat(" \n -------------------------------- \n Check convergence and divergences \n -------------------------------- \n")
fit.target.pars <- c('logit_y_mix_0','gpscale','lscale[1]','lscale[2]',"lp__")

po <- model_fit$draws(inc_warmup = TRUE,
                      #format = 'draws_df',
                      variables = fit.target.pars
)
su <- as.data.table(posterior::summarise_draws(po))
tmp <- paste0(outfile.base,"-convergence.csv")
write.csv(su, file = tmp, row.names = TRUE)
head(su)
su[,min(ess_bulk)]
su[,max(rhat)]

# Trace plots ----

color_scheme_set("mix-blue-red")
p <- bayesplot:::mcmc_trace(po,
                            pars = fit.target.pars,
                            n_warmup = 5e2,
                            facet_args = list(ncol = 1, labeller = label_parsed))+
  bayesplot:::legend_move("bottom")+ggplot2::labs(x="Iteration")+ggplot2::theme(text = element_text(size = 16))+
  xaxis_text(size = 16)+facet_text(size = 16)+legend_text(size = 16)+yaxis_text(size = 16)+xaxis_ticks(size = 14)+
  yaxis_ticks(size = 14)
p
ggsave(file = paste0(outfile.base,'-trace_allpars.pdf'), p, w = 10, h = length(tmp)*10)

tmp <- su$variable[which.min(su$ess_bulk)]
po <- model_fit$draws(inc_warmup = TRUE,
                      #format = 'draws_df',
                      variables = tmp
)
p <- bayesplot:::mcmc_trace(po,
                            pars = tmp,
                            n_warmup = 5e2,
                            facet_args = list(ncol = 1, labeller = label_parsed))+
  bayesplot:::legend_move("bottom")+ggplot2::labs(x="Iteration")+ggplot2::theme(text = element_text(size = 16))+
  xaxis_text(size = 16)+facet_text(size = 16)+legend_text(size = 16)+yaxis_text(size = 16)+xaxis_ticks(size = 14)+
  yaxis_ticks(size = 14)+coord_cartesian(ylim=c(500,1500))

p
ggsave(file = paste0(outfile.base,'-trace_lwstneff.pdf'), p, w = 10, h = 5)

# Pairs plots ----
#
cat("\n ----------- make pairs plots: start ----------- \n")
pd <- model_fit$draws(inc_warmup = FALSE,
                      variables = c(fit.target.pars))
bayesplot:::color_scheme_set("mix-blue-pink")
p <- bayesplot:::mcmc_pairs(pd,
                            pars = c(fit.target.pars),
                            diag_fun = "dens",
                            off_diag_fun = "hex"
)
ggsave(p, file = paste0(outfile.base, "-HMC-pairs_transmission_pars.pdf"), w=length(fit.target.pars)*2, h=length(fit.target.pars)*2)
cat("\n ----------- make pairs plots: end ----------- \n")


