
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
    source_dir = '~/Documents/GitHub/source.attr.with.infection.time',
    indir = '~/Box\ Sync/Roadmap/source_attribution',
    outdir = '/Users/alexb/Documents/GitHub/source.attr.with.infection.time/out_Amsterdam/mm_bgUnif_piGP_221027b-agegps_sensanalysis_210216_MSM-618873',
    stanModelFile = 'mm_bgUnif_piGP_221027b',
    job_tag = 'agegps_sensanalysis_210216_MSM'
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

  args <- list()
  args[['source_dir']] <- args_line[[2]]
  args[['stanModelFile']] <- args_line[[4]]
  args[['outdir']] <- args_line[[6]]
  args[['job_tag']] <- args_line[[8]]
}
args
replicate <- args$rep
cat(paste0("rep ", replicate))
## load functions
source(file.path(args$source_dir, 'R', 'functions_simulation_scenarios.R'))

cat(" \n --------------------------------  with arguments -------------------------------- \n")

## read stanin
cat('\nReading Stan input data...')
infile.stanin <- list.files(args$outdir, pattern=paste0('rep_1','_stanin.RData$'), recursive=TRUE)[1]
stopifnot(length(infile.stanin)>0)
stopifnot(length(infile.stanin)<=1)
tmp <- load(file.path(args$outdir, infile.stanin))
stopifnot(c('args','stan_data')%in%tmp)

do <- do[order(PAIR_ID),]

tmp <- paste0(outfile.base,'-rep_',replicate,'-fitted_stan_model.rds')
cat("\n Read fitted dat ", tmp , "\n")
model_fit <- readRDS(file = tmp)

cat(" \n -------------------------------- \n Check convergence and divergences \n -------------------------------- \n")
if(grepl('Vanilla',args$stanModelFile)){
  fit.target.pars <- c('logit_y_mix_0','y_mix','log_alpha1_pair[1]','log_phi_pair[1]',"lp__")
}else if(grepl('Reg',args$stanModelFile)){
  fit.target.pars <- c('y_mix[1]','logit_y_mix_0','beta_age_src[1]','beta_age_src[2]','beta_age_src[3]','beta_age_src[4]','beta_age_src[5]',
                       'beta_age_rcp[1]','beta_age_rcp[2]','beta_age_rcp[3]','beta_age_rcp[4]','beta_age_rcp[5]',"lp__")
}else if(grepl('GP',args$stanModelFile)){
  fit.target.pars <- c('logit_y_mix_0','gpscale','lscale[1]','lscale[2]',"lp__")
}
if(grepl('1DGP',args$stanModelFile)){
  if(grepl('_agesrcrec_',args$job_tag)){
      fit.target.pars <- c('logit_y_mix_0',"lp__","gpscale","lscale")
  }else{
      fit.target.pars <- c('logit_y_mix_0',"lp__","gpscale","lscale")
  }
}

po <- model_fit$draws(inc_warmup = TRUE,
                      #format = 'draws_df',
                      variables = fit.target.pars
)
su <- as.data.table(posterior::summarise_draws(po))
tmp <- paste0(outfile.base,'-gp_scenario',i,"-convergence.csv")
write.csv(su, file = tmp, row.names = TRUE)
head(su)
su[,min(ess_bulk)]
su[,max(rhat)]

args$file_stanModelGQs <- file.path(args$source_dir, 'stan_model_files',paste0(args$stanModelFile,'_GQs','.stan'))
gq_program <- args$file_stanModelGQs
mod_gq <- cmdstan_model(gq_program)
fit_gq <- mod_gq$generate_quantities(model_fit, data = stan_data, seed = 123)
log_lik <- fit_gq$draws("log_lik", format = "matrix")
LOO <- loo::loo(log_lik)
print(LOO)
saveRDS(LOO, file = paste0(outfile.base, "-LOO.rds"))
saveRDS(LOO, file = paste0(outfile.base, "-LOO.csv"))

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


## GP diagnostics----
fit.target.pars = c("lscale[1]","lscale[2]","gpscale")
if(grepl('1DGP',args$stanModelFile)) fit.target.pars <- c('logit_y_mix_0',"lp__","gpscale","lscale")
po <- model_fit$draws(inc_warmup = FALSE,
                      format = 'draws_df',
                      variables = fit.target.pars
)
po <- as.data.table(po)
setnames(po, colnames(po), gsub('^\\.','',colnames(po)))
po <- melt(po, id.vars = c('chain','iteration','draw'))
po <- po[,
         list( q = quantile(value, probs = c(0.5) ),
               stat = c('M')
         ),
         by = c('variable')
]
po <- dcast.data.table(po, variable~stat, value.var = 'q')

diagnostic <- function(l,l_hat) l_hat + 0.01 >= l
m_QE <- function(c,l,S) ceiling(1.75 * c / (l/S))
l_QE <- function(c,m,S) round(S * 1.75 * c / m, 3)

l_hat1 <- po[1,2]
l_hat2 <- po[2,2]
c(l_hat1, l_hat2)
check1 <- diagnostic(2, l_hat1)
check2 <- diagnostic(2, l_hat2)

print(paste0("Check for lengthscale 1: ", check1))
print(paste0("Check for lengthscale 2: ", check2))


