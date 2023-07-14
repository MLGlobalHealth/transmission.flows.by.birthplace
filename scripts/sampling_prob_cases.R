
sp <- readRDS(file=paste0(outfile.base,'-sampling_prob_cases','.RDS'))
sp[, N_inf_CL:=  N/(1-p0.025)]
sp[, N_inf_UL:=  N/(1-p0.975)]
sp[, psi_CL:=  N_seq/N_inf_CL]
sp[, psi_CU:=  N_seq/N_inf_UL]

samples <- readRDS(file=file.path(outdir, paste0('samples_',job_tag,"_",args$trsm,'.rds')))
shape_msm <- data.table(reshape::melt(samples$wb_shape_grp))
setnames(shape_msm,c('iterations','Var.2'),c('iter','mg'))
shape_msm[, trsm:='MSM']
shape_msm[, par:='shape']

scale_msm <- data.table(reshape::melt(samples$wb_scale_grp))
setnames(scale_msm,c('iterations','Var.2'),c('iter','mg'))
scale_msm[, trsm:='MSM']
scale_msm[, par:='scale']

ds <- rbind(shape_msm,scale_msm)
ds <- dcast(ds,trsm+mg+iter~par,value.var="value")

dat <- tidyr::crossing(year=seq(2010,2021,1),month=seq(1,12,1))
#dat <- tidyr::crossing(year=seq(2010,2021,1))
ds <- merge(dat,ds,all=T)
ds <- data.table(ds)
ds[, time:=(2022+(1/12))-(year + (month/12))] # +(1/12) to add 1 month so we count until end Dec 2018/start of Jan 2019. but I think it makes sense someone infected Dec 2018 has 0 prob of being diagnosed by end of 2018 due to time taken to detect virus
ds[, p:=1 - pweibull(time,shape=shape,scale=scale)]

mean_y <- ds[, list(p=quantile(p,prob=c(0.025,0.5,0.975)),
                    qlabel=c('p0.025','p0.5','p0.975')),
             by=c('trsm','mg','year')] # summarise quantiles for each year
mean_y <- merge(mean_y,dmap,by.x='mg',by.y='mgid')
mean_y <- dcast(mean_y,trsm+mwmb+year~qlabel,value.var=c("p"))


spy <- dinf[YEAR_OF_INF_EST>=2010, list(N=length(TO_SEQUENCE_ID),
                                           N_seq = length(TO_SEQUENCE_ID[SEQ==T])),
             by=c('LOC_BIRTH_POS','YEAR_OF_INF_EST')]
# add rows for years with missing values for each birth region? would need to make the N=0.0001 or similar..
spy <- merge(spy,mean_y,by.x=c('LOC_BIRTH_POS','YEAR_OF_INF_EST'),by.y=c('mwmb','year'),all.x=T)
spy[, N_inf:= N/(1-p0.5)]
spy[, N_inf_CL:=  N/(1-p0.025)]
spy[, N_inf_UL:=  N/(1-p0.975)]

#spy <- spy[, list(N=sum(N),N_inf=sum(N_inf),N_seq=sum(N_seq)),by=c('LOC_BIRTH_POS')]
# calculate sampling prob among incident cases 
spy[, psi:= N_seq/N_inf]
spy[, psi_CL:=  N_seq/N_inf_CL]
spy[, psi_CU:=  N_seq/N_inf_UL]

sp[, YEAR_OF_INF_EST:= 'Average 2010-2021']
sp[, LOC_BIRTH_POS:= factor(LOC_BIRTH_POS,
                          levels=c('Netherlands','W.Europe,\nN.America,Oceania','Suriname &\nDutch Caribbean',
                                   'S. America &\n Caribbean','E. & C. Europe','MENA','Other'))]

g1 <- ggplot(sp) + geom_bar(aes(x=YEAR_OF_INF_EST,y=psi,fill=LOC_BIRTH_POS),stat='identity',position=position_dodge(width=0.9)) +
  geom_errorbar(aes(x=YEAR_OF_INF_EST,ymin=psi_CL, ymax=psi_CU,fill=LOC_BIRTH_POS),position=position_dodge(width=0.9),width=0.5, colour="black")	+
  scale_fill_npg() +
  labs(x='',fill='Birthplace of incident case', y='Sequence sampling fraction') +
  theme_bw() +
  theme(legend.pos='bottom')+
        #axis.title.x = element_blank(),
        #axis.text.x = element_text(angle=60, vjust = 0.95,hjust = 0.9)) +#+ #,
  coord_cartesian(ylim = c(0,1)) +
  scale_y_continuous(labels = scales::label_percent(accuracy = 1L),breaks=seq(0,1,0.2))

g2 <- ggplot(subset(spy,YEAR_OF_INF_EST<2022)) + geom_bar(aes(x=YEAR_OF_INF_EST,y=psi,fill=LOC_BIRTH_POS),stat='identity',position=position_dodge(width=0.9)) +
  geom_errorbar(aes(x=YEAR_OF_INF_EST,ymin=psi_CL, ymax=psi_CU,fill=LOC_BIRTH_POS),position=position_dodge(width=0.9),width=0.5, colour="black")	+
  scale_fill_npg() +
  labs(x='',fill='Birthplace of incident case', y='Sequence sampling fraction') +
  theme_bw() +
  theme(legend.pos='bottom')+
  #axis.title.x = element_blank(),
  #axis.text.x = element_text(angle=60, vjust = 0.95,hjust = 0.9)) +#+ #,
  coord_cartesian(ylim = c(0,1)) +
  scale_y_continuous(labels = scales::label_percent(accuracy = 1L),breaks=seq(0,1,0.2))

g <- ggarrange(g1,g2,ncol=2,align='h',widths=c(0.3,0.7),common.legend=T,legend='bottom')

ggsave(file = paste0(outfile.base,'-sequence_sampling_fraction_cases.png'),
       g, w = 9, h = 6)



g1 <- ggplot(sp) + geom_bar(aes(x=LOC_BIRTH_POS,y=psi,fill=LOC_BIRTH_POS),stat='identity',position=position_dodge(width=0.9)) +
  geom_errorbar(aes(x=LOC_BIRTH_POS,ymin=psi_CL, ymax=psi_CU,fill=LOC_BIRTH_POS),position=position_dodge(width=0.9),width=0.5, colour="black")	+
  scale_fill_npg() +
  labs(x='',fill='Birthplace of incident case', y='Sequence sampling fraction') +
  theme_bw() +
  theme(legend.pos='bottom')+
  #axis.title.x = element_blank(),
  #axis.text.x = element_text(angle=60, vjust = 0.95,hjust = 0.9)) +#+ #,
  coord_cartesian(ylim = c(0,1)) +
  scale_y_continuous(labels = scales::label_percent(accuracy = 1L),breaks=seq(0,1,0.2))

spy[, YEAR:= factor(YEAR_OF_INF_EST,levels=c('2010','2011','2012','2013','2014','2015','2016','2017','2018','2019','2020','2021'))]
g2 <- ggplot(subset(spy,YEAR_OF_INF_EST<2022)) + geom_bar(aes(x=LOC_BIRTH_POS,y=psi,fill=YEAR),stat='identity',position=position_dodge(width=0.9)) +
  geom_errorbar(aes(x=LOC_BIRTH_POS,ymin=psi_CL, ymax=psi_CU,fill=YEAR),position=position_dodge(width=0.9),width=0.5, colour="black")	+
  scale_fill_npg() +
  labs(x='',fill='Birthplace of incident case', y='Sequence sampling fraction') +
  theme_bw() +
  theme(legend.pos='bottom')+
  #axis.title.x = element_blank(),
  #axis.text.x = element_text(angle=60, vjust = 0.95,hjust = 0.9)) +#+ #,
  coord_cartesian(ylim = c(0,1)) +
  scale_y_continuous(labels = scales::label_percent(accuracy = 1L),breaks=seq(0,1,0.2))

g <- ggarrange(g1,g2,ncol=2,align='h',widths=c(0.3,0.7),common.legend=F,legend='bottom')

ggsave(file = paste0(outfile.base,'-sequence_sampling_fraction_cases_year.png'),
       g, w = 9, h = 6)


# facets
g1 <- ggplot(sp) + geom_bar(aes(x=YEAR_OF_INF_EST,y=psi,fill=LOC_BIRTH_POS),stat='identity',position=position_dodge(width=0.9)) +
  geom_errorbar(aes(x=YEAR_OF_INF_EST,ymin=psi_CL, ymax=psi_CU,fill=LOC_BIRTH_POS),position=position_dodge(width=0.9),width=0.5, colour="black")	+
  scale_fill_npg() +
  labs(x='',fill='Birthplace of incident case', y='Sequence sampling fraction') +
  theme_bw() +
  theme(legend.pos='bottom',
        axis.text.x = element_blank())+
  #axis.title.x = element_blank(),
  #axis.text.x = element_text(angle=60, vjust = 0.95,hjust = 0.9)) +#+ #,
  coord_cartesian(ylim = c(0,1)) +
  scale_y_continuous(labels = scales::label_percent(accuracy = 1L),breaks=seq(0,1,0.2))

spy[, LOC_BIRTH_POS:= factor(LOC_BIRTH_POS,
                            levels=c('Netherlands','W.Europe,\nN.America,Oceania','Suriname &\nDutch Caribbean',
                                     'S. America &\n Caribbean','E. & C. Europe','MENA','Other'))]
g2 <- ggplot(subset(spy,YEAR_OF_INF_EST<2022)) + geom_bar(aes(x=YEAR_OF_INF_EST,y=psi,fill=LOC_BIRTH_POS),stat='identity',position=position_dodge(width=0.9)) +
  geom_errorbar(aes(x=YEAR_OF_INF_EST,ymin=psi_CL, ymax=psi_CU,fill=LOC_BIRTH_POS),position=position_dodge(width=0.9),width=0.5, colour="black")	+
  scale_fill_npg() +
  facet_wrap(.~LOC_BIRTH_POS) +
  labs(x='',fill='Birthplace of incident case', y='Sequence sampling fraction') +
  theme_bw() +
  theme(legend.pos='bottom',
        strip.background = element_blank())+
  #axis.title.x = element_blank(),
  #axis.text.x = element_text(angle=60, vjust = 0.95,hjust = 0.9)) +#+ #,
  coord_cartesian(ylim = c(0,1)) +
  scale_y_continuous(labels = scales::label_percent(accuracy = 1L),breaks=seq(0,1,0.2))

g <- ggarrange(g1,g2,ncol=2,align='h',widths=c(0.3,0.7),common.legend=T,legend='bottom')

ggsave(file = paste0(outfile.base,'-sequence_sampling_fraction_cases_year_facets.png'),
       g, w = 9, h = 6)
