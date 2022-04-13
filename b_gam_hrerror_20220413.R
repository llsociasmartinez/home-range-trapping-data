Sys.setenv(LANG = "en")
memory.limit(size=100000000)
source(paste0(getwd(),"/scripts/a_functions_sim_20220413.R"))
source(paste0(getwd(),"/scripts/b_functions_gam_hrerror_20220413.R"))

#load data
pathmike<-yourpath#
filenames<-list.files(pathmike,pattern="results_all_long")
load(paste0(pathmike,filenames[[1]]))

#create dataset for gam---------
data.gam<-results%>% 
  mutate(id=as.factor(id))%>% mutate_if(is.character,as.factor)

data.gam<-data.gam %>% 
  mutate(intsem=interaction(sampling,estimator,model),
         ointsem=ordered(intsem)) %>% select(-models_diag) 

#filter extreme outliers and create weights for eventual wighted regresion
data.gam2<-data.gam %>% filter(area<10,area>-10) %>%
  group_by(intsem,area.c,trap.d,nbobs) %>%
  group_map(~{
    dt<-.x %>%
      mutate(sdev=sd(area,na.rm=T),
             ssize=length(which(!is.na(area))))
    return(dt)},.keep=T) %>% bind_rows %>% ungroup

data.gam3<-data.gam2 %>% na.omit %>% 
  group_by(intsem,area.c,trap.d) %>% 
  group_map(~{
    dt<-.x %>% 
      mutate(
        ssizee=ssize/mean(ssize,na.rm=T),
        weight=(1/(sdev+0.01)),
        logweight=log10(weight),
        difmin=abs(min(logweight)),
        mlw=((logweight+difmin)+0.01)*ssizee,
        mmlw=mlw/mean(mlw))
    return(dt)},.keep=T) %>% bind_rows %>% ungroup


data.gam4<-data.gam3%>% 
  group_by(id) %>% group_split

dtm<-data.gam4 %>% bind_rows
dtm<-dtm %>% select(id,area,intsem,ointsem,area.c,trap.d,nbobs,mmlw)

#save datasets
save(data.gam3,dtm,file=paste0(path,"datagam_",
                               dtm$id %>% unique %>% length,".Rdata"))

#place the knots of nbobs towards the lower end
xlim<-dtm$nbobs %>% range
knbobs<-15#number of knots
b<-10#base logarithm
knots = list(`log10(nbobs)` = (seq(log(xlim[1],base = b), log(xlim[2],base=b), length.out = knbobs)))

knots
#type of spline
spl<-c("cr","tp","tp")
t<-Sys.time()
gmodel <- bam(form.ti,
              data = dtm,
              discrete=T,
              nthreads=4,
              method="fREML",
              knots=knots#,
              #family=scat(link="identity"),#for eventual correction
              #weights=mmlw #weighted regresion
)
print(Sys.time()-t)
gmodel
save(gmodel,file=paste0(getwd(),"/files/gam_gauss_weight_ti_ids",dtm$id %>% unique %>% length,".Rdata"))

load(paste0(getwd(),"/files/gam_ti_ids",length(data.gam4),".Rdata"))

#checkings for gam fit
diags <-gmodel %>% appraise
diags
kscheck<-gmodel %>% gam.check
kscheck
