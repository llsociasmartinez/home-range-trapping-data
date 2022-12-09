Sys.setenv(LANG = "en")
memory.limit(size=100000000)
source(paste0(getwd(),"/Ecography/scripts/b_functions_sim_20220918.R"))
source(paste0(getwd(),"/Ecography/scripts/c_functions_gam_hrerror_20220918.R"))


pathcomplete<-paste0(getwd(),"/Ecography/files/simulations/complete/") #for complete datasets to be used in further scripts

#load data
filenames<-list.files(pathcomplete,pattern="results_all_long")
load(paste0(pathcomplete,filenames[[1]]))

#create dataset for gam---------
data.gam<-results%>% mutate(id=as.factor(id))%>% mutate_if(is.character,as.factor)

data.gam<-data.gam %>% 
  mutate(approach=interaction(sampling,estimator,model),
         oapproach=ordered(approach),
         lg10nbobs=log10(nbobs)) %>% select(-models_diag) 

#filter extreme outliers and create weights for eventual weighted regresion
data.gam2<-data.gam %>% filter(area<10,area>-10) %>%
  group_by(approach,area.c,trap.d,nbobs) %>%
  group_map(~{
    dt<-.x %>%
      mutate(ssize=length(which(!is.na(area))))#how many datapoints on which the gam fit will be based
    return(dt)},.keep=T) %>% bind_rows %>% ungroup

data.gam3<-data.gam2 %>% na.omit %>% 
  group_by(approach,area.c,trap.d) %>% 
  group_map(~{
    dt<-.x %>% mutate(
      ssizee=ssize/mean(ssize,na.rm=T),#how does the ssize compare to the mean across all nbobs in this approach
      weightt=lg10nbobs*ssizee)#weight based on how much ids produced data and how many number of obs per id
    return(dt)},.keep=T) %>% bind_rows %>% ungroup

datagam<-data.gam3 %>% select(id,area,sampling,model,estimator,approach,oapproach,area.c,trap.d,nbobs,lg10nbobs,ssizee,weightt)

#save datasets
save(datagam,file=paste0(pathcomplete,"datagam_", datagam$id %>% unique %>% length,".Rdata"))

#place the knots of nbobs towards the lower end
xlim<-datagam$nbobs %>% range
knbobs<-15#number of knots
b<-10#base logarithm
knots <- list(lg10nbobs = (seq(log(xlim[1],base = b), log(xlim[2],base=b), length.out = knbobs)))

#type of spline
spl<-c("cr","tp","tp")
#runs 4H as it is now on laptop, 12H with weights + not converge
t<-Sys.time()
gmodel <- bam(form.ti,
              data = datagam,
              discrete=T,
              nthreads=4,
              method="fREML",
              knots=knots,
              #weights=weightt #weighted regresion does not converge
)
print(Sys.time()-t) 

save(gmodel,file=paste0(pathcomplete,"gam_gauss_ti_ids",datagam$id %>% unique %>% length,".Rdata"))
#save(gmodel,file=paste0(pathcomplete,"gam_gauss_weight_ti_ids",datagam$id %>% unique %>% length,".Rdata"))

#END------------