source(paste0(getwd(),"/scripts/a_functions_sim_20220413.R"))
Sys.setenv(LANG = "en")
memory.limit(size=100000000)

#set parameters---------
mac.linux<-FALSE #set to FALSE if using windows
path<-paste0(getwd(),"/files/sim_iid/")
pathmike<-paste0(getwd(),"/files/Mike_files/") #not working


allids.equal<-TRUE #all analyses applied to all individuals
nbids<-50 #number of individuals to simulate

paralel<-FALSE #should computation be carried in paralel
nbcores<-2 #number of cores to use
chunk<-1 #how many individuals per chunk to each core
save<-TRUE #should results be saved

#set the 95% area for simulation
true.area=100*10000 %#% 'm^2'
#the radius of attraction for traps in m
attract=200

#create the needed objects to carry the simulation more efficiently-----------
#p.intervals sets the variables values
presimul<-set_variables_process(allids.equal=T,
                                nbids=nbids,
                                true.area,
                                p.intervals=list(auto.c=list(from=c(0),to=c(0),by=c(1)),
                                                 area.c=list(from=c(0.2),to=c(1),by=c(0.2)),
                                                 trap.d=list(from=c(2),to=c(42),by=c(10)),
                                                 essize=list(from=c(2000),to=c(2000),by=c(1))),
                                bufferzone=0,
                                timewind=1 %#% 'hours',each=10 %#% 'minutes',
                                save=save,pathway=path)

#load presimulated datasets
#load(file=paste0(path,"presimul_allidsequal",".Rdata"))

#simulate the data--------
ids<-presimul$cmbn[["id"]] %>% unique
chks<-split(ids, ceiling(seq_along(ids)/chunk))
f<-1;last<-length(chks);todo<-f:last
if(paralel==T){plan(multisession, workers = nbcores)}
st<-Sys.time()
simulations<-future_map_dfr(.progress=T,.options = furrr_options(seed = TRUE),
                            .x=todo,~sim_capt(allids.equal,
                                              .x,
                                              presimul$cmbn[which(presimul$cmbn$id%in%chks[[.x]]),],
                                              presimul$cmbn.grid,
                                              presimul$cmbn.mdltmp,
                                              true.area,attract,save=save,pathway=path))
future:::ClusterRegistry("stop")

print(Sys.time()-st)




#check that all ids have the desired nb of observations-----------
#and correct for those that doesnt
simulations<-resim_capt_enough(ids.sim=simulations,cmbn.grid=presimul$cmbn.grid,cmbn.mdltmp=presimul$cmbn.mdltmp,
                                   attract,timewind=3 %#% 'hours',each=10 %#% 'minutes',buffer=1.4,
                                   save=save,pathway=path,paralel=T,nbcores=nbcores,chunk=10,mac.linux)
#if already simulated------
# folds<-c("sim_iid",paste0("sim_iid_",1:5))
# filenames<-map(folds,~list.files(paste0(getwd(),"/files/Mike_files/",.x,"/"),pattern="sim_total_2500"))
# #filenames<-map(filenames,~.x[grepl("all",.x)==F])
# folder<-paste0("/files/Mike_files/",folds,"/")
# simulations <-map2(filenames,folder,~{
#   all<-lapply(.x, function(.file){  #load.Rdata(filename, objname)
#     #print(.file)
#     load.Rdata(paste0(getwd(),.y,.file),"x")
#     dat<-x
#     #print(colnames(dat))
#     #print(dat)
# 
#     dat    # return the dataframe
#   })
#   all<-do.call(rbind,all)
#   all}
# )

# nbids_simulated<-1:(map_dbl(simulations,~.x$id %>% unique %>% length) %>% sum)
# ncombs<-simulations[[1]] %>% group_by(id) %>% summarise(nb=n()) %>% slice(1) %>% pull(nb)
# simulations<-simulations %>% bind_rows
# simulations$id<-map(nbids_simulated,~rep(.x,ncombs)) %>% unlist

colnames(simulations)[which(colnames(simulations)=="captures")]<-"data.ctmm"
save(simulations,file=paste0(path,"simulations_all_",
                            min(simulations$id %>% unique),"_",
                            max(simulations$id %>% unique),".Rdata"))

#load simulations-------------
#load(paste0(path,"simulations_all_1_600",".Rdata"))


#fit the estimators and the models------------
#set characteristics of analyses
#estimators to calculate
estimators<-c("akima_bic","ctmm_akde","amt_kde","amt_mcp","amt_locoh")
#asymp models to calculate
asympmodels<-c("micmen","monmol")


sampr<-list(normal=T,farthest=T,boot=F,boottimes=3) #sampling procedures, bootstrapping is impresively slow
recursive.model<-TRUE #should models be recursive (for each datapoint per id)
iv<-NULL #initial value for asymptotic model

paralel<-T #paralelised
nbcores<-2
chunk<-1 #nb ids done per partial saving (nb ids iterated each time)
save<-T #save results

#individuals to do from 1:600 #################
#I will do the first 20 myself until we get some way of computing the rest
idstodo<-c(1:nbids)#[c(1:600)%!in%c(1:4,20,310)]#600#5:6#1:4#:8#

#estimators----------
#p.intervals sets the points on which to conduct the estimations
st<-Sys.time()
estimators_df<-fit_estimators(allids.equal=T,
                               simulations[which(simulations$id%in%idstodo),],
                               presimul$cmbn.grid,
                               p.intervals=list(essize=list(from=c(2,10,20,30,60,100,500),
                                                            to=c(10,20,30,60,100,500,2000),
                                                            by=c(1,2,5,10,20,100,500))),
                               estimators,akima.res=1000,
                               sampr=sampr,
                               pathway=path,save=save,paralel=paralel,nbcores=nbcores,chunk=chunk,mac.linux)
print(Sys.time()-st)
# #str(estimators_df,max.level=1)

#if already fitted
# filenames<-list.files(paste0(getwd(),"/files/sim_iid/"),pattern="estimators_df_all")
# filenames<-list.files(paste0(getwd(),"/files/sim_iid/"),pattern="estimators_id")
# filenames<-list.files(pathmike,pattern="estimators_id")
# estimators_df <-map_dfr(filenames,~{
#   all<-lapply(.x, function(.file){  #load.Rdata(filename, objname)
#     #print(.file)
#     load.Rdata(paste0(pathmike,.file),"x")
#     dat<-x
#     #print(colnames(dat))
#     #print(dat)
# 
#     dat    # return the dataframe
#   })
#   all}
# )
# estimators_df[1:10,]
# load(paste0(getwd(),"/files/sim_iid/",filenames[[1]]))
# # results %>% head()
#str(estimators_df$raws[[50]],max.level=3)

#asymptotic models----------------
nbcores<-4
chunk<-2
paralel<-TRUE
save<-TRUE

results<-fit_asympmodels(estimators_df,estimators,asympmodels,recursiv.model=T,iv=NULL,
                         pathway=pathmike,save=save,paralel=paralel,nbcores=nbcores,chunk=chunk,mac.linux)
#pathmike
#if already fitted----------
#filenames<-list.files(paste0(getwd(),"/files/sim_iid/"),pattern="results_id")
# filenames<-list.files(pathmike,pattern="results_id")
# results <-map_dfr(filenames,~{
#   all<-lapply(.x, function(.file){  #load.Rdata(filename, objname)
#     #print(.file)
#     #load.Rdata(paste0(getwd(),"/files/sim_iid/",.file),"x")
#     load.Rdata(paste0(pathmike,.file),"x")
#     dat<-x
#     #print(colnames(dat))
#     #print(dat)
# 
#     dat    # return the dataframe
#   })
#   all}
# )

#str(results,max.level=3)
#results<-results %>% select(-data.amt,-data.ctmm,-rowstouse,-mtchgr,-idx,-obssample)
#colnames(results)

#convert results to longformat (for GAM modelling and plotting)------------------
results<-convert_long(results %>% select(-data.amt,-data.ctmm,-rowstouse,-mtchgr,-idx),estimators,true.area,pathway=pathmike,save=save)
save(results,file=paste0(path,"results_all_long_",
                             min(results$id %>% unique),"_",
                             max(results$id %>% unique),".Rdata"))



