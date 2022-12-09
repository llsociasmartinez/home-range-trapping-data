#case study-----
Sys.setenv(LANG = "en")
memory.limit(size=100000000)
source(paste0(getwd(),"/Ecography/scripts/b_functions_sim_20220918.R"))
source(paste0(getwd(),"/Ecography/scripts/c_functions_gam_hrerror_20220918.R"))
library(sf)


#set parameters------
mac.linux=F
#estimators to use
estimators<-c("akima_bic","ctmm_akde","amt_kde","amt_mcp","amt_locoh")
estimators.lbs<-c("AKDE","KDE","MCP","LoCoH","BicubIt")
#asymp models to fit
asympmodels<-c("raw","micmen","monmol")
asympmodels.lbs<-c("Raw","MicMen","MonMol")
#sampling regimes
sampregs<-c("normal","fps") 
sampregs.lbs<-c("Time-ordered","Distance-ordered")
sampr<-list(normal=T,farthest=T,boot=F,boottimes=3) #sampling procedures, bootstrapping is impresively slow
recursive.model<-TRUE #should models be recursive (for each datapoint per id)
iv<-NULL #initial value for asymptotic model
paralel<-T #paralelised
nbcores<-2
chunk<-1 #nb ids done per partial saving (nb ids iterated each time)
save<-T #save results

#plots
pathplots<-paste0(getwd(),"/Ecography/plots/")
dina4.600.h<-3508*2#pixels at 600 dpi
dina4.600.w<-2480*2
singlecolumn<-c(945,945*2)#1st for 300 dpi, 2nd 600
singlehalfcolumn<-c(1476,1476*2)
doublecolumn<-c(1961,1961*2)
plotsaveformat<-c("tiff","pdf")

#A. DEER-------------
#Data published by:
#Chandler, R. B. et al. 2021. Modeling abundance, distribution, movement and space use with camera and telemetry data. - Ecology n/a: e03583.
pathfilesdeer<-paste0(getwd(),"/Ecography/files/deer/")
pathplotsdeer<-paste0(getwd(),"/Ecography/plots/deer/")
pathfilesdeerours<-paste0(pathfilesdeer,"Socias-Martinez etal/")
pathcomplete<-paste0(getwd(),"/Ecography/files/simulations/complete/") 

#1.True home range----
#preprocesing to avoid unnecessary files-----
# load(paste0(pathfilesdeer,"R/deer_scr_telem.RData"))
#explore the telemetry data
# telem.obj<-grep("telem",names(.GlobalEnv),value=TRUE)
# telem.obj<-telem.obj[telem.obj%!in%grep("cam|dtime",names(.GlobalEnv),value=T)]
# funst<-paste0(telem.obj,"%>% as_tibble %>% group_by(deerID) %>% dplyr::summarise(",telem.obj,"=n())")
# reduce(map(funst,~{x<-eval(parse(text=.x));x}),left_join)

# save(telem,cam.telem.sub,file=paste0(pathfilesdeerours,"deer_data.Rdata"))
#load(paste0(pathfilesdeerours,"deer_data.Rdata"))

#transform to movebank data type and then to telemetry object
# telem.ctmm<-telem %>% rename(individual.local.identifier=deerID,
#                  timestamp=datetime,
#                  location.long=Longitude,
#                  location.lat=Latitude) %>% ctmm::as.telemetry(.,keep=T)
# 
# telem.ctmm$gps_area_akde<-map_dbl(telem.ctmm,~{
#   GUESS <- ctmm.guess(.x,interactive=FALSE) # automated model guess
#   M.OUF <- ctmm.fit(.x,GUESS) # in general, use ctmm.select instead
#   AKDE <- akde(.x,M.OUF, weights = TRUE) # AKDE
#   
#   area<-summary(AKDE,units=F)$CI[,2]
#   area
# })
# dt<-data.frame(ids=names(telem.ctmm)[-length(names(telem.ctmm))],area=telem.ctmm$gps_area_akde,stringsAsFactors = F)
# 
# ggplot(data=dt)+
#   geom_violin(aes(x=1,y=area))+
#   geom_point(aes(x=1,y=area),alpha=0.5,size=2)
# 
# trial<-telem.ctmm$`103`
# plot(trial,UD=AKDE)
# dt.plot(trial)
# 
# #mean,median--------------
# #median or mean area
# #arround 6 km^2----
# dt$area %>% median
# m.hrarea<-dt$area %>% mean

#save(telem,cam.telem.sub,telem.ctmm,m.hrarea,file=paste0(pathout,"deer_data.Rdata"))
load(paste0(pathfilesdeerours,"deer_data.Rdata"))

#2.trap data-----------
#2.1 formating to use our functions-----
#sample with different ordering procedures, estimators and asymptotic models
#use functions already created, adapt if necessary
#use GAM to predict error

#cameras positions-----------
traps.locs<-cam.telem.sub$camera.locs
traps.locs.df<-tibble(individual.local.identifier=traps.locs %>% rownames,
                      UTME=traps.locs[,"UTME"] %>% unname,
                      UTMN=traps.locs[,"UTMN"] %>% unname) %>% arrange(individual.local.identifier)
longslats <- st_as_sf(traps.locs.df, coords = c("UTME", "UTMN"),crs = "epsg:32617")
longslats <- st_transform(longslats, crs= "+proj=longlat +lat_1=26.1947403815134 +lon_1=-81.2922921013087 +lat_2=26.2215124508171 +lon_2=-81.1817418675517 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")#cbind(telem.ctmm[[2]]$longitude,telem.ctmm[[2]]$latitude)
longslats <- longslats %>% mutate(location.long = unlist(map(longslats$geometry,1)),location.lat = unlist(map(longslats$geometry,2)))
traps.locs.df<-traps.locs.df %>% #select(individual.local.identifier) %>%  
  mutate(location.long=longslats$location.long, location.lat=longslats$location.lat,
         timestamp=paste0(sample(seq(as.Date('1999/01/01'), as.Date('2000/01/01'), by="day"), n())," ","12:00:00"))#add random time

#individual locations-----
ids.traps.t<-cam.telem.sub$camera.data
days<-1:dim(ids.traps.t)[3]
daysnames<-dimnames(ids.traps.t)[[3]]
ids.traps.t.df<-map2_dfr(days,daysnames,~{
  dta<-ids.traps.t[,,.x]
  dta_df<-as_tibble(dta) %>% mutate(timestamp=.y) %>% rename(individual.local.identifier=deerID)
  dta_df
})
#counted as 1 if detected (even if several times) within 1H
ids.traps.t.df<-ids.traps.t.df %>% filter(n>=1) %>% select(-n)

#add the cameras longitud and latitude
ids.traps.t.df<-ids.traps.t.df %>% 
  mutate(location.long=traps.locs.df$location.long[match(cameraID,traps.locs.df$individual.local.identifier)],
         location.lat=traps.locs.df$location.lat[match(cameraID,traps.locs.df$individual.local.identifier)]) %>% 
  select(-cameraID)
ids.traps.t.df

#convert to ctmm and obtain x y projection-------
#combine camera positions and individual locations so all are projected in the same manner by ctmm
ids.traps.locs<-bind_rows(ids.traps.t.df %>% mutate(type="id"),
                          traps.locs.df %>% mutate(type="trap"))
ids.traps.locs.ctmm<-ids.traps.locs %>% ctmm::as.telemetry(.,keep=c("individual.local.identifier","type"))
ids.traps.locs.xyt<-map_dfr(ids.traps.locs.ctmm,~{.x@.Data %>% set_names(.x@names)}) %>% select(type,individual.local.identifier,x,y,t,longitude,latitude,timestamp)
ids.traps.locs.xyt[,c("x","y","t")]<-lapply(ids.traps.locs.xyt[,c("x","y","t")],as.numeric)

#extract x y projections for ids locations-----
individuals<-ids.traps.t.df$individual.local.identifier %>% unique %>% sort
ids.traps.ctmm<-ids.traps.locs.ctmm[individuals]
ids.traps.xyt<-ids.traps.locs.xyt %>% filter(type=="id")

#extract x y projections for camera trap locations-----
cameras<-traps.locs.df$cameraID %>% unique %>% sort
traps.locs.ctmm<-ids.traps.locs.ctmm[cameras]
traps.locs.xyt<-ids.traps.locs.xyt %>% filter(type=="trap")


#visualize gps vs traps positions--------
ids.gps<-telem %>% rename(individual.local.identifier=deerID,longitude=Longitude,latitude=Latitude)
both<-names(ids.traps.ctmm)[which(names(ids.traps.ctmm) %in% ids.gps$individual.local.identifier==T)]
#red dots are traps where individuals with collars were photographed
ggplot()+
  geom_point(data=ids.gps,aes(x=longitude,y=latitude,color=individual.local.identifier %>% as.factor),alpha=0.5)+
  geom_point(data=traps.locs.xyt,aes(x=longitude,y=latitude),color="black") +
  geom_point(data=ids.traps.xyt[ids.traps.xyt$individual.local.identifier%in%both,],
             aes(x=longitude,y=latitude,shape=individual.local.identifier),color="red")+
  theme_classic()

#2.2 Estimate trap density----------
#for each trap find the CLOSEST TRAP
#do so for the x and y dimensions
#what is the mean distance of the closest trap----
m.dists<-calc_meandist_xy(traps.locs.xyt)
m.dist<-(m.dists$dist.x+m.dists$dist.y)/2

#what is the radius of the mean home range gps-----
r<-sqrt(m.hrarea/pi)
#derive trap density----
trap.d<-r/m.dist
#27 1st method
#4 2nd


#2.3 calculate home range area from trapping data-----
ids.traps.tofit<-tibble(id=1:length(ids.traps.ctmm),individual.local.identifier=names(ids.traps.ctmm),data.ctmm=ids.traps.ctmm)

#individuals to do from 1:600 #################
#I will do the first 20 myself until we get some way of computing the rest
nbids<-max(ids.traps.tofit$id)
idstodo<-c(1:nbids)#[c(1:600)%!in%c(1:4,20,310)]#600#5:6#1:4#:8#
#traps.locs.xyt2 %>% View

capture.grid<-map(idstodo,~traps.locs.xyt)
#do any loc in trial correspond to any loc in grid
#trial$data.ctmm$fake$loc %in% paste0(capture.grid[[1]]$x," ",capture.grid[[1]]$y)
#estimators----------
#p.intervals sets the points on which to conduct the estimations
p.intervals<-map(ids.traps.ctmm,~1:length(.x@.Data[[1]])) %>% unname

nbcores<-4
# st<-Sys.time()#takes a while
# estimators_df<-fit_estimators(allids.equal=F,
#                               ids.traps.tofit,
#                               p.intervals=p.intervals,
#                               capture.grid=capture.grid,
#                               estimators=estimators,
#                               akima.res=1000,
#                               regular=F,#the capture grid is not regular
#                               extrapol=F,#should akima::interp() extrapolate outside point triangulation? if not then all points outside=0
#                               varslocfe=c("x","y"),
#                               sampr=sampr,
#                               pathwaypartial=pathfilesdeerours,pathwaycomplete=pathfilesdeerours,
#                               save=save,paralel=paralel,nbcores=nbcores,chunk=chunk,mac.linux)
# print(Sys.time()-st)
#if already fitted
filenames<-list.files(pathfilesdeerours,pattern="estimators_df_all")
load(paste0(pathfilesdeerours,filenames[[1]]))


#asympmodels
results<-fit_asympmodels(estimators_df,
                         estimators,
                         asympmodels[asympmodels!="raw"],recursiv.model=T,iv=NULL,
                         pathwaypartial=pathfilesdeerours,pathwaycomplete=pathfilesdeerours,save=save,paralel=paralel,nbcores=nbcores,chunk=chunk,mac.linux)
future:::ClusterRegistry("stop")

#long format
results<-convert_long(results %>% select(-idx),estimators,true.area=m.hrarea,pathwaycomplete=pathfilesdeerours,save=save)
save(results,file=paste0(pathfilesdeerours,"results_all_deer_long_",
                         min(results$id %>% unique),"_",
                         max(results$id %>% unique),".Rdata"))

#2.4 compare with simulations through GAM-----
#load GAM
filenames<-list.files(pathcomplete,pattern="gam_gauss_ti")
load(paste0(pathcomplete,filenames[[1]]))

#predict using this imput as the newdata for predict with gam----
results<-results %>% mutate(trap.d=trap.d,
                            approach=interaction(sampling,estimator,model),
                            lg10nbobs=log10(nbobs)) %>% select(-models_diag)
#prepare for predictions
area.cs<-c(1,.7,.5,.2)

# topred<-map(area.cs,~results %>% mutate(area.c=.x)) %>% set_names(paste0("deer_area",area.cs))
# preds.deer<-pred_gam(gmodel,topred,specials=c("[()]","log10"),#log10 is for the case when log10(nbobs) was done in the model formula
#                   varint=list(c("approach","lg10nbobs","area.c","trap.d")),
#                   nms=names(topred))
#save(preds.deer,file=paste0(pathfilesdeerours,"preds_deer",".Rdata"))
load(file=paste0(pathfilesdeerours,"preds_deer",".Rdata"))
dtat_deer<-map_dfr(preds.deer,"pred")%>% 
  mutate(cip=(upperfitili-lowerfitili)/fitili,
         model=factor(model,levels=asympmodels,labels=asympmodels.lbs),
         estimator=factor(estimator,levels=estimators,labels=estimators.lbs),
         sampling=factor(sampling,levels=sampregs,labels=sampregs.lbs),
         nbobs=10^lg10nbobs %>% round)

#B. JAGUAR-----
pathfilesjaguar<-paste0(getwd(),"/Ecography/files/jaguar/")
pathplotsjaguar<-paste0(getwd(),"/Ecography/plots/jaguar/")
pathfilesjaguarours<-paste0(pathfilesjaguar,"Socias-Martinez etal/")
pathcomplete<-paste0(getwd(),"/Ecography/files/simulations/complete/") 
#data from
#Harmsen, B. J. et al. 2010. Chapter 18 - The ecology of jaguars in the cockscomb basin wildlife sanctuary, belize. - 
#In: The biology and conservation of wild felids. Oxford University Press Oxford, United Kingdom, pp. 403–416.
jag.traps.locs<-read_csv(paste0(pathfilesjaguar,"JagsTraps.csv"))
jag.ids.traps<-read_csv(paste0(pathfilesjaguar,"JagsCaptures.csv"))

#mean of 12 trappings per id and in 3 different traps
samplingsize<-jag.ids.traps %>% group_by(ID) %>% 
  summarise(n=n(),#nbobs per id
            ntraps=Station %>% unique %>% length) %>% ungroup %>% summarise(mn=mean(n),mntraps=mean(ntraps))#nbobs in diff traps

#obtain cameras positions in longlat-----------
jag.ids.traps<- jag.ids.traps %>% rename(UTME=X,UTMN=Y)%>% arrange(UTME,UTMN) %>% mutate(timestamp=paste0(as.Date(Date,format="%d-%b-%y")," ",Time)) %>% 
  rename(individual.local.identifier=ID) %>%select(-Date,-Time)
jag.traps.locs<-jag.traps.locs %>% rename(UTME=X,UTMN=Y)%>% arrange(UTME,UTMN) %>% 
  mutate(timestamp=jag.ids.traps$timestamp[1]) %>% rename(individual.local.identifier=`Station ID`) %>% select(individual.local.identifier,UTME,UTMN,timestamp)
jags.ids.traps.locs<-bind_rows(jag.ids.traps,jag.traps.locs) 
longslats <- st_as_sf(jags.ids.traps.locs, coords = c("UTME", "UTMN"),crs ="+proj=utm +zone=16N +datum=WGS84")#"EPSG:28992"
longslats <- st_transform(longslats, crs= "+proj=longlat")#cbind(telem.ctmm[[2]]$longitude,telem.ctmm[[2]]$latitude) "+init=epsg:2419"
longslats <- longslats %>% mutate(location.long = unlist(map(longslats$geometry,1)),location.lat = unlist(map(longslats$geometry,2)))
jags.ids.traps.locs<-jags.ids.traps.locs %>% mutate(location.long=longslats$location.long, location.lat=longslats$location.lat)

#visualize grid
ggplot()+
  geom_point(data=jags.ids.traps.locs %>% filter(is.na(Sex)==T),aes(x=location.long,y=location.lat),color="black")


#convert to ctmm
jags.ids.traps.locs.ctmm<-jags.ids.traps.locs %>%  select(individual.local.identifier,timestamp,location.long,location.lat) %>% ctmm::as.telemetry(.)
ids.traps.idx<-jags.ids.traps.locs$individual.local.identifier  %>% unique %>% sort #need to arrange like this because ctmm::as.telemetry() arranges by name
jg<-jags.ids.traps.locs %>% group_by(individual.local.identifier) %>% summarise(sex=Sex %>% unique,
                                                                                idx=which(ids.traps.idx==individual.local.identifier),
                                                                                .groups="keep") %>% ungroup %>% arrange(idx)
jags.traps.locs.ctmm<-jags.ids.traps.locs.ctmm[which(is.na(jg$sex)==T)]
jags.ids.traps.ctmm<-jags.ids.traps.locs.ctmm[which(is.na(jg$sex)==F)]

#export xyt
jags.traps.locs.xyt<-map2_dfr(jags.traps.locs.ctmm,names(jags.traps.locs.ctmm),~.x@.Data  %>%  set_names(.x@names) %>%bind_cols%>% mutate(individual.local.identifier=.y))
jags.ids.traps.xyt<-map2_dfr(jags.ids.traps.ctmm,names(jags.ids.traps.ctmm),~.x@.Data  %>%  set_names(.x@names) %>% bind_cols %>% mutate(individual.local.identifier=.y))

#fit our approaches
ids.traps.tofit<-tibble(id=1:length(jags.ids.traps.ctmm),individual.local.identifier=names(jags.ids.traps.ctmm),data.ctmm=jags.ids.traps.ctmm)

#intervals to do the recursive analyses
p.intervals<-map(jags.ids.traps.ctmm,~1:length(.x@.Data[[1]])) %>% unname
#grid
capture.grid<-map(idstodo,~jags.traps.locs.xyt)

nbids<-max(ids.traps.tofit$id);idstodo<-c(1:nbids)


nbcores<-4
# st<-Sys.time()
# estimators_df<-fit_estimators(allids.equal=F,
#                               ids.traps.tofit,
#                               p.intervals=p.intervals,
#                               capture.grid=capture.grid,
#                               estimators=estimators,
#                               akima.res=1000,
#                               regular=F,#the capture grid is not regular
#                               extrapol=F,#should akima::interp() extrapolate outside point triangulation? if not then all points outside=0
#                               varslocfe=c("x","y"),
#                               sampr=sampr,
#                               pathwaypartial=pathfilesjaguarours,pathwaycomplete=pathfilesjaguarours,
#                               save=save,paralel=paralel,nbcores=nbcores,chunk=chunk,mac.linux)
# 
# print(Sys.time()-st)

#if already fitted
filenames<-list.files(pathfilesjaguarours,pattern="estimators_df_all")
load(paste0(pathfilesjaguarours,filenames[[1]]))

#asympmodels
results<-fit_asympmodels(estimators_df,
                         estimators,
                         asympmodels[asympmodels!="raw"],recursiv.model=T,iv=NULL,
                         pathwaypartial=pathfilesjaguarours,pathwaycomplete=pathfilesjaguarours,
                         save=save,paralel=paralel,nbcores=nbcores,chunk=chunk,mac.linux)

results<-results %>% mutate(sex=jags.ids.traps.locs$Sex[match(individual.local.identifier,jags.ids.traps.locs$individual.local.identifier)])
results.fm<-results %>% filter(sex!="U") %>% group_by(sex,.drop=T) %>% group_split
m.hrarea.fm<-c(10500000,33400000)# in m2 from Rabinowitz and Nottingam (1986)
#here need to take into account sex on home range size before continuing

results<-map2_dfr(results.fm,m.hrarea.fm,~convert_long(.x %>% select(-idx),estimators,true.area=.y,pathwaycomplete=pathfilesjaguarours,save=T))
save(results,file=paste0(pathfilesjaguarours,"results_all_jaguar_long_",
                         min(results$id %>% unique),"_",
                         max(results$id %>% unique),".Rdata"))

#2.trap.d------
m.dists<-calc_meandist_xy(jags.traps.locs.xyt)
m.dist<-(m.dists$dist.x+m.dists$dist.y)/2

#what is the radius of the mean home range gps-----
r.fm<-(m.hrarea.fm/pi) %>% sqrt

#derive trap density----
trap.d.fm<-r.fm/m.dist

#predict using this imput as the newdata for predict with gam,
results$trap.d<-NA
results$trap.d[which(results$sex=="F")]<-trap.d.fm[1]
results$trap.d[which(results$sex=="M")]<-trap.d.fm[2]
results<-results %>% mutate(approach=interaction(sampling,estimator,model),
                            lg10nbobs=log10(nbobs)) %>% select(-models_diag)
#prepare for predictions

#prepare for predictions
area.cs<-c(1,.7,.5,.2)
topred<-map(area.cs,~results %>% mutate(area.c=.x)) %>% set_names(paste0("jaguar_area",area.cs))
#predict----------
# preds.jaguar<-pred_gam(gmodel,topred,specials=c("[()]","log10"),#log10 is for the case when log10(nbobs) was done in the model formula
#                      varint=list(c("approach","lg10nbobs","area.c","trap.d")),
#                      nms=names(topred))
# save(preds.jaguar,file=paste0(pathfilesjaguarours,"preds_jaguar.Rdata"))
load(file=paste0(pathfilesjaguarours,"preds_jaguar.Rdata"))

dtat_jaguar<-map_dfr(preds.jaguar,"pred")%>% 
  mutate(cip=(upperfitili-lowerfitili)/fitili,
         model=factor(model,levels=asympmodels,labels=asympmodels.lbs),
         estimator=factor(estimator,levels=estimators,labels=estimators.lbs),
         sampling=factor(sampling,levels=sampregs,labels=sampregs.lbs),
         nbobs=10^lg10nbobs %>% round) 

dtattrim_jaguar<-dtat_jaguar
dtattrim_jaguar$area[which(dtat_jaguar$area>5)]<-NA
#max(dtat$area,na.rm=T)

yvar<-"fit"
vars<-c("nbobs","trap.d","area.c")
vars.fix<-c("Number observations","Trap density","Area covered")
xlbs<-c("Number of observations","Trap density","Area covered")
varvary<-"nbobs"

mrgins.plot<-c(.1,.1,.3,.3)
mrgin.axis.title.x<-c(.5, 0, 0, 0)
mrgin.axis.title.y<-c(0, .3, 0, 0)
mrgin.legend<-c(0,0,.5, 0)
mrgin.box<-c(0, 0, 0, 0)
tag.position<-c(0.85, 0.1)

thlemon.fmm<-c(
  "theme(
  plot.margin = unit(mrgins.plot, 'cm'),
      axis.title.x = element_text(size = 8, margin=margin(mrgin.axis.title.x,unit='cm')),
      axis.title.y = element_text(size = 8, margin=margin(mrgin.axis.title.y,unit='cm')),
      axis.text.x=element_blank(),
      axis.text.y=element_blank(),

      axis.line = element_line(colour = 'black', size = 0.2),
      axis.ticks = element_line(colour = 'black', size = 0.2),
      axis.ticks.length=unit(0.5, 'mm'),
      
      legend.title=element_text(size = 4),#,angle = 90, vjust = 0.5),
      legend.text=element_text(size=4),#,angle = 90, vjust = 0.5),
      legend.key.height= unit(0.1, 'cm'),
      legend.key.width= unit(0.1, 'cm'),
      legend.position='bottom',
      legend.justification='right',
      legend.margin=margin(mrgin.legend,unit='cm'),
      legend.box.margin=margin(mrgin.box,unit='cm'),

      strip.text = element_text(size = 4),
      
      panel.spacing.x=unit(0.1, 'lines'),
      panel.spacing.y=unit(0.1, 'lines'),
      plot.tag = element_text(size = 4),
      plot.tag.position = tag.position,
      plot.background = element_rect(colour ='black'))")

thlemon2.fmm<-paste0("theme(
axis.text.x = element_text(size=5,angle = 65, vjust = 0.5,margin = margin(t = 2)),
axis.text.y = element_text(size=6,margin = margin(r = 2)))")

n<-8
mycolors<-viridis(n)[seq(n,1,length.out = 4)]
pls2d_jagu_deer<-map(1,~{
  dtats<-bind_rows(dtat_deer %>% mutate(sp="deer",sex="M"),
                   dtat_jaguar %>% mutate(sp="jagu")) %>%
    filter(sampling!="Distance-ordered") %>% mutate(sex=factor(sex,levels=c("M","F")),
                                                    sp=factor(sp,levels=c("deer","jagu"),labels=c("Deer","Jaguar")))
  dtattrim<-bind_rows(dtat_deer %>% mutate(sp="deer",sex="M"),
                      dtattrim_jaguar%>% mutate(sp="jagu"))%>%
    filter(sampling!="Distance-ordered") %>% mutate(sex=factor(sex,levels=c("M","F")),
                                                    sp=factor(sp,levels=c("deer","jagu"),labels=c("Deer","Jaguar")))
  
  dtams<-dtattrim %>% group_by(sp,sex,model,estimator,sampling,area.c,nbobs) %>% summarise(marea=mean(area,na.rm=T),.groups="keep") %>% ungroup
  
  gg1<-paste0("geom_line(data=dtats,aes(x =",varvary,", y = ",yvar,",color=factor(area.c,levels=c('1','0.7','0.5','0.2'))))")
  gg2<-"scale_x_log10(breaks=c(2,5,10,20,40,80,160,320,600,1000,2000))"
  gg3<-paste0("geom_ribbon(data=dtats,aes(x =",varvary,", y = ",yvar,",ymin=lowerfit, ymax=upperfit,fill=factor(area.c,levels=c('1','0.7','0.5','0.2'))), alpha=0.4)")
  gg1sim<-paste0("geom_point(data=dtats,aes(x =",varvary,", y = area","),size=1.8,shape=21,fill='black',color='transparent',alpha=0.3,stroke=0,position = position_jitter(width=.05,seed = 062020))")
  gg12sim<-paste0("geom_point(data=dtams,aes(x =",varvary,", y = marea","),size=0.75,shape=23,fill='red',color='red',stroke=0,alpha=0.3)")
  
  trapds<-dtattrim %>% group_by(sp,sex) %>% summarise(trap.d=unique(trap.d,na.rm=T),.groups="keep") %>% pull(trap.d)
  true.areas<-c(m.hrarea,m.hrarea.fm[2],m.hrarea.fm[1])/1000000#get Km2
  
  xlb<-xlbs[1];fixv1<-vars[2];fixv2<-vars[3];fixval1<-dtats[[vars[2]]] %>% unique;fixval2<-dtats[[vars[3]]] %>% unique
  fixlb<-paste0(xlbs[2]," = ",dtats[[vars[2]]] %>% unique %>% round(digits=1))
  
  pl2d<-ggplot()+
    geom_hline(yintercept=1,color="black",size=.2,alpha=.1)+#true area
    eval(parse(text=gg1))+
    eval(parse(text=gg2))+
    eval(parse(text=gg3))+
    geom_line(data=dtats,aes(x =nbobs, y = area,group=individual.local.identifier),alpha=0.4,linetype="dotted")+
    eval(parse(text=gg1sim))+
    eval(parse(text=gg12sim))+
    scale_color_manual(values = mycolors,labels=c(paste0(c(100,70,50,20),"%")))+
    scale_fill_manual(values = mycolors,labels=c(paste0(c(100,70,50,20),"%")))+
    facet_rep_grid(estimator~sp+sex+model+sampling)+
    #labeller = labeller(sp = function(x) {rep("", length(x))},sex = function(x) {rep("", length(x))}))+
    coord_cartesian(ylim=c(0,3))+
    theme_classic()+
    labs(fill="Area covered",
         color="Area covered",
         y = "Estimated area / True home range area",
         x=xlb,
         #tag = fixlb,
         size= 4)
  pl2d<-lemonworkarround(pl2d,theme1=thlemon.fmm,theme2=thlemon2.fmm,nr=5,ncl=3*3)
  library(grid)
  # Create tags
  tags<-paste0("True area = ",true.areas %>% round(digits=1)," Km2\n","Trap density = ",trapds%>% round(digits=1))
  posx<-c(0.08,0.38,0.70)
  grobtags <- map2(tags,posx,~grobTree(textGrob(.x, x=.y,  y=0.045, hjust=0,gp=gpar(col="black", fontsize=4))))
  pl2d<-pl2d + annotation_custom(grobtags[[1]]) + annotation_custom(grobtags[[2]]) + annotation_custom(grobtags[[3]])
  
  pl2d
  
})

ggsave(file=paste0(pathplots,"figure_6",".","pdf"),
       device="pdf",
       plot=pls2d_jagu_deer[[1]],
       height=(doublecolumn[2]/12)*9,width=doublecolumn[2],
       units="px",dpi=600)

ggsave(file=paste0(pathplots,"figure_6",".","tiff"),
       device = grDevices::tiff,
       plot=pls2d_jagu_deer[[1]],
       height=(doublecolumn[2]/12)*9,width=doublecolumn[2],
       units="px",dpi=600)

#save plots----
#for some reason R session aborts when ggsave with tiff, appears to be a bug
# walk(plotsaveformat,~{
#   #print(.x)
#   ggsave(file=paste0(pathplots,"figure_7",".",.x),
#          device=.x,
#          plot=pls2d_jagu_deer[[1]],
#          height=(doublecolumn[2]/12)*9,width=doublecolumn[2],
#          units="px",dpi=600)
# })

#Figure SI-------
#texts for figures--------
#load figure captions etc & convert to grob for plotting with graphs and save as pdf
filenames<-list.files(pathplots,pattern="texts SI")
texts.si<-readLines(paste0(pathplots,filenames[[1]]), encoding = "UTF-8")
texts.si<-gsub(texts.si,pattern="\r\n",replacement="")
#texts.si<-strsplit(texts.si,split="newnewnew")[[1]]
texts.grob<-map(texts.si,~{
  tit.bod<-strsplit(.x,split="]] ")[[1]]; 
  tit<-ggparagraph(text = tit.bod[1], face = "bold", size = 11, color = "black")  
  bod<-ggparagraph(text = tit.bod[2], face = "plain", size = 8, color = "black")
  tit.bod.grob<-ggplotGrob(ggarrange(NULL,NULL,NULL,
                                     NULL,tit,NULL,
                                     NULL,bod,NULL,
                                     NULL,NULL,NULL,nrow = 4, ncol =3,
                                     heights = c(0.01,0.05,0.15,0.01),
                                     widths = c(0.01,0.15,0.01)))
  #list(tit.bod.grob,tit.bod.grob)
  tit.bod.grob
}) #%>% unlist(recursive=F)


pls2d_jagu_deer_complete<-map(1,~{
  dtats<-bind_rows(dtat_deer %>% mutate(sp="deer",sex="M"),
                   dtat_jaguar %>% mutate(sp="jagu")) %>% 
    mutate(sex=factor(sex,levels=c("M","F")),
           sp=factor(sp,levels=c("deer","jagu"),labels=c("Deer","Jaguar")))
  dtattrim<-bind_rows(dtat_deer %>% mutate(sp="deer",sex="M"),
                      dtattrim_jaguar%>% mutate(sp="jagu"))%>% 
    mutate(sex=factor(sex,levels=c("M","F")),
           sp=factor(sp,levels=c("deer","jagu"),labels=c("Deer","Jaguar")))
  
  dtams<-dtattrim %>% group_by(sp,sex,model,estimator,sampling,area.c,nbobs) %>% summarise(marea=mean(area,na.rm=T),.groups="keep") %>% ungroup
  
  gg1<-paste0("geom_line(data=dtats,aes(x =",varvary,", y = ",yvar,",color=factor(area.c,levels=c('1','0.7','0.5','0.2'))))")
  gg2<-"scale_x_log10(breaks=c(2,5,10,20,40,80,160,320,600,1000,2000))"
  gg3<-paste0("geom_ribbon(data=dtats,aes(x =",varvary,", y = ",yvar,",ymin=lowerfit, ymax=upperfit,fill=factor(area.c,levels=c('1','0.7','0.5','0.2'))), alpha=0.4)")
  gg1sim<-paste0("geom_point(data=dtats,aes(x =",varvary,", y = area","),size=1.8,shape=21,fill='black',color='transparent',alpha=0.3,stroke=0,position = position_jitter(width=.05,seed = 062020))")
  gg12sim<-paste0("geom_point(data=dtams,aes(x =",varvary,", y = marea","),size=0.75,shape=23,fill='red',color='red',stroke=0,alpha=0.3)")
  
  trapds<-dtattrim %>% group_by(sp,sex) %>% summarise(trap.d=unique(trap.d,na.rm=T),.groups="keep") %>% pull(trap.d)
  true.areas<-c(m.hrarea,m.hrarea.fm[2],m.hrarea.fm[1])/1000000#get Km2
  
  xlb<-xlbs[1];fixv1<-vars[2];fixv2<-vars[3];fixval1<-dtats[[vars[2]]] %>% unique;fixval2<-dtats[[vars[3]]] %>% unique
  fixlb<-paste0(xlbs[2]," = ",dtats[[vars[2]]] %>% unique %>% round(digits=1))
  
  pl2d<-ggplot()+
    geom_hline(yintercept=1,color="black",size=.2,alpha=.1)+#true area
    eval(parse(text=gg1))+
    eval(parse(text=gg2))+
    eval(parse(text=gg3))+
    geom_line(data=dtats,aes(x =nbobs, y = area,group=individual.local.identifier),alpha=0.4,linetype="dotted")+
    eval(parse(text=gg1sim))+
    eval(parse(text=gg12sim))+
    scale_color_manual(values = mycolors,labels=c(paste0(c(100,70,50,20),"%")))+
    scale_fill_manual(values = mycolors,labels=c(paste0(c(100,70,50,20),"%")))+
    facet_rep_grid(estimator~sp+sex+model+sampling)+
    #labeller = labeller(sp = function(x) {rep("", length(x))},sex = function(x) {rep("", length(x))}))+
    coord_cartesian(ylim=c(0,3))+
    theme_classic()+
    labs(fill="Area covered",
         color="Area covered",
         y = "Estimated area / True home range area",
         x=xlb,
         #tag = fixlb,
         size= 4)
  pl2d<-lemonworkarround(pl2d,theme1=thlemon.fmm,theme2=thlemon2.fmm,nr=5,ncl=3*5)
  library(grid)
  # Create tags
  tags<-paste0("True area = ",true.areas %>% round(digits=1)," Km2\n","Trap density = ",trapds%>% round(digits=1))
  posx<-c(0.08,0.38,0.70)
  grobtags <- map2(tags,posx,~grobTree(textGrob(.x, x=.y,  y=0.045, hjust=0,gp=gpar(col="black", fontsize=4))))
  pl2d<-pl2d + annotation_custom(grobtags[[1]]) + annotation_custom(grobtags[[2]]) + annotation_custom(grobtags[[3]])
  
  pl2d
  
})

#intertwin texts and graphs
pls.si.8<-pmap(list(pl=pls2d_jagu_deer_complete,txt=texts.grob[7]),
               function(pl,txt){list(pl,txt)}) %>% unlist(.,recursive=FALSE)


#save plots----
ggsave(plot=marrangeGrob(pls.si.8, nrow = 2, ncol = 1,heights=c(.6,.4),widths=c(1),top=NULL),
       file=paste0(pathplots,"SI/","si_figure_8",".","pdf"),
       device="pdf",
       height=dina4.600.h,width=dina4.600.w,
       units="px",dpi=600)
ggsave(plot=marrangeGrob(pls.si.8, nrow = 2, ncol = 1,heights=c(.6,.4),widths=c(1),top=NULL),
       file=paste0(pathplots,"SI/","si_figure_8",".","tiff"),
       device=grDevices::tiff,
       height=dina4.600.h,width=dina4.600.w,
       units="px",dpi=600)

# walk(plotsaveformat,~{
#   ggsave(plot=marrangeGrob(pls.si.7, nrow = 2, ncol = 1,heights=c(.6,.4),widths=c(1),top=NULL),
#          file=paste0(pathplots,"SI/","si_figure_6",".",.x),
#          device=.x,
#          height=dina4.600.h,width=dina4.600.w,
#          units="px",dpi=600)
# })

#C. Map to add to figure 6---------
#data from Sentinel 2 Scihub Copernicus, manually downloaded for the region of interest 
#see package sen2r for automatization
#parts of the code inspired from the following ressources
#https://gis.stackexchange.com/questions/211157/open-jpeg2000-sentinel-2-in-r
# see Miguel Castro Gómez tutorial : https://www.youtube.com/watch?v=1dAhrc-kw8o
#some material available at https://rus-copernicus.eu/portal/the-rus-library/train-with-rus/
pkgs<-c("tidyverse","rgdal","raster","leaflet","rasterVis","gridExtra","RColorBrewer","plotly","sf","reshape2","RStoolbox","caret","viridis","ggthemes")
new_pkgs<-pkgs[!pkgs %in% installed.packages()[,"Package"]]
install.packages(new_pkgs)
sapply(pkgs,require, character.only=T)

#plot the satellite data
pathmaps<-c(paste0(getwd(),"/Ecography/files/copernicus_florida/"),
            paste0(getwd(),"/Ecography/files/copernicus_belize/"))
s2_stacks<-map(pathmaps,~{
  s2<-list.files(.x, recursive = T, full.names = T,pattern="B0[2348]_10m.jp2$")
  s2<-lapply(s2,function(x)raster(x))
  #s2[1] %>% print
  #T17RMJ Area
  #options(repr.plot.width=41,repr.plot.height=20)
  m<-rbind(c(1,2));layout(m)
  s2_stack<-stack(s2)
  
  plotRGB(s2_stack,r=3,g=2,b=1,scale=maxValue(s2[[2]]),stretch='hist')
  plotRGB(s2_stack,r=4,g=3,b=2,scale=maxValue(s2[[2]]),stretch='hist')
  return(s2_stack)
}) %>% set_names(c("deer","jaguar"))

#crop the area of interest----------
#define the study area with the trapping grids
traps.locs2<-cam.telem.sub$camera.locs
deer.traps.locs2<-tibble(individual.local.identifier=traps.locs2 %>% rownames,
                         UTME=traps.locs2[,"UTME"] %>% unname,
                         UTMN=traps.locs2[,"UTMN"] %>% unname) %>% 
  arrange(individual.local.identifier) 
jags.traps.locs2<-jag.traps.locs%>% select(-timestamp)
traps.locs2<-list(deer.traps.locs2,jags.traps.locs2)

cls<-c("UTME","UTMN")
buf<-c(2000,5000)
rngs<-map(traps.locs2,~{lapply(cls,function(cl)range(.x[[cl]])) %>% set_names(cls)}) %>% set_names(c("deer","jaguar"))
rngs<-rngs %>% melt %>% rename(coord=L2,dataset=L1) %>% 
  mutate(var=rep(c("min","max"),nrow(.)/2),
         valuebuf=case_when(var=="min"&dataset=="deer"~value-buf[1],
                            var=="max"&dataset=="deer"~value+buf[1],
                            var=="min"&dataset=="jaguar"~value-buf[2],
                            var=="max"&dataset=="jaguar"~value+buf[2]))
#the four points
frpts<-expand_grid(UTME=c("min","max"),UTMN=c("min","max"),.name_repair = "minimal")
frpts<-rngs %>% group_by(dataset) %>% group_split %>% map(.,~{
  frpts %>% rowwise %>%  
    mutate(val.UTME=rngs$valuebuf[which(rngs$coord=="UTME" & rngs$var==UTME & rngs$dataset==.x$dataset %>% unique)],
           val.UTMN=rngs$valuebuf[which(rngs$coord=="UTMN" & rngs$var==UTMN & rngs$dataset==.x$dataset %>% unique)])
})

polys<-map2(frpts,s2_stacks,~{
  pol = st_polygon(list(cbind(.x$val.UTME[c(1,3,3,1,1)],.x$val.UTMN[c(1,1,2,2,1)])))
  polc = st_sfc(pol, crs=.y@crs %>% as.character)
  # Generate empty raster layer and rasterize points
  rpolc<- raster(crs = .y@crs %>% as.character, 
                 vals = 0, 
                 ext = extent(c(.x$val.UTME[c(1)], .x$val.UTME[c(3)],
                                .x$val.UTMN[c(1)], .x$val.UTMN[c(2)])),
                 resolution=c(10,10)) %>%
    rasterize(as(polc,"Spatial"), .)
}
)

#crop
s2_stacks_crop<-map2(s2_stacks,polys,~{raster::crop(.x, .y)})
#visualise crop
walk2(s2_stacks_crop,pathmaps,~{
  s2<-list.files(.y, recursive = T, full.names = T,pattern="B0[2348]_10m.jp2$")
  s2<-lapply(s2,function(x)raster(x))
  plotRGB(.x,r=3,g=2,b=1,scale=maxValue(s2[[2]]),stretch='hist')
})


#calculate polygons for home ranges and overlay them to the plots
ids.traps.t<-cam.telem.sub$camera.data
days<-1:dim(ids.traps.t)[3];daysnames<-dimnames(ids.traps.t)[[3]]
ids.traps.t.df<-map2_dfr(days,daysnames,~{
  dta<-ids.traps.t[,,.x]
  dta_df<-as_tibble(dta) %>% mutate(timestamp=.y) %>% rename(individual.local.identifier=deerID)
  dta_df
})
#counted as 1 if detected (even if several times) within 1H
ids.traps.t.df<-ids.traps.t.df %>% filter(n>=1) %>% select(-n)
#add the cameras longitud and latitude
ids.traps.t.df<-ids.traps.t.df %>% 
  mutate(Sex="M",
         UTME=traps.locs.df$UTME[match(cameraID,traps.locs.df$individual.local.identifier)],
         UTMN=traps.locs.df$UTMN[match(cameraID,traps.locs.df$individual.local.identifier)]) %>% 
  select(-cameraID)

#calculate KDE to plot
ids.traps2<-list(ids.traps.t.df,jag.ids.traps)

#calculate the KDE home ranges and extract polygons
hr_amt_kde<-map2(ids.traps2,s2_stacks,~{
  #.x<-ids.traps2[[2]];.y<-s2_stacks[[2]]
  #.x$Sex %>% unique
  sbs<-.x %>% filter(Sex!="U") %>% group_by(Sex) %>% group_split %>%  
    map_dfr(.,~{.x<-.x %>% group_by(Sex,individual.local.identifier) %>% 
      summarise(n=n(),.groups="keep") %>% arrange(-n)
    .x[1:3,]})
  res<-.x  %>% filter(individual.local.identifier%in%sbs$individual.local.identifier) %>% 
    mutate(timestamp=as.POSIXct(timestamp,tryFormats = c("%Y-%m-%d %H:%M:%OS"))) %>% 
    make_track(.x=UTME,.y=UTMN,.t=timestamp,
               id = individual.local.identifier,
               sex=Sex,
               crs=st2@crs %>% as.character) %>% nest_by(sex,id) %>% ungroup
  res<-res %>% mutate(hr=map(data,function(idx){
    bandw<- tryCatch({hr_kde_ref_scaled(idx,num.of.parts=1,levels=0.95)},error=function(e){NA})
    if(is.na(bandw)==F){ 
      amt_kde<-tryCatch({
        amt::hr_kde(idx,levels=0.95,h=bandw,rand_buffer=1e-05,keep.data=F,units=T) %>% list}, error=function(e){NA %>% list})
    }else{amt_kde<-NA %>% list}
    amt_kde
  }))
  
})

hr_amt_kde<-hr_amt_kde %>%map(.,~.x %>% mutate(isop=map(hr,function(h)if(length(h[[1]])>1){hr_isopleths(h[[1]])}else{NA})))
#extract sf objects for plotting
hr_amt_kde_pl<-map(hr_amt_kde,~{.x$isop  %>% bind_rows %>% mutate(ids=.x$id,sex=.x$sex)})

#plot all together
n<-8
nn<-3
mycolors<-viridis(n)[seq(n,1,length.out = nn)]
#divide by spexies and sex
maps<-pmap(list(.x=list(s2_stacks_crop[[1]],s2_stacks_crop[[2]],s2_stacks_crop[[2]]),
                .y=list(traps.locs2[[1]],traps.locs2[[2]],traps.locs2[[2]]),
                .w=list(hr_amt_kde_pl[[1]],hr_amt_kde_pl[[2]] %>% filter(sex=="M"),hr_amt_kde_pl[[2]] %>% filter(sex=="F")),
                .z=list(pathmaps[[1]],pathmaps[[2]],pathmaps[[2]])),function(.x,.y,.w,.z){
                  s2<-list.files(.z, recursive = T, full.names = T,pattern="B0[2348]_10m.jp2$")
                  s2<-lapply(s2,function(x)raster(x))
                  pl <- RStoolbox::ggRGB(.x, r=3, g=2, b=1, scale=maxValue(s2[[2]]), stretch = 'hist') +
                    geom_sf(data=.w,aes(fill=ids),colour = NA,alpha=0.8)+
                    scale_fill_manual(values = mycolors,guide="none")+
                    geom_point(data=.y,aes(x=UTME,y=UTMN),shape=21,color="yellow", fill="black",size=0.75)+
                    facet_wrap(~sex)+
                    theme_bw()+
                    theme(text = element_text(size=4),
                          axis.title=element_blank(),
                          axis.text.y = element_text(angle=90, hjust=1),
                          strip.background = element_blank(),
                          strip.text.x = element_blank())
                  pl
                })

#save the maps plots
plsmaps.save <- lapply(maps, ggplotGrob)
walk(plotsaveformat,~{
  ggsave(plot=marrangeGrob(plsmaps.save, nrow = 1, ncol = 3, top=NULL),
         file=paste0(pathplots,"maps_trial",".",.x),
         device=.x,
         height=(doublecolumn[2]/12)*9,width=doublecolumn[2],
         units="px",dpi=600)
})

#END------------
