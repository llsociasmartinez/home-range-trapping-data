Sys.setenv(LANG = "en")
rm()
gc()
library(ctmm)
library(amt)
library(tidyverse)
library(furrr)
library(minpack.lm)
library(mvoutlier)
library(miceadds)
library(rdist)
library(ggpubr)
library(hrbrthemes)
library(plotly)

"%!in%"<-Negate("%in%")

#area covered by trapping area
xlimit<-function(area.c,r,true.area){
  #print(list(area.c,r,true.area))
  areacovered=area.c*true.area
  find.line<-function (h) ((r^2)*acos(h/r))-(h*sqrt((r^2)-(h^2)))-areacovered
  h<-uniroot(find.line,c(-r,r),tol = 0.000001)
  h<-h$root
  return(h)}

#create time slots to simulate data depending on autocorrelation, desired observations etc
findtmp<-function(tobs,twind,ech,auto.c){
  #return a sequence of times separated
  ts<-map_dbl(0:tobs,~.x*auto.c+ .x*twind)
  #print(ts)
  tmp<-map(ts,~seq(.x,
                   (.x+twind),
                   ech)) %>% 
    unlist %>% unique
  #print(tmp[1:100])
  return(tmp)
}


#create grid from telemetry objects with d meters spacing
#or a predefined grid if pre=list() is provided and data=NULL
captgrid<-function(data,dist,ex,pre){
  #gather coordinates of generated movement in the data for x and y
  if(is.null(data)==F){
  x<-lapply(data,function(z) c(max(z$x),min(z$x))) %>% unlist
  y<-lapply(data,function(z) c(max(z$y),min(z$y))) %>% unlist
  #create grid every d*d m
  gx<-seq(min(x)-(ex*dist%#%"m"),max(x)+(ex*dist%#%"m"),dist%#%"m")
  gy<-seq(min(y)-(ex*dist%#%"m"),max(y)+(ex*dist%#%"m"),dist%#%"m")
  }else{
    #create trap locations on a grid every d*d m
    gx <- seq(pre[["min_x"]],pre[["max_x"]], dist%#%"m")
    gy <- seq(pre[["min_y"]], pre[["max_y"]], dist%#%"m")
  }
  #create intersections (i.e. trap locations)
  gridxy<-expand.grid(gx,gy) %>% round(.,0)%>% 
    "colnames<-"(c("x","y"))%>% arrange(x,y)
  return(gridxy)
}

#capture the animals passing near=d from traps
capt<-function(data, gxy, d){
  #print(gxy)
  #Fun calculate the Euclidean distance of a location from all trap locations 
  sld <- function(i,x, y, gxy){
    #which gxy are within a d*d square
    idx<-which(abs(gxy$x-x)<d)
    idy<-which(abs(gxy$y-y)<d)
    idxy<-intersect(idx,idy)

    #if any
    if(length(idxy)!=0){
      #calculate euclidean distance to such gxy
      dist<-map_dbl(.x=idxy,~sqrt((x - gxy$x[.x])^2 + (y - gxy$y[.x])^2))
      #find which are lower than d and if several, which is lowest
      caught<-idxy[which(dist == min(dist))]
      #print(list(i=i,idx=idx,idy=idy,idxy=idxy,dist=dist,caught=caught))
      #this should return one value, either the trap index or NA
      return(caught)
      }else{return(NA)}
  }
  #for each movement location find closest trap in d distance
  #vector nrow movement data, for each movement NA or the index of the trap at d distance
  st <- sapply(1:nrow(data), function(i) sld(i,data$x[i], data$y[i], gxy))
  # gxy[st[is.na(st)==F],c("x","y")] %>% nrow
  # data[is.na(st)==F,c("x","y")] %>% nrow
  if(any(is.na(st)==F)==T){
    data[is.na(st)==F,c("x","y")]<-gxy[st[is.na(st)==F],c("x","y")]
    data$loc<-paste(data$x,data$y,sep=" ")
    
    sctraps<-data[is.na(st)==F,]
    row.names(sctraps)<-1:nrow(sctraps)

    return(sctraps)
  }else{return("no data")}

}


#filter autocorrelated trapping data out
findautocor<-function(data,tt){
  #if there is any data
  if(data[[1]][1] =="no data"){return(data)}
  #with more than 1 row
  if(nrow(data)>1){
  data$t.diff<-c(diff(data$t),0)
  t.diff<-data$t.diff
  
  #fun find for each data point which others are acorrelated (i.e. difference time<tt)
  noacor<-function(y){
    #cumsum dfference in time, which are bigger than tt
    noacor<-cumsum(t.diff[y:length(t.diff)]) %>% 
      ">"(tt) %>% which
    if(length(noacor)!=0){
      #find the first that is not autocorrelated,
      #& correct for row of start to search (focal row)
      noacor<-noacor %>% min %>% "+"(y)
    }else{
      noacor<-y
    }
    if(y==length(t.diff)){noacor<-y}
    return(noacor)
  }
  data$noacor<-sapply(1:nrow(data),function(y)noacor(y))
  #for each line, record and go to next trap noautocorrelated, repeat
  ac<-c(1)
  for(i in 1:(nrow(data))){
    if(i %in% ac){
      ac[length(ac)+1]<-data$noacor[i]
    }else{next}
  }
  ac<-ac[-length(ac)]
  #remove rows not in ac
  #remove data$t.diff & $noacor
  data<-data[ac,-c(ncol(data)-1,ncol(data))]
  data@row.names<-1:nrow(data)
  }
  return(data)
}

#apply akima.bicub(), draw contour and compute area
bic<-function(data,xtraps,ytraps,z1,resol){
  #trap locations & times the id was in each
  tb<-table(data$loc) %>% as.data.frame(.) %>% mutate_if(is.factor,as.character)
  subs.mat<-function(x2,y2,tb){
    #create a list with vectors that represent the coordinates (paste(x,y)) columnwise
    cidxy<-map(x2,~paste(.x,y2,sep=" "))
    #nms<-as.character(1:length(cidxy))
    #map it and pass it to match, convert to matrix
    zz<-map_dfc(cidxy,~tb[match(.x,tb$Var1),"Freq"]) %>% as.matrix
    #substitute nas by 0s
    zz[is.na(zz)==T]<-0
    return(zz)
  }
  capts<-suppressMessages(subs.mat(xtraps,ytraps,tb)) %>% t
  #gives message about names of columns that is irrelevant

  #apply bicubic interpolation
  akima.bic<-akima::bicubic.grid(x=xtraps,
                                 y=ytraps,
                                 z=capts,
                                 nx=resol,
                                 ny=resol)
    #Next convert the interpolated distribution into a probability density that sums to 1
    PDF <- akima.bic$z/sum(akima.bic$z)
    #Calculate the CDF (needed for contours and calculating the 95% area)
    CDF <- ctmm:::pmf2cdf(PDF) 
    #95% home range area
    dV <- prod(diff(akima.bic$x)[1],diff(akima.bic$y)[1]) #Area of a grid cell
    area <- sum(CDF <= 0.95) * dV #Area of all grids within the 95% contour
    return(area)
}

bic_n<-function(data,gridxy,obssample,resol){
  #if there is any data
  if(data[[1]][1]!="no data"){
    library(akima)
    library(ctmm)
    x1<-gridxy$x %>% unique #%>% as.character
    y1<-gridxy$y %>% unique #%>% as.character
    #to mimic the x y space in the movement model, 
    #y should correspond to the rows and be positive to negative
    #x the columns and the dimnames follow along
    z1<-matrix(0,
              nrow=length(y1),
              ncol=length(x1))
    
    area<-map_dbl(.x=obssample,~bic(data[1:.x,],x1,y1,z1,resol))
    return(area)}else{return("no data")}
}

#create the sampling regime (i.e. ordering of observations for further calculation of incremental home range analyses)
#farthest = apply farthest point sampling, boot= bootstrapping but needs modification for implementation
samplingregime<-function(data,obssamp,allids.equal=T,normal=T,farthest=T,boot=F,boottimes){
  if(data[[1]][1]=="no data"){return(data)}
  #nrow(data)
    indxs<-tibble(normal=obssamp)
    if(farthest==TRUE){
      #needs to be reworked to be ok when ids are different (irrelevant current project)
      #Farthest point sampling 
      # if(allids.equal==F){
      #   coords<-cbind(data$x,data$y)
      #   fps<-farthest_point_sampling(coords, metric = "euclidean")
      #   indxs$farthest<-fps
      # }
      #if all ids are same then need to do the fps at each line
      fps<-map(obssamp,~{
        coords<-cbind(data$x[1:.x],data$y[1:.x])
        fps<-farthest_point_sampling(coords, metric = "euclidean")
        fps
        })
      indxs$fps<-fps
      #print(indxs)
      #print(bind_cols(indxs))
     }
    if(boot==T){
      #for bootstrapping, just create new combinations of rows 
      nrows<-1:nrow(data)
      boots<-map_dfc(1:boottimes,~map(obssamp,function(obssmp) {sample(nrows[1:obssmp],replace=T)})) %>% set_names(paste0("boot",c(1:boottimes)))
      indxs$boots<-boots
    }
    
    indxs<-bind_cols(indxs)
    
    # if(farthest==TRUE & allids.equal==TRUE){
    #   indxs$farthest<-fps
    # }
    #print(indxs)
    return(indxs)
}

#set all combinations, grids, simulations to do.
set_variables_process<-function(allids.equal=T,
                             nbids=4,
                            true.area=100*10000 %#% 'm^2',
                            p.intervals=list(auto.c=list(from=c(0),to=c(0),by=c(1)),
                                             area.c=list(from=c(0.2),to=c(1),by=c(0.2)),
                                             trap.d=list(from=c(2),to=c(42),by=c(10)),
                                             essize=list(from=c(2000),to=c(2000),by=c(1))),
                            bufferzone=0,
                            timewind=1 %#% 'hours',each=10 %#% 'minutes',
                            save,pathway){

  
  r <- sqrt(true.area / pi)
  diam<-r*2
  #The spatial autocorrelation based on true hr area desired
  find.sig <- function(sig,area)(-2*log(0.05)*pi*sig)-area
  sig <- uniroot(find.sig, c(-1E10, 1E10), tol = 0.0001, area = true.area)$root
  #The velocity autocorrelation timescale (directional persitence)
  tau_v <- 5 %#% 'min'
  #the parameters
  params<-map(p.intervals,~{pmap(list(f=.x$from,t=.x$to,b=.x$by),function(f,t,b){seq(f,t,by=b)}) %>% 
      unlist %>% unique}) 
  area.c.coord<-map_dbl(params$area.c,~xlimit(area.c=.x,r=r,true.area=true.area)) %>% round
  
  if(allids.equal==T){
    #this option will use only the max value for essize/nbobs to simulate a maximum nb of observations
    #and then trap them in different grids and areas covered 
    params$essize<-max(params$essize)
    ids<-list(1:nbids)
    params<-prepend(params,ids) %>% set_names(c("id",names(p.intervals)))

    #if params with only one value, no interval get one value back
    cmbn<-expand.grid(params) %>% set_names(c("id",names(p.intervals))) %>% as_tibble
    #print(head(cmbn))
    #the model/s
    #there are as many models as params$auto.c
    #as many tmp as essize
    cmbn.mdltmp<-tibble(auto.c=params$auto.c,essize=max(params$essize))
    cmbn.mdltmp<-cmbn.mdltmp %>% 
      mutate(model=map(auto.c,~{if(.x==0){.x<-NA}else{.x<-.x %#% 'day'}
        ctmm(tau=switch(is.na(.x)==T, NULL, c(.x, tau_v)),isotropic=TRUE,sigma=sig)}),
        tmp=map2(essize,auto.c,~findtmp(tobs=.x,twind=timewind,ech=each,auto.c=.y)))
    #print(head(cmbn.mdltmp))

    #the trapping grid/s 
    #as many as trap.d*area.c
    interval.m<-r/params$trap.d %>% round
    diam<-sqrt(true.area/pi)*2
    pre <- map(.x=area.c.coord,~list(.x-bufferzone,.x+diam+bufferzone,-(diam/2)-bufferzone,(diam/2)+bufferzone) %>% 
                 "names<-"(c("min_x","max_x","min_y","max_y")))
    #for all combinations of area.c and trap.d, create gridxy
    cmbn.grid<-expand.grid(interval=1:length(interval.m),predef=1:length(pre))
    cmbn.grid<-cmbn.grid %>% mutate(trap.d=params$trap.d[interval],
                                    area.c=params$area.c[predef],
                                    gridxy=map2(.x=interval,.y=predef,
                                     ~captgrid(data=NULL,dist=interval.m[.x],pre=pre[[.y]]))) %>% 
      select(-interval,-predef)
    #print(head(cmbn.grid))
    #save the params combns in case
    presimul<-list(cmbn=cmbn,cmbn.grid=cmbn.grid,cmbn.mdltmp=cmbn.mdltmp)
    if(save==T){save(presimul,file=paste0(pathway,"/presimul_allidsequal",".Rdata"))}
    
    
  }else{
    
     cmbn<-map_dfc(params,~if(length(.x)>1){sample(.x,nbids,replace=T)}else{rep(.x,nbids)})
     cmbn$area.c.coord<-area.c.coord[match(cmbn$area.c,params$area.c)]
     id<-c(1:nbids)
     cmbn<-cbind(id,cmbn) %>% as_tibble
     
     #values of variables for combinations
     auto.c.lvl<-cmbn$auto.c %>% unique %>% sort
     essize.lvl<-cmbn$essize %>% unique %>% sort
     trap.d.lvl<-cmbn$trap.d %>% unique %>% sort
     area.c.lvl<-cmbn$area.c %>% unique %>% sort
     area.c.coord.lvl<-area.c.coord[match(area.c.lvl,params$area.c)]
     #the model/s
     #there are as many models as params$auto.c
     #as many tmp as essize

     cmbn.mdltmp<-expand.grid(auto.c.lvl,essize.lvl) %>% set_names("auto.c","essize")
     model<-map(.x=auto.c.lvl,~{
       if(.x==0){.x<-NA}else{.x<-.x %#% 'day'}
       ctmm(tau=switch(is.na(.x)==T, NULL, c(.x, tau_v)),isotropic=TRUE,sigma=sig)})%>% 
       "names<-"(c(auto.c.lvl))

     cmbn.mdltmp<-cmbn.mdltmp %>% 
       mutate(model=model[match(auto.c,auto.c.lvl)],
              tmp=map2(essize,auto.c,~findtmp(tobs=.x,twind=timewind,ech=each,auto.c=.y)))
     
     #the trapping grid/s 
     #as many as trap.d*area.c
     interval.m<-r/trap.d.lvl %>% round
     pre <- map(.x=area.c.coord.lvl,~list(.x-bufferzone,.x+diam+bufferzone,-(diam/2)-bufferzone,(diam/2)+bufferzone) %>% 
                  "names<-"(c("min_x","max_x","min_y","max_y")))
     #for all combinations of area.c and trap.d, create gridxy
     cmbn.grid<-expand.grid(interval=1:length(interval.m),predef=1:length(pre))
     cmbn.grid<-cmbn.grid %>% mutate(trap.d=trap.d.lvl[interval],
                                     area.c=area.c.lvl[predef],
                                     gridxy=map2(.x=interval,.y=predef,
                                                 ~captgrid(data=NULL,dist=interval.m[.x],pre=pre[[.y]]))) %>% 
       select(-interval,-predef)
     presimul<-list(cmbn=cmbn,cmbn.grid=cmbn.grid)
     if(save==T){save(cmbn,file=paste0(pathway,"/presimul",".Rdata"))}
    }
  return(presimul)
}

#simulate movement and capture
sim_capt<-function(allids.equal,chkid,cmbn,cmbn.grid,cmbn.mdltmp,tr.area,attract,save,pathway){
  if (allids.equal == T){
    #for each individual simulate all auto.c with maximum essize
    #here each id is simulated with as many models (=nb of different autocorrelations)
    # print(cmbn)
    # print(cmbn$id)
    ids<-cmbn$id %>% unique
    ids.sim<-imap_dfr(1:length(ids),
                      ~{x<-cmbn.mdltmp %>%
                        mutate(id=ids[.y],
                               sim=map2(model,
                                        tmp,function(mdl,tm){ctmm::simulate(mdl,t = tm)})) %>%
                        select(-model,-tmp)
                      x})
    #print("simulation finished")
    
    #2.b "Capture" the animal in grid system------
    #for each sim, capture by as many cmbn.grid
    
    st<-Sys.time()
    ids.sim <- ids.sim %>%
      mutate(trap=map(sim,
                      ~cmbn.grid %>% mutate(captures=map2(gridxy,area.c,function(grd,areac){
                        #reduce the number of simuls for trapping if the area.c is bigger than 0.5
                        if(areac>=0.5){nr<-(nrow(.x)*0.25) %>% round}else{nr<-nrow(.x)}
                        cpt<-capt(data=.x[1:nr,],gxy=grd, attract)
                        cpt})) %>%
                        select(-gridxy)))
    
    #remove sim, unnest and remove captures occurred same time
    ids.sim<-ids.sim %>% #mutate(nrsim=map_dbl(sim,nrow)) %>% 
      select(-sim)%>%
      unnest(trap) %>% mutate(captures=map(captures,~.x[!duplicated(.x$t),]))
    
    #print("trapping finished")
    # print(Sys.time()-st)
    
    #Remove autocorrelated data traps-----
    ids.sim<-ids.sim %>% mutate(captures=map2(.x=captures,.y=auto.c,~if(.y!=0){findautocor(data=.x,tt=.y %#% "day")}else{.x}))
    #print("rm auto.c finished")
    
  }else{
    #matches per id
    mtchmt<-match(paste(cmbn$auto.c,cmbn$essize),paste(cmbn.mdltmp$auto.c,cmbn.mdltmp$essize))
    mtchgr<-match(paste(cmbn$trap.d,cmbn$area.c),paste(cmbn.grid$trap.d,cmbn.grid$area.c))
    #simulate
    ids.sim<-cmbn %>% mutate(sim=map2(cmbn.mdltmp$model[mtchmt],cmbn.mdltmp$tmp[mtchmt],~ctmm::simulate(.x,t = .y)))
    #trap
    ids.sim <- ids.sim %>% mutate(captures=pmap(list(sm=sim,grd=cmbn.grid$gridxy[mtchgr],areac=cmbn.grid$area.c[mtchgr]),
                                                function(sm,grd,areac){
                                                  #reduce the number of simuls for trapping if the area.c is bigger than 0.5
                                                  if(areac>=0.5){nr<-(nrow(sm)*0.25) %>% round}else{nr<-nrow(sm)}
                                                  cpt<-capt(data=sm[1:nr,],gxy=grd, attract)
                                                  cpt}))
    #Remove autocorrelated data traps-----
    ids.sim<-ids.sim %>% mutate(captures=map2(.x=captures,.y=auto.c,~if(.y!=0){findautocor(data=.x,tt=.y %#% "day")}else{.x}))
    
  }
  
  if(save==T){save(ids.sim,file=paste0(pathway,"/sim_",chkid,".Rdata"))}
  return(ids.sim)
  
}

#resimulate if any id was not captured enough
resim_capt_enough<-function(ids.sim,cmbn.grid,cmbn.mdltmp,attract,timewind=6 %#% 'hours',each=10 %#% 'minutes',buffer=1.5,save,pathway,paralel=F,nbcores=2,chunk,mac.linux=F){
  #desired ssize attained?
  ids.sim<-ids.sim %>% mutate(desi=map2_dbl(captures,essize,~nrow(.x)-(.y))) %>% select(id,everything())
  nope<-which(ids.sim$desi<0)#which have less traps than should
  if(length(nope)==0){
    ids.sim<-ids.sim %>% mutate(captures=map2(captures,essize,~.x[1:.y,])) %>% select(-desi)
    if(save==T){save(ids.sim,file=paste0(pathway,"/sim_total_",nrow(ids.sim),".Rdata"))}
    return(ids.sim)
  }
  if(save==T){save(ids.sim[nope,],file=paste0(pathway,"/sim_nope_",nrow(ids.sim[nope,]),".Rdata"))}
  print(list("not enough data in",nope))

  #save what is already ok
  idx<-1:nrow(ids.sim)
  ids.sim.ok<-ids.sim[idx%!in%nope,]
  ids.sim.ok<-ids.sim.ok %>% mutate(captures=map2(captures,essize,~.x[1:.y,]))
  if(save==T){save(ids.sim.ok,file=paste0(pathway,"/sim_ok_partial_",nrow(ids.sim.ok),".Rdata"))}
  
  chks<-split(nope, ceiling(seq_along(nope)/chunk))
  f<-1;last<-length(chks);todo<-f:last
  #now I have the ids with their parameters, the models, grids and
  if(paralel==T & mac.linux==F){plan(multisession, workers = nbcores)}
  if(paralel==T & mac.linux==T){plan(multicore, workers = nbcores)}
  ids.sim.nope.ok<-future_map_dfr(.progress=T,.options = furrr_options(seed = TRUE),
                                  .x=todo,~{
                                    data<-ids.sim[chks[[.x]],]
                                    notenough<-1:nrow(data)
                                    #match for using the models and grids
                                    mtchmt<-match(paste(data$auto.c,data$essize),paste(cmbn.mdltmp$auto.c,cmbn.mdltmp$essize))
                                    mtchgr<-match(paste(data$trap.d,data$area.c),paste(cmbn.grid$trap.d,cmbn.grid$area.c))
                            
                                    while(length(notenough)>0){
                                      ae<-data[notenough,c("auto.c","essize")]
                                      ae<-ae[!duplicated(ae),]
                                      #need to resimulate but with more tmp
                                      #how far are we
                                      howfar<-abs(data$desi[notenough])/data$essize[notenough]
                                      newessize<-(data$essize[notenough]/howfar)*buffer
                                      ids.sim.notenough<-data[notenough,] %>% 
                                        mutate(tmp=map2(newessize,auto.c,~findtmp(tobs=.x,twind=timewind,ech=each,auto.c=.y)))
                                      ids.sim.notenough<-ids.sim.notenough %>% mutate(sim=map2(cmbn.mdltmp$model[mtchmt[notenough]],tmp,~ctmm::simulate(.x,t = .y)))
                                      
                                      ids.sim.notenough <- ids.sim.notenough %>% 
                                        mutate(captures=map2(sim,cmbn.grid$gridxy[mtchgr[notenough]],~capt(data=.x,traps=.y, attract)))
                                      ids.sim.notenough<-ids.sim.notenough %>% select(-tmp,-sim)
                                      #Remove autocorrelated data traps-----
                                      ids.sim.notenough<-ids.sim.notenough %>% 
                                        mutate(captures=map2(.x=captures,.y=auto.c,~if(.y!=0){findautocor(data=.x,tt=.y %#% "day")}else{.x}))
                                      #resubstitute to partial data
                                      data[notenough,]<-ids.sim.notenough
                                      #still not enough?
                                      data<-data %>% mutate(desi=map2_dbl(captures,essize,~nrow(.x)-(.y)))
                                      notenough<-which(data$desi<0)
                                      print(list("still not enough data in",notenough))}
                                    
                                    return(data)})
  future:::ClusterRegistry("stop")
  
  ids.sim[nope,]<-ids.sim.nope.ok 
  ids.sim<-ids.sim %>% mutate(captures=map2(captures,essize,~.x[1:.y,])) %>% select(-desi)
  if(save==T){save(ids.sim,file=paste0(pathway,"/sim_total_",nrow(ids.sim.ok),".Rdata"))}
  return(ids.sim)
  
}

#fit estimators to a given dataset
estimator_fit<-function(estimators,data1.ctmm,data1.amt,x1,y1,z1,akima.res){
  res_estim<-list()
  nr<-nrow(data1.amt)

  if("akima_bic" %in% estimators){res_estim$akima_bic<-bic(data1.amt,x1,y1,z1,resol=akima.res)}

  if("ctmm_akde" %in% estimators){
    ctmm.model <- tryCatch({ctmm.guess(data1.ctmm,interactive=FALSE)},error=function(e){NA})
    ctmm.model <-tryCatch({ctmm.fit(data1.ctmm,ctmm.model)},error=function(e){NA})
    if(is.na(ctmm.model[[1]])==F){
      res_estim$ctmm_akde<-tryCatch({ctmm::akde(data1.ctmm, CTMM = ctmm.model, weights = TRUE)},error=function(e){NA})
      if(is.na(res_estim$ctmm_akde[[1]][1])==F){
        res_estim$ctmm_akde<-summary(res_estim$ctmm_akde, units = FALSE)$CI[,2]}
      }else{res_estim$ctmm_akde<-NA }}

  if("amt_kde" %in% estimators){
    bandw<- tryCatch({hr_kde_ref_scaled(data1.amt,num.of.parts=1,levels=0.95)},error=function(e){NA})
    if(is.na(bandw)==F){ 
      res_estim$amt_kde<-tryCatch({amt::hr_kde(data1.amt,levels=0.95,h=bandw,rand_buffer=1e-05) %>% list}, error=function(e){NA})
    }else{res_estim$amt_kde<-NA}}
  
  if("amt_mcp" %in% estimators){res_estim$amt_mcp<-if(nr>2){
    tryCatch({amt::hr_mcp(data1.amt,levels=0.95,rand_buffer=1e-05) %>% list},error=function(e){NA %>% list})}else{NA %>% list}}
  
  if("amt_locoh" %in% estimators){res_estim$amt_locoh<-if(nr>2){
    tryCatch({amt::hr_locoh(data1.amt,n = ceiling(sqrt(nr)),rand_buffer=1e-05) %>% list},error=function(e){NA %>% list})}else{NA %>% list}}
  
  # if("amt_kde_modctmm" %in% estimators){res_estim$amt_kde_modctmm<-tryCatch({ amt::hr_akde(data1.amt,levels=0.95,model=ctmm.model) %>% list},
  #                                                                           error=function(e){NA})}
  return(res_estim)
  }

#fit michaelis or monomolecular asymptotic model to a given set or areas and times
#times in this analyses correspond to the number of observations since autocorrelation is 0
model_fit<-function(areas,times,initval=NULL,type="micmen",pl.data=F){
  res<-list(asym=NA,mconfest=NA,munstab=NA,gfit=NA)
  #print(list(area=areas,time=times))
  datat<-data.frame(area=areas,
                    time=times)
  
  datat<-datat[complete.cases(datat),]
  
  if(pl.data==T){
    plot(y=datat$area,x=datat$time)
  }
  
  if(type=="micmen"){
    form<-formula(area~(asym*time)/(k+time))
    if(is.null(initval)==T){ 
      initval<-tryCatch({stats::getInitial(area ~ SSmicmen(time,asym,k),data=datat)}, 
                        error=function(e){NA})
      
    }else{initval<-list(asym=initval["asym"],k=2*log(2)/initval["asym"])}}
  
  if(type=="monmol"){
    form<-formula(area~asym*(1-exp(-exp(lrc) * time)))
    #a numeric vector of the same length as input. 
    #It is the value of the expression Asym*(1 - exp(-exp(lrc)*input)). 
    #If all of the arguments Asym and lrc are names of objects, 
    #the gradient matrix with respect to these names is attached as 
    #an attribute named gradient.
    if(is.null(initval)==T){
      initval<-tryCatch({stats::getInitial(area ~ SSasympOrig(time,asym,lrc), data=datat)}, 
                        error=function(e){NA})
    }else{initval<-list(asym=initval["asym"],lrc=initval["lrc"])}}
  #print(initval)
  if(is.na(initval[1])==F){#if initval is not NA
    library(minpack.lm)
    m<-tryCatch({nlsLM(form,start=initval,data=datat,
                       control = nls.lm.control(maxiter = 500))}, 
                error=function(e){NA})
    
    #print(list(form=form,init=initval,len.m=length(m),m=m))
    # print(nlsLM(form,start=initval,data=datat,
    #             control = nls.lm.control(maxiter = 500)))
    if(length(m)>1){#if m is not NA
      coeffs<-coefficients(m) 
      asym<-coeffs[1]
      res$asym<-unname(asym)
      #goodness of fit
      res$gfit<-cor(datat[["area"]],predict(m))
      
      # check model stability
      #removing one datavalue at a time
      all.coeffs=matrix(NA,nrow=nrow(datat),ncol=length(coefficients(m)))
      colnames(all.coeffs)=names(coefficients(m))
      
      all.coeffs<-map(1:nrow(datat),~{
        xx<-tryCatch({nlsLM(form,start=initval,data=datat[-.x,],
                            control = nls.lm.control(maxiter = 500))}, 
                     error=function(e){NA})
        if(length(xx)>1){return(coefficients(xx))}else{return(rep(NA,length(coeffs)))}
      })%>% do.call(rbind,.)
      
      munstab<-apply(all.coeffs,MARGIN=2,FUN=function(x)max(x,na.rm=T)-min(x,na.rm=T))[1]
      res$munstab<-unname(munstab)
      
      #confidence interval estimate
      mcf <- function(mdl){
        tryCatch({
          #need to do it twice otherwise trycatch return empty object
          confint(mdl);mc<-confint(mdl);return(mc)},
          error=function(e){NA})}
      
      mconfest<-suppressMessages(mcf(m))
      
      if(length(mconfest)>1){
        res$mconfest<-if(any(is.na(mconfest[1,])==T)){NA}else{max(mconfest[1,])-min(mconfest[1,])} 
      }
      #print(list(asym=asym,mconfest=mconfest,munstab=munstab,gfit=gfit))
      #return(list(asym,mconfest,munstab,gfit))
      return(res)
    }else{#print("NO MODEL")
      return(res)}
  }else{#print("NO STARTING VALUES")
    return(res)
  }}


#fit estimators and asymptotic models to the simulation data
fit_estimators<-function(allids.equal,captures,cmbn.grid,
                                p.intervals=list(essize=list(from=c(2,10,20,30,60,100,500),to=c(10,20,30,60,100,500,2000),by=c(1,2,5,10,20,100,500))),
                                estimators,akima.res=1000,
                                sampr=list(normal=T,farthest=F,boot=F,boottimes=20),
                                pathway,save,paralel=F,nbcores=2,chunk=1,mac.linux=F){

  #nb observations to sample-----
  params<-map(p.intervals,~{pmap(list(f=.x$from,t=.x$to,b=.x$by),function(f,t,b){seq(f,t,by=b)}) %>%
      unlist %>% unique})
  captures<-captures %>% mutate(obssample=map(essize,~params$essize[which(params$essize<=.x)]))

  #transform ctmm telemetry to trackxy amt
  dat.amt<-pmap(list(data=captures$data.ctmm,id=captures$id,idx=1:length(captures$data.ctmm)),function(data,id,idx){
    dat<-data@.Data %>% "names<-"(colnames(data)) %>% do.call(cbind,.) %>% as.data.frame(stringsAsFactors=F)
    dat[,c("t","x","y")]<-lapply(dat[,c("t","x","y")],function(x)as.numeric(x))
    dat$id<-id
    dat$idx<-idx
    dat }) %>% bind_rows %>% make_track(x, y, t, loc=loc, id = id, idx = idx)
  captures<-captures %>% mutate(data.amt=dat.amt %>% group_by(idx) %>% group_split)
  
  #add sampling regimes & grid matches for bicubic-----
  captures<-captures %>%
    mutate(idx=1:nrow(captures),
           mtchgr=match(paste(trap.d,area.c),paste(cmbn.grid$trap.d,cmbn.grid$area.c)))

  #divide the dataset to run either seq or in paralel by individual
  if(allids.equal==T){
    ids <- captures$id %>% unique
    chks<-split(ids, ceiling(seq_along(ids)/chunk))
  }else{
    idxs<-captures$idx
    chks<-split(idxs, ceiling(seq_along(idxs)/chunk))
  }

  f<-1;last<-length(chks);todo<-f:last

  if(paralel==T & mac.linux==F){plan(multisession, workers = nbcores)}
  if(paralel==T & mac.linux==T){plan(multicore, workers = nbcores)}
  #print(todo)
  estimators_df<-future_map_dfr(.progress=T,.options = furrr_options(seed = TRUE),
                  todo,function(td){
                               #estimate areas
                               if(allids.equal==T){
                                 estimators_id<-captures[which(captures$id%in%chks[[td]]),]}else{
                                 estimators_id<-captures[chks[[td]],]
                                 }

                           #create sampling regimes
                           estimators_id<-estimators_id %>% mutate(rowstouse=map2(data.ctmm,obssample,function(dat,obssam){samplingregime(dat,obssam,allids.equal,normal=sampr$normal,farthest=sampr$farthest,boot=sampr$boot,sampr$boottimes)}))

                           estimators_id<-estimators_id %>%
                                 mutate(raws=pmap(list(idxx=idx,d.ctmm=data.ctmm,d.amt=data.amt,rwstuse=rowstouse,grididx=mtchgr,obssamp=obssample),
                                                  function(idxx,d.ctmm,d.amt,rwstuse,grididx,obssamp){
                                                    #akima
                                                    grid<-cmbn.grid$gridxy[[grididx]]
                                                    x1<-grid$x %>% unique
                                                    y1<-grid$y %>% unique
                                                    z1<-matrix(0,nrow=length(y1),ncol=length(x1))
                                                    
                                                    #ctmm
                                                    dd.ctmm<-d.ctmm[1:nrow(d.ctmm),]
                                                    dd.ctmm$t<-d.ctmm$t
                                                    dd.ctmm$vx <- NULL
                                                    dd.ctmm$vy <- NULL
                                                    
                                                    #calculate area per estimator each subrows and each sampling
                                                    estim_sampl<-map2_dfr(rwstuse,names(rwstuse),
                                                                         function(rowst,samp){
                                                                           #how many samples per each id an treatment
                                                                           nbobssamp<-1:length(rowst)
                                                                           
                                                                           #calculate area at each nb of observations
                                                                           estim_nbobs<-map_dfr(nbobssamp,function(nbobssam){
                                                                             
                                                                             # FPS-------------
                                                                             #a bit more complicated, need to recalculate for each new sample the areas until there
                                                                             #calculated with fps until there also
                                                                             if(samp=="fps"){#| grepl("boot", samp, fixed = TRUE)==TRUE
                                                                               #extract list of obs reordered with fps for the given sample
                                                                               subrows<-rowst[[nbobssam]]
                                                                               #follow alternative curve by fps until current nbobs by
                                                                               #mapping over obssample untill current
                                                                               res_estim<-map_dfr(obssamp[1:nbobssam],function(obssamptill){
                                                                                 #use all subrows until there (sort of recursive version)
                                                                                 subrowss<-subrows[1:obssamptill]
                                                                                 res_estims<-suppressMessages(estimator_fit(estimators,dd.ctmm[subrowss,],d.amt[subrowss,],x1,y1,z1,akima.res))
                                                                                 res_estims
                                                                               })
                                                                               
                                                                               #NORMAL----
                                                                               }else{
                                                                               subrows<-1:rowst[nbobssam]
                                                                               res_estim<-suppressMessages(estimator_fit(estimators,dd.ctmm[subrows,],d.amt[subrows,],x1,y1,z1,akima.res))
                                                                               }
                                                                             #collect all info into 1 row
                                                                             res_estim<-res_estim %>% as_tibble
                                                                             res_estim[]<-lapply(res_estim,function(eachcol){eachcol[1]<-list(eachcol);eachcol})
                                                                             res_estim<-res_estim[1,]
                                                                             
                                                                             res_estim
                                                                           })
                                                                           
                                                                           #add nb of observations used for each row and the type of sampling
                                                                           estim_nbobs$nbobs<-obssamp
                                                                           estim_nbobs$sampling<-samp
                                                                           
                                                                           estim_nbobs
                                                                           })
                                                    
                                                    #retrieve areas from amt estimators
                                                    estim_sampl<-estim_sampl %>% 
                                                      mutate(across(.names="{col}",all_of(contains("amt")),
                                                                    ~map(., function(el){
                                                                      #recursive to account for fps collection of areas until each sample (the alternative courve)
                                                                      ar<-map_dbl(el,function(ella){
                                                                        if(is.na(ella)==F){ars<-ella %>% amt::hr_area(.) %>% pull(area);ars}else{NA}})
                                                                      return(ar)})),
                                                             model="raw")
                                                    estim_sampl
                                                    })
                                        )
                          
                          save(estimators_id,file=paste0(pathway,"/estimators_id_",unique(estimators_id$id),".Rdata"))

                          estimators_id<-estimators_id %>% select(-all_of(c("data.ctmm","data.amt","rowstouse")))
                          return(estimators_id)
                          })

  #stop cluster avoid running afterwards in vain
  future:::ClusterRegistry("stop")
  #SAVE whole info
  if(save==T){save(estimators_df,file=paste0(pathway,"/estimators_df_all_",
                                          min(c(map(estimators_df,~.x$id %>% unique) %>% unlist)),"_",
                                          max(c(map(estimators_df,~.x$id %>% unique) %>% unlist)),".Rdata"))}
  return(estimators_df)
  }


  
fit_asympmodels<-function(estimators_df,
                          estimators,asympmodels,recursiv.model=T,iv=NULL,
                          pathway,save,paralel,nbcores,chunk,mac.linux=F){
  ids <- estimators_df$id %>% unique
  chks<-split(ids, ceiling(seq_along(ids)/chunk))
  
  f<-1;last<-length(chks);todo<-f:last
  
  if(paralel==T & mac.linux==F){plan(multisession, workers = nbcores)}
  if(paralel==T & mac.linux==T){plan(multicore, workers = nbcores)}
  
  #print(todo)
  results<-future_map_dfr(.progress=T,.options = furrr_options(seed = TRUE), #future_
                      todo,function(td){

                        results_id<-estimators_df[which(estimators_df$id%in%chks[[td]]),]

                        #asymptotic models-------
                        #for each obssample, do the modelling
                        results_id <- results_id %>%
                          mutate(raws=pmap(list(rws=raws,obssam=obssample),
                                           function(rws,obssam){
                                             #print(list(rawsp=rws,obssamplep=obssam))
                                             nrows<-1:nrow(rws);initval<-NULL

                                             area_asymp<-map_dfr(asympmodels,function(asympmod){
                                               if(recursiv.model==T){
                                                 #use group_by sampling
                                                 #print(list(asympmod=asympmod))
                                                 model_res<-rws %>% group_by(sampling) %>%
                                                   mutate(across(.names="{col}",all_of(c(estimators)),
                                                                 ~{areas.col<-.;
                                                                 mfit<-map2(areas.col,1:n(),function(areas,nrws){

                                                                   #if structure corresponds to fps collection
                                                                   if(length(areas)>1){ areass<-areas %>% unlist
                                                                   #if normal
                                                                   }else{ areass<-areas.col[1:nrws] %>% unlist(.,recursive=TRUE)}
                                                                   #print(list(areas=areas,areass=areass))
                                                                   x<-model_fit(areass,obssam[1:nrws],initval,type=asympmod);
                                                                   return(x)
                                                                 })
                                                                 mfit
                                                                 }),
                                                          model=asympmod)
                                               }else{
                                                 model_res<-.x %>%
                                                   summarise(across(.names="{col}",all_of(c(estimators)),~{x<-model_fit(.,obssam,initval,type=asympmod);x[[1]]})) %>%
                                                   mutate(model=asympmod)
                                               }
                                               model_res})
                                             rws<-rbind(rws,area_asymp)
                                             rws
                                           }))
                        results_id<-results_id %>% select(-all_of(c("obssample")))
                        #leave out raw data for fps at a given point, its equivalent to normal
                        results_id$raws<-map(results_id$raws,~{on<-which(.x$model!="raw");tw<-which(.x$sampling!="fps");.x[unique(c(on,tw) %>% sort),]})
                        #add to big table
                        
                         #print(idssave)
                        if(save==T){save(results_id,file=paste0(pathway,"/results_id_",
                                                                min(results_id$id %>% unique),"_",
                                                                max(results_id$id %>% unique),".Rdata"))}
                        return(results_id)
                      })
  
  #stop cluster avoid running afterwards in vain
  future:::ClusterRegistry("stop")
  
  results<-estimators_df
  #SAVE whole info
  if(save==T){save(results,file=paste0(pathway,"/results_all_",
                                       min(results$id %>% unique),"_",
                                       max(results$id %>% unique),".Rdata"))}
  return(results)
}



#convert the estimators and models to a long format,
#also returns the area as a percentate of true area
convert_long<-function(captures,estimators,true.area,
                       pathway,save){
    #unnest areas per individual and treatment
    results<-unnest(captures,cols = c(raws))
    #tranform columns of estimator into rows
    results<-results %>% pivot_longer(cols=all_of(estimators),values_to="area",names_to="estimator")
    #collect asymptotic model diagnostics and the asymptote
    results<-results %>% mutate(area2=map(area,"asym"),
                            area3=case_when(model=="raw"~area,
                                            TRUE~area2),
                            models_diag=map(area,~c(map(.x[2:3],function(md) md/true.area)),.x[4]),
                            area=map_dbl(area3,~if(is.na(.x)==T){NA}else{if(.x==0){0}else{.x/true.area}}),
                            error=map_dbl(area3,~if(is.na(.x)==T){NA}else{(.x-true.area)/true.area})) %>% select(-area2,-area3)
    
    results<-results %>% mutate(model=factor(model,levels=c("raw","micmen","monmol")),
                                sampling=factor(sampling,levels=c("normal","fps")),
                                estimator=factor(estimator,levels=c("ctmm_akde","amt_kde","amt_mcp","amt_locoh","akima_bic")))
    
  #save
  if(save==T){save(results,file=paste0(pathway,"/results_all_long_",
                                       min(results$id %>% unique),"_",
                                       max(results$id %>% unique),".Rdata"))}
  return(results)
}


#compile functions-----
funct<-lsf.str() %>% as.vector
library(compiler)
map(.x=funct,~cmpfun(get(.x)))
