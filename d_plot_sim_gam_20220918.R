Sys.setenv(LANG = "en")
memory.limit(size=100000000)
source(paste0(getwd(),"/Ecography/scripts/b_functions_sim_20220918.R"))
source(paste0(getwd(),"/Ecography/scripts/c_functions_gam_hrerror_20220918.R"))

pathcomplete<-paste0(getwd(),"/Ecography/files/simulations/complete/") #for complete datasets to be used in further scripts
pathplots<-paste0(getwd(),"/Ecography/plots/")
#load GAM
filenames<-list.files(pathcomplete,pattern="gam_gauss_ti")
load(paste0(pathcomplete,filenames[[1]]))

#load the data
filenames<-list.files(pathcomplete,pattern="datagam")
load(paste0(pathcomplete,filenames[[1]]))
#datagam<-dtm


#set parameters----
true.area=100*10000 %#% 'm^2'
estimators<-c("ctmm_akde","amt_kde","amt_mcp","amt_locoh","akima_bic")# 
asympmodels<-c("raw","micmen","monmol")
sampregs<-c("normal","fps")


#A.Manuscript figures----------
#indications by Ecography in pixels for 300dpi
singlecolumn<-c(945,945*2)#1st for 300 dpi, 2nd 600
singlehalfcolumn<-c(1476,1476*2)
doublecolumn<-c(1961,1961*2)
plotsaveformat<-c("tiff","pdf")


#A.1 Published articles trapping home range---------
#Figure 1
articles<-read.csv(file=paste0(getwd(),"/Ecography/files/table_articles_telemetry.csv"))
ws<-articles$words_search %>% unique
articles<-articles %>% mutate(search=factor(words_search,levels=c(ws[2],ws[c(1,3)]),labels=c("a. trap NOT home range","b. trap AND home range","c. telemetry AND home range")))
pl_articles<-ggplot(articles)+
  geom_line(aes(x=year,y=nb_studies,linetype=search))+
  scale_linetype_manual(values=c(1,3,2)) + 
  theme_classic()+
  labs(linetype = "Type of study",x="Year",y="Number of studies")+
  theme(axis.title = element_text(size = 10),
        axis.text = element_text(size=8),
        axis.text.x = element_text(angle = 25, vjust = 0.5),
        legend.title=element_text(size = 10),
        legend.text=element_text(size=8),
        legend.position = c(0.3, 0.8))
  
pl_articles

#ggsave(paste0(pathplots,"figure_1",".pdf"),height=(singlehalfcolumn[2]/7)*5,width=singlehalfcolumn[2],units="px",dpi=600)
walk(plotsaveformat,~{
  ggsave(plot=pl_articles,
         file=paste0(pathplots,"figure_1",".",.x),
         device=.x,
         height=(singlehalfcolumn[2]/7)*5,width=singlehalfcolumn[2],
         units="px",dpi=600)
})


#figure 2 and SI1 made with ppoint
#A.2 GAM results-------
#Figures 3 to 5-----------
#accuracy and reliability of estimations in the different conditions
nms<-c(paste0("nbobs_trapd_area.c_",c("1","0.6","0.2")),
       paste0("nbobs_areac_trap.d_",c("4","22","42")),
       paste0("trapd_areac_nbobs_",c("7","100","1000")))
topred<-map(nms,~{
  nmsspl<-str_split(.x,pattern="_",simplify = T)
  var<-nmsspl[1,3]
  val<-nmsspl[1,4] %>% as.numeric
  tocr<-
    list(
      approach=levels(gmodel$model$approach),
      nbobs=10^seq(log10(min(datagam$nbobs)),log10(2000),length.out=50),
      trap.d=seq(min(datagam$trap.d),max(datagam$trap.d),by=10),
      area.c=seq(min(datagam$area.c),max(datagam$area.c),by=0.2))
  tocr[[which(names(tocr)==var)]]<-val
  tocr<-tocr %>% cross_df(.filter = NULL) %>% mutate(lg10nbobs=log10(nbobs))
  return(tocr)
})

#predict with adapted pred.gam function (takes ages with so many observations)
# preds3d<-pred_gam(gmodel,topred,specials=c("[()]","log10"),#log10 is for the case when log10(nbobs) was done in the model formula
#                   varint=list(c("approach","lg10nbobs","area.c","trap.d")),
#                   nms=names(topred))
# save(preds3d,file=paste0(pathcomplete,"/preds3d.Rdata"))
load(file=paste0(pathcomplete,"/preds3d.Rdata"))

#extract data
dtat<-map(preds3d,"pred")
#reconf
dtat<-map(dtat,~{
  .x$approach<-str_replace_all(.x$approach,"\\.","--")
  .x<-bind_cols(.x,colsplit(.x$approach, "--", names = c("sampling", "estimator", "model"))) %>% 
    mutate(cip=(upperfitili-lowerfitili)/fitili,
           model=factor(model,levels=asympmodels,labels=c("Raw","MicMen","MonMol")),
           estimator=factor(estimator,levels=estimators,labels=c("AKDE","KDE","MCP","LoCoH","BicubIt")),
           sampling=factor(sampling,levels=sampregs,labels=c("Time-ordered","Distance-ordered")))
  .x})

#trim extreme values to be able to observe interesting variation with color
dtattrim<-map(dtat,~{.x[which(.x$fitili<0),"fitili"]<-0;.x[which(.x$fitili>2),"fitili"]<-2;.x})

#params for contour plots defining accuracy and reliability niches
sensib<-0.1
brks <- c(1-sensib, 1+sensib)
brks2<-c(1+(sensib*1.5))
brks3<-c(1-(sensib*1.5))
brksci1<-0.4
brksci2<-0.6
brksci3<-1

#create a list of variables etc to do all plots in the same function using map
#pairs of variables that change
vars<-list(c("nbobs","trap.d"),c("nbobs","trap.d"),c("nbobs","trap.d"),
           c("nbobs","area.c"),c("nbobs","area.c"),c("nbobs","area.c"),
           c("trap.d","area.c"),c("trap.d","area.c"),c("trap.d","area.c"))
#labs for the plots
labss<-list(c("Number of observations","Trap density"),c("Number of observations","Trap density"),c("Number of observations","Trap density"),
            c("Number of observations","Area covered"),c("Number of observations","Area covered"),c("Number of observations","Area covered"),
            c("Trap density","Area covered"),c("Trap density","Area covered"),c("Trap density","Area covered"))

#themes for plots because of bug in lemon package
thlemon1<-"theme(axis.title.x = element_text(size = 10, margin = margin(t = 18)),
      axis.title.y = element_text(size = 10, margin = margin(r = 15)),
      axis.text.y=element_blank(),
      axis.text.x=element_blank(),
      axis.line = element_line(colour = 'black', size = 0.2),
      axis.ticks = element_line(colour = 'black', size = 0.2),
      axis.ticks.length=unit(0.5, 'mm'),
      legend.key.height= unit(0.8, 'cm'),
      legend.key.width= unit(0.2, 'cm'),
      legend.title=element_text(size = 6),
      legend.text=element_text(size=6),
      strip.text = element_text(size = 6),
      panel.spacing.x=unit(0.1, 'lines'),
      panel.spacing.y=unit(0.1, 'lines'),
      plot.tag = element_text(size = 6),
      plot.tag.position = c(0.85, 0.015),
      plot.background = element_rect(colour ='black'))"
thlemon2<-paste0("theme(axis.text.x = element_text(size=5,angle = 65, vjust = 0.5,margin = margin(t = 2)),
                 axis.text.y = element_text(size=6,margin = margin(r = 2)))")
#plot
pls3d<-map(1:length(dtattrim),~{
  #.x<-1
  dta<-dtattrim[[.x]] 
  gg1<-paste0("ggplot(data=dta,aes(x =",vars[[.x]][1],", y =",vars[[.x]][2],",z=fitili))")
  #print(gg1)
  gg2<-"NULL"
  gg3<-"NULL"
  tvars<-unique(vars %>% unlist)
  idxx<-which(tvars%!in%vars[[.x]])
  vars.fix<-c("Number observations","Trap density","Area covered")
  fix<-vars.fix[idxx]
  fixv<-dta[[tvars[idxx]]] %>% unique
  if(grepl("nbobs",vars[[.x]][1])==T){
    gg2<-"scale_x_log10(breaks=c(2,5,10,20,40,80,160,320,600,1000,2000))"}
  if(grepl("area.c",vars[[.x]][2])==T){
    gg3<-"scale_y_continuous(breaks=c(0.2,0.4,0.6,0.8,1),labels=c('20','40','60','80','100'))"
  }else{
    fixv<-paste0(fixv*100,"%")}
  
  pl3d<-eval(parse(text=gg1))+
    eval(parse(text=gg2))+
    eval(parse(text=gg3))+
    geom_raster(aes(fill = fitili), interpolate=TRUE) +
    stat_contour(breaks=brks,color="black")+
    geom_dl(aes(z=fitili,label=..level..), method=list("bottom.pieces",cex=.45, hjust=1),
            stat="contour",breaks = brks) +
    stat_contour(aes(z=cip), breaks=brksci1,color="black",linetype="dotted")+
    stat_contour(aes(z=cip), breaks=brksci2,color="black",linetype="dotdash")+
    stat_contour(aes(z=cip), breaks=brksci3,color="black",linetype="dashed")+
    scale_fill_gradient2(low = hcl.colors(20,palette="Zissou 1")[1],mid="yellow",high = hcl.colors(20,palette="Zissou 1")[20],
                         midpoint=1,labels=c(expression(""<=0),'50','100','150',expression("">=200)),limits=c(0,2))+
    facet_rep_grid(estimator~model+sampling)+
    lemon::coord_capped_cart(bottom='both') +
    theme_classic()+
    labs(fill = "Estimated \narea (%)",x=labss[[.x]][1],y=labss[[.x]][2],tag = paste0(fix," = ",fixv),size= 4)
  
    pl3d<-lemonworkarround(pl3d,theme1=thlemon1,theme2=thlemon2,nr=5,ncl=5)
    
  pl3d
  
})

pls3d.save <- lapply(pls3d, ggplotGrob)


# ggsave(paste0(pathplots,"area_nbobs_trapd_areac_3d_","all",".pdf"), 
#        marrangeGrob(pls3d.save, nrow = 1, ncol = 1, top=NULL),
#        height=(doublecolumn[2]/12)*9,width=doublecolumn[2],units="px",dpi=600)
#tiff only saves last plot
# walk(plotsaveformat,~{
#   ggsave(plot=marrangeGrob(pls3d.save, nrow = 1, ncol = 1, top=NULL),
#          file=paste0(pathplots,"area_nbobs_trapd_areac_3d_","all",".",.x),
#          device=.x,
#          height=(doublecolumn[2]/12)*9,width=doublecolumn[2],
#          units="px",dpi=600)
# })

figs_345<-pls3d.save[c(1,4,7)]


# ggsave(paste0(pathplots,"figures_3_4_5",".pdf"), 
#        marrangeGrob(figs_456, nrow = 1, ncol = 1, top=NULL),
#        height=(doublecolumn[2]/12)*9,width=doublecolumn[2],units="px",dpi=600)

walk(plotsaveformat,function(format){
  walk2(figs_345,1:length(figs_345),~{
    ggsave(plot=.x,#marrangeGrob(figs_345, nrow = 1, ncol = 1, top=NULL),
           file=paste0(pathplots,"figure_",(.y+2),".",format),
           device=format,
           height=(doublecolumn[2]/12)*9,width=doublecolumn[2],
           units="px",dpi=600)
    
  })

})

#B.Supplementary information-------------
dina4.600.h<-3508*2#pixels at 600 dpi
dina4.600.w<-2480*2
#texts for figures--------
#load figure captions etc & convert to grob for plotting with graphs and save as pdf
filenames<-list.files(paste0(pathplots,"SI/"),pattern="texts")
texts.si<-readLines(paste0(pathplots,"SI/",filenames[[1]]), encoding = "UTF-8")
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

#B.1 interactions between nbobs, trap.d and area.c (additional) ---------
si_figs_2_to_7<-pls3d.save[!c(1:length(pls3d.save))%in%c(1,4,7)]
#intertwin texts and graphs
pls.save<-pmap(list(pl=si_figs_2_to_7,txt=texts.grob[1:6]),
               function(pl,txt){list(pl,txt)}) %>% unlist(.,recursive=FALSE)
#save
# ggsave(paste0(pathplots,"si_figures_1_to_6",".pdf"), 
#        marrangeGrob(pls.save, nrow = 2, ncol = 1,heights=c(.6,.4),widths=c(1),top=NULL),
#        height=dina4.600.h,width=dina4.600.w,units="px",dpi=600)

walk(plotsaveformat,~{
  ggsave(plot=marrangeGrob(pls.save, nrow = 2, ncol = 1,heights=c(.6,.4),widths=c(1),top=NULL),
         file=paste0(pathplots,"SI/","si_figures_2_to_7",".",.x),
         device=.x,
         height=dina4.600.h,width=dina4.600.w,
         units="px",dpi=600)
})

#B.2. GAM diagnostics---------
#B.2.a Basis dimensions sufficiency------
kcheck<- gam.checking(gmodel,k.sample = 5000, k.rep = 200)
rncheck<-rownames(kcheck);rownames(kcheck) <- NULL
kcheck.df<-cbind(term=rncheck,kcheck %>% as.data.frame(stringsAsFactors=F)) 
kcheck.tb <- reporter::create_table(kcheck.df) %>% reporter::titles("SI Table 1: Basis dimensions sufficiency tests for GAM")
kcheck.rpt <- reporter::create_report(orientation = "portrait", font="Arial",font_size=8,paste0(pathplots,"si_table_1.pdf"),output_type = "PDF") %>% 
  reporter::add_content(kcheck.tb)
reporter::write_report(kcheck.rpt)

#B.2.b Residuals----------- 
restrial<-residuals.gam(gmodel)
fittrial<-fitted(gmodel)
dtm<-datagam %>% mutate(residuals=restrial,
                    #residuals.cor=residuals/(lg10nbobs*ssizee),#correct for weighted regression
                    fitted=fittrial,
                    model=factor(model,levels=asympmodels,labels=c("Raw","MicMen","MonMol")),
                    estimator=factor(estimator,levels=estimators,labels=c("AKDE","KDE","MCP","LoCoH","BicubIt")),
                    sampling=factor(sampling,levels=sampregs,labels=c("Time-ordered","Distance-ordered")))
#reduce datapoints for tractability
set.seed(0620)
redusamp<-4*2
srws<-sample(size=nrow(dtm)/redusamp, x=1:nrow(dtm), replace=FALSE) %>% sort

#create all residual plots in one go
brksx<-c(2,5,10,20,40,80,160,320,600,1000,2000)
vars.resi<-c("residuals")#,"residuals.cor")
vars.other<-c("fitted","nbobs","trap.d","area.c","qq","density")
vr.labs<-c("Residuals")#,"Weighted residuals")
vo.labs<-c("Fitted","Number of observations","Trap density","Area covered","Theoretical quantiles","Density")
varss<-expand.grid(vars.resi,vars.other) %>% set_names(c("resi","other"))
labss<-expand.grid(vr.labs,vo.labs) %>% set_names(c("resi","other"))

pls.resid<-map(1:nrow(varss),function(.x){
  gg0<-"ggplot(data=dtm[srws,])";gg1<-gg2<-gg3<-"NULL"
  if(varss[.x,"other"]%in%vars.other[1:4]){gg1<-paste0("geom_point(aes(x=",varss[.x,"other"],",y=",varss[.x,"resi"],"),size=.5,alpha=0.2)");gg3<-"geom_hline(yintercept=0,size=.2)"}
  if(varss[.x,"other"]=="nbobs"){gg2<-"scale_x_log10(breaks=brksx)"}
  if(varss[.x,"other"]=="qq"){gg0<-paste0("ggplot(data=dtm[srws,],aes(sample=",varss[.x,"resi"],"))");gg2<-"stat_qq(size=.5)";gg3<-"stat_qq_line()"}
  if(varss[.x,"other"]=="density"){gg0<-paste0("ggplot(data=dtm[srws,],aes(x=",varss[.x,"resi"],"))");gg2<-"geom_density()";gg3<-"geom_vline(xintercept=0,size=.2)"}
  pl<-eval(parse(text=gg0))+ eval(parse(text=gg1))+ eval(parse(text=gg2))+ eval(parse(text=gg3))+
    lemon::facet_rep_grid(estimator~model+sampling,scales="free")+
    theme_classic()+ 
    #theme(axis.text.x = element_text(size=7,angle = 25, vjust = 0.5))+
    labs(x=labss[.x,"other"],y=labss[.x,"resi"])
  
  pl<-lemonworkarround(pl,theme1=thlemon1,theme2=thlemon2,nr=5,ncl=5)
  pl
  
})

#convert graphs to grobs
pls.resid.grob <- lapply(pls.resid, ggplotGrob)

#intertwin texts and graphs
pls.save<-pmap(list(pl=pls.resid.grob,txt=texts.grob[8:13]),
               function(pl,txt){list(pl,txt)}) %>% unlist(.,recursive=FALSE)
# pls.save<-pmap(list(pl=pls.resid.grob,txt=texts.grob[7:10],idx=1:length(pls.resid)),
#                function(pl,txt,idx){if((idx %% 2) == 0){list(pl,grid::nullGrob())}else{list(pl,txt,grid::nullGrob(),grid::nullGrob())}}) %>% unlist(.,recursive=FALSE)
#save
# ggsave(paste0(pathplots,"si_figures_7_to_12",".pdf"), 
#        marrangeGrob(pls.save, nrow = 2, ncol = 1,heights=c(.6,.4),widths=c(1),top=NULL),
#        height=dina4.600.h,width=dina4.600.w,units="px",dpi=600)

walk(plotsaveformat,~{
  ggsave(plot=marrangeGrob(pls.save, nrow = 2, ncol = 1,heights=c(.6,.4),widths=c(1),top=NULL),
         file=paste0(pathplots,"SI/","si_figures_9_to_14",".",.x),
         device=.x,
         height=dina4.600.h,width=dina4.600.w,
         units="px",dpi=600)
})

#B.3 Asymptotic models' fit on data------------
#si fig 13
#load raw data
filenames<-list.files(pathcomplete,pattern="estimators_df_all")
load(paste0(pathcomplete,filenames[[1]]))

#what were the sampling ticks for nbobs

nbobsvals<-datagam$nbobs %>% unique %>% sort

#subset the data for a representative plot
estimators_ac1td12<-estimators_df %>% filter(area.c==1,trap.d==12)
estimators_ac1td12<-map2_dfr(estimators_ac1td12$raws,estimators_ac1td12$id,~.x %>% mutate(id=.y))

#long format
estimators_ac1td12_l<-estimators_ac1td12 %>% pivot_longer(cols=all_of(c(estimators)),values_to="area",names_to="estimator") 
#adapt to create a group variable in the plots for apropriate line drawing & convert to long long format
estimators_ac1td12_l<-estimators_ac1td12_l %>% 
  mutate(area=pmap(list(.x=area,.y=nbobs,.z=sampling),function(.x,.y,.z){
    nbobsh<-case_when(.z=="fps"~nbobsvals[which(nbobsvals<=.y)],TRUE~.y)
    dt<-tibble(area=.x %>% unlist,nbobs=nbobsh,nbobstot=.y)
    dt})) %>% select(-nbobs) %>% unnest(area)

#calculate mean areas pr condition
estimators_ac1td12_ls<-estimators_ac1td12_l %>% 
  mutate(model=factor(model,levels=asympmodels,labels=c("Raw","MicMen","MonMol")),
         estimator=factor(estimator,levels=estimators,labels=c("AKDE","KDE","MCP","LoCoH","BicubIt")),
         sampling=factor(sampling,levels=sampregs,labels=c("Time-ordered","Distance-ordered"))) %>% 
  group_by(sampling,estimator,model,nbobstot,nbobs) %>% 
  summarise(sd=sd(area,na.rm=T)/true.area,area=mean(area,na.rm=T)/true.area,.groups="keep") %>% ungroup 

#Fit the asymptotic models to the mean areas and recover the predictions (not only asymptotes)
estimators_ac1td12_am<-estimators_ac1td12_ls %>% 
  group_by(sampling,model,estimator) %>% 
  group_map(~{
    smp<-unique(.x$sampling) %>% as.character
    if(smp=="Distance-ordered"){
      #fit the asympmodels for each different curve generated with farthest point sampling
      .x<-.x %>% group_by(nbobstot) %>% group_split
      .x<-map_dfr(.x,function(grp){grp<-fit_asympmodels2(data=grp,nbobstotval=grp[nrow(grp),][["nbobs"]]);grp}) %>% ungroup
      .x
    }else{
      .x<-map_dfr(1:nrow(.x),function(rw){dt<-fit_asympmodels2(data=.x[1:rw,],nbobstotval=.x[rw,][["nbobs"]]);dt})
      .x}
    .x},.keep=T) %>% bind_rows

#add model as a variable to facet the plot according to it
estimators_ac1td12_am_l<-estimators_ac1td12_am %>% 
  pivot_longer(cols=all_of(c("micmen","monmol")),names_to="asympmodel",values_to="asympmodelval") %>% 
  mutate(asympmodel=factor(asympmodel,levels=c("micmen","monmol"),labels=c("MicMen","MonMol")))


#recover the asymptotes from simulations 
#(already available from simulations unlike the asymptotic model predictions calculated above)
filenames<-list.files(pathcomplete,pattern="results_all_long")
load(paste0(pathcomplete,filenames[[1]]))
results_subs<-results %>% filter(area.c==1,trap.d==12) %>% filter(model!="Raw") %>% 
  group_by(sampling,estimator,model,nbobs) %>% 
  summarise(sd=sd(area,na.rm=T),area=median(area,na.rm=T),.groups="keep") %>% ungroup 

#plot together
pl.asympmodels.fit<-ggplot()+
  geom_hline(yintercept=1,color="black",size=.2)+#true area
  #time-ordered data
  geom_point(data=estimators_ac1td12_am_l%>% filter(sampling=="Time-ordered"),
             aes(x=nbobs,y=area,color=nbobs),alpha=0.3,size=1.5)+
  geom_line(data=estimators_ac1td12_am_l%>% filter(sampling=="Time-ordered"),
            aes(x=nbobs,y=asympmodelval,color=nbobstot,group=nbobstot),alpha=0.7,size=.5)+
  #distance-ordered
  geom_point(data=estimators_ac1td12_am_l%>% filter(sampling=="Distance-ordered"),
             aes(x=nbobs,y=area,group=nbobstot,color=nbobstot),alpha=0.3,size=1.5)+
  geom_line(data=estimators_ac1td12_am_l %>% filter(sampling=="Distance-ordered"),
            aes(x=nbobs,y=asympmodelval,group=nbobstot,color=nbobstot),alpha=0.7,size=.5)+
  paletteer::scale_color_paletteer_c("ggthemes::Classic Orange-Blue",direction=-1,trans="log10")+
  scale_x_log10(breaks=c(2,5,10,20,40,80,160,320,600,1000,2000))+
  facet_rep_grid(estimator~asympmodel+sampling)+#
  coord_cartesian(ylim=c(0,2.5))+
  theme_classic()+
  labs(color = "Partial \ndataset \nlength",
       x="Number of observations",
       y="Estimated area / True home range area",
       tag = paste0("Trap density = 12 | Area covered = 1"),
       size= 4)
#themes for plots because of bug in lemon package

pl.asympmodels.fit<-lemonworkarround(pl.asympmodels.fit,theme1=thlemon1,theme2=thlemon2,nr=5,ncl=4)

#convert graphs to grobs
pls.asympmodels.fit.grob <- lapply(pl.asympmodels.fit %>% list, ggplotGrob)

#intertwin texts and graphs
pls.save<-pmap(list(pl=pls.asympmodels.fit.grob,txt=texts.grob[14]),function(pl,txt){list(pl,txt)}) %>% unlist(.,recursive=FALSE)

#save
# ggsave(paste0(pathplots,"si_figure_15",".pdf"), 
#        marrangeGrob(pls.save, nrow = 2, ncol = 1,heights=c(.6,.4),widths=c(1),top=NULL),
#        height=dina4.600.h,width=dina4.600.w,units="px",dpi=600)
walk(plotsaveformat,~{
  ggsave(plot=marrangeGrob(pls.save, nrow = 2, ncol = 1,heights=c(.6,.4),widths=c(1),top=NULL),
         file=paste0(pathplots,"SI/","si_figure_15",".",.x),
         device=.x,
         height=dina4.600.h,width=dina4.600.w,
         units="px",dpi=600)
})

#C. GAM predictions fit the data---------
#C.1 Predict values---------
#create the combinations needed
vars<-c("lg10nbobs","trap.d","area.c")
fixedvars<-combn(vars, 2, simplify = T) %>% t %>% as_tibble %>% 
  rowwise %>% mutate(V3=vars[vars%!in%c(V1,V2)]) %>% ungroup
fixedvals<-list(lg10nbobs=c(7,100,1000) %>% log10 %>% as.character,
                trap.d=c("2","22","42"),
                area.c=c("1","0.6","0.2"))
varyingvals<-list(
  lg10nbobs=seq(log10(2),log10(2000),length.out=15),
  trap.d=seq(min(datagam$trap.d),max(datagam$trap.d),length.out=4) %>% round,
  area.c=seq(min(datagam$area.c),max(datagam$area.c),length.out=4) %>% round(digits=2))


#create the combinations needed
fixedvars<-fixedvars %>% 
  mutate(idfix1=1:nrow(.),
         V=pmap(list(idfix1=idfix1,x=V1,y=V2,z=V3),
                function(idfix1,x,y,z){
  res<-expand.grid(fixedvals[[paste0(x)]],
                   fixedvals[[paste0(y)]],stringsAsFactors = F) %>% 
    set_names(c(paste0(x),paste0(y))) 
  res<-res %>% mutate(varvary=z,idfix1=idfix1,idfix2=1:nrow(res))
  res[[paste0(z)]]<-list(varyingvals[[paste0(z)]])
  res[["approach"]]<-list(levels(gmodel$model$approach))#levels(gmodel$model$approach)-------
  res<-res  %>% unnest(paste0(z)) %>% unnest(approach)#approach------
  res<-res %>% select(varvary,idfix1,idfix2,approach,all_of(vars)) #approach------%>% as.data.frame
  res<-cbind(map_dfc(res %>% select(-approach,-varvary),as.numeric),varvary=res$varvary,approach=res$approach)#reduce(list(map_dfc(res %>% select(-approach),as.numeric),approach=res$approach,idsim=res$idsim),cbind)
  res}))

topred<-fixedvars$V %>% bind_rows %>% group_by(idfix1,idfix2) %>% group_split

#predict with adapted pred.gam function (takes ages with so much data)
# preds2d<-pred_gam(gmodel,topred,specials=c("[()]","log10"),#log10 is for the case when log10(nbobs) was done in the model formula
#                   varint=list(c("approach","lg10nbobs","area.c","trap.d")),#list(c("approach","lg10nbobs","area.c","trap.d")),
#                   nms=names(topred))
#save(preds2d,file=paste0(pathcomplete,"preds2d.Rdata"))
load(file=paste0(pathcomplete,"preds2d.Rdata"))

#extract data
dtat<-map(preds2d,"pred")
#reconf
dtat<-map(dtat,~{
  .x$approach<-str_replace_all(.x$approach,"\\.","--")
  .x<-bind_cols(.x,colsplit(.x$approach, "--", names = c("sampling", "estimator", "model"))) %>% 
    mutate(cip=(upperfitili-lowerfitili)/fitili,
           model=factor(model,levels=asympmodels,labels=c("Raw","MicMen","MonMol")),
           estimator=factor(estimator,levels=estimators,labels=c("AKDE","KDE","MCP","LoCoH","BicubIt")),
           sampling=factor(sampling,levels=sampregs,labels=c("Time-ordered","Distance-ordered")),
           nbobs=10^lg10nbobs %>% round)
  .x})

#trim extreme values to be able to observe interesting variation with color
dtattrim<-map(dtat,~{
  .x[which(.x$fitili<0),"fitili"]<-0
  .x[which(.x$fitili>2),"fitili"]<-2
  .x})

#need to know for each group which variable is the one to plot
varvary<-map_chr(dtat,~.x$varvary %>% unique)
varvary[which(varvary=="lg10nbobs")]<-"nbobs"
#create a list of variables etc to do all plots in the same function using map
#pairs of variables that change
vars<-list(c("nbobs","trap.d"),c("nbobs","trap.d"),c("nbobs","trap.d"),
           c("nbobs","area.c"),c("nbobs","area.c"),c("nbobs","area.c"),
           c("trap.d","area.c"),c("trap.d","area.c"),c("trap.d","area.c"))
#labs for the plots
labss<-list(c("Number of observations","Trap density"),c("Number of observations","Trap density"),c("Number of observations","Trap density"),
            c("Number of observations","Area covered"),c("Number of observations","Area covered"),c("Number of observations","Area covered"),
            c("Trap density","Area covered"),c("Trap density","Area covered"),c("Trap density","Area covered"))
yvar<-"fit"
vars<-c("nbobs","trap.d","area.c")
vars.fix<-c("Number observations","Trap density","Area covered")
xlbs<-c("Number of observations","Trap density","Area covered")


pls2d<-map(1:length(dtattrim),~{#length(dtattrim)
  dta<-dtattrim[[.x]]
  gg1<-paste0("ggplot(data=dta,aes(x =",varvary[.x],", y = ",yvar,"))")
  gg2<-"NULL"
  gg1sim<-paste0("geom_point(data=dtasim,aes(x =jitter(",varvary[.x],",amount=0.1)",", y = area","),alpha=0.2)")
  
  if(grepl("nbobs",varvary[.x])==T){
    gg2<-"scale_x_log10(breaks=c(2,5,10,20,40,80,160,320,600,1000,2000))"
    xlb<-xlbs[1]
    fixv1<-vars[2];fixv2<-vars[3]
    fixval1<-dta[[vars[2]]] %>% unique;fixval2<-dta[[vars[3]]] %>% unique
    fixlb<-paste0(xlbs[2]," = ",dta[[vars[2]]] %>% unique, " | ",
                  xlbs[3]," = ",dta[[vars[3]]] %>% unique)}
  if(grepl("trap.d",varvary[.x])==T){
    gg2<-"scale_x_continuous(breaks=c(2,12,22,32,42),labels=c('2','12','22','32','42'))"
    xlb<-xlbs[2]
    fixv1<-vars[1];fixv2<-vars[3]
    fixval1<-dta[[vars[1]]] %>% unique;fixval2<-dta[[vars[3]]] %>% unique
    fixlb<-paste0(xlbs[1]," = ",dta[[vars[1]]] %>% unique, " | ",
                  xlbs[3]," = ",dta[[vars[3]]] %>% unique)}
  if(grepl("area.c",varvary[.x])==T){
    gg2<-"scale_x_continuous(breaks=c(0.2,0.4,0.6,0.8,1),labels=c('20','40','60','80','100'))"
    xlb<-xlbs[3]
    fixv1<-vars[1];fixv2<-vars[2]
    fixval1<-dta[[vars[1]]] %>% unique;fixval2<-dta[[vars[2]]] %>% unique
    fixlb<-paste0(xlbs[1]," = ",dta[[vars[1]]] %>% unique, " | ",
                  xlbs[2]," = ",dta[[vars[2]]] %>% unique)}
  #as.character needed for point float issues with 0.6 area.covered
  dtasimeval<-paste0("datagam[which(as.character(datagam$",
                     fixv1,")==as.character(",fixval1,") & as.character(datagam$",
                     fixv2,")==as.character(",fixval2,")),]")
  
  dtasim<-eval(parse(text=dtasimeval)) %>% mutate(model=factor(model,levels=asympmodels,labels=c("Raw","MicMen","MonMol")),
                                                 estimator=factor(estimator,levels=estimators,labels=c("AKDE","KDE","MCP","LoCoH","BicubIt")),
                                                 sampling=factor(sampling,levels=sampregs,labels=c("Time-ordered","Distance-ordered")))
  
  
  pl2d<-eval(parse(text=gg1))+
    geom_hline(yintercept=1,color="red",size=.2,alpha=.5)+
    eval(parse(text=gg2))+
    geom_line() +
    geom_ribbon(aes(ymin=lowerfit, ymax=upperfit), alpha=0.2)+
    eval(parse(text=gg1sim))+
    facet_rep_grid(estimator~model+sampling)+
    coord_cartesian(ylim=c(0,2.5))+
    #coord_capped_cart(bottom='both') +
    theme_classic()+
    labs(y = "Estimated area / True home range area",
         x=xlb,
         tag = fixlb,
         size= 4)
  pl2d<-lemonworkarround(pl2d,theme1=thlemon1,theme2=thlemon2,nr=5,ncl=5)
  pl2d
  
})

pls2d.grob <- lapply(pls2d, ggplotGrob)

#intertwin texts and graphs
pls2d.save<-pmap(list(pl=pls2d.grob,
                      txt=rep(texts.grob[15],length(pls2d))),
               function(pl,txt){list(pl,txt)}) %>% unlist(.,recursive=FALSE)

#save
# ggsave(paste0(pathplots,"si_figures_16_to_44",".pdf"), 
#        marrangeGrob(pls2d.save, nrow = 2, ncol = 1,heights=c(.6,.4),widths=c(1),top=NULL),
#        height=dina4.600.h,width=dina4.600.w,units="px",dpi=600)
walk(plotsaveformat,~{
  ggsave(plot=marrangeGrob(pls2d.save, nrow = 2, ncol = 1,heights=c(.6,.4),widths=c(1),top=NULL),
         file=paste0(pathplots,"SI/","si_figures_16_to_44",".",.x),
         device=.x,
         height=dina4.600.h,width=dina4.600.w,
         units="px",dpi=600)
})
#D NAs in home range size calculation--------
#si fig 43
filenames<-list.files(pathcomplete,pattern="results_all_long")
load(paste0(pathcomplete,filenames[[1]]))

results<-results %>% mutate(approach=interaction(sampling,estimator,model),
                            model=factor(model,levels=asympmodels,labels=c("Raw","MicMen","MonMol")),
                            estimator=factor(estimator,levels=estimators,labels=c("AKDE","KDE","MCP","LoCoH","BicubIt")),
                            sampling=factor(sampling,levels=sampregs,labels=c("Time-ordered","Distance-ordered")))


results.s<-results%>% group_by(approach) %>% 
  summarise(na=sum(is.na(area)==T),
            nona=sum(is.na(area)==F),
            pna=na/(na+nona)) %>% dplyr::arrange(pna)

#D.1 % NAs overall---------------
results<-results %>% mutate(nas=is.na(area) %>% as.numeric) 
pna.all<-(which(results$nas==1) %>% length)/(results$nas %>% length)

#D.2 % NAs per type nbobs*trapd, all area.c confounded--------------------
#si figure 18
results.mna<-results %>% 
  group_by(sampling,estimator,model,nbobs,trap.d) %>% summarise(mnas=mean(nas),.groups="keep") %>% ungroup

pl.NA<-ggplot(data=results.mna,aes(x=nbobs,y=trap.d,z=mnas)) + 
  scale_x_log10(breaks=brksx)+
  geom_raster(aes(fill = mnas), interpolate=TRUE) +
  facet_rep_grid(estimator~model+sampling)+
  coord_capped_cart(bottom='both') +
  theme_classic()+
  xlab("Number of observations")+
  labs(fill = "Propotion of NA",
       x="Number of observations",
       y="Trap density",
       size= 4)
pl.NA<-lemonworkarround(pl.NA,theme1=thlemon1,theme2=thlemon2,nr=5,ncl=5)
pl.NA.grob <- lapply(pl.NA %>% list, ggplotGrob)

#intertwin texts and graphs
pl.NA.save<-pmap(list(pl=pl.NA.grob,
                      txt=texts.grob[16]),
                 function(pl,txt){list(pl,txt)}) %>% unlist(.,recursive=FALSE)

#save
# ggsave(paste0(pathplots,"si_figure_45",".pdf"), 
#        marrangeGrob(pl.NA.save, nrow = 2, ncol = 1,heights=c(.6,.4),widths=c(1),top=NULL),
#        height=dina4.600.h,width=dina4.600.w,units="px",dpi=600)

walk(plotsaveformat,~{
  ggsave(plot=marrangeGrob(pl.NA.save, nrow = 2, ncol = 1,heights=c(.6,.4),widths=c(1),top=NULL),
         file=paste0(pathplots,"SI/","si_figure_45",".",.x),
         device=.x,
         height=dina4.600.h,width=dina4.600.w,
         units="px",dpi=600)
})
#END-----