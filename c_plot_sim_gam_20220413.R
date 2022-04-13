Sys.setenv(LANG = "en")
memory.limit(size=100000000)
source(paste0(getwd(),"/scripts/a_functions_sim_20220413.R"))
source(paste0(getwd(),"/scripts/b_functions_gam_hrerror_20220413.R"))

library(stringr)
library(reshape2)
library(directlabels)
library(lemon)
library(gridExtra)


#load GAM
filenames<-list.files(path,pattern="gam_ti")
load(paste0(path,filenames[[1]]))

#load the data
filenames<-list.files(path,pattern="datagam")
load(paste0(path,filenames[[1]]))


#Figure 1------------
#plots articles number---------
articles<-read.csv(file=paste0(getwd(),"/files/table_articles_telemetry.csv"))
articles<-articles %>% mutate(search=factor(words_search,labels=c("a. trap AND home range",
                                                                  "b. trap NOT home range",
                                                                  "c. telemetry AND home range")))
pl_articles<-ggplot(articles)+
  geom_line(aes(x=year,y=nb_studies,linetype=search))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 25, vjust = 0.5),
        legend.position = c(0.2, 0.8))+
  labs(linetype = "Type of study",
       x="Year",
       y="Number of studies")
ggsave(paste0(getwd(),"/plots/articles_type",".pdf"),
       height=5,width=7)

#Figures 4 to 6 + supmat-----------
#accuracy and reliability of estimations in the different conditions

datagam<-dtm
nms<-c(paste0("nbobs_trapd_area.c_",c("1","0.6","0.2")),
       paste0("nbobs_areac_trap.d_",c("4","22","42")),
       paste0("trapd_areac_nbobs_",c("7","100","1000")))


topred<-map(nms,~{
  nmsspl<-str_split(.x,pattern="_")
  var<-nmsspl[[3]]
  val<-nmsspl[[4]]
  tocr<-
    list(
      intsem=levels(gmodel$model$intsem)
      nbobs=10^seq(log10(min(datagam$nbobs)),log10(2000),length.out=50)
      trap.d=seq(min(datagam$trap.d),max(datagam$trap.d),by=10)
      area.c=seq(min(datagam$area.c),max(datagam$area.c),by=0.2))
  tocr[[which(tocr==var)]]<-val
  tocr<-tocr %>% cross_df(.filter = NULL)
  return(tocr)
})

#predict with adapted pred.gam function (takes ages with so much data)
preds3d<-pred_gam(gmodel,topred,specials=c("[()]","log10"),
                  varint=list(c("intsem","nbobs","area.c","trap.d")),
                  nms=names(topred))
save(preds3d,file=paste0(getwd(),"/files/preds3d.Rdata"))
load(file=paste0(getwd(),"/files/preds3d.Rdata"))

#extract data
dtat<-map(preds3d,"pred")
#reconf
dtat<-map(dtat,~{
  .x$intsem<-str_replace_all(.x$intsem,"\\.","--")
  .x<-bind_cols(.x,colsplit(.x$intsem, "--", names = c("sampling", "estimator", "model")))
  .x<-.x %>% 
    mutate(cip=(upperfitili-lowerfitili)/fitili,
           sampling=factor(sampling,levels=levels(datagam$sampling)),
           model=factor(model,levels=levels(datagam$model)),
           estimator=factor(estimator,levels=levels(datagam$estimator)))
  .x})

#trim extreme values to be able to observe interesting variation with color
dtattrim<-map(dtat,~{
  .x[which(.x$fitili<0),"fitili"]<-0
  .x[which(.x$fitili>2),"fitili"]<-2
  .x})

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
#plot
pls3d<-map(1:length(dtattrim),~{
  dta<-dtattrim[[.x]] %>% mutate(model=factor(model,labels=c("Raw","MicMen","MonMol")),
                                 estimator=factor(estimator,labels=c("AKDE","KDE","MCP","LoCoH","BicubIt")),
                                 sampling=factor(sampling,labels=c("Time-ordered","Distance-ordered")))
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
    geom_dl(aes(z=fitili,label=..level..), method=list("bottom.pieces",cex=.5, hjust=1),
            stat="contour",breaks = brks) +
    stat_contour(aes(z=cip), breaks=brksci1,
                 color="black",linetype="dotted",)+
    stat_contour(aes(z=cip), breaks=brksci2,
                 color="black",linetype="dotdash",)+
    stat_contour(aes(z=cip), breaks=brksci3,
                 color="black",linetype="dashed",)+
    scale_fill_gradient2(low = hcl.colors(20,palette="Zissou 1")[1],
                         mid="yellow",
                         high = hcl.colors(20,palette="Zissou 1")[20],
                         midpoint=1,
                         labels=c(expression(""<=0),'50','100','150',expression("">=200)))+
    facet_rep_grid(estimator~model+sampling)+
    coord_capped_cart(bottom='both') +
    theme_classic()+
    labs(fill = "Estimated \narea (%)",
         x=labss[[.x]][1],
         y=labss[[.x]][2],
         tag = paste0(fix,"=",fixv),
         size= 4)+
    theme(axis.text.x = element_text(size=7,angle = 25, vjust = 0.5),
          plot.background = element_rect(colour ='black'),
          plot.tag = element_text(size = rel(1)),
          plot.tag.position = c(0.85, 0.005))
  pl3d
  
})
#pls3d[[1]]
pls3d.save <- lapply(pls3d, ggplotGrob)

ggsave(paste0(getwd(),"/plots/trial3d_","all",".pdf"), 
       marrangeGrob(pls3d.save, nrow = 1, ncol = 1),height=9,width=12)



#Suplementary information figures (gam diagnostics)------------
#residual vs. fitted plot
restrial<-residuals.gam(gmodel)
fittrial<-fitted(gmodel)
dtm<-dtm %>% mutate(residuals=restrial,
                    fitted=fittrial,
                    sampling=factor(sampling,levels=levels(datagam$sampling)),
                    model=factor(model,levels=levels(datagam$model)),
                    estimator=factor(estimator,levels=levels(datagam$estimator)))%>% 
  mutate(model=factor(model,labels=c("Raw","MicMen","MonMol")),
         estimator=factor(estimator,labels=c("AKDE","KDE","MCP","LoCoH","BicubIt")),
         sampling=factor(sampling,labels=c("Time-ordered","Distance-ordered")))

pl.residuals<-ggplot(data=dtm[sample(size=nrow(dtm),
                                     x=1:nrow(dtm),
                                     replace=FALSE),])+
  geom_point(aes(x=fitted,y=residuals),size=.5)+
  geom_hline(yintercept=0)+
  facet_wrap(estimator~model+sampling,scales="free")+
  theme_classic()+
  theme(strip.text.x = element_blank())+
  labs(x="Fitted",
       y="Residual")
ggsave(paste0(getwd(),"/plots/residuals_fitted",".pdf"),
       height=7,width=10)

pl.nbobsresiduals<-ggplot(data=dtm[sample(size=nrow(dtm),
                                     x=1:nrow(dtm),
                                     replace=FALSE),])+
  geom_point(aes(x=nbobs,y=residuals),size=.5)+
  scale_x_log10()+
  geom_hline(yintercept=0)+
  facet_wrap(estimator~model+sampling,scales="free")+
  theme_classic()+
  theme(strip.text.x = element_blank())+
  labs(x="Number of observations",
       y="Residual")
ggsave(paste0(getwd(),"/plots/residuals_nbobs",".pdf"),
       height=7,width=10)

pl.trapdresiduals<-ggplot(data=dtm[sample(size=nrow(dtm)/10,
                                          x=1:nrow(dtm),
                                          replace=FALSE),])+
  geom_point(aes(x=trap.d,y=residuals),size=.5)+
  #scale_x_log10()+
  facet_wrap(estimator~model+sampling,scales="free")+
  theme_classic()+
  theme(strip.text.x = element_blank())+
  labs(x="Trap density",
       y="Residual")
pl.trapdresiduals

pl.areacresiduals<-ggplot(data=dtm[sample(size=nrow(dtm)/10,
                                          x=1:nrow(dtm),
                                          replace=FALSE),])+
  geom_point(aes(x=area.c,y=residuals),size=.5)+
  #scale_x_log10()+
  facet_wrap(estimator~model+sampling,scales="free")+
  theme_classic()+
  theme(strip.text.x = element_blank())+
  labs(x="Area covered",
       y="Residual")
pl.areacresiduals


#residual plot #################
pl.residualsqq<-ggplot(data=dtm[sample(size=nrow(dtm)/10,
                                     x=1:nrow(dtm),
                                     replace=FALSE),],
                       aes(sample=residuals))+
  stat_qq(size=.5) + stat_qq_line()+
  theme_classic()+
  facet_wrap(estimator~model+sampling,scales="free")+
  theme_classic()+
  theme(strip.text.x = element_blank())+
  labs(x="Theoretical quantiles",
       y="Residual")
#pl.residualsqq
ggsave(paste0(getwd(),"/plots/residuals_qqplot",".pdf"),
       height=9,width=12)

#density plot #############################
pl.residualdensity<-ggplot(data=dtm[sample(size=nrow(dtm)/10,
                                       x=1:nrow(dtm),
                                       replace=FALSE),],
                       aes(x=residuals))+
  geom_density()+
  theme_classic()+
  facet_wrap(estimator~model+sampling,scales="free")+
  theme_classic()+
  theme(strip.text.x = element_blank())+
  labs(y="Density",
       x="Residual")
#pl.residualsqq
ggsave(paste0(getwd(),"/plots/residuals_density",".pdf"),
       height=9,width=12)

