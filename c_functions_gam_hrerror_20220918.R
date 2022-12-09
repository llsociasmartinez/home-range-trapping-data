library(mgcv)
library(gratia)
library(reporter)#export table gam.check


"%!in%"<-Negate("%in%")
#formulas for GAM---------
form.ti<-as.formula(area~
                      approach
                    +s(lg10nbobs,k=knbobs,bs=spl[1],by=approach)
                    +s(trap.d,k=5,bs=spl[2],by=approach)
                    +s(area.c,k=5,bs=spl[3],by=approach)
                    +ti(lg10nbobs,trap.d,k=c(knbobs,5),bs=spl[c(1,2)],by=approach)
                    +ti(lg10nbobs,area.c,k=c(knbobs,5),bs=spl[c(1,3)],by=approach)
                    +ti(trap.d,area.c,k=c(5,5),bs=spl[c(2,3)],by=approach)
                    +ti(lg10nbobs,area.c,trap.d,k=c(knbobs,5,5),bs=spl[c(1,2,3)],by=approach)
)


#function from ricardo in 
#https://stackoverflow.com/questions/5577221/how-can-i-load-an-object-into-a-variable-name-that-i-specify-from-an-r-data-file
loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  mget(ls()[ls() != "fileName"])
}

#prepare the datasets for predicting with pred_gam
prep_data<-function(data,varsfac,varscont,varint){
  allvars<-c(varsfac,varscont)
  
  #create list to expand
  var.sq<-map(allvars,~{data[[.x]] %>% unique %>% sort}) %>% set_names(allvars)
  #list of means
  var.m<-map(varscont,~{data[[.x]]  %>% mean(na.rm=T)}) 
  var.mm<-var.sq
  lv<-length(allvars)
  var.mm[(length(varsfac)+1):lv]<-var.m
  
  #which are of no interest?
  #substitute with mean
  var.sq[which(allvars%!in%varint)]<-var.mm[which(allvars%!in%varint)]
  
  p.data<-expand.grid(var.sq)
  
}

#predict with new data
#this function was inspired by Gavin Simpson's functions to predict using the lp matrix
#in https://fromthebottomoftheheap.net/
pred_gam <- function(modell, newdata, varint, specials, nms,alpha = 0.05,
                  unconditional = FALSE,exclude=NA) {
  vars<-varint %>% unlist %>% unique
  preds<-map2(newdata,varint,function(data,vint){
    #idxmdl<-1
    # modell<-gmodel
    # data<-trial[[1]] %>%
    #   mutate(approach=interaction(sampling,estimator,model))#,
    #          #`log10(nbobs)`=log10(nbobs))
    # data<-data[which(data$approach%in%modell$model$approach),]
    # data$approach<-factor(data$approach,levels=levels(modell$model$approach))
    # 
    # nms<-c("area.c","trap.d","nbobs")
    # varint<-list(c("approach","area.c"),c("approach","trap.d"),c("approach","nbobs"))
    # vint<-varint[[3]]
    # specials<-c("[()]","log10")
    #nm<-names(model)[idxmdl]
    #modell<-model[[idxmdl]]
   
    "%!in%"<-Negate("%in%")
    trms<-modell$terms
    terms<-attr(trms,"term.labels")
    for(i in specials){
      terms<-gsub(i,"",terms)
    }
    terms
    termsn<-paste(terms,collapse="|")
    nointerms<-terms[which(terms%!in%vint)] 
    nointerms
    #nointerms
    xp <- predict.gam(object=modell, newdata=data, type = "lpmatrix")
    
    #find columns not interest
    cno<-NULL
    if(length(nointerms)>0){
      cno<-grep(paste(nointerms,collapse="|"), 
                colnames(xp),value=TRUE) 
      
    }
    cno
    #modell$model$approach %>% unique
    library(gtools)
    #create combinations for interactions
    combos<-map(1:length(vint),~combinations(n = length(vint), r = .x, v = vint, repeats.allowed = FALSE))
    
    colsint<-map(1:length(combos),
                 ~map(1:nrow(combos[[.x]]),function(nr){
                   #.x=2
                   #nr=1
                   #combos[[.x]][nr,]
                   # .x=1
                   # nr=2
                   matches<-c(combos[[.x]][nr,])
                   
                   notmatch<-which(vars%!in%matches)
                   tomatch<-which(vars%in%matches)
                   #need to find those columns that contain the notmatch
                   tf<-lapply(vars, grepl, colnames(xp))%>% set_names(vars)
                   tf<-bind_cols(tf) 
                   #those that match need to be all in at the same time while those that do not match might be in 
                   #match those cols where tomatch is all True and notmatch is all false
                   #cno#those that we dont want
                   cls<-intersect(which(rowSums(tf[,tomatch]) == length(tomatch)),which(rowSums(tf[,notmatch])==0L))
                   cls<-cls[cls%!in%cno]
                   #cls
                   #xp[,cls] %>% colnames
                   return(cls)
                 }) 
    ) %>% unlist %>% unique %>% sort
    colnames(xp[,colsint])
    #colnames(xp)
    ## zero out cols of X related to splines for other variables
    #do not exclude intercept
    xp[, generics::setdiff(2:ncol(xp),colsint)] <- 0
    
    ## zero out the parametric cols
    fit <- xp %*% coef(modell)
    se <- sqrt(rowSums((xp %*% vcov(modell, unconditional = unconditional)) * xp))
    crit <- qt(alpha/2, df.residual(modell), lower.tail = FALSE)
    upr <- fit + (crit * se)
    lwr <- fit - (crit * se)
    
    fitili<-modell$family$linkinv(fit)
    uprili <- modell$family$linkinv(upr)
    lwrili <- modell$family$linkinv(lwr)
    pred <- reduce(list(data,data.frame(fit=fit,upperfit=upr,lowerfit=lwr,fitili=fitili,upperfitili=uprili,lowerfitili=lwrili)),
                   bind_cols) 
    #pred %>% View
    return(list(pred=pred,xp=xp))})
  
  #print(str(preds,max.level=1))
  preds<-preds%>% set_names(nms)
  return(preds)
  
}

#from greengrass62 answer
#https://stackoverflow.com/questions/42042822/how-to-save-edf-from-mgcvgam-check-and-skip-the-plots
gam.checking <- function(b, k.sample = 5000, k.rep = 200) {
  mgcv:::k.check(b, subsample = k.sample, n.rep = k.rep)
}

#END------------
