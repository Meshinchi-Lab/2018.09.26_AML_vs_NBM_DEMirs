#Jenny Smith


#original Author: Emilia Lim 

#Puprpose: DE analysis of miRNA

setwd(file.path(TARGET, "RNA/miRNAseq/analysis/2018.09.26_AML_vs_NBM_DEMirs/"))

calcWilcoxTest = function(df,libsA,libsB,log=T,paired=F,aname="A",bname="B",dcores=F){
  #df = your df of expression values
  #libsA/libsB = vector of strings of colnames in either group A or group B
  #dcores = number of cores to use for analysis. "F" will use only one. 
  libsA = intersect(colnames(df),libsA)
  libsB = intersect(colnames(df),libsB)
  df = df[,c(libsA,libsB)]
  if(log){
    df[df<1] = 1
    df = log2(df)
  }
  libAidx = c(1:length(libsA))
  libBidx = c((length(libsA)+1):(dim(df)[2]))
  
  if(dcores){
    library(doMC)
    library(foreach)
    registerDoMC(dcores)
    wt_res = foreach(i = 1:nrow(df))%dopar%{
      values = as.numeric(df[i,])
      x = values[libAidx]
      y = values[libBidx]
      wt = wilcox.test(x,y,paired=paired)
      if(log){
        fc = mean(x)-mean(y)
        c(mean(values),mean(x),mean(y),fc,wt$p.value)
      }else{
        fc = log2(mean(x))-log2(mean(y))
        c(log2(mean(values)),log2(mean(x)),log2(mean(y)),fc,wt$p.value)
      }
    }
    wt_res = do.call(rbind,wt_res)
    rownames(wt_res) = rownames(df)
  }else{
    wt_res = apply(df,1,function(values){
      x = values[libAidx]
      y = values[libBidx]
      wt = wilcox.test(x,y,paired=paired)
      if(log){
        fc = mean(x)-mean(y)
        c(mean(values),mean(x),mean(y),fc,wt$p.value)
      }else{
        fc = log2(mean(x))-log2(mean(y))
        c(log2(mean(values)),log2(mean(x)),log2(mean(y)),fc,wt$p.value)
      }
    })
    wt_res = t(wt_res)
  }
  colnames(wt_res) = c("log2_base_mean",paste("log2_mean_",aname,"_n",length(libsA),sep=""), 
                       paste("log2_mean_",bname,"_n",length(libsB),sep=""), 
                       paste("log2_fold_change_",aname,"..",bname,sep=""),"p_val")
  wt_res = data.frame(wt_res)
  adj_p_val = p.adjust(wt_res$p_val,method="BH")
  wt_res = cbind(wt_res,adj_p_val)
}
