#BD project 2 fuctions
#Christopher Hall 1aug23
install.packages('BiocManager')
install.packages('devtools')
BiocManager::install('flowCore')
devtools::install_github('hally166/flowSpectrum')

library(flowCore)
library(flowSpectrum)

similarityANDplotter<-function(fSet,maxchnl,negCTRL){
  #select pos events
  gate <- openCyto::gate_mindensity(fSet, channel = maxchnl) # fr is a flowFrame
  
  #gate pos
  samp<-exprs(fSet)
  samp2<-samp[samp[, maxchnl] > gate@min,]
  
  #calculate medians
  medians<-apply(samp2, 2, median)
  neg_median<-apply(exprs(negCTRL),2,median)
  #remove background
  medians2<-medians-neg_median
  medians2<-t(medians2)
  medians3<-medians2[,c(2:17,20:35 ,40:71)]
  #normalise
  normed<-medians3/max(medians3)
  normed
}

#function to wrap it all up
spectral_differences<-function(fluor,fSet,maxchnl,negCTRL,rangex){
  normalised_spectra<-fsApply(fSet, function(x) similarityANDplotter(x,maxchnl=maxchnl,negCTRL=negCTRL)) 
  
  #calculate similarity
  similarity_matrix<-proxy::simil(normalised_spectra,diag=TRUE)
  #similarity_matrix<-round(similarity_matrix,2)
  
  #save file
  png(filename = paste0(fluorophore,".png"), width = 2000, height = 2000)
  par(mfrow = c(2, 1))
  # create base scatter plot
  plot(factor(colnames(normalised_spectra),levels = as.list(colnames(normalised_spectra))), normalised_spectra[1,],type="p",pch=1,las=2,ylab='',xlab='')
  
  # overlay line plot 
  lines(factor(colnames(normalised_spectra),levels = as.list(colnames(normalised_spectra))), normalised_spectra[1,], col='green', lwd=1)
  lines(factor(colnames(normalised_spectra),levels = as.list(colnames(normalised_spectra))), normalised_spectra[2,], col='yellow', lwd=1)
  lines(factor(colnames(normalised_spectra),levels = as.list(colnames(normalised_spectra))), normalised_spectra[3,], col='red', lwd=1)
  lines(factor(colnames(normalised_spectra),levels = as.list(colnames(normalised_spectra))), normalised_spectra[4,], col='purple', lty=1)
  lines(factor(colnames(normalised_spectra),levels = as.list(colnames(normalised_spectra))), normalised_spectra[5,], col='orange', lty=1)
  lines(factor(colnames(normalised_spectra),levels = as.list(colnames(normalised_spectra))), normalised_spectra[6,], col='pink', lty=1)
  
  legend("topright", legend=c(rownames(normalised_spectra)),
         col=c("green", "yellow", "red", "purple", "orange","pink"), lty=1, cex=0.8)
  title(paste(fluor,"Normalised Spectra"))
  
  #plot differences
  differences<-sweep(normalised_spectra, 2, normalised_spectra[1,], "-")
  plot(factor(colnames(differences),levels = as.list(colnames(differences))), differences[1,],type="l",pch=1,las=2,ylab='',xlab='',ylim=c(-rangex,rangex)) #you could play with round() here
  lines(factor(colnames(differences),levels = as.list(colnames(differences))), differences[1,], col='green', lwd=1)
  lines(factor(colnames(differences),levels = as.list(colnames(differences))), differences[2,], col='yellow', lwd=1)
  lines(factor(colnames(differences),levels = as.list(colnames(differences))), differences[3,], col='red', lwd=1)
  lines(factor(colnames(differences),levels = as.list(colnames(differences))), differences[4,], col='purple', lty=1)
  lines(factor(colnames(differences),levels = as.list(colnames(differences))), differences[5,], col='orange', lty=1)
  lines(factor(colnames(differences),levels = as.list(colnames(differences))), differences[6,], col='pink', lty=1)
  
  legend("topright", legend=c(rownames(differences)),
         col=c("green", "yellow", "red", "purple", "orange","pink"), lty=1, cex=0.8)
  title(paste(fluor, "Spectral Differences"))
  dev.off()
  MASS::write.matrix(as.matrix(similarity_matrix),file = paste0(fluorophore,"_similarity.csv"),sep=",")
}

#christopher Hall
#whole analysis workflow

setwd("C:\\fcsfiles\\") #set this to your analysis folder

fluorophore<-"BUV496"

#data
fs<-read.flowSet(list.files("C:\\fcsfiles\BUV496", full.names = TRUE,pattern = ".fcs")) #play with REGEX
neg<-read.FCS("C:\\fcsfiles\\A7 Well_091.fcs")

#work out max channel
flowSpectrum::spectralplot(fs[[1]])
maxchnl="UV7-A"

#run the script
spectral_differences(fluor=fluorophore,fSet=fs,maxchnl=maxchnl,negCTRL=neg,rangex=0.2)



