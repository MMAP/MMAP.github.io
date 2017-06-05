
args<-commandArgs(TRUE)
# command line args
y<-read.table(args[1],stringsAsFactors=F,sep="\t")
str(y)
# make a filename for the saved object
w<-vector("list",length(rownames(y)))
#q()

rname<-rownames(y)
names(w)<-rname
for( i in 1:length(rownames(y))){
#num_loci,n,sey,markers,maf,scores,cov
#  w[[i]] <-as.list(y[i,])
  w[[i]] <-list( y[i,2],y[i,3],y[i,5],y[i,6],y[i,7])
  names(w[[i]])<-c( "n","sey","maf","scores","cov")
#  print(names(w[[i]]))
# get marker names
  markers<-unlist(strsplit(y[i,]$markers,"/"))
#  cat("\n***** markers\n") 
#  print(markers)

# -------------------------------------------------------
# maf
#  print(w[[i]]$maf)
  w[[i]]$maf<-as.numeric(unlist(strsplit(w[[i]]$maf,"/")))
  names(w[[i]]$maf)<-markers
#  cat("\n********** MAF w[i]\n") 

# -------------------------------------------------------
# scores
  w[[i]]$scores<-as.numeric(unlist(strsplit(w[[i]]$scores,"/")))
  names(w[[i]]$scores)<-markers
#  cat("\n********** SCORES w[i]\n") 
#  print(w[[i]]$scores)

# -------------------------------------------------------
#cov
#  cat("\n********** COV w[i]\n") 
  w[[i]]$cov<-matrix(as.numeric(unlist(strsplit(w[[i]]$cov,"/"))),nrow=y[i,]$num_loci,ncol=y[i,]$num_loci);
  rownames(w[[i]]$cov)<-markers
  colnames(w[[i]]$cov)<-markers
#  cat("\n**** dim \n") 
#  print(w[[i]]$num_loci)
#  cat("\n**** dim \n") 
#  print(dim(w[[i]]$cov))
#  print(w[[i]]$cov)

}
  attr(w,"family") <- "gaussian"
  class(w) <- "skatCohort"
newname<-args[2]
assign(newname,w)
print(newname)

#filename <- paste(args[3], ".csv", sep="") 
filename <- paste(args[3], ".Rdata",sep="") 
#save(w,file="oldshool.rdata")
save(list=newname,file=filename)
q()



