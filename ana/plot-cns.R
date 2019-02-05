library(copynumber)
library(reshape)
library(tools)

# load("../../lpWGS phylogenies/bin_locations_4401.Rdata")
load("./bin_locations_4401.Rdata")

args = commandArgs(trailingOnly=TRUE)
data.dir <- args[1]

if(0){
   d <- read.table("cnprofile-2.txt",header=FALSE)
   names(d) <- c("sample","chromosome","index","cn")
   nsample <- length(unique(d$sample))

   md <- melt(d, id=c("sample","chromosome","index"))
   data <- cast(md,chromosome+index~sample)
   data <- cbind(bins_4401, data)
   data <- data[,-c(3,4,5)]

   #par(ask=F)
   plotGenome(data=data, ylab="copy number", layout=c(2,3), sample=c(1:nsample), equalRange=FALSE,q=0, col="red" )
}

#dir <- "../sim-data/"
dir <- data.dir
files <- list.files(dir,"^sim\\-data\\-\\d+\\-cn")
for(f in files){
   stub <- file_path_sans_ext(f)
   fout <- paste(dir,"plot-",stub,".pdf",sep="")

   d <- read.table(paste(dir,f,sep=""),header=FALSE)
   names(d) <- c("sample","chromosome","index","cn")
   nsample <- length(unique(d$sample))

   md <- melt(d, id=c("sample","chromosome","index"))
   data <- cast(md,chromosome+index~sample)
   data <- cbind(bins_4401, data)
   data <- data[,-c(3,4,5)]

   #par(ask=F)
   pdf(fout, height=10,width=20)
   plotGenome(data=data, ylab="copy number", layout=c(2,3), sample=c(1:nsample), equalRange=FALSE,q=0, col="red" )
   dev.off()

}
