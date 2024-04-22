# Takes a samtools coverage file called "[sample]_coverage.txt" and a sample name from the command line and produces a coverage plot 
# Developed by Sophia Cameron-Christie and Edana Lord 22-11-2016
# Usage: Rscript make_coverage_plots1.2.R <sample> <chr to plot> <coordinate 1 (optional) <coordinate 2 (optional) >
# eg: Rscript make_coverage_plots1.2.R N1001 MT-human
# eg: Rscript make_coverage_plots2.0.R MS10560 "gi|251831106|ref|NC_012920.1|" 2000 3000


#get arguments from the command line
args<-commandArgs(TRUE)
if(length(args)==0){
	print("No sample specified; cannot make plot")
	stopifnot(length(args)>0)
} else if(length(args)==1){
	samp<-as.character(args[1])
	print("No chr specified; defaulting to MT")
	plotchr<-"MT"
} else if(length(args)==2){
	samp<-as.character(args[1])
	plotchr<-as.character(args[2])
} else if(length(args)==3){
	samp<-as.character(args[1])
	print("No chr specified; defaulting to MT")
	plotchr<-"MT"
	print("two coordinates required to plot specific region!")
} else if(length(args)==4){
	samp<-as.character(args[1])
	plotchr<-as.character(args[2])
	coors1<-as.numeric(args[3])
	coors2<-as.numeric(args[4])
}  else {
	print("Too many arguments given from the command line")
}

print(paste0("plotting sample ",samp))
print(paste0("plotting chromosome ",plotchr))
if(length(args)==4){
print(paste0("start coordinate ", coors1))
print(paste0("end coordinate ", coors2))
}

#get the coverage file

samp.filename<-paste0(samp,"_coverage.txt")
print(samp.filename)
sampcov <- read.table(samp.filename, header = FALSE)
colnames(sampcov)<-c("chr","phys","coverage")
sampcov<-sampcov[sampcov[,"chr"]==plotchr,]
if(length(args)==4){
	sampcov<-sampcov[(sampcov[,"phys"]>=coors1)&(sampcov[,"phys"]<=coors2),]
}

#get the percent and mean coverage

total.cov <- max(sampcov$phys)
percent.cov<-((total.cov-sum(sampcov$coverage==0, na.rm=T))/total.cov)*100
meancov<-mean(sampcov$coverage)

#make the plot

if(length(args)==4){
	plot.filename<-paste0(samp,"_coverage_percent-",round(percent.cov,digits=1),"_average-",round(meancov,digits=1),"_bp",coors1,"-",coors2,".pdf")
} else {
	plot.filename<-paste0(samp,"_coverage_percent-",round(percent.cov,digits=1),"_average-",round(meancov,digits=1),".pdf")
}
pdf(plot.filename,width=9,height=3)
plot(sampcov$phys,sampcov$coverage,type="h", main="Coverage Distribution",xlab="Genome Position",ylab="Read Depth")
#,ylim=c(coors1,coors2)
abline(h=meancov,col="red",lwd=1)
dev.off()

