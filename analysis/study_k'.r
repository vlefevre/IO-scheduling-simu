###################### BEGIN R STUFF ##########################
library(reshape2)
library(data.table)
library(ggplot2)
library(scales) 
library(grid)
library(plyr)
###################### END R STUFF ##########################

#########################
## Parameters that can be modified for the plot.
#########################


## We plot between n_kp_min and n_kp
n_kp_min =1 
n_kp =15 

n_expes = 10 #number of expes

img="evol_k'.pdf" #location where image is saved

#########################
## Plotting impact of k'
#########################


#N_XX is for XX normalized with respect to the best one attainable (when K' is maximum).
DF <- data.frame(Kprime=rep(NA,0),id=rep(NA,0),N_SysEff=rep(NA, 0),N_Dilation=rep(NA, 0),stringsAsFactors=FALSE)

for (j in 1:n_expes) {
	x <-paste("study_kp/simu_all_set",j,sep = "")
	dat = read.table(x, header = TRUE)
	#We compute K'. To do this we take only three significant numbers (otherwise there are roundup errors, even with this, but we filter them later).
	setDT(dat)[]

	#Next step is to compute the max as a function of Kprime
	df <- data.frame(
			SysEff=rep(NA, nrow(dat)),
			Dilation=rep(NA, nrow(dat)),
			Kprime=rep(NA, nrow(dat)),
			id=rep(NA,nrow(dat)),
			stringsAsFactors=FALSE) 

	for (i in 1:nrow(dat)) {
		if(i == 1) {
			df[1,] <- c(dat[1,SysEff],dat[1,Dilation],dat[1,Kprime],j)
		} else {
			if(dat[i,SysEff] > df[i-1,"SysEff"]) {
				df[i,] <- c(dat[i,SysEff],dat[i,Dilation],dat[i,Kprime],j)
			} 
			if(dat[i,SysEff] < df[i-1,"SysEff"]) {
				df[i,] <- c(df[i-1,"SysEff"],df[i-1,"Dilation"],dat[i,Kprime],j)
			}
			if(dat[i,SysEff] == df[i-1,"SysEff"]) {
				df[i,] <- c(dat[i,SysEff],min(df[i-1,"Dilation"],dat[i,Dilation]),dat[i,Kprime],j)
			}
		}
	}
	setDT(df)[ , N_SysEff := SysEff/SysEff[nrow(dat)]]
	setDT(df)[ , N_Dilation := Dilation/Dilation[nrow(dat)]]

	#Removing columns useless for this plot (syseff+Dilation, we only need the normalized versions).
	dat <- subset( df, select = -c(SysEff, Dilation) )
	DF <-rbind(DF,dat)
}
#DF now contains everything :)


#We want to do averages over all sets for all values of K'.

dat_S <-  subset( DF, select = c(Kprime, N_SysEff) )
dat_D <-  subset( DF, select = c(Kprime, N_Dilation) )

# Run the functions length, mean, and sd on the value of "N_SysEff/N_Dilation" for each Kprime
tgc_S <- ddply(dat_S, c("Kprime"), summarise,
               N    = length(N_SysEff),
               mean = mean(N_SysEff),
               sd   = sd(N_SysEff),
               se   = sd / sqrt(N)
)
tgc_D <- ddply(dat_D, c("Kprime"), summarise,
               N    = length(N_Dilation),
               mean = mean(N_Dilation),
               sd   = sd(N_Dilation),
               se   = sd / sqrt(N)
)

# We now filter roundup errors: we only take the elements when we have all expe values (n_expes)
tgc_S_Filtered=subset(tgc_S, (N==n_expes & Kprime <n_kp & Kprime >= n_kp_min))
tgc_D_Filtered=subset(tgc_D, (N==n_expes & Kprime <n_kp & Kprime >= n_kp_min))




###NOW WE PLOT !!

pSysEff = ggplot(data=tgc_S_Filtered, aes(x=Kprime,y=mean))+
    geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1,colour="grey30",alpha=0.3) +#theme_bw()+
	xlab("")+scale_x_continuous(breaks=n_kp_min:n_kp)+
	ylab("Normalized SysEff")+
	geom_line(size=1, colour="#000099")+ geom_point(colour="#000099")




pDil = ggplot(data=tgc_D_Filtered, aes(x=Kprime,y=mean))+
    geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1,colour="grey30",alpha=0.3) +#theme_bw()+
	xlab("K' = Tmax / Tmin")+scale_x_continuous(breaks=n_kp_min:n_kp)+
	ylab("Normalized Dilation")+
	geom_line(size=1, colour="#FF9999")+ geom_point(colour="#FF9999")


#This is to print both figures as one, sharing xlab. 
#To print individually, remove all commented text above.
cairo_pdf(file =img, width=6, height=7)
pushViewport(viewport(layout = grid.layout(2, 1)))
print(pSysEff, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(pDil, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
dev.off()
