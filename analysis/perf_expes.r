library(reshape2)
library(data.table)
library(ggplot2)

#########################
##GOING THROUGH SYSTEM EFFICIENCY RESULTS
#########################

dat_syseff = read.table("expes_syseff", header = TRUE)
syseff_melt <- melt(dat_syseff,
    id.vars=c("Set"),
    variable.name="schedule",
    value.name="measurement"
)

#For System Efficiency, we want to normalize with respect to the Upper Limit.
setDT(syseff_melt)[ , perf := measurement/measurement[(schedule=="syseff_UL")], by=.(Set)]

#Then we will plot only the results and not the Upper Limit
syseff_clean <- syseff_melt[ which(schedule!='syseff_UL'),] 

#renaming algo names
SE_ori_names=c("syseff_periodic_expes","syseff_periodic_simu","syseff_online_expes","syseff_online_simu","syseff_cong_expes")
SE_legends=c("Periodic (expe)", "Periodic (simu)", "Online (expe)", "Online (simu)","Congestion")
#SE_colors=c("#999999", "#E69F00", "#56B4E9", "#009E73", "#CC79A7")

#Careful, colors and lines are in the order of the adta in the database!
SE_colors=c("#CC79A7","#E69F00","#CC79A7","#E69F00","#999999")
SE_lines=c("solid","solid", "dashed", "dashed","solid")



pl_cong<-ggplot(syseff_clean, aes(x = Set, y = perf, colour = schedule,linetype = schedule)) +
	xlab("Set")+scale_x_continuous(breaks=1:10)+
	ylab("System Efficiency / Upper bound")+
	theme(legend.title=element_blank())+
	geom_line(size=1)+ geom_point(aes(shape=schedule),size = 4)

pl_cong+
	scale_shape_discrete(name="",breaks=SE_ori_names,labels=SE_legends)+
	scale_colour_manual(name="",values=SE_colors,breaks=SE_ori_names,labels=SE_legends)+ 
	scale_linetype_manual(name="",values=SE_lines,breaks=SE_ori_names,labels=SE_legends)+ 
#legend inside plotting area
	theme(legend.justification=c(1,0), legend.position=c(1,0))
                         

ggsave("../fig/perf_comparison.pdf", width=8, height=5,device=cairo_pdf)

#########################
## NOW FOR DILATION
#########################

dat_dil = read.table("expes_dil", header = TRUE)
dil_melt <- melt(dat_dil,
    id.vars=c("Set"),
    variable.name="schedule",
    value.name="perf"
)


DIL_ori_names=c("dil_periodic_expes","dil_periodic_simu","dil_online_expes","dil_online_simu")
DIL_legends=c("Periodic (expe)", "Periodic (simu)", "Online (expe)", "Online (simu)")
#DIL_colors=c("#999999", "#E69F00", "#56B4E9", "#009E73")

#Careful, colors and lines are in the order of the adta in the database!
DIL_colors=c("#CC79A7","#E69F00","#CC79A7","#E69F00")
DIL_lines=c("solid","solid", "dashed", "dashed")


pl_dil<-ggplot(dil_melt, aes(x = Set, y = perf, colour = schedule,linetype = schedule)) +
	xlab("Set")+scale_x_continuous(breaks=1:10)+
	ylab("Dilation")+
	theme(legend.title=element_blank())+
	geom_line(size=1)+ geom_point(aes(shape=schedule),size = 4)
pl_dil+
	scale_shape_discrete(name="",breaks=DIL_ori_names,labels=DIL_legends)+
	scale_colour_manual(name="",values=DIL_colors,breaks=DIL_ori_names,labels=DIL_legends)+ 
	scale_linetype_manual(name="",values=DIL_lines,breaks=DIL_ori_names,labels=DIL_legends)+ 
#legend inside plotting area
	theme(legend.justification=c(1,1), legend.position=c(1,1))
                         
ggsave("../fig/dilation_comparison.pdf", width=8, height=5,device=cairo_pdf)

