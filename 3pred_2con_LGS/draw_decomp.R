
# Create graph of mechanistic decomposition
rm(list=ls())
require(ggfortify)
require(ggplot2)
require(GGally)
require(scales)
require(cowplot)
require(ggrepel)
require(reshape2)
require(ggforce)
#read them in
C1_final_mechanisms_3con <- read.csv(file = "C:/Users/aiteu/Documents/eco_working_group/C1_final_mechanisms_3con_sig0_3.csv")
C2_final_mechanisms_3con <- read.csv(file = "C:/Users/aiteu/Documents/eco_working_group/C2_final_mechanisms_3con_sig0_3.csv")
C3_final_mechanisms_3con <- read.csv(file = "C:/Users/aiteu/Documents/eco_working_group/C3_final_mechanisms_3con_sig0_3.csv")

C1_final_mechanisms_3con_2pred <- read.csv(file = "C:/Users/aiteu/Documents/eco_working_group/C1_final_mechanisms_3con_sig0_2.csv")
C2_final_mechanisms_3con_2pred <- read.csv(file = "C:/Users/aiteu/Documents/eco_working_group/C2_final_mechanisms_3con_sig0_2.csv")
C3_final_mechanisms_3con_2pred <- read.csv(file = "C:/Users/aiteu/Documents/eco_working_group/C3_final_mechanisms_3con_sig0_2.csv")


min_y<-min(C1_final_mechanisms_3con,C2_final_mechanisms_3con ,C3_final_mechanisms_3con,C1_final_mechanisms_3con_2pred,C2_final_mechanisms_3con_2pred,C3_final_mechanisms_3con_2pred)  
max_y<-max(C1_final_mechanisms_3con,C2_final_mechanisms_3con ,C3_final_mechanisms_3con,C1_final_mechanisms_3con_2pred,C2_final_mechanisms_3con_2pred,C3_final_mechanisms_3con_2pred )  

min_y<- max(min_y,max_y)*-1
max_y <- max(min_y,max_y)

C1_pos<-melt(C1_final_mechanisms_3con)
C2_pos<-melt(C2_final_mechanisms_3con)
C3_pos<-melt(C3_final_mechanisms_3con)

C1_pos$Competitor<-rep("Competitor 1",100)
C2_pos$Competitor<-rep("Competitor 2",100)
C3_pos$Competitor<-rep("Competitor 3",100)

pos<-rbind(C1_pos,C2_pos,C3_pos)
pos$delta<-substring(pos$variable,4)

levels<-c("C1_r_bar","C2_r_bar","C3_r_bar","C1_delta_0","C2_delta_0","C3_delta_0","C1_delta_P","C2_delta_P","C3_delta_P","C1_delta_E","C2_delta_E","C3_delta_E","C1_delta_EP","C2_delta_EP","C3_delta_EP")

pos$variable=factor(pos$variable,levels = levels)

pos$delta<-factor(pos$delta, levels=unique(pos$delta))

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


cols<-c(1,6,3,8,2)


cbPalette <- cbbPalette[cols]
 

fad<-cbPalette
xlab=c(expression(bar("r")[i]-bar("r")[r]) ,
       expression(bar(Delta)[i]^0),
       expression(bar(Delta)[i]^P),
       expression(bar(Delta)[i]^E),
       expression(bar(Delta)[i]^{E*P}))


p1<-pos[which(pos$Competitor== "Competitor 1"),]
g1_bar<-ggplot(data=p1, aes(x=variable, y=value,  fill=variable)) + 
  geom_bar(aes(x=delta),position = "dodge",stat = "summary", fun.y = "mean") +
  stat_summary(aes(x=delta),position = "dodge",fun.data = mean_se, geom = "errorbar")+
  scale_fill_manual(values=fad)+scale_x_discrete(labels=xlab)+xlab("mechanistic partitioning")+ylab("growth rate\nwhen rare")+ylim(min_y,max_y)+
  theme(legend.position = "none")+ggtitle("Competitor 1") +
  theme(plot.title = element_text(hjust = 0.5))


g1_bar

p2<-pos[which(pos$Competitor== "Competitor 2"),]
g2_bar<-ggplot(data=p2, aes(x=variable, y=value,  fill=variable)) + 
  geom_bar(aes(x=delta),position = "dodge",stat = "summary", fun.y = "mean") +
  stat_summary(aes(x=delta),position = "dodge",fun.data = mean_se, geom = "errorbar")+
  scale_fill_manual(values=fad)+scale_x_discrete(labels=xlab)+xlab("mechanistic partitioning")+ylab("growth rate\nwhen rare")+ylim(min_y,max_y)+
  theme(legend.position = "none")+ggtitle("Competitor 2") +
  theme(plot.title = element_text(hjust = 0.5))


g2_bar



p3<-pos[which(pos$Competitor== "Competitor 3"),]
g3_bar<-ggplot(data=p3, aes(x=variable, y=value,  fill=variable)) + 
  geom_bar(aes(x=delta),position = "dodge",stat = "summary", fun.y = "mean") +
  stat_summary(aes(x=delta),position = "dodge",fun.data = mean_se, geom = "errorbar")+
  scale_fill_manual(values=fad)+scale_x_discrete(labels=xlab)+xlab("mechanistic partitioning")+ylab("growth rate\nwhen rare")+ylim(min_y,max_y)+
  theme(legend.position = "none")+ggtitle("Competitor 3") +
  theme(plot.title = element_text(hjust = 0.5))


g3_bar

C1_no<-melt(C1_final_mechanisms_3con_2pred)
C2_no<-melt(C2_final_mechanisms_3con_2pred)
C3_no<-melt(C3_final_mechanisms_3con_2pred)

C1_no$Competitor<-rep("Competitor 1",100)
C2_no$Competitor<-rep("Competitor 2",100)
C3_no$Competitor<-rep("Competitor 3",100)

no<-rbind(C1_no,C2_no,C3_no)
no$delta<-substring(no$variable,4)

levels<-c("C1_r_bar","C2_r_bar","C3_r_bar","C1_delta_0","C2_delta_0","C3_delta_0","C1_delta_P","C2_delta_P","C3_delta_P","C1_delta_E","C2_delta_E","C3_delta_E","C1_delta_EP","C2_delta_EP","C3_delta_EP")

no$variable=factor(no$variable,levels = levels)

no$delta<-factor(no$delta, levels=unique(no$delta))


p1<-no[which(no$Competitor== "Competitor 1"),]
g1_bar2<-ggplot(data=p1, aes(x=variable, y=value,  fill=variable)) + 
  geom_bar(aes(x=delta),position = "dodge",stat = "summary", fun.y = "mean") +
  stat_summary(aes(x=delta),position = "dodge",fun.data = mean_se, geom = "errorbar")+
  scale_fill_manual(values=fad)+scale_x_discrete(labels=xlab)+xlab("mechanistic partitioning")+ylab("growth rate\nwhen rare")+ylim(min_y,max_y)+
  theme(legend.position = "none")+ggtitle("Competitor 1") +
  theme(plot.title = element_text(hjust = 0.5))


g1_bar2

p22<-no[which(no$Competitor== "Competitor 2"),]
g2_bar2<-ggplot(data=p22, aes(x=variable, y=value,  fill=variable)) + 
  geom_bar(aes(x=delta),position = "dodge",stat = "summary", fun.y = "mean") +
  stat_summary(aes(x=delta),position = "dodge",fun.data = mean_se, geom = "errorbar")+
  scale_fill_manual(values=fad)+scale_x_discrete(labels=xlab)+xlab("mechanistic partitioning")+ylab("growth rate\nwhen rare")+ylim(min_y,max_y)+
  theme(legend.position = "none")+ggtitle("Competitor 2") +
  theme(plot.title = element_text(hjust = 0.5))


p33<-no[which(no$Competitor== "Competitor 3"),]
g3_bar2<-ggplot(data=p33, aes(x=variable, y=value,  fill=variable)) + 
  geom_bar(aes(x=delta),position = "dodge",stat = "summary", fun.y = "mean") +
  stat_summary(aes(x=delta),position = "dodge",fun.data = mean_se, geom = "errorbar")+
  scale_fill_manual(values=fad)+scale_x_discrete(labels=xlab)+xlab("mechanistic partitioning")+ylab("growth rate\nwhen rare")+ylim(min_y,max_y)+
  theme(legend.position = "none")+ggtitle("Competitor 3") +
  theme(plot.title = element_text(hjust = 0.5))


pp<- plot_grid(g1_bar,g2_bar,g3_bar,g1_bar2,g2_bar2,g3_bar2,
               labels=c("A","B","C","D","E","F"),
               align = 'vh',
               hjust = -1,
               nrow = 2)
pp
save_plot("3con_2parms.pdf",plot=pp,base_height=12,units = "in")
