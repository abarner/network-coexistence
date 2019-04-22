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
C1_final_mechanisms_pos <- read.csv(file = "E:/eco_workinggroup/C1_final_mechanisms_3con_pos_cor.csv")
C2_final_mechanisms_pos <- read.csv(file = "E:/eco_workinggroup/C2_final_mechanisms_3con_pos_cor.csv")
C3_final_mechanisms_pos <- read.csv(file = "E:/eco_workinggroup/C3_final_mechanisms_3con_pos_cor.csv")

C1_final_mechanisms_no <- read.csv(file = "E:/eco_workinggroup/C1_final_mechanisms_3con_no_cor.csv")
C2_final_mechanisms_no <- read.csv(file = "E:/eco_workinggroup/C2_final_mechanisms_3con_no_cor.csv")
C3_final_mechanisms_no <- read.csv(file = "E:/eco_workinggroup/C3_final_mechanisms_3con_no_cor.csv")


C1_final_mechanisms_neg <- read.csv(file = "E:/eco_workinggroup/C1_final_mechanisms_3con_negative_cor.csv")
C2_final_mechanisms_neg <- read.csv(file = "E:/eco_workinggroup/C2_final_mechanisms_3con_negative.csv")
C3_final_mechanisms_neg <- read.csv(file = "E:/eco_workinggroup/C3_final_mechanisms_3con_negative.csv")


C1_pos<-melt(C1_final_mechanisms_pos)
C2_pos<-melt(C2_final_mechanisms_pos)
C3_pos<-melt(C3_final_mechanisms_pos)

C1_pos$Competitor<-rep("Competitor 1",100)
C2_pos$Competitor<-rep("Competitor 2",100)
C3_pos$Competitor<-rep("Competitor 3",100)

pos<-rbind(C1_pos,C2_pos,C3_pos)
pos$delta<-substring(pos$variable,4)

levels<-c("C1_r_bar","C2_r_bar","C3_r_bar","C1_delta_0","C2_delta_0","C3_delta_0","C1_delta_P","C2_delta_P","C3_delta_P","C1_delta_E","C2_delta_E","C3_delta_E","C1_delta_EP","C2_delta_EP","C3_delta_EP")

pos$variable=factor(pos$variable,levels = levels)

pos$delta<-factor(pos$delta, levels=unique(pos$delta))






cbPalette <- c("#E69F00", "#009E73","#56B4E9")
ll_Palette<-NULL
for(i in 1:length(cbPalette)){
  ll_Palette<-c(ll_Palette,colorspace::lighten(cbPalette[i],amount=0.5))
  
}  

fad<-c(cbPalette,rep(ll_Palette,4))
xlab=c(expression(bar("r")[i]-bar("r")[r]) ,
                                         expression(bar(Delta)[i]^0),
                                        expression(bar(Delta)[i]^P),
                                        expression(bar(Delta)[i]^E),
                                         expression(bar(Delta)[i]^{E*P}))

g_bar<-ggplot(data=pos, aes(x=variable, y=value,  fill=Competitor )) + 
  geom_bar(aes(x=delta),position = "dodge",stat = "summary", fun.y = "mean") +
  stat_summary(aes(x=delta),position = "dodge",fun.data = mean_se, geom = "errorbar")+
  scale_fill_manual(values=fad) +scale_x_discrete(labels=xlab)+xlab("mechanistic partitioning")+ylab("growth rate\nwhen rare")+ 
  theme(legend.position = c(0.85, 0.9))+ theme(legend.title = element_blank())
g_bar
                                                                         

C1_no<-melt(C1_final_mechanisms_no)
C2_no<-melt(C2_final_mechanisms_no)
C3_no<-melt(C3_final_mechanisms_no)

C1_no$Competitor<-rep("Competitor 1",100)
C2_no$Competitor<-rep("Competitor 2",100)
C3_no$Competitor<-rep("Competitor 3",100)

no<-rbind(C1_no,C2_no,C3_no)
no$delta<-substring(no$variable,4)

levels<-c("C1_r_bar","C2_r_bar","C3_r_bar","C1_delta_0","C2_delta_0","C3_delta_0","C1_delta_P","C2_delta_P","C3_delta_P","C1_delta_E","C2_delta_E","C3_delta_E","C1_delta_EP","C2_delta_EP","C3_delta_EP")

no$variable=factor(no$variable,levels = levels)

no$delta<-factor(no$delta, levels=unique(no$delta))

g_bar1<-ggplot(data=no, aes(x=variable, y=value,  fill=Competitor )) + 
  geom_bar(aes(x=delta),position = "dodge",stat = "summary", fun.y = "mean") +
  stat_summary(aes(x=delta),position = "dodge",fun.data = mean_se, geom = "errorbar")+
  scale_fill_manual(values=fad) +scale_x_discrete(labels=xlab)+xlab("mechanistic partitioning")+ylab("growth rate\nwhen rare")+ 
  theme(legend.position = "none")+ theme(legend.title = element_blank())
g_bar1



C1_neg<-melt(C1_final_mechanisms_neg)
C2_neg<-melt(C2_final_mechanisms_neg)
C3_neg<-melt(C3_final_mechanisms_neg)

C1_neg$Competitor<-rep("Competitor 1",100)
C2_neg$Competitor<-rep("Competitor 2",100)
C3_neg$Competitor<-rep("Competitor 3",100)

neg<-rbind(C1_neg,C2_neg,C3_neg)
neg$delta<-substring(neg$variable,4)

levels<-c("C1_r_bar","C2_r_bar","C3_r_bar","C1_delta_0","C2_delta_0","C3_delta_0","C1_delta_P","C2_delta_P","C3_delta_P","C1_delta_E","C2_delta_E","C3_delta_E","C1_delta_EP","C2_delta_EP","C3_delta_EP")

neg$variable=factor(neg$variable,levels = levels)

neg$delta<-factor(neg$delta, levels=unique(neg$delta))

g_bar2<-ggplot(data=neg, aes(x=variable, y=value,  fill=Competitor )) + 
  geom_bar(aes(x=delta),position = "dodge",stat = "summary", fun.y = "mean") +
  stat_summary(aes(x=delta),position = "dodge",fun.data = mean_se, geom = "errorbar")+
  scale_fill_manual(values=fad) +scale_x_discrete(labels=xlab)+xlab("mechanistic partitioning")+ylab("growth rate\nwhen rare")+ 
  theme(legend.position = "none")+ theme(legend.title = element_blank())
g_bar2


pp<- plot_grid(g_bar,g_bar1,g_bar2,
                labels=c("A","B","C"),
                align = 'vh',
                hjust = -1,
                nrow = 3)

ggplot2::ggsave("C:/Users/User/Documents/3con.pdf",plot=pp,width = 12, height = 12,units = "in")
