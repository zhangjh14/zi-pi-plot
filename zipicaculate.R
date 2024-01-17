#this function is used to calculate the Among-module connectivity Pi and Within-module connectivity Zi
#the input is igraph created by graph.adjacency() from adjacency matrices
#the out is a data.frame catains the degree, module, pi, zi of each node in the graph
zipiCalculate<-function(net){
  E(net)$weight = NA
  fc.net<-cluster_fast_greedy(net,weights =NULL)
  modularity.net<-modularity(net,membership(fc.net))
  comps.group4<-membership(fc.net)
  
  seqdeg<-igraph::degree(net,v = V(net))
  Nnodes=length(seqdeg)
  Z=seqdeg
  Z[]=0
  P=Z
  
  # sp_g_Im<-read.csv("sp_g_1m.csv",header = T)
  # fc.sp_g_I$membership<-sp_g_Im[,2]
  Membership=membership(fc.net)
  Seq=seq(1:Nnodes)
  for(i in 1:Nnodes){
    L=Membership==Membership[i]
    neighbs=neighbors(net,i)
    Kis=sum(L[neighbs])
    SUM=0
    SUMsq=0	
    SUMP=0
    Miv=Seq[L]
    for(j in 1:sum(L)){
      neighbsj=neighbors(net,Miv[j])
      Kjs=sum(L[neighbsj])
      SUM=SUM+Kjs
      SUMsq=SUMsq+Kjs^2
    }
    Z[i]=(Kis-SUM/sum(L))/sqrt(SUMsq/sum(L)-(SUM/sum(L))^2)
    if(Kis-SUM/sum(L)==0){Z[i]=0}
    for(k in 1:max(Membership)){
      Lp=Membership==k
      Kisp=sum(Lp[neighbs])
      SUMP=SUMP+(Kisp/seqdeg[i])^2
      }
    P[i]=1-SUMP
    
  }
  result<-data.frame(degree=seqdeg,module=Membership,Pi=P,Zi=Z)
  return(result)
}


#zi-pi plot
attribute<- zipiCalculate(net)
ggplot(attribute)+
  geom_point(aes(x=Pi,y=Zi))+
  geom_vline(xintercept = 0.62,linetype = "dashed", color = "blue")+
  geom_hline(yintercept = 2.5,linetype = "dashed", color = "blue")+
  xlab("Among-module connectivity Pi")+
  ylab("Within-module connectivity Zi")+
  geom_text(data = zipilabel,aes(x,y,label=label))+
  theme(axis.line = element_line(linetype = "solid"), 
        axis.ticks = element_line(colour = "black"),
        axis.ticks.length=unit(-0.25, "cm"),
        panel.border = element_rect(fill = NA, 
                                    colour = "grey20"),
        panel.grid.major = element_line(colour = NA), 
        panel.grid.minor = element_line(colour = NA), 
        axis.title = element_text(size = 14, 
                                  face = "bold"), 
        axis.text = element_text(size = 12, 
                                 face = "bold", 
                                 colour = "black"), 
        axis.text.x = element_text(size = 12,
                                   vjust = -2), 
        axis.text.y = element_text(size = 12,
                                   hjust = -0.25), 
        plot.title = element_text(size = 15,
                                  face = "bold"), 
        panel.background = element_rect(fill = "white"), 
        legend.position = "none",
        plot.background = element_rect(colour = NA)) 

