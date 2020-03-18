#Making phylo figures
library(ape)
library(phytools)
library(ggtree)
library(ggplot2)
library(cowplot)

#taxa<-read.csv("C:/Users/Brian/Google Drive/DNH_scale/L48_taxa.csv")
status<-read.csv("C:/Users/Brian/Google Drive/DNH_scale/status4.csv",stringsAsFactors = F)
status<-status[which(status$status != "NA"),]

#need figures of: 
  #(pruned) base trees (gbmb, gbotb)
  #genus-only trees (gbmb,gbotb)
  #genus+fam trees (gbmb,gbotb)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#base trees

#gbotb
gbotb_base <- read.tree("smith_and_brown_2018_trees/v0.1/v0.1/GBOTB.tre")
gbotb_base <- drop.tip(phy =gbotb_base,tip = setdiff(x = gbotb_base$tip.label,y = status$binomial))

#Get taxa
angios <- getMRCA(phy = gbotb_base, tip = intersect(x = status$binomial[which(status$Group=="Angiosperms")],y = gbotb_base$tip.label) )
gymnos <- getMRCA(phy = gbotb_base, tip = intersect(x = status$binomial[which(status$Group=="Gymnosperms")],y = gbotb_base$tip.label) )


#sort(table(status$Family),decreasing = T)[1:10]
asteraceae <- getMRCA(phy = gbotb_base, tip = intersect(x = status$binomial[which(status$Family=="Asteraceae")],y = gbotb_base$tip.label) )
fabaceae <- getMRCA(phy = gbotb_base, tip = intersect(x = status$binomial[which(status$Family=="Fabaceae")],y = gbotb_base$tip.label) )
poaceae <- getMRCA(phy = gbotb_base, tip = intersect(x = status$binomial[which(status$Family=="Poaceae")],y = gbotb_base$tip.label) )
cyperaceae <- getMRCA(phy = gbotb_base, tip = intersect(x = status$binomial[which(status$Family=="Cyperaceae")],y = gbotb_base$tip.label) )
scrophulariaceae  <- getMRCA(phy = gbotb_base, tip = intersect(x = status$binomial[which(status$Family=="Scrophulariaceae")],y = gbotb_base$tip.label) )

l1 <- get_legend(ggtree(tr = gbotb_base,layout = "circular") %<+% status+geom_treescale(y = length(gbotb_base$tip.label)-1000,offset = 100,width = 100) + 
                   geom_tippoint(aes(subset=status %in% c('I')), color='red', alpha=.8,position = ggplot2::position_nudge(x = 10))+
                   geom_tippoint(aes(subset=status %in% c('N')), color='blue', alpha=.5  )+
                   geom_cladelabel(node = angios, label="Angiosperms", offset = 20,barsize = 3,color="grey")+
                   geom_tippoint(aes(color=status),show.legend = T)+theme(legend.position = "right")+labs(color = "Native Status\n")+
                   scale_color_manual(labels = c("Native", "Non-native"), values = c("blue", "red")))
#plot(l1)

t1 <- ggtree(tr = gbotb_base,layout = "circular") %<+% status+geom_treescale(y = length(gbotb_base$tip.label)-1000,offset = 100,width = 100) + 
  geom_tippoint(aes(subset=status %in% c('I')), color='red', alpha=.8,size=2,position = ggplot2::position_nudge(x = 25))+
  geom_tippoint(aes(subset=status %in% c('N')), color='blue', alpha=.5,size=2, position = ggplot2::position_nudge(x = 10)  )+
  geom_cladelabel(node = asteraceae, label="Asteraceae", offset = 35,angle = 90,offset.text = 50,barsize = 3,color="grey")+
  geom_cladelabel(node = fabaceae, label="Fabaceae", offset = 35,angle = 90,offset.text = 50,barsize = 3,color="brown")+
  geom_cladelabel(node = poaceae, label="Poaceae", offset = 35,angle = 90,offset.text = 50,barsize = 3,color="tan")+
  geom_cladelabel(node = cyperaceae, label="Cyperaceae", offset = 35,angle = 90,offset.text = 50,barsize = 3,color="dark grey")+
  geom_cladelabel(node = scrophulariaceae, label="Scrophulariaceae", offset = 35,angle = 0,offset.text = 50,barsize = 3,color="#CC6600")+
  ggtitle(label = "GBOTB")


#plot(t1)

combined_gbotb_base <- ggdraw(t1) +
  draw_plot(plot = l1,x = .8,y = .8,height = .2,width = .2)

combined_gbotb_base
ggsave(filename = "DNH_scale_GBOT_base.pdf",plot = combined_gbotb_base,device = "pdf",path = "C:/Users/Brian/Desktop/",width = 20,height = 20,units = "cm")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#gbmb

gbmb_base <- read.tree("smith_and_brown_2018_trees/v0.1/v0.1/GBMB.tre")
gbmb_base <- drop.tip(phy =gbmb_base,tip = setdiff(x = gbmb_base$tip.label,y = status$binomial))

#Get taxa
angios <- getMRCA(phy = gbmb_base, tip = intersect(x = status$binomial[which(status$Group=="Angiosperms")],y = gbmb_base$tip.label) )
gymnos <- getMRCA(phy = gbmb_base, tip = intersect(x = status$binomial[which(status$Group=="Gymnosperms")],y = gbmb_base$tip.label) )


#sort(table(status$Family),decreasing = T)[1:3]
asteraceae <- getMRCA(phy = gbmb_base, tip = intersect(x = status$binomial[which(status$Family=="Asteraceae")],y = gbmb_base$tip.label) )
fabaceae <- getMRCA(phy = gbmb_base, tip = intersect(x = status$binomial[which(status$Family=="Fabaceae")],y = gbmb_base$tip.label) )
poaceae <- getMRCA(phy = gbmb_base, tip = intersect(x = status$binomial[which(status$Family=="Poaceae")],y = gbmb_base$tip.label) )



l1 <- get_legend(ggtree(tr = gbmb_base,layout = "circular") %<+% status+geom_treescale(y = length(gbmb_base$tip.label)-1000,offset = 100,width = 100) + 
                   geom_tippoint(aes(subset=status %in% c('I')), color='red', alpha=.8,position = ggplot2::position_nudge(x = 10))+
                   geom_tippoint(aes(subset=status %in% c('N')), color='blue', alpha=.5  )+
                   geom_cladelabel(node = angios, label="Angiosperms", offset = 20,barsize = 3,color="grey")+
                   geom_tippoint(aes(color=status),show.legend = T)+theme(legend.position = "right")+labs(color = "Native Status\n")+
                   scale_color_manual(labels = c("Native", "Non-native"), values = c("blue", "red")))
#plot(l1)

t1 <- ggtree(tr = gbmb_base,layout = "circular") %<+% status+geom_treescale(y = length(gbmb_base$tip.label)-1000,offset = 100,width = 100) + 
  geom_tippoint(aes(subset=status %in% c('I')), color='red', alpha=.8,size=2,position = ggplot2::position_nudge(x = 25))+
  geom_tippoint(aes(subset=status %in% c('N')), color='blue', alpha=.5,size=2, position = ggplot2::position_nudge(x = 10)  )+
  geom_cladelabel(node = asteraceae, label="Asteraceae", offset = 35,angle = 90,offset.text = 50,barsize = 3,color="grey")+
  geom_cladelabel(node = fabaceae, label="Fabaceae", offset = 35,angle = 90,offset.text = 50,barsize = 3,color="brown")+
  geom_cladelabel(node = poaceae, label="Poaceae", offset = 35,angle = 90,offset.text = 50,barsize = 3,color="tan")+ggtitle(label = "GBMB")


#plot(t1)

combined_gbmb_base <- ggdraw(t1) +
  draw_plot(plot = l1,x = .8,y = .8,height = .2,width = .2)

combined_gbmb_base




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#gbmb_plus_congeners

gbmb_plus_congeners <- read.tree("sunplin_trees/genus_addition_only/gbmb_genus_additions_only/gbmb_genus_addition_only616.tre") #616 = my area code
gbmb_plus_congeners <- drop.tip(phy =gbmb_plus_congeners,tip = setdiff(x = gbmb_plus_congeners$tip.label,y = status$binomial))

#Get taxa
angios <- getMRCA(phy = gbmb_plus_congeners, tip = intersect(x = status$binomial[which(status$Group=="Angiosperms")],y = gbmb_plus_congeners$tip.label) )
gymnos <- getMRCA(phy = gbmb_plus_congeners, tip = intersect(x = status$binomial[which(status$Group=="Gymnosperms")],y = gbmb_plus_congeners$tip.label) )


#sort(table(status$Family),decreasing = T)[1:3]
asteraceae <- getMRCA(phy = gbmb_plus_congeners, tip = intersect(x = status$binomial[which(status$Family=="Asteraceae")],y = gbmb_plus_congeners$tip.label) )
fabaceae <- getMRCA(phy = gbmb_plus_congeners, tip = intersect(x = status$binomial[which(status$Family=="Fabaceae")],y = gbmb_plus_congeners$tip.label) )
poaceae <- getMRCA(phy = gbmb_plus_congeners, tip = intersect(x = status$binomial[which(status$Family=="Poaceae")],y = gbmb_plus_congeners$tip.label) )



l1 <- get_legend(ggtree(tr = gbmb_plus_congeners,layout = "circular") %<+% status+geom_treescale(y = length(gbmb_plus_congeners$tip.label)-1000,offset = 100,width = 100) + 
                   geom_tippoint(aes(subset=status %in% c('I')), color='red', alpha=.8,position = ggplot2::position_nudge(x = 10))+
                   geom_tippoint(aes(subset=status %in% c('N')), color='blue', alpha=.5  )+
                   geom_cladelabel(node = angios, label="Angiosperms", offset = 20,barsize = 3,color="grey")+
                   geom_tippoint(aes(color=status),show.legend = T)+theme(legend.position = "right")+labs(color = "Native Status\n")+
                   scale_color_manual(labels = c("Native", "Non-native"), values = c("blue", "red")))
#plot(l1)

t1 <- ggtree(tr = gbmb_plus_congeners,layout = "circular") %<+% status+geom_treescale(y = length(gbmb_plus_congeners$tip.label)-1000,offset = 100,width = 100) + 
  geom_tippoint(aes(subset=status %in% c('I')), color='red', alpha=.8,size=2,position = ggplot2::position_nudge(x = 25))+
  geom_tippoint(aes(subset=status %in% c('N')), color='blue', alpha=.5,size=2, position = ggplot2::position_nudge(x = 10)  )+
  geom_cladelabel(node = asteraceae, label="Asteraceae", offset = 35,angle = 90,offset.text = 50,barsize = 3,color="grey")+
  geom_cladelabel(node = fabaceae, label="Fabaceae", offset = 35,angle = 90,offset.text = 50,barsize = 3,color="brown")+
  geom_cladelabel(node = poaceae, label="Poaceae", offset = 35,angle = 90,offset.text = 50,barsize = 3,color="tan")+ggtitle(label = "GBMB plus congeners")


#plot(t1)

combined_gbmb_plus_congeners <- ggdraw(t1) +
  draw_plot(plot = l1,x = .8,y = .8,height = .2,width = .2)

combined_gbmb_plus_congeners

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#gbotb_plus_congeners

gbotb_plus_congeners <- read.tree("sunplin_trees/genus_addition_only/gbotb_genus_additions_only/gbotb_genus_addition_only616.tre") #616 = my area code
gbotb_plus_congeners <- drop.tip(phy =gbotb_plus_congeners,tip = setdiff(x = gbotb_plus_congeners$tip.label,y = status$binomial))

#Get taxa
angios <- getMRCA(phy = gbotb_plus_congeners, tip = intersect(x = status$binomial[which(status$Group=="Angiosperms")],y = gbotb_plus_congeners$tip.label) )
gymnos <- getMRCA(phy = gbotb_plus_congeners, tip = intersect(x = status$binomial[which(status$Group=="Gymnosperms")],y = gbotb_plus_congeners$tip.label) )


#sort(table(status$Family),decreasing = T)[1:3]
asteraceae <- getMRCA(phy = gbotb_plus_congeners, tip = intersect(x = status$binomial[which(status$Family=="Asteraceae")],y = gbotb_plus_congeners$tip.label) )
fabaceae <- getMRCA(phy = gbotb_plus_congeners, tip = intersect(x = status$binomial[which(status$Family=="Fabaceae")],y = gbotb_plus_congeners$tip.label) )
poaceae <- getMRCA(phy = gbotb_plus_congeners, tip = intersect(x = status$binomial[which(status$Family=="Poaceae")],y = gbotb_plus_congeners$tip.label) )
cyperaceae <- getMRCA(phy = gbotb_plus_congeners, tip = intersect(x = status$binomial[which(status$Family=="Cyperaceae")],y = gbotb_plus_congeners$tip.label) )
scrophulariaceae  <- getMRCA(phy = gbotb_plus_congeners, tip = intersect(x = status$binomial[which(status$Family=="Scrophulariaceae")],y = gbotb_plus_congeners$tip.label) )



l1 <- get_legend(ggtree(tr = gbotb_plus_congeners,layout = "circular") %<+% status+geom_treescale(y = length(gbotb_plus_congeners$tip.label)-1000,offset = 100,width = 100) + 
                   geom_tippoint(aes(subset=status %in% c('I')), color='red', alpha=.8,position = ggplot2::position_nudge(x = 10))+
                   geom_tippoint(aes(subset=status %in% c('N')), color='blue', alpha=.5  )+
                   geom_cladelabel(node = angios, label="Angiosperms", offset = 20,barsize = 3,color="grey")+
                   geom_tippoint(aes(color=status),show.legend = T)+theme(legend.position = "right")+labs(color = "Native Status\n")+
                   scale_color_manual(labels = c("Native", "Non-native"), values = c("blue", "red")))
#plot(l1)

t1 <- ggtree(tr = gbotb_plus_congeners,layout = "circular") %<+% status+geom_treescale(y = length(gbotb_plus_congeners$tip.label)-1000,offset = 100,width = 100) + 
  geom_tippoint(aes(subset=status %in% c('I')), color='red', alpha=.8,size=2,position = ggplot2::position_nudge(x = 25))+
  geom_tippoint(aes(subset=status %in% c('N')), color='blue', alpha=.5,size=2, position = ggplot2::position_nudge(x = 10)  )+
  geom_cladelabel(node = asteraceae, label="Asteraceae", offset = 35,angle = 90,offset.text = 50,barsize = 3,color="grey")+
  geom_cladelabel(node = fabaceae, label="Fabaceae", offset = 35,angle = 90,offset.text = 50,barsize = 3,color="brown")+
  geom_cladelabel(node = poaceae, label="Poaceae", offset = 35,angle = 90,offset.text = 50,barsize = 3,color="tan")+
  geom_cladelabel(node = cyperaceae, label="Cyperaceae", offset = 35,angle = 90,offset.text = 50,barsize = 3,color="dark grey")+
  geom_cladelabel(node = scrophulariaceae, label="Scrophulariaceae", offset = 35,angle = 0,offset.text = 50,barsize = 3,color="#CC6600")+
  ggtitle(label = "GBOTB plus congeners")


#plot(t1)

combined_gbotb_plus_congeners <- ggdraw(t1) +
  draw_plot(plot = l1,x = .8,y = .8,height = .2,width = .2)

combined_gbotb_plus_congeners
ggsave(filename = "DNH_scale_GBOTB_plus_congeners.pdf",plot = combined_gbotb_plus_congeners,device = "pdf",path = "C:/Users/Brian/Desktop/",width = 20,height = 20,units = "cm")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#gbmb_plus_congeners_and_confamilials

gbmb_plus_congeners_and_confamilials <- read.tree("sunplin_trees/genus_and_family_additions/gbmb_genus_and_family_additions/gbmb_and_family_additions616.tre") #616 = my area code
gbmb_plus_congeners_and_confamilials <- drop.tip(phy =gbmb_plus_congeners_and_confamilials,tip = setdiff(x = gbmb_plus_congeners_and_confamilials$tip.label,y = status$binomial))

#Get taxa
angios <- getMRCA(phy = gbmb_plus_congeners_and_confamilials, tip = intersect(x = status$binomial[which(status$Group=="Angiosperms")],y = gbmb_plus_congeners_and_confamilials$tip.label) )
gymnos <- getMRCA(phy = gbmb_plus_congeners_and_confamilials, tip = intersect(x = status$binomial[which(status$Group=="Gymnosperms")],y = gbmb_plus_congeners_and_confamilials$tip.label) )


#sort(table(status$Family),decreasing = T)[1:3]
asteraceae <- getMRCA(phy = gbmb_plus_congeners_and_confamilials, tip = intersect(x = status$binomial[which(status$Family=="Asteraceae")],y = gbmb_plus_congeners_and_confamilials$tip.label) )
fabaceae <- getMRCA(phy = gbmb_plus_congeners_and_confamilials, tip = intersect(x = status$binomial[which(status$Family=="Fabaceae")],y = gbmb_plus_congeners_and_confamilials$tip.label) )
poaceae <- getMRCA(phy = gbmb_plus_congeners_and_confamilials, tip = intersect(x = status$binomial[which(status$Family=="Poaceae")],y = gbmb_plus_congeners_and_confamilials$tip.label) )



l1 <- get_legend(ggtree(tr = gbmb_plus_congeners_and_confamilials,layout = "circular") %<+% status+geom_treescale(y = length(gbmb_plus_congeners_and_confamilials$tip.label)-1000,offset = 100,width = 100) + 
                   geom_tippoint(aes(subset=status %in% c('I')), color='red', alpha=.8,position = ggplot2::position_nudge(x = 10))+
                   geom_tippoint(aes(subset=status %in% c('N')), color='blue', alpha=.5  )+
                   geom_cladelabel(node = angios, label="Angiosperms", offset = 20,barsize = 3,color="grey")+
                   geom_tippoint(aes(color=status),show.legend = T)+theme(legend.position = "right")+labs(color = "Native Status\n")+
                   scale_color_manual(labels = c("Native", "Non-native"), values = c("blue", "red")))
#plot(l1)

t1 <- ggtree(tr = gbmb_plus_congeners_and_confamilials,layout = "circular") %<+% status+geom_treescale(y = length(gbmb_plus_congeners_and_confamilials$tip.label)-1000,offset = 100,width = 100) + 
  geom_tippoint(aes(subset=status %in% c('I')), color='red', alpha=.8,size=2,position = ggplot2::position_nudge(x = 25))+
  geom_tippoint(aes(subset=status %in% c('N')), color='blue', alpha=.5,size=2, position = ggplot2::position_nudge(x = 10)  )+
  geom_cladelabel(node = asteraceae, label="Asteraceae", offset = 35,angle = 90,offset.text = 50,barsize = 3,color="grey")+
  geom_cladelabel(node = fabaceae, label="Fabaceae", offset = 35,angle = 90,offset.text = 50,barsize = 3,color="brown")+
  geom_cladelabel(node = poaceae, label="Poaceae", offset = 35,angle = 90,offset.text = 50,barsize = 3,color="tan")+ggtitle(label = "GBMB plus congeners and confamilials")


#plot(t1)

combined_gbmb_plus_congeners_and_confamilials <- ggdraw(t1) +
  draw_plot(plot = l1,x = .8,y = .8,height = .2,width = .2)

combined_gbmb_plus_congeners_and_confamilials

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#gbotb_plus_congeners_and_confamilials

gbotb_plus_congeners_and_confamilials <- read.tree("sunplin_trees/genus_and_family_additions/gbotb_genus_and_family_additions/gbotb_genus_and_family_additions616.tre") #616 = my area code
gbotb_plus_congeners_and_confamilials <- drop.tip(phy =gbotb_plus_congeners_and_confamilials,tip = setdiff(x = gbotb_plus_congeners_and_confamilials$tip.label,y = status$binomial))

#Get taxa
angios <- getMRCA(phy = gbotb_plus_congeners_and_confamilials, tip = intersect(x = status$binomial[which(status$Group=="Angiosperms")],y = gbotb_plus_congeners_and_confamilials$tip.label) )
gymnos <- getMRCA(phy = gbotb_plus_congeners_and_confamilials, tip = intersect(x = status$binomial[which(status$Group=="Gymnosperms")],y = gbotb_plus_congeners_and_confamilials$tip.label) )


#sort(table(status$Family),decreasing = T)[1:3]
asteraceae <- getMRCA(phy = gbotb_plus_congeners_and_confamilials, tip = intersect(x = status$binomial[which(status$Family=="Asteraceae")],y = gbotb_plus_congeners_and_confamilials$tip.label) )
fabaceae <- getMRCA(phy = gbotb_plus_congeners_and_confamilials, tip = intersect(x = status$binomial[which(status$Family=="Fabaceae")],y = gbotb_plus_congeners_and_confamilials$tip.label) )
poaceae <- getMRCA(phy = gbotb_plus_congeners_and_confamilials, tip = intersect(x = status$binomial[which(status$Family=="Poaceae")],y = gbotb_plus_congeners_and_confamilials$tip.label) )
cyperaceae <- getMRCA(phy = gbotb_plus_congeners_and_confamilials, tip = intersect(x = status$binomial[which(status$Family=="Cyperaceae")],y = gbotb_plus_congeners_and_confamilials$tip.label) )
scrophulariaceae  <- getMRCA(phy = gbotb_plus_congeners_and_confamilials, tip = intersect(x = status$binomial[which(status$Family=="Scrophulariaceae")],y = gbotb_plus_congeners_and_confamilials$tip.label) )



l1 <- get_legend(ggtree(tr = gbotb_plus_congeners_and_confamilials,layout = "circular") %<+% status+geom_treescale(y = length(gbotb_plus_congeners_and_confamilials$tip.label)-1000,offset = 100,width = 100) + 
                   geom_tippoint(aes(subset=status %in% c('I')), color='red', alpha=.8,position = ggplot2::position_nudge(x = 10))+
                   geom_tippoint(aes(subset=status %in% c('N')), color='blue', alpha=.5  )+
                   geom_cladelabel(node = angios, label="Angiosperms", offset = 20,barsize = 3,color="grey")+
                   geom_tippoint(aes(color=status),show.legend = T)+theme(legend.position = "right")+labs(color = "Native Status\n")+
                   scale_color_manual(labels = c("Native", "Non-native"), values = c("blue", "red")))
#plot(l1)

t1 <- ggtree(tr = gbotb_plus_congeners_and_confamilials,layout = "circular") %<+% status+geom_treescale(y = length(gbotb_plus_congeners_and_confamilials$tip.label)-1000,offset = 100,width = 100) + 
  geom_tippoint(aes(subset=status %in% c('I')), color='red', alpha=.8,size=2,position = ggplot2::position_nudge(x = 25))+
  geom_tippoint(aes(subset=status %in% c('N')), color='blue', alpha=.5,size=2, position = ggplot2::position_nudge(x = 10)  )+
  geom_cladelabel(node = asteraceae, label="Asteraceae", offset = 35,angle = 90,offset.text = 50,barsize = 3,color="grey")+
  geom_cladelabel(node = fabaceae, label="Fabaceae", offset = 35,angle = 90,offset.text = 50,barsize = 3,color="brown")+
  geom_cladelabel(node = poaceae, label="Poaceae", offset = 35,angle = 90,offset.text = 50,barsize = 3,color="tan")+
  geom_cladelabel(node = cyperaceae, label="Cyperaceae", offset = 35,angle = 90,offset.text = 50,barsize = 3,color="dark grey")+
  geom_cladelabel(node = scrophulariaceae, label="Scrophulariaceae", offset = 35,angle = 0,offset.text = 50,barsize = 3,color="#CC6600")+
  ggtitle(label = "GBOTB plus congeners and confamilials")


#plot(t1)

combined_gbotb_plus_congeners_and_confamilials <- ggdraw(t1) +
  draw_plot(plot = l1,x = .8,y = .8,height = .2,width = .2)

combined_gbotb_plus_congeners_and_confamilials

ggsave(filename = "DNH_scale_GBOTB_plus_congeners_and_confamilials.pdf",plot = combined_gbotb_plus_congeners_and_confamilials,device = "pdf",path = "C:/Users/Brian/Desktop/",width = 20,height = 20,units = "cm")

