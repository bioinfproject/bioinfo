#
#
# Definition of ranges and positions for functions
#
# Transform edges to graph object
link_tbl=as_tbl_graph(links)
#link_tbl
#
# Total nodes
total.nodes=vcount(link_tbl)
#total.nodes
#
# Total terms
total.terms = nrow(drop_na(as.data.frame(nodes$Term)))
#total.terms
#
# Where start the entry names
start.entry = total.terms + 1 
#start.entry
#
# Only names of entry
only.entrys = nodes$Entry[start.entry:total.nodes]
#only.entrys
#
# xxx
id = rep("",total.terms)
#id
#
# Label names without terms
nodes.label = c(id,only.entrys)
#nodes.label
#
# Colors configuration
#
colors=c(brewer.pal(12,"Set3")[c(c(1:8),c(11:12))],
         brewer.pal(8,"Set2"),
         brewer.pal(12,"Paired"))
many.colors=rep(colors,12)
#plot(x=30:1, y=1:30, pch=19, cex=10, col=many.colors)
#dev.off()
proteins.color=rep("gray60",1000)#9D9D9D
colors.set=c(many.colors[nodes$num[1:total.terms]],
             proteins.color[nodes$num[start.entry:total.nodes]])
#
#
# Edges color by cluster
#
num.edges.colors=as.numeric(nodes$Freq[1:total.terms]*100)
col.edges.cluster=rep(colors[1:total.terms],num.edges.colors)
######################################
############# Plot 1 circle
######################################
#
# Control label size: loop for label size of geom_node_text(gnt)
if ((ecount(link_tbl) > 0) && (ecount(link_tbl) <= 50)) {
  gnt.label.size=4
  #print(gnt.label.size)
} else if ((ecount(link_tbl) >= 50.1) && (ecount(link_tbl) <= 75)) {
  gnt.label.size=3.5
  #print(gnt.label.size)
} else if ((ecount(link_tbl) >= 75.1) && (ecount(link_tbl) <= 100)) {
  gnt.label.size=3
  #print(gnt.label.size)
} else if (ecount(link_tbl) >= 100.1) {
  gnt.label.size=2.5
  #print(gnt.label.size)
}
#
# Control node size: loop for node size of geom_node_point(gnp)
if ((ecount(link_tbl) > 0) && (ecount(link_tbl) <= 50)) {
  gnp.node.size=5
  #print(gnp.node.size)
} else if ((ecount(link_tbl) >= 50.1) && (ecount(link_tbl) <= 75)) {
  gnp.node.size=3.7
  #print(gnp.node.size)
} else if ((ecount(link_tbl) >= 75.1) && (ecount(link_tbl) <= 100)) {
  gnp.node.size=2.5
  #print(gnp.node.size)
} else if (ecount(link_tbl) >= 100.1) {
  gnp.node.size=1.5
  #print(gnp.node.size)
}
#
{style.1="linear"
graph1=ggraph(link_tbl,
              layout = style.1,
              circular = TRUE) + 
  geom_edge_arc(colour=col.edges.cluster,
                width = 1,
                alpha=1,
                position = "identity",
                start_cap = circle(1.2, 'mm'),
                end_cap = circle(1.2, 'mm'))+
  geom_node_point(colour= colors.set,
                  size=gnp.node.size,
                  show.legend=T ) +
  geom_node_text(label = nodes$Entry,
                 position = "identity",
                 size = gnt.label.size,
                 repel = T,
                 fontface = "bold") +
  ggtitle('KEGG Pathways')+
  coord_fixed() +
  ggforce::theme_no_axes(base.theme = theme_minimal())
#
graph1=graph1 + theme(plot.title=element_text(hjust=0,
                                              size=10,face="bold"))
#
plot.1=ggdraw() + draw_plot(graph1, -0.17, 0, 1, 1)
#
# Save Plot 1
#
{freq=rev(as.numeric(nodes$Freq[1:total.terms]))
  png(file="./plots/Plot1_circle.png",
      width = 20,
      height = 10,
      units = "in",
      res=600,bg="white")
  par(mfrow=c(1,2),mar=c(20,0,3,20))
  plot(plot.1)
  pie(freq,
      col=rev(colors.set[1:total.terms]),
      border= NA,
      labels=as.character(as.character(freq)),
      clockwise = T,
      radius=0.4,
      bg="white",
      cex=1.2)
  plot(plot.1)
  pie(freq,
      col=rev(colors.set[1:total.terms]),
      border= NA,
      labels=as.character(as.character(freq)),
      clockwise = T,
      radius=0.4,
      bg="white",
      cex=1.2)
  par(mar=c(1, 3.8, 15, 0.2))
  legend(x=-.45,
         y=-.75,
         legend=unique(nodes$Term[c(1:total.terms)]),
         cex=1.2,
         bg="transparent",
         border = NA,
         bty = "n",
         pt.bg = colors,
         x.intersp=1.2,
         y.intersp=1.5,
         col=NA,
         text.font = 2,
         pch=22,
         pt.cex=2.5,
         ncol=1)
  legend(x=-.45,
         y=-.57,
         legend = "Pathways",
         pt.cex = NA,
         col=NA,
         pch=NA,
         pt.bg = NA,
         bty = "n",
         cex=1.5,
         ncol = 1,
         text.font = 2)
dev.off()
dev.off()}}
#
# 
######################################
############# Plot 2 random
######################################
#
# Control node size: loop for node size of geom_node_point(gnp)
if ((ecount(link_tbl) > 0) && (ecount(link_tbl) <= 75)) {
  style.2="kk"
  graph2=ggraph(link_tbl,
                layout = style.2) + 
    geom_edge_link(colour= col.edges.cluster,
                   width = 1,
                   position = "identity") +
    geom_node_point(colour=colors.set,
                    size=gnp.node.size) +
    geom_node_text(label = nodes$Entry,
                   position = "identity",
                   size = gnt.label.size,
                   repel = T,
                   fontface = "bold") +
    theme_graph(title_size = 11,
                title_face = "bold") +
    ggtitle("KEGG Pathways")+
    coord_fixed() +
    ggforce::theme_no_axes(base.theme = theme_minimal())
  #
  graph2=graph2 + theme(plot.title=element_text(hjust=0,
                                                size=10,face="bold"))
  #
  plot.2=ggdraw() + draw_plot(graph2, -0.19, 0, 1, 1)
  
  {freq=rev(as.numeric(nodes$Freq[1:total.terms]))
    png(file="./plots/Plot2_random.png",
        width = 20,
        height = 10,
        units = "in",
        res=600,
        bg="white")
    par(mfrow=c(1,2),mar=c(20,0,3,20))
    plot(plot.2)
    pie(freq,
        col=rev(colors.set[1:total.terms]),
        border= NA,
        labels=as.character(as.character(freq)),
        clockwise = T,
        radius=0.4,
        bg="white",
        cex=1.2)
    plot(plot.2)
    pie(freq,
        col=rev(colors.set[1:total.terms]),
        border= NA,
        labels=as.character(as.character(freq)),
        clockwise = T,
        radius=0.4,
        bg="white",
        cex=1.2)
    par(mar=c(1, 3.8, 15, 0.2))
    legend(x=-.45,
           y=-.75,
           legend=unique(nodes$Term[c(1:total.terms)]),
           cex=1.2,
           bg="transparent",
           border = NA,
           bty = "n",
           pt.bg = colors,
           x.intersp=1.2,
           y.intersp=1.5,
           col=NA,
           text.font = 2,
           pch=22,
           pt.cex=2.5,
           ncol=1)
    legend(x=-.45,
           y=-.57,
           legend = "Pathways",
           pt.cex = NA,
           col=NA,
           pch=NA,
           pt.bg = NA,
           bty = "n",
           cex=1.5,
           ncol = 1,
           text.font = 2)
    dev.off()
    dev.off()}
}else if ((ecount(link_tbl) > 75.1)){
  print("There are many nodes for random style")
}
#
#
######################################
############# Plot 3 circle, If fold change is included
######################################
#
#loop for data with fold change
#
# Color for fold change
#
#
if (colnames(nodes[,2]) == "Exp") {
  rbPal <- colorRampPalette(c('green','black','red'))
  xcol <- rbPal(100)[as.numeric(cut(nodes$Exp,breaks = 100))]
  xcol=xcol[!is.na(xcol)]
  colors.fold.change=c(colors[nodes$num[1:total.terms]],xcol)
  #
  #
  graph3=ggraph(link_tbl,layout = style.1,
                circular = TRUE) + 
    geom_edge_arc(colour = col.edges.cluster,
                  width = 1,
                  label_alpha = 2,
                  position = "identity",
                  start_cap = circle(1.2, 'mm'),
                  end_cap = circle(1.2, 'mm'))+
    geom_node_point(aes(colour=nodes$Exp),
                    size=gnp.node.size,
                    position = "identity") +
    scale_color_gradient(low = c("green","black"),high = "red")+
    geom_node_text(label = nodes$Entry,
                   position = "identity",
                   size = gnt.label.size,
                   repel = T,
                   fontface = "bold",
                   hjust='outward') +
    ggtitle("KEGG Pathways")+
    labs(colour="Fold change",size = "Log10 Total proteins")+
    theme_graph(title_size = 11, title_face = "bold") + 
    scale_fill_manual(name = "Source")+
    coord_fixed() +
    ggforce::theme_no_axes(base.theme = theme_minimal())
  #
  graph3 = graph3 + geom_node_point(colour=colors.fold.change,
                                    size=gnp.node.size)
  #
  graph3 = graph3 + theme(legend.key.width=unit(1,"cm"),
                          legend.key.height=unit(.7,"cm"), 
                          legend.position = c(1.11, .9),
                          legend.title=element_text(size=16,face="bold"),
                          legend.text=element_text(size=13,face="bold"),
                          plot.title=element_text(hjust=.1,
                                                  size=10,
                                                  margin=margin(b =0,unit="pt"),
                                                  face="bold"))
  #
  plot.3=ggdraw() + draw_plot(graph3, -0.17, 0, 1, 1)
  #
  {freq=rev(as.numeric(nodes$Freq[1:total.terms]))
  png(file="./plots/Plot3_circle_foldchange.png",
      width = 20,
      height = 10,
      units = "in",
      res=600,bg="white")
  par(mfrow=c(1,2),mar=c(20,0,3,20))
  plot(plot.3)
  pie(freq,
      col=rev(colors.set[1:total.terms]),
      border= NA,
      labels=as.character(as.character(freq)),
      clockwise = T,
      radius=0.4,
      bg="white",
      cex=1.2)
  plot(plot.3)
  pie(freq,
      col=rev(colors.set[1:total.terms]),
      border= NA,
      labels=as.character(as.character(freq)),
      clockwise = T,
      radius=0.4,
      bg="white",
      cex=1.2)
  par(mar=c(1, 3.8, 15, 0.2))
  legend(x=-.45,
         y=-.75,
         legend=unique(nodes$Term[c(1:total.terms)]),
         cex=1.2,
         bg="transparent",
         border = NA,
         bty = "n",
         pt.bg = colors,
         x.intersp=1.2,
         y.intersp=1.5,
         col=NA,
         text.font = 2,
         pch=22,
         pt.cex=2.5,
         ncol=1)
  legend(x=-.45,
         y=-.57,
         legend = "Pathways",
         pt.cex = NA,
         col=NA,
         pch=NA,
         pt.bg = NA,
         bty = "n",
         cex=1.5,
         ncol = 1,
         text.font = 2)
  dev.off()
  dev.off()}
} else {
  print("There are no fold change values")
} 
#
#
######################################
############# Plot 4 random, If fold change is included
######################################
#
#
if ((ecount(link_tbl) > 0) && (ecount(link_tbl) <= 75)) {
  if (colnames(nodes[,2]) == "Exp") {
    print("There are fold change values")
    graph4 = ggraph(link_tbl,layout = style.2) + 
      scale_color_gradient(low = c("green","black"), high = "red")+
      geom_edge_link(colour= col.edges.cluster,
                     width = 1,
                     position = "identity") +
      geom_node_point(aes(colour= nodes$Exp),
                      size=gnp.node.size,
                      position = "identity") +
      theme_graph(title_size = 11,
                  title_face = "bold") + 
      ggtitle("KEGG Pathways")+
      labs(colour="Fold change",size = "Log10 Total proteins")+
      coord_fixed() +
      ggforce::theme_no_axes(base.theme = theme_minimal())
    #
    graph4 = graph4 + geom_node_point(colour=colors.fold.change,
                                      size=gnp.node.size)+
      geom_node_text(label = nodes$Entry,
                     position ="identity",
                     size = gnt.label.size,
                     repel = T,
                     fontface = "bold")
    #
    graph4 = graph4 + theme(legend.key.width=unit(1,"cm"),
                            legend.key.height=unit(0.7,"cm"),
                            legend.position = c(1.11,.9),
                            legend.title=element_text(size=16,face="bold"),
                            legend.text=element_text(size=13,face="bold"),
                            plot.title=element_text(hjust=.1,
                                                    size=10,
                                                    margin=margin(b =0,unit="pt"),
                                                    face="bold"))
    #
    plot.4=ggdraw() + draw_plot(graph4, -0.19, 0, 1, 1)
    #
    {freq=rev(as.numeric(nodes$Freq[1:total.terms]))
      png(file="./plots/Plot4_random_foldchange.png",
          width = 20,
          height = 10,
          units = "in",
          res=600,bg="white")
      par(mfrow=c(1,2),mar=c(20,0,3,20))
      plot(plot.4)
      pie(freq,
          col=rev(colors.set[1:total.terms]),
          border= NA,
          labels=as.character(as.character(freq)),
          clockwise = T,radius=0.4,
          bg="white",cex=1.2)
      plot(plot.4)
      pie(freq,
          col=rev(colors.set[1:total.terms]),
          border= NA,
          labels=as.character(as.character(freq)),
          clockwise = T,radius=0.4,
          bg="white",cex=1.2)
      par(mar=c(1, 3.8, 15, 0.2))
      legend(x=-.45,y=-.75,
             legend=unique(nodes$Term[c(1:total.terms)]),
             cex=1.2,bg="transparent",border = NA,
             bty = "n",pt.bg = colors,x.intersp=1.2,
             y.intersp=1.5,col=NA,
             text.font = 2,pch=22,
             pt.cex=2.5,ncol=1)
      legend(x=-.45,y=-.57,legend = "Pathways",
             pt.cex = NA,col=NA,
             pch=NA,pt.bg = NA,
             bty = "n",cex=1.5,
             ncol = 1,text.font = 2)
      dev.off()
      dev.off()}
  } else {
    print("There are no fold change values")
  } 
} else {
  print("There are many nodes for random style")
}
#
#
######################################
############# Plot 5 Chord Plot, without fold change
######################################
#
# Control label size: loop for chord label size
if ((ecount(link_tbl) > 0) && (ecount(link_tbl) <= 50)) {
  chord.node.size=1
  #print(chord.node.size)
} else if ((ecount(link_tbl) >= 50.1) && (ecount(link_tbl) <= 75)) {
  chord.node.size=.85
  #print(chord.node.size)
} else if ((ecount(link_tbl) >= 75.1) && (ecount(link_tbl) <= 100)) {
  chord.node.size=.68
  #print(chord.node.size)
} else if (ecount(link_tbl) >= 100.1) {
  chord.node.size=0.5
  #print(chord.node.size)
}
#
{png(file="./plots/Plot5_ChordPlot.png",width = 25,height = 10,
     units = "in",
     res=600,bg="white")
  simple.chord = function() {
    chordDiagram(links,grid.col = colors.set,
                 transparency = .3, directional = 0,
                 symmetric = T,
                 annotationTrack = "grid",
                 preAllocateTracks = 1.5)
    circos.trackPlotRegion(track.index = 1,
                           panel.fun = function(x, y){
                             xlim = get.cell.meta.data("xlim")
                             ylim = c(-.1, 1)
                             sector.name = get.cell.meta.data("sector.index")
                             circos.text(mean(xlim), ylim[1] + .1,
                                         sector.name, facing = "clockwise", 
                                         niceFacing = TRUE,
                                         adj = c(0, .5), col = "black",
                                         cex =chord.node.size,
                                         font = par("font"))
                             circos.par(track.margin=c(1,5))
                           }, bg.border = NA)
    legend(x=1.03,y=.42, 
           legend=unique(paste(nodes$Term[1:total.terms]," _ ",
                               "(",nodes$Freq[1:total.terms],")")),
           cex=1,bg="transparent",border = NA,pt.bg=colors,bty = "n",
           x.intersp=1.5,y.intersp=1.5,col=NA,text.font = 2,pch=22,
           pt.cex=2)
   text(1.05, .45, "Pathway", adj = c(0,0), cex = 1.5,font = 2)
    circos.clear()}
  simple.chord()
  dev.off()}

#
#
######################################
############# Plot 6 Chord Plot, If fold change is included
######################################
#
#
if (colnames(nodes[,2]) == "Exp") {
  print("si hay datos de expresion")
  #
  min.value.exp=round(min(nodes$Exp[!is.na(nodes$Exp)])) ## obtain negative value
  #min.value.exp
  max.value.exp=round(min.value.exp*-1) ## obtain positive value
  #max.value.exp
  #
  legend_image <- as.raster(matrix(rbPal(100), ncol=1))
  #
  png(file="./plots/Plot6_ChordPlot_FoldChange.png",width = 25,height = 10,
       units = "in",
       res=600,bg="white")
    simple.chord = function() {
      chordDiagram(links,grid.col = colors.fold.change,
                   transparency = .3, directional = 0,
                   symmetric = T,
                   annotationTrack = "grid",
                   preAllocateTracks = 1.5)
      circos.trackPlotRegion(track.index = 1,
                             panel.fun = function(x, y){
                               xlim = get.cell.meta.data("xlim")
                               ylim = c(-.1, 1)
                               sector.name = get.cell.meta.data("sector.index")
                               circos.text(mean(xlim), ylim[1] + .1,
                                           sector.name, facing = "clockwise", 
                                           niceFacing = TRUE,
                                           adj = c(0, .5), col = "black",
                                           cex =chord.node.size,
                                           font = par("font"))
                               circos.par(track.margin=c(1,5))
                             }, bg.border = NA)
      legend(x=1.03,y=0, 
             legend=unique(paste(nodes$Term[1:total.terms]," _ ",
                                 "(",nodes$Freq[1:total.terms],")")),
             cex=.9,bg="transparent",border = NA,pt.bg=colors,bty = "n",
             x.intersp=1.5,y.intersp=1.5,col=NA,text.font = 2,pch=22,
             pt.cex=2)
      ##  positions             a=x, b= +y,c= x,d= -y
      rasterImage(legend_image, 1.09, .4, 1.15,.16) 
      text(x=1.19, y = seq(0.4,.16,l=5),
           labels = seq(max.value.exp,
                        min.value.exp,l=5),cex=.8,font = 2)
      text(1.05, .45, "Fold change", adj = c(0,0), cex = 1.2,font = 2)
      text(1.05, .03, "Pathway", adj = c(0,0), cex = 1.2,font = 2)
      circos.clear()}
    simple.chord()
    dev.off()
} else {
  print("There are no fold change values")
}
#
#
######################################
############# Plot 7 UpSet Plot, without fold change
######################################
#
#
{val=length(only.entrys) # nrow number (total genes)
#
ups=as.data.frame.matrix(t(table(links)))
ups[["log.Pval"]]=c(-log10(nodes$P[1:total.terms]),rep(NA,val-total.terms))
ups[["term"]]=c(nodes$Entry[1:total.terms],rep(NA,(val-total.terms)))
#
max(ups$log.Pval[1:total.terms]) + 5
#
# Bar plot with P value
#ggplot(ups) +
#  geom_bar(aes(x = term, weight = log.Pval),fill="wheat3") +
#  geom_line(aes(x = c(num[1:total.terms],rep(NA,35)), y = -log10(0.05),color=-log10(0.05)),color = "black", size = .5 )+
#  geom_point(aes(x = c(num[1:total.terms],rep(NA,35)), y = -log10(0.05),color=-log10(0.05)),
#             position = position_dodge(0.3),color = "black", size = 5)+
#  theme(axis.text.x = element_text(face="bold",size=8, angle=90),
#        axis.text.y = element_text(face="bold",size=8),
#        axis.title.y = element_text(size=11,face="bold"),
#        axis.ticks.x = element_blank(),
#        axis.line.x = element_blank(),
#        axis.title.x = element_blank(),
#        legend.title = element_text(colour = "black",size=5,face="bold"),
#        legend.text = element_text(colour = "black",size =5))+#coord_fixed(ratio=5)+
#  scale_y_continuous("-log10(P-value)",limit = c(0,round(max(ups$log.Pval[0:total.terms]+2))),
#                     breaks=seq(0,round(max(ups$log.Pval[0:total.terms]+3),digits = 0),3)) +
#  scale_x_discrete("",limit=nodes$Entry[1:total.terms])+
#  annotate("text", x = length(ups$log.Pval[1:total.terms]), y = min(ups$log.Pval[1:total.terms])+3, label = "P = 0.05 (-log10)",colour = "black", size = 3.5)+
#  annotate("pointrange", x = length(ups$log.Pval[1:total.terms]), y =  min(ups$log.Pval[1:total.terms])+2.3, ymin = 0, ymax = 0,
#           colour = "black", size = .5)
#
# 
if ((length(nodes$Entry[1:total.terms]) > 0) && (length(nodes$Entry[1:total.terms]) <= 9)) {
  text.annotate.size=3
} else if (length(nodes$Entry[1:total.terms]) >= 9.1) {
  text.annotate.size=2
}
# Function for bar plot
myplot2 <- function(mydata, term, log.Pval) {
  plot <- (ggplot(ups) +
             geom_bar(aes(x = term, weight = log.Pval),fill="darkgoldenrod2") +
             geom_line(aes(x = c(nodes$num[1:total.terms],rep(NA,length(rownames(ups))-length(nodes$Term[1:total.terms]))), y = -log10(0.05),color=-log10(0.05)),color = "blue4", size = .3 )+
             geom_point(aes(x = c(nodes$num[1:total.terms],rep(NA,length(rownames(ups))-length(nodes$Term[1:total.terms]))), y = -log10(0.05),color=-log10(0.05)),
                        position = position_dodge(0.3),color = "blue4", size = 1.5)+
             theme(axis.text.x = element_text(face="bold",size=8, angle=90),
                   axis.text.y = element_text(face="bold",size=8),
                   axis.title.y = element_text(size=11,face="bold"),
                   axis.ticks.x = element_blank(),
                   axis.line.x = element_blank(),
                   axis.title.x = element_blank(),
                   legend.title = element_text(colour = "black",size=5,face="bold"),
                   legend.text = element_text(colour = "black",size =5))+#coord_fixed(ratio=5)+
             scale_y_continuous("-log10(P-value)",limit = c(0,round(max(ups$log.Pval[0:total.terms]+2))),
                                breaks=seq(0,round(max(ups$log.Pval[0:total.terms]+3),digits = 0),3)) +
             scale_x_discrete("",limit=nodes$Entry[1:total.terms])+
             annotate("text", x = length(ups$log.Pval[1:total.terms])-1, 
                      y = min(ups$log.Pval[1]),
                      label = 'bold("P = 0.05 (-log10)")',colour = "black", size = text.annotate.size, parse = TRUE)+
             annotate("pointrange", x = length(ups$log.Pval[1:total.terms])-1, y =  min(ups$log.Pval[1])-1.5, ymin = 0, ymax = 0,
                      colour = "blue4", size = .3))
}
#
if ((length(unique(links$Entry)) > 0) && (length(unique(links$Entry)) <= 6)) {
  set.name.size=1.5
  #print(gnp.node.size)
} else if ((length(unique(links$Entry)) >= 6.1) && (length(unique(links$Entry)) <= 9)) {
  set.name.size=1.25
  #print(gnp.node.size)
} else if ((length(unique(links$Entry)) >= 9.1) && (length(unique(links$Entry)) <= 13)) {
  set.name.size=1
  #print(gnp.node.size)
} else if (length(unique(links$Entry)) >= 13.1) {
  set.name.size=0.8
  #print(gnp.node.size)
}
#
{png(file="./plots/Polt7_UpSetR.png",width = 10,height = 7,units = "in",res=700,bg="white")
upset(ups,sets=nodes$Entry[1:total.terms],
      sets.bar.color = "darkgoldenrod2",order.by ="freq",empty.intersections = NULL,
      point.size=2,mainbar.y.label="KEGG Pathways Functional Annotation",sets.x.label = "",
      main.bar.color="black",matrix.color="black",shade.color="wheat3",
      line.size=0.5,show.numbers = "yes",group.by = "degree",
      matrix.dot.alpha = 0.5,mb.ratio = c(0.7, 0.3),
      #c(intersection size title, intersection size tick labels, set size title, set size tick labels, set names, numbers above bars)
      text.scale = c(1.5,1.5,0,1.5,set.name.size,1.5), 
      scale.sets = "identity",scale.intersections = "identity",
      shade.alpha = 0.35,set_size.angles = 0,
      attribute.plots = list(gridrows = 45, plots = list(list(plot = myplot2,x="term",y = "log.Pval", queries = F)),ncols=2))
dev.off()}
}
#
#
######################################
############# Plot 8 UpSet Plot, If fold change is included
######################################
#
if (colnames(nodes[,2]) == "Exp") {
  print("si hay datos de expresion")
  val=length(only.entrys) # nrow number (total genes)
  #
  ups=as.data.frame.matrix(t(table(links)))
  ups[["num"]]=c(1:val)
  ups[["exp"]]=nodes$Exp[start.entry:length(nodes$Entry)]
  ups[["log.Pval"]]=c(-log10(nodes$P[1:total.terms]),rep(NA,val-total.terms))
  ups[["term"]]=c(nodes$Entry[1:total.terms],rep(NA,(val-total.terms)))
  
  #
  # Scatterplot, fold change data
  #ggplot(ups, aes(x=c(1:val),exp))+ geom_point(aes(color=exp),size=4) + 
  # scale_color_gradient(low = c("green","black"), high = "red")+
  #coord_fixed(ratio=5)+ theme(axis.text.x = element_text(size = 8,face="bold"),
  #                           axis.text.y = element_text(size = 8,face="bold"),
  #                          axis.title.x = element_text(size=11,face="bold"),
  #                         axis.title.y = element_text(size=11,face="bold"),
  #                        legend.text = element_text(size=8))+
  #scale_y_continuous("Fold change",limit = c(-5,5)) + scale_x_continuous("Proteins",breaks=seq(1,val,20)) + labs(color = "")
  #
  # Function for configuration plot
  myplot <- function(mydata, x, y) {
    plot <- (ggplot(data = mydata, aes_string(x = c(1:val), y = y)) + 
               geom_point(aes(color=exp),size=2) +
               scale_color_gradient(low = c("green","black"), high = "red")+
               theme(axis.text.x = element_text(size = 8,face="bold"),
                     axis.text.y = element_text(size = 8,face="bold"),
                     axis.title.x = element_text(size=11,face="bold"),
                     axis.title.y = element_text(size=11,face="bold"),
                     legend.text = element_text(size=8,face="bold"))+
               coord_fixed(ratio=7)+
               scale_y_continuous("Fold change",limit = c(-5,5)) +
               labs(color = "") +
               scale_x_continuous("Protein",breaks=seq(1,val,20)))
  }
  #
  #
if ((length(nodes$Entry[1:total.terms]) > 0) && (length(nodes$Entry[1:total.terms]) <= 9)) {
  text.annotate.size=3
} else if (length(nodes$Entry[1:total.terms]) >= 9.1) {
  text.annotate.size=2
}
  # Bar plot with P value
  # Function for bar plot
  myplot2 <- function(mydata, term, log.Pval) {
    plot <- (ggplot(ups) +
               geom_bar(aes(x = term, weight = log.Pval),fill="darkgoldenrod2") +
               geom_line(aes(x = c(nodes$num[1:total.terms],rep(NA,length(rownames(ups))-length(nodes$Term[1:total.terms]))), y = -log10(0.05),color=-log10(0.05)),color = "blue4", size = .3 )+
               geom_point(aes(x = c(nodes$num[1:total.terms],rep(NA,length(rownames(ups))-length(nodes$Term[1:total.terms]))), y = -log10(0.05),color=-log10(0.05)),
                          position = position_dodge(0.3),color = "blue4", size = 1.5)+
               theme(axis.text.x = element_text(face="bold",size=8, angle=90),
                     axis.text.y = element_text(face="bold",size=8),
                     axis.title.y = element_text(size=11,face="bold"),
                     axis.ticks.x = element_blank(),
                     axis.line.x = element_blank(),
                     axis.title.x = element_blank(),
                     legend.title = element_text(colour = "black",size=5,face="bold"),
                     legend.text = element_text(colour = "black",size =5))+#coord_fixed(ratio=5)+
               scale_y_continuous("-log10(P-value)",limit = c(0,round(max(ups$log.Pval[0:total.terms]+2))),
                                  breaks=seq(0,round(max(ups$log.Pval[0:total.terms]+3),digits = 0),3)) +
               scale_x_discrete("",limit=nodes$Entry[1:total.terms])+
               annotate("text", x = length(ups$log.Pval[1:total.terms])-1, 
                        y = min(ups$log.Pval[1]),
                        label = 'bold("P = 0.05 (-log10)")',colour = "black", size = text.annotate.size, parse = TRUE)+
               annotate("pointrange", x = length(ups$log.Pval[1:total.terms])-1, y =  min(ups$log.Pval[1])-1.5, ymin = 0, ymax = 0,
                        colour = "blue4", size = .5))
  }
  #
  if ((length(unique(links$Entry)) > 0) && (length(unique(links$Entry)) <= 6)) {
    set.name.size=1.5
    #print(gnp.node.size)
  } else if ((length(unique(links$Entry)) >= 6.1) && (length(unique(links$Entry)) <= 9)) {
    set.name.size=1.25
    #print(gnp.node.size)
  } else if ((length(unique(links$Entry)) >= 9.1) && (length(unique(links$Entry)) <= 13)) {
    set.name.size=1
    #print(gnp.node.size)
  } else if (length(unique(links$Entry)) >= 13.1) {
    set.name.size=0.8
    #print(gnp.node.size)
  }
  #
  png(file="./plots/Polt8_UpSetR_FoldChange.png",width = 10,height = 7,units = "in",res=700,bg="white")
    upset(ups,sets=nodes$Entry[1:total.terms],
          sets.bar.color = "darkgoldenrod2",order.by ="freq",empty.intersections = NULL,
          point.size=2,mainbar.y.label="KEGG Pathways Functional Annotation",sets.x.label = "",
          main.bar.color="black",matrix.color="black",shade.color="wheat3",
          line.size=0.5,show.numbers = "yes",group.by = "degree",
          matrix.dot.alpha = 0.5,mb.ratio = c(0.7, 0.3),
          #c(intersection size title, intersection size tick labels, set size title, set size tick labels, set names, numbers above bars)
          text.scale = c(1.5,1.5,0,1.5,set.name.size,1.5), 
          scale.sets = "identity",scale.intersections = "identity",
          shade.alpha = 0.35,set_size.angles = 0,
          attribute.plots = list(gridrows = 45, plots = list(list(plot = myplot2,x = "term", y = "log.Pval", queries = T),
                                                             list(plot = myplot ,x="num",y = "exp", queries = T)),ncols=2))
    dev.off()

} else {
  print("There are no fold change values")
} 
dev.off()
############

