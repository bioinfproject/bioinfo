#
####### Packages
#
{
  library(tidyverse)
  library(tidygraph)
  library(viridis)
  library(ggraph)
  library(circlize)
  library(RColorBrewer)
  library(igraph)
  library(cowplot)
  library(grid)
  library(networkD3)
  library(UpSetR)
  library(ComplexHeatmap)
  library(gridBase)
  library(reshape2)
}
#
#
nodes = read_csv("qwertyuiop")
links = read_csv("asdfghjkl")
#
# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#
# |||||||||     Definition of ranges and positions for functions
#
# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#
{
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
  colors=c(brewer.pal(12,"Set3")[c(c(3:8),c(11:12))],
           brewer.pal(8,"Set2"),
           brewer.pal(12,"Paired"))
  many.colors=rep(colors,12)
  #plot(x=30:1, y=1:30, pch=19, cex=10, col=many.colors)
  #dev.off()
  #proteins.color=rep("blue",1000)#9D9D9D #D7BC9C
  #col.set=c(many.colors[nodes$num[1:total.terms]],
  #             proteins.color[nodes$num[start.entry:total.nodes]])
  #
  #
  # Edges color by cluster
  #
  num.edges.colors=as.numeric(nodes$Freq[1:total.terms]*100)
  col.edges.cluster=rep(many.colors[1:total.terms],num.edges.colors)
  #
  ## Add attributes tu network
  E(link_tbl)$GO = links$GO
  V(link_tbl)$Label = c(as.character(nodes$Term[1:total.terms]),
                        as.character(nodes$Entry[start.entry:total.nodes]))
  V(link_tbl)$Entry = c(as.character(rep(NA,total.terms)),
                        as.character(nodes$Entry[start.entry:total.nodes]))
  #
  #
  # ///////////////////////////////////////////
  correl.width = function (x) {
    if (x <= 0) {
      x/x * 0.5
    } else {
      x
    }
  }
  #
  correl.point = function (x) {
    if (x <= 0) {
      x/x * 3
    } else {
      x
    }
  }
  #
  correl.text = function (x) {
    if (x <= 0) {
      x/x
    } else {
      x
    }
  }
  # ///////////////////////////////////////////
  #
  #
}
#
#
if (ncol(nodes) == 7) {
  print("there are expression values")
  #
  V(link_tbl)$logFC = nodes$Exp
  #
  rbPal <- colorRampPalette(c('green','black','red'))
  legend_image <- as.raster(matrix(rbPal(100), ncol=1))
  xcol <- rbPal(100)[as.numeric(cut(nodes$Exp,breaks = 100))]
  xcol=xcol[!is.na(xcol)]
  colors.fold.change=c(many.colors[nodes$num[1:total.terms]],xcol)
  #
  min.value.exp=round(min(nodes$Exp[!is.na(nodes$Exp)]),digits=0)
  max.value.exp=round(max(nodes$Exp[!is.na(nodes$Exp)]),digits=0)
  #___________________________________________________________________________________________ Simple plots
  #
  # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  #
  # |||||||||     Simple Plots
  #
  # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  layouts=c("nicely","kk","dh")
  #
  if ((total.terms > 1) && (total.terms <= 20)) {
    width.png=20
    print(paste("el ancho de la figura es:",width.png))
  } else if ((total.terms >= 21) && (total.terms <= 40)) {
    width.png=25
    print(paste("el ancho de la figura es:",width.png))
  } else if (total.terms > 41) {
    width.png=30
  }
  #
  for (i in layouts){
    graf1 = ggraph(link_tbl,layout = "igraph", algorithm = i) +
      geom_edge_link(aes(colour = GO),alpha = 1,width = correl.width(-0.01296*length(only.entrys)+1.798))+ #<---------------------
    geom_node_point(aes(colour= logFC),
                    size= correl.point(-0.04259*length(only.entrys)+5.980),show.legend = NA, #<---------------------
                    position = "identity")+
      scale_edge_colour_manual(values=colors.fold.change,name = "ASPECT")+
      geom_node_text(aes(label=Label),
                     position ="identity",
                     size = correl.text(-0.02778*length(only.entrys)+4.639), #<---------------------
                     repel = T,
                     #vjust=3,
                     #hjust=.3,
                     fontface = "bold")+
      theme_classic() +
      theme(legend.key.width=unit(.4,"cm"),
            legend.key.height=unit(.4,"cm"),
            legend.title=element_text(size=10,face="bold"),
            legend.text=element_text(size=8),
            axis.line=element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.title = element_blank())
    graf1 = graf1 + scale_color_gradient(low = c("green","black"), high = "red") +
      labs(caption = paste("NeVOmics:",sep=" ",format(Sys.time(),"%a(%d)/%b/%Y   %X")))
    {
      png(file=paste0("zxcvbnmNeVOmics_Plot_Network_Des_Layout_",i,".png"),
          width = width.png,  #<---------------------
          height = 15, 
          units = "cm",
          res=900,bg="white")
      plot(graf1)
      dev.off()}
  }
  #
  for (i in layouts){
    graf1 = ggraph(link_tbl,layout = "igraph", algorithm = i) +
      geom_edge_link(aes(colour = GO),alpha = 1,width = correl.width(-0.01296*length(only.entrys)+1.798))+ #<---------------------
    geom_node_point(aes(colour= logFC),
                    size= correl.point(-0.04259*length(only.entrys)+5.980),show.legend = NA, #<---------------------
                    position = "identity")+
      scale_edge_colour_manual(values=colors.fold.change,name = "ASPECT")+
      geom_node_text(aes(label=name),
                     position ="identity",
                     size = correl.text(-0.02778*length(only.entrys)+4.639), #<---------------------
                     repel = T,
                     #vjust=3,
                     #hjust=.3,
                     fontface = "bold")+
      theme_classic() +
      theme(legend.key.width=unit(.4,"cm"),
            legend.key.height=unit(.4,"cm"),
            legend.title=element_text(size=10,face="bold"),
            legend.text=element_text(size=8),
            axis.line=element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.title = element_blank())
    graf1 = graf1 + scale_color_gradient(low = c("green","black"), high = "red") +
      labs(caption = paste("NeVOmics:",sep=" ",format(Sys.time(),"%a(%d)/%b/%Y   %X")))
    {
      png(file=paste0("zxcvbnmNeVOmics_Plot_Network_Id_Layout_",i,".png"),
          width = width.png,  #<---------------------
          height = 15, 
          units = "cm",
          res=900,bg="white")
      plot(graf1)
      dev.off()}
  }
  #
  for (i in layouts){
    graf1 = ggraph(link_tbl,layout = "igraph", algorithm = i) +
      geom_edge_link(aes(colour = GO),alpha = 1,width = correl.width(-0.01296*length(only.entrys)+1.798))+ #<---------------------
    geom_node_point(aes(colour= logFC),
                    size= correl.point(-0.04259*length(only.entrys)+5.980),show.legend = NA, #<---------------------
                    position = "identity")+
      scale_edge_colour_manual(values=colors.fold.change,name = "ASPECT")+
      geom_node_text(aes(label=Entry),
                     position ="identity",
                     size = correl.text(-0.02778*length(only.entrys)+4.639), #<---------------------
                     repel = T,
                     #vjust=3,
                     #hjust=.3,
                     fontface = "bold")+
      theme_classic() +
      theme(legend.key.width=unit(.4,"cm"),
            legend.key.height=unit(.4,"cm"),
            legend.title=element_text(size=10,face="bold"),
            legend.text=element_text(size=8),
            axis.line=element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.title = element_blank())
    graf1 = graf1 + scale_color_gradient(low = c("green","black"), high = "red") +
      labs(caption = paste("NeVOmics:",sep=" ",format(Sys.time(),"%a(%d)/%b/%Y   %X")))
    {
      png(file=paste0("zxcvbnmNeVOmics_Plot_Network_Layout_",i,".png"),
          width = width.png,  #<---------------------
          height = 15, 
          units = "cm",
          res=900,bg="white")
      plot(graf1)
      dev.off()}
  }
  #
  #
  #
  #\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
  if (length(only.entrys) > 80) {
    modelos = c('nicely','kk','dh', 'drl')
    V(link_tbl)$Etiqueta = c(as.character(nodes$Term[1:total.terms]),
                             as.character(rep(NA,total.nodes - total.terms)))
    for (i in modelos){
      graf1 = ggraph(link_tbl,layout = "igraph", algorithm = i) +
        geom_edge_link(aes(colour = GO),alpha = 1,width = correl.width(-0.01296*length(only.entrys)+1.798))+ #<---------------------
      geom_node_point(aes(colour= logFC),
                      size= 3,#correl.point(-0.04259*length(only.entrys)+5.980),show.legend = NA, #<---------------------
                      position = "identity")+
        scale_edge_colour_manual(values=colors.fold.change,name = "ASPECT")+
        geom_node_text(aes(label=Etiqueta),
                       position ="identity",
                       size = 3, #<---------------------
                       repel = T,
                       #vjust=3,
                       #hjust=.3,
                       fontface = "bold")+
        theme_classic() +
        theme(legend.key.width=unit(.4,"cm"),
              legend.key.height=unit(.4,"cm"),
              legend.title=element_text(size=10,face="bold"),
              legend.text=element_text(size=8),
              axis.line=element_blank(),
              axis.text = element_blank(),
              axis.ticks = element_blank(),
              axis.title = element_blank())
      graf1 = graf1 + scale_color_gradient(low = c("green","black"), high = "red") +
        labs(caption = paste("NeVOmics:",sep=" ",format(Sys.time(),"%a(%d)/%b/%Y   %X")))
      {
        png(file=paste0("zxcvbnmNeVOmics_Plot_Network_Des_Layout_XL_",i,".png"),
            width = width.png,  #<---------------------
            height = 15, 
            units = "cm",
            res=900,bg="white")
        plot(graf1)
        dev.off()}
    }
  }
  #\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
  #___________________________________________________________________________________________ Complex plots
  #
  if (total.terms <= 30) {
    # --------------------------------------------------------------------- one bar configuration
    #
    # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    #
    # |||||||||      Bar plot
    #
    # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    #
    {
      for.bar = data.frame(c(data.frame("Entry"=c(nodes$Entry[1:total.terms],rep(NA,31-total.terms))),
                             data.frame("Freq"=c(nodes$Freq[1:total.terms],rep(NA,31-total.terms))),
                             data.frame("Fdr"=c(nodes$FDR[1:total.terms],rep(NA,31-total.terms))),
                             data.frame("Term"=c(nodes$Term[1:total.terms],rep(NA,31-total.terms))),
                             data.frame("P"=c(nodes$P[1:total.terms],rep(NA,31-total.terms)))))
      for.bar["log.pval"] = -log10(for.bar$P)
      for.bar["log.fdr"] = -log10(for.bar$Fdr)
    }
    ###
    #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    {
      barplot.enrichment = ggplot(for.bar) +
        geom_bar(aes(x = Entry, weight = log.pval),fill=rev(many.colors[1:(total.terms+1)]),width=.7) +
        geom_line(aes(x = c(1:31), y = -log10(0.05),
                      color=-log10(0.05)),
                  color = c(rep("transparent",31-total.terms),
                            rep("black",total.terms)), size = .2 )+
        geom_point(aes(x = c(1:31), y = -log10(0.05),color=-log10(0.05)),
                   position = position_dodge(0.3),
                   color = c(rep("transparent",31-total.terms),
                             rep("black",total.terms)), size = 1)+
        geom_line(aes(x = c(1:31), y = rev(log.fdr),
                      color=log.fdr),
                  color = c(rep("transparent",31-total.terms),
                            rep("blue",total.terms)), size = .2 )+
        geom_point(aes(x = c(1:31), y = rev(log.fdr),color=log.fdr),
                   position = position_dodge(0.3),shape=17,
                   color = c(rep("transparent",31-total.terms),
                             rep("blue",total.terms)), size = 1)+
        #scale_y_continuous(position = "right")+
        coord_flip()+
        scale_x_discrete(limits=rev(for.bar$Entry))+
        geom_text(aes(x=Entry,y=log.pval+.35,label=Freq),size=2.5,fontface ="bold")+ #<---------------------
      #ggtitle("-log(P-value) and\n -log(Corrected P-value)")+
      theme(axis.title.x =element_blank(),# element_text(colour="black",size=27,face="bold",hjust = 0,vjust = 10),
            axis.title.y = element_blank(),
            #plot.title = element_text(vjust = 4,hjust=.08,size = 12),
            axis.text.y = element_text(colour=c("white",rep("black",30)),size=8,
                                       face="bold",hjust = 1.3,vjust = .5), #element_text(size = 8,color="black"),#tama?o y color de texto de eje x
            axis.text.x = element_text(colour="black",size=7,face="bold",vjust = 0),
            axis.ticks.length =unit(.1, "cm"), #reduce la longitud de las marcas en la linea de escala de los ejes x y
            axis.ticks = element_line(color = "black",size = 1),
            axis.ticks.y.left = element_line(colour = "black",size=.5,unit(0, "cm")),
            #axis.ticks.x.top = element_line(colour = "black",size=.5,unit(3, "cm")),
            axis.line.y = element_blank(), axis.title.x.top =  element_blank(), 
            plot.margin = margin(1.5,1,0.5,0, "cm"))+
        scale_y_continuous(position = "right",
                           limit = c(0,round(max(for.bar$log.pval[1:total.terms]+1.5))),
                           breaks=seq(0,round(max(for.bar$log.pval[1:total.terms]+1.5),digits = 0),1))+
        labs(title = expression(bold("ASPECT")),
             subtitle = expression(bold("-log10("*italic("P")*"-value & adj."*italic("P")*"-value)")))+
        theme(plot.title = element_text(vjust = -.7,hjust=1.5,size = 12),
              plot.subtitle = element_text(vjust = -.7,hjust=.6,size = 8))
    }
    ###
    #
    # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    #
    # |||||||||      Complex Plots
    #
    # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    #
    layouts=c("nicely","kk","dh")
    #
    ###
    for (i in layouts){
      graf2 = ggraph(link_tbl,layout = "igraph",algorithm = i) +
        geom_edge_link(colour = col.edges.cluster,alpha = 1,
                       width = correl.width(-0.01296*length(only.entrys)+1.798))+ #<---------------------
      geom_node_point(aes(colour= logFC),
                      size= correl.point(-0.04259*length(only.entrys)+5.980), #<---------------------
                      position = "identity")+
        geom_node_text(aes(label=Label),
                       position ="identity",
                       size = correl.text(-0.02778*length(only.entrys)+4.639), #<---------------------
                       repel = T,
                       #vjust=3,
                       #hjust=.3,
                       fontface = "bold")+
        theme_classic()
      graf2 = graf2 +  theme(legend.key.width=unit(.4,"cm"),
                             legend.key.height=unit(.4,"cm"),
                             legend.position = c(1.06, .9),
                             legend.title=element_text(size=10,face="bold"),
                             legend.text=element_text(size=8),
                             axis.line=element_blank(),
                             axis.text = element_blank(),
                             axis.ticks = element_blank(),
                             axis.title = element_blank())
      graf2 = graf2 + scale_color_gradient(low = c("green","black"), high = "red") +
        labs(caption = paste("NeVOmics:",sep=" ",format(Sys.time(),"%a(%d)/%b/%Y   %X")))
      {
        # combine plots
        plot.with.inset = ggdraw() +
          draw_plot(graf2,width = .725, height = 1) +
          #                         x , y , width , height
          draw_plot(barplot.enrichment, .71,.01, .25, .87)
        ## Save plot
        png(file=paste("zxcvbnmNeVOmics_Plot_Network_Bar_Des_Layout_",i,".png"),
            width = 25,
            height = 15,
            units = "cm",
            res=900,bg="white")
        plot(plot.with.inset)
        dev.off()}
    }
    #
    for (i in layouts){
      graf2 = ggraph(link_tbl,layout = "igraph",algorithm = i) +
        geom_edge_link(colour = col.edges.cluster,alpha = 1,
                       width = correl.width(-0.01296*length(only.entrys)+1.798))+ #<---------------------
      geom_node_point(aes(colour= logFC),
                      size= correl.point(-0.04259*length(only.entrys)+5.980), #<---------------------
                      position = "identity")+
        geom_node_text(aes(label=name),
                       position ="identity",
                       size = correl.text(-0.02778*length(only.entrys)+4.639), #<---------------------
                       repel = T,
                       #vjust=3,
                       #hjust=.3,
                       fontface = "bold")+
        theme_classic()
      graf2 = graf2 +  theme(legend.key.width=unit(.4,"cm"),
                             legend.key.height=unit(.4,"cm"),
                             legend.position = c(1.06, .9),
                             legend.title=element_text(size=10,face="bold"),
                             legend.text=element_text(size=8),
                             axis.line=element_blank(),
                             axis.text = element_blank(),
                             axis.ticks = element_blank(),
                             axis.title = element_blank())
      graf2 = graf2 + scale_color_gradient(low = c("green","black"), high = "red") +
        labs(caption = paste("NeVOmics:",sep=" ",format(Sys.time(),"%a(%d)/%b/%Y   %X")))
      {
        # combine plots
        plot.with.inset = ggdraw() +
          draw_plot(graf2,width = .725, height = 1) +
          #                         x , y , width , height
          draw_plot(barplot.enrichment, .71,.01, .25, .87)
        ## Save plot
        png(file=paste("zxcvbnmNeVOmics_Plot_Network_Bar_Id_Layout_",i,".png"),
            width = 25,
            height = 15,
            units = "cm",
            res=900,bg="white")
        plot(plot.with.inset)
        dev.off()}
    }
    #
    for (i in layouts){
      graf2 = ggraph(link_tbl,layout = "igraph",algorithm = i) +
        geom_edge_link(colour = col.edges.cluster,alpha = 1,
                       width = correl.width(-0.01296*length(only.entrys)+1.798))+ #<---------------------
      geom_node_point(aes(colour= logFC),
                      size= correl.point(-0.04259*length(only.entrys)+5.980), #<---------------------
                      position = "identity")+
        geom_node_text(aes(label=Entry),
                       position ="identity",
                       size = correl.text(-0.02778*length(only.entrys)+4.639), #<---------------------
                       repel = T,
                       #vjust=3,
                       #hjust=.3,
                       fontface = "bold")+
        theme_classic()
      graf2 = graf2 +  theme(legend.key.width=unit(.4,"cm"),
                             legend.key.height=unit(.4,"cm"),
                             legend.position = c(1.06, .9),
                             legend.title=element_text(size=10,face="bold"),
                             legend.text=element_text(size=8),
                             axis.line=element_blank(),
                             axis.text = element_blank(),
                             axis.ticks = element_blank(),
                             axis.title = element_blank())
      graf2 = graf2 + scale_color_gradient(low = c("green","black"), high = "red") +
        labs(caption = paste("NeVOmics:",sep=" ",format(Sys.time(),"%a(%d)/%b/%Y   %X")))
      {
        # combine plots
        plot.with.inset = ggdraw() +
          draw_plot(graf2,width = .725, height = 1) +
          #                         x , y , width , height
          draw_plot(barplot.enrichment, .71,.01, .25, .87)
        ## Save plot
        png(file=paste("zxcvbnmNeVOmics_Plot_Network_Bar_Layout_",i,".png"),
            width = 25,
            height = 15,
            units = "cm",
            res=900,bg="white")
        plot(plot.with.inset)
        dev.off()}
    }
    #
    #
    #\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    if (length(only.entrys) > 80) {
      modelos = c('nicely','kk','dh', 'drl')
      V(link_tbl)$Etiqueta = c(as.character(nodes$Term[1:total.terms]),
                               as.character(rep(NA,total.nodes - total.terms)))
      for (i in modelos){
        graf2 = ggraph(link_tbl,layout = "igraph",algorithm = i) +
          geom_edge_link(colour = col.edges.cluster,alpha = 1,
                         width = correl.width(-0.01296*length(only.entrys)+1.798))+ #<---------------------
        geom_node_point(aes(colour= logFC),
                        size= 3, #correl.point(-0.04259*length(only.entrys)+5.980), #<---------------------
                        position = "identity")+
          geom_node_text(aes(label=Etiqueta),
                         position ="identity",
                         size = 3, #<---------------------
                         repel = T,
                         #vjust=3,
                         #hjust=.3,
                         fontface = "bold")+
          theme_classic()
        graf2 = graf2 +  theme(legend.key.width=unit(.4,"cm"),
                               legend.key.height=unit(.4,"cm"),
                               legend.position = c(1.06, .9),
                               legend.title=element_text(size=10,face="bold"),
                               legend.text=element_text(size=8),
                               axis.line=element_blank(),
                               axis.text = element_blank(),
                               axis.ticks = element_blank(),
                               axis.title = element_blank())
        graf2 = graf2 + scale_color_gradient(low = c("green","black"), high = "red") +
          labs(caption = paste("NeVOmics:",sep=" ",format(Sys.time(),"%a(%d)/%b/%Y   %X")))
        {
          # combine plots
          plot.with.inset = ggdraw() +
            draw_plot(graf2,width = .725, height = 1) +
            #                         x , y , width , height
            draw_plot(barplot.enrichment, .71,.01, .25, .87)
          ## Save plot
          png(file=paste("zxcvbnmNeVOmics_Plot_Network_Bar_Id_Layout_XL_",i,".png"),
              width = 25,
              height = 15,
              units = "cm",
              res=900,bg="white")
          plot(plot.with.inset)
          dev.off()}
      }
    }
    #\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    #
    #
    #
  } else if ((total.terms >= 31) && (total.terms <= 60)) {
    # ------------------------------------------------------------------------  two bars configuration
    #
    #
    # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    #
    # |||||||||      Bar plots
    #
    # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    #
    nodes.1_30 = nodes[1:30,]
    #
    {
      for.bar.1_30 = data.frame(c(data.frame("Entry"=c(nodes.1_30$Entry,rep(NA,31-nrow(nodes.1_30)))),
                                  data.frame("Freq"=c(nodes.1_30$Freq,rep(NA,31-nrow(nodes.1_30)))),
                                  data.frame("Fdr"=c(nodes.1_30$FDR,rep(NA,31-nrow(nodes.1_30)))),
                                  data.frame("Term"=c(nodes.1_30$Term,rep(NA,31-nrow(nodes.1_30)))),
                                  data.frame("P"=c(nodes.1_30$P,rep(NA,31-nrow(nodes.1_30))))))
      for.bar.1_30["log.pval"] = -log10(for.bar.1_30$P)
      for.bar.1_30["log.fdr"] = -log10(for.bar.1_30$Fdr)
    }
    #
    barplot.p.val.1_30 = ggplot(for.bar.1_30) +
      geom_bar(aes(x = Entry, weight = log.pval),fill=rev(many.colors[1:(nrow(nodes.1_30)+1)]),width=.7) +
      geom_line(aes(x = c(1:31), y = -log10(0.05),
                    color=-log10(0.05)),
                color = c(rep("transparent",31-nrow(nodes.1_30)),
                          rep("black",nrow(nodes.1_30))), size = .2 )+
      geom_point(aes(x = c(1:31), y = -log10(0.05),color=-log10(0.05)),
                 position = position_dodge(0.3),
                 color = c(rep("transparent",31-nrow(nodes.1_30)),
                           rep("black",nrow(nodes.1_30))), size = 1)+
      geom_line(aes(x = c(1:31), y = rev(log.fdr),
                    color=log.fdr),
                color = c(rep("transparent",31-nrow(nodes.1_30)),
                          rep("blue",nrow(nodes.1_30))), size = .2 )+
      geom_point(aes(x = c(1:31), y = rev(log.fdr),color=log.fdr),
                 position = position_dodge(0.3),shape=17,
                 color = c(rep("transparent",31-nrow(nodes.1_30)),
                           rep("blue",nrow(nodes.1_30))), size = 1)+
      #scale_y_continuous(position = "right")+
      coord_flip()+
      scale_x_discrete(limits=rev(for.bar.1_30$Entry))+
      geom_text(aes(x=Entry,y=log.pval+.35,label=Freq),size=2.5,fontface ="bold")+ #<---------------------
    #ggtitle("-log(P-value) and\n -log(Corrected P-value)")+
    theme(plot.title = element_blank(),
          axis.title.x = element_blank(),# element_text(colour="black",size=27,face="bold",hjust = 0,vjust = 10),
          axis.title.y = element_blank(),
          axis.text.y = element_text(colour=c("white",rep("black",30)),size=8,
                                     face="bold",hjust = 1.3,vjust = .5), #element_text(size = 8,color="black"),#tama?o y color de texto de eje x
          axis.text.x = element_text(colour="black",size=7,face="bold",vjust = 0),
          axis.ticks.length =unit(.1, "cm"),
          axis.ticks = element_line(color = "black",size = 1),
          axis.ticks.y.left = element_line(colour = "black",size=.5,unit(0, "cm")),
          #axis.ticks.x.top = element_line(colour = "black",size=.5,unit(3, "cm")),
          axis.line.y = element_blank(), axis.title.x.top =  element_blank(), 
          plot.margin = margin(1.5,1,0.5,0, "cm"))+
      scale_y_continuous(position = "right",
                         limit = c(0,round(max(for.bar.1_30$log.pval[1:nrow(nodes.1_30)]+1.5))),
                         breaks=seq(0,round(max(for.bar.1_30$log.pval[1:nrow(nodes.1_30)]+1.5),digits = 0),1))+
      labs(title = expression(bold("ASPECT")),
           subtitle = expression(bold("-log10("*italic("P")*"-value & adj."*italic("P")*"-value)")))+
      theme(plot.title = element_text(vjust = 2,hjust=-.3,size = 15),
            plot.subtitle = element_text(vjust = -.7,hjust=.6,size = 8))
    #
    # Second part of the data
    #
    nodes.31_final = nodes[31:total.terms,]
    #
    {
      for.bar.31_final = data.frame(c(data.frame("Entry"=c(nodes.31_final$Entry,rep(NA,31-nrow(nodes.31_final)))),
                                      data.frame("Freq"=c(nodes.31_final$Freq,rep(NA,31-nrow(nodes.31_final)))),
                                      data.frame("Fdr"=c(nodes.31_final$FDR,rep(NA,31-nrow(nodes.31_final)))),
                                      data.frame("Term"=c(nodes.31_final$Term,rep(NA,31-nrow(nodes.31_final)))),
                                      data.frame("P"=c(nodes.31_final$P,rep(NA,31-nrow(nodes.31_final))))))
      for.bar.31_final["log.pval"] = -log10(for.bar.31_final$P)
      for.bar.31_final["log.fdr"] = -log10(for.bar.31_final$Fdr)
    }
    
    barplot.p.val.31_final = ggplot(for.bar.31_final) +
      geom_bar(aes(x = Entry, weight = log.pval),fill=rev(many.colors[(nrow(nodes.1_30)+1):(total.terms+1)]),width=.7) +
      geom_line(aes(x = c(1:31), y = -log10(0.05),
                    color=-log10(0.05)),
                color = c(rep("transparent",31-nrow(nodes.31_final)),
                          rep("black",nrow(nodes.31_final))), size = .2 )+
      geom_point(aes(x = c(1:31), y = -log10(0.05),color=-log10(0.05)),
                 position = position_dodge(0.3),
                 color = c(rep("transparent",31-nrow(nodes.31_final)),
                           rep("black",nrow(nodes.31_final))), size = 1)+
      geom_line(aes(x = c(1:31), y = rev(log.fdr),
                    color=log.fdr),
                color = c(rep("transparent",31-nrow(nodes.31_final)),
                          rep("blue",nrow(nodes.31_final))), size = .2 )+
      geom_point(aes(x = c(1:31), y = rev(log.fdr),color=log.fdr),
                 position = position_dodge(0.3),shape=17,
                 color = c(rep("transparent",31-nrow(nodes.31_final)),
                           rep("blue",nrow(nodes.31_final))), size = 1)+
      #scale_y_continuous(position = "right")+
      coord_flip()+
      scale_x_discrete(limits=rev(for.bar.31_final$Entry))+
      geom_text(aes(x=Entry,y=log.pval+.35,label=Freq),size=2.5,fontface ="bold")+ #<---------------------
    #ggtitle("-log(P-value) and\n -log(Corrected P-value)")+
    theme(plot.title = element_text(colour="white"),
          axis.title.x =element_blank(),# element_text(colour="black",size=27,face="bold",hjust = 0,vjust = 10),
          axis.title.y = element_blank(),
          axis.text.y = element_text(colour=c("white",rep("black",30)),size=8,
                                     face="bold",hjust = 1.3,vjust = .5), #element_text(size = 8,color="black"),#tama?o y color de texto de eje x
          axis.text.x = element_text(colour="black",size=7,face="bold",vjust = 0),
          axis.ticks.length =unit(.1, "cm"),
          axis.ticks = element_line(color = "black",size = 1),
          axis.ticks.y.left = element_line(colour = "black",size=.5,unit(0, "cm")),
          #axis.ticks.x.top = element_line(colour = "black",size=3,unit(3, "cm")),
          axis.line.y = element_blank(), axis.title.x.top =  element_blank(), 
          plot.margin = margin(1.5,1,0.5,0, "cm"))+
      scale_y_continuous(position = "right",
                         limit = c(0,round(max(for.bar.1_30$log.pval[1:nrow(nodes.1_30)]+1.5))),
                         breaks=seq(0,round(max(for.bar.1_30$log.pval[1:nrow(nodes.1_30)]+1.5),digits = 0),1))+
      labs(title = expression(bold("xxxxxxxx")),
           subtitle = expression(bold("-log10("*italic("P")*"-value & adj."*italic("P")*"-value)")))+
      theme(plot.title = element_text(vjust = 0.5,hjust=1.5,size = 12),
            plot.subtitle = element_text(vjust = -.7,hjust=.6,size = 8))
    #
    #
    #
    # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    #
    # |||||||||      Complex Plots
    #
    # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    #
    #
    graf2 = ggraph(link_tbl,layout = "igraph",algorithm = "kk") +
      geom_edge_link(colour = col.edges.cluster,alpha = 1,
                     width = correl.width(-0.01296*length(only.entrys)+1.798))+ #<---------------------
    geom_node_point(aes(colour= logFC),
                    size= correl.point(-0.04259*length(only.entrys)+5.980), #<---------------------
                    position = "identity")+
      geom_node_text(aes(label=Label),
                     position ="identity",
                     size = correl.text(-0.02778*length(only.entrys)+4.639), #<---------------------
                     repel = T,
                     #vjust=3,
                     #hjust=.3,
                     fontface = "bold")+
      theme_classic()
    graf2 = graf2 +  theme(legend.key.width=unit(.4,"cm"),
                           legend.key.height=unit(.4,"cm"),
                           legend.position = c(1.06, .9),
                           legend.title=element_text(size=10,face="bold"),
                           legend.text=element_text(size=8),
                           axis.line=element_blank(),
                           axis.text = element_blank(),
                           axis.ticks = element_blank(),
                           axis.title = element_blank())
    graf2 = graf2 + scale_color_gradient(low = c("green","black"), high = "red") +
      labs(caption = paste("NeVOmics:",sep=" ",format(Sys.time(),"%a(%d)/%b/%Y   %X")))
    #
    ###
    for (i in layouts){
      graf2 = ggraph(link_tbl,layout = "igraph",algorithm = i) +
        geom_edge_link(colour = col.edges.cluster,alpha = 1,
                       width = correl.width(-0.01296*length(only.entrys)+1.798))+ #<---------------------
      geom_node_point(aes(colour= logFC),
                      size= correl.point(-0.04259*length(only.entrys)+5.980), #<---------------------
                      position = "identity")+
        geom_node_text(aes(label=Label),
                       position ="identity",
                       size = correl.text(-0.02778*length(only.entrys)+4.639), #<---------------------
                       repel = T,
                       #vjust=3,
                       #hjust=.3,
                       fontface = "bold")+
        theme_classic()
      graf2 = graf2 +  theme(legend.key.width=unit(.4,"cm"),
                             legend.key.height=unit(.4,"cm"),
                             legend.position = c(1.06, .9),
                             legend.title=element_text(size=10,face="bold"),
                             legend.text=element_text(size=8),
                             axis.line=element_blank(),
                             axis.text = element_blank(),
                             axis.ticks = element_blank(),
                             axis.title = element_blank())
      graf2 = graf2 + scale_color_gradient(low = c("green","black"), high = "red") +
        labs(caption = paste("NeVOmics:",sep=" ",format(Sys.time(),"%a(%d)/%b/%Y   %X")))
      {
        plots_combined  = ggdraw() +
          draw_plot(barplot.p.val.1_30, x = 0, y = 0, width = .5, height = 1) +
          draw_plot(barplot.p.val.31_final, x = .45, y = 0, width = .5, height = 1)
        # combine plots
        plot.with.inset = ggdraw() +
          draw_plot(graf2,  width = .54, height = 1) +
          #                         x , y , width , height
          draw_plot(plots_combined, .53,.01, .45, .85) # 2 barplots
        
        ## Save plot
        png(file=paste("zxcvbnmNeVOmics_Plot_Network_Bar_Des_Layout_",i,".png"),
            width = 35,
            height = 15,
            units = "cm",
            res=900,bg="white")
        plot(plot.with.inset)
        dev.off()}
    }
    #
    for (i in layouts){
      graf2 = ggraph(link_tbl,layout = "igraph",algorithm = i) +
        geom_edge_link(colour = col.edges.cluster,alpha = 1,
                       width = correl.width(-0.01296*length(only.entrys)+1.798))+ #<---------------------
      geom_node_point(aes(colour= logFC),
                      size= correl.point(-0.04259*length(only.entrys)+5.980), #<---------------------
                      position = "identity")+
        geom_node_text(aes(label=name),
                       position ="identity",
                       size = correl.text(-0.02778*length(only.entrys)+4.639), #<---------------------
                       repel = T,
                       #vjust=3,
                       #hjust=.3,
                       fontface = "bold")+
        theme_classic()
      graf2 = graf2 +  theme(legend.key.width=unit(.4,"cm"),
                             legend.key.height=unit(.4,"cm"),
                             legend.position = c(1.06, .9),
                             legend.title=element_text(size=10,face="bold"),
                             legend.text=element_text(size=8),
                             axis.line=element_blank(),
                             axis.text = element_blank(),
                             axis.ticks = element_blank(),
                             axis.title = element_blank())
      graf2 = graf2 + scale_color_gradient(low = c("green","black"), high = "red") +
        labs(caption = paste("NeVOmics:",sep=" ",format(Sys.time(),"%a(%d)/%b/%Y   %X")))
      {
        plots_combined  = ggdraw() +
          draw_plot(barplot.p.val.1_30, x = 0, y = 0, width = .5, height = 1) +
          draw_plot(barplot.p.val.31_final, x = .45, y = 0, width = .5, height = 1)
        # combine plots
        plot.with.inset = ggdraw() +
          draw_plot(graf2,  width = .54, height = 1) +
          #                         x , y , width , height
          draw_plot(plots_combined, .53,.01, .45, .85) # 2 barplots
        
        ## Save plot
        png(file=paste("zxcvbnmNeVOmics_Plot_Network_Bar_Id_Layout_",i,".png"),
            width = 35,
            height = 15,
            units = "cm",
            res=900,bg="white")
        plot(plot.with.inset)
        dev.off()}
    }
    #
    for (i in layouts){
      graf2 = ggraph(link_tbl,layout = "igraph",algorithm = i) +
        geom_edge_link(colour = col.edges.cluster,alpha = 1,
                       width = correl.width(-0.01296*length(only.entrys)+1.798))+ #<---------------------
      geom_node_point(aes(colour= logFC),
                      size= correl.point(-0.04259*length(only.entrys)+5.980), #<---------------------
                      position = "identity")+
        geom_node_text(aes(label=Entry),
                       position ="identity",
                       size = correl.text(-0.02778*length(only.entrys)+4.639), #<---------------------
                       repel = T,
                       #vjust=3,
                       #hjust=.3,
                       fontface = "bold")+
        theme_classic()
      graf2 = graf2 +  theme(legend.key.width=unit(.4,"cm"),
                             legend.key.height=unit(.4,"cm"),
                             legend.position = c(1.06, .9),
                             legend.title=element_text(size=10,face="bold"),
                             legend.text=element_text(size=8),
                             axis.line=element_blank(),
                             axis.text = element_blank(),
                             axis.ticks = element_blank(),
                             axis.title = element_blank())
      graf2 = graf2 + scale_color_gradient(low = c("green","black"), high = "red") +
        labs(caption = paste("NeVOmics:",sep=" ",format(Sys.time(),"%a(%d)/%b/%Y   %X")))
      {
        plots_combined  = ggdraw() +
          draw_plot(barplot.p.val.1_30, x = 0, y = 0, width = .5, height = 1) +
          draw_plot(barplot.p.val.31_final, x = .45, y = 0, width = .5, height = 1)
        # combine plots
        plot.with.inset = ggdraw() +
          draw_plot(graf2,  width = .54, height = 1) +
          #                         x , y , width , height
          draw_plot(plots_combined, .53,.01, .45, .85) # 2 barplots
        
        ## Save plot
        png(file=paste("zxcvbnmNeVOmics_Plot_Network_Bar_Layout_",i,".png"),
            width = 35,
            height = 15,
            units = "cm",
            res=900,bg="white")
        plot(plot.with.inset)
        dev.off()}
    }
    #
    #
    #\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    if (length(only.entrys) > 80) {
      modelos = c('nicely','kk','dh', 'drl')
      V(link_tbl)$Etiqueta = c(as.character(nodes$Term[1:total.terms]),
                               as.character(rep(NA,total.nodes - total.terms)))
      for (i in modelos){
        graf2 = ggraph(link_tbl,layout = "igraph",algorithm = i) +
          geom_edge_link(colour = col.edges.cluster,alpha = 1,
                         width = correl.width(-0.01296*length(only.entrys)+1.798))+ #<---------------------
        geom_node_point(aes(colour= logFC),
                        size= 3, #correl.point(-0.04259*length(only.entrys)+5.980), #<---------------------
                        position = "identity")+
          geom_node_text(aes(label=Etiqueta),
                         position ="identity",
                         size = 3, #<---------------------
                         repel = T,
                         #vjust=3,
                         #hjust=.3,
                         fontface = "bold")+
          theme_classic()
        graf2 = graf2 +  theme(legend.key.width=unit(.4,"cm"),
                               legend.key.height=unit(.4,"cm"),
                               legend.position = c(1.06, .9),
                               legend.title=element_text(size=10,face="bold"),
                               legend.text=element_text(size=8),
                               axis.line=element_blank(),
                               axis.text = element_blank(),
                               axis.ticks = element_blank(),
                               axis.title = element_blank())
        graf2 = graf2 + scale_color_gradient(low = c("green","black"), high = "red") +
          labs(caption = paste("NeVOmics:",sep=" ",format(Sys.time(),"%a(%d)/%b/%Y   %X")))
        {
          plots_combined  = ggdraw() +
            draw_plot(barplot.p.val.1_30, x = 0, y = 0, width = .5, height = 1) +
            draw_plot(barplot.p.val.31_final, x = .45, y = 0, width = .5, height = 1)
          # combine plots
          plot.with.inset = ggdraw() +
            draw_plot(graf2,  width = .54, height = 1) +
            #                         x , y , width , height
            draw_plot(plots_combined, .53,.01, .45, .85) # 2 barplots
          
          ## Save plot
          png(file=paste("zxcvbnmNeVOmics_Plot_Network_Bar_Des_Layout_XL_",i,".png"),
              width = 35,
              height = 15,
              units = "cm",
              res=900,bg="white")
          plot(plot.with.inset)
          dev.off()}
      }
    }
    #\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    #
    #
    #
    #>
  }
  #
  #___________________________________________________________________________________________ Chrod plots
  #
  #
  # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  #
  # |||||||||     Chord plots (with expression values)
  #
  # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  #
  if (length(only.entrys) <= 70) {
    if ((total.terms > 1) && (total.terms <= 15)) {
      print("grafico 1 y 2")
      for (i in -1:0) {
        #
        #chord 1
        #
        go.term = drop_na(nodes[,c("Entry","Term")])
        colnames(go.term)[1] <- "GO"
        #
        go.entry.term = merge(links,go.term,by="GO")#[c(3,2)]
        #
        go.paste.term=paste0(go.entry.term$GO," ~ ",go.entry.term$Term)
        #
        term.entry = data.frame("Term"=go.paste.term,"Entry"=go.entry.term$Entry)
        term.entry=as.data.frame.matrix(table(term.entry)) 
        term.entry.mat = as.matrix(term.entry) 
        #########
        row_sum = sum(rowSums(abs(term.entry.mat)))
        col_sum = sum(colSums(abs(term.entry.mat)))
        small_gap = 1
        big_gap = 15
        nr = nrow(term.entry.mat)
        nc = ncol(term.entry.mat)
        n_sector = nr + nc
        row_sector_degree = (360 - small_gap*(n_sector - 2) - big_gap*2) * (row_sum/(row_sum + col_sum)) +
          small_gap*(nr-1)
        ###
        start_degree = 90 - (180 - row_sector_degree)/2
        #
        circos.clear()
        # gap entre cada cuerda
        circos.par(gap.after = c(rep(1, nrow(term.entry.mat)-1), 15, rep(1, ncol(term.entry.mat)-1),15),
                   #circos.par(gap.after = c(rep(0.5, nrow(a)-1), 20, rep(0.5, ncol(a)-1),20),
                   start.degree = start_degree)
        #circos.par(start.degree = 80, clock.wise = T)
        circos.par(track.margin = c(0,0))
        chord.plot.1 = function () {
          chordDiagram(term.entry.mat, order = c(paste0(go.term$GO," ~ ",go.term$Term),only.entrys),
                       grid.col = colors.fold.change,
                       transparency = 0,
                       directional = i,
                       annotationTrack = "grid",
                       #    c(% del datio del circulo, alto del track)
                       annotationTrackHeight = c(0.08, .01),
                       preAllocateTracks = 2)
          circos.track(track.index = 2,ylim = c(0, 1), track.height = uh(5, "cm"),
                       track.margin = c(5, uh(5, "cm")),
                       panel.fun = function(x, y) {
                         xcenter = get.cell.meta.data("xcenter")
                         circos.lines(c(xcenter, xcenter), c(0, uy(.3, "cm")),
                                      col = 1:4)
                       },bg.border =  NA)
          circos.track(track.index = 2, panel.fun = for(si in get.all.sector.index()[1:total.terms]) {
            xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 2)
            ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 2)
            circos.text(mean(xlim), mean(ylim), si, sector.index = si, track.index = 2,
                        facing = "downward",
                        niceFacing = TRUE,
                        adj = c(0,0.5),
                        cex = .5, 
                        font = 2,
                        col = NA)
            
          }, bg.border =  NA)
          circos.track(track.index = 2, panel.fun = for(si in get.all.sector.index()[start.entry:total.nodes]) {
            xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 2)
            ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 2)
            circos.text(mean(xlim), mean(ylim), si, sector.index = si, track.index = 2,
                        facing = "clockwise",
                        niceFacing = TRUE,
                        adj = c(0,0.5),
                        cex = -0.003039*length(only.entrys) + 0.4913, #<---------------------
                        font = 2,
                        col = "black")
          }, bg.border = NA)
        }
        #
        #
        {
          png(file=paste0("zxcvbnmNeVOmics_Plot_Chord_1_",i,".png"),
              width = 30,
              height = 10,
              units = "cm",
              res=1024,bg="white")
          chord.plot.1()
          rasterImage(legend_image, -1, .8, -1.08,.53)
          text(x=-1.03, y = seq(0.8,.53,l=5),
               labels = round(seq(max.value.exp,min.value.exp,length.out = 5),digits=0),font = 1,
               pos = 4, cex=.5) #<---------------------
          text(x=-1,y= .83, "logFC", adj = c(0,0) ,font = 2, pos = 3,
               cex = .6) #<---------------------
          text(x=.4, y=.9, "ASPECT", adj = c(0,0), 
               cex = .9,font = 2) #<---------------------
          dev.off()}
        #
        #chord 2
        #
        go.term = drop_na(nodes[,c("Entry","Term")])
        colnames(go.term)[1] <- "GO"
        
        go.entry.term = merge(links,go.term,by="GO")#[c(3,2)]
        
        ##
        
        go.entry = data.frame("GO"=go.entry.term$GO,"Entry"=go.entry.term$Entry)
        go.entry=as.data.frame.matrix(table(go.entry)) 
        go.entry.mat = as.matrix(go.entry) 
        
        #########
        circos.clear()
        row_sum = sum(rowSums(abs(go.entry.mat)))
        col_sum = sum(colSums(abs(go.entry.mat)))
        small_gap = 1
        big_gap = 5
        nr = nrow(go.entry.mat)
        nc = ncol(go.entry.mat)
        n_sector = nr + nc
        row_sector_degree = (360 - small_gap*(n_sector - 2) - big_gap*2) * (row_sum/(row_sum + col_sum)) +
          small_gap*(nr-1)
        ###
        start_degree = 90 - (180 - row_sector_degree)/2
        #
        # gap entre cada cuerda
        circos.par(gap.after = c(rep(.5, nrow(go.entry.mat)-1), 5, rep(.5, ncol(go.entry.mat)-1),5),
                   #circos.par(gap.after = c(rep(0.5, nrow(a)-1), 20, rep(0.5, ncol(a)-1),20),
                   start.degree = start_degree)
        #circos.par(start.degree = 80, clock.wise = T)
        circos.par(track.margin = c(0,0))
        chord.plot.2 = function () {
          chordDiagram(go.entry.mat, order = nodes$Entry, # c(nodes$Term[1:total.terms],only.entrys)
                       grid.col = colors.fold.change,
                       transparency = 0,
                       directional = i,
                       annotationTrack = "grid",
                       #    c(% del datio del circulo, alto del track)
                       annotationTrackHeight = c(0.08, .01),
                       preAllocateTracks = 2)
          circos.track(track.index = 2,ylim = c(0, 1), track.height = uh(5, "cm"),
                       track.margin = c(5, uh(5, "cm")),
                       panel.fun = function(x, y) {
                         xcenter = get.cell.meta.data("xcenter")
                         circos.lines(c(xcenter, xcenter), c(0, uy(.2, "cm")),
                                      col = 1:4)
                       },bg.border =  NA)
          circos.track(track.index = 2, panel.fun = for(si in get.all.sector.index()[1:total.terms]) {
            xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 2)
            ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 2)
            circos.text(mean(xlim), mean(ylim), si, sector.index = si, track.index = 2,
                        facing = "clockwise",
                        niceFacing = TRUE,
                        adj = c(0.15,0.5),
                        cex = .5, #<---------------------
                        font = 2,
                        col = "black")
            
          }, bg.border =  NA)
          circos.track(track.index = 2, panel.fun = for(si in get.all.sector.index()[start.entry:total.nodes]) {
            xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 2)
            ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 2)
            circos.text(mean(xlim), mean(ylim), si, sector.index = si, track.index = 2,
                        facing = "clockwise",
                        niceFacing = TRUE,
                        adj = c(0.15,0.5),
                        cex = -0.003039*length(only.entrys) + 0.4913, #<---------------------
                        font = 2,
                        col = "black")
          }, bg.border = NA)
        }
        #
        {
          png(file=paste0("zxcvbnmNeVOmics_Plot_Chord_2_",i,".png"),
              width = 30,
              height = 10,
              units = "cm",
              res=900,bg="white")
          lgd_points = Legend(at = nodes$Term[1:total.terms], type = "grid",
                              size = unit(0, "mm"),
                              labels_gp = gpar(col = "white"),
                              legend_gp = gpar(col = "white"), 
                              title_position = "topleft",
                              title = "")
          lgd_list_vertical = packLegend(lgd_points)
          plot.new()
          circle_size = unit(1, "snpc") # snpc unit gives you a square region
          pushViewport(viewport(x = 0, y = 0.5, width = circle_size, height = circle_size,
                                just = c("left", "center")))
          par(omi = gridOMI(), new = TRUE)
          chord.plot.2()
          upViewport()
          pushViewport(viewport(x = circle_size, y = 0.5, width = grobWidth(lgd_list_vertical),
                                height = grobHeight(lgd_list_vertical), just = c("left", "center")))
          grid.draw(lgd_list_vertical)
          upViewport()
          rasterImage(legend_image, -1, .9, -1.08,.63)
          text(x=-1.03, y = seq(0.9,.63,l=5),
               labels = round(seq(max.value.exp,min.value.exp,length.out = 5),digits=0),font = 1,
               pos = 4, cex=.5) #<---------------------
          text(x=-1,y= .93, "logFC", adj = c(0,0) ,font = 2, pos = 3,
               cex = .6) #<---------------------
          text(x=1, y=1, "ASPECT", adj = c(0,0), 
               cex = .9,font = 2) #<---------------------
          legend(x=1.1,y=.95, 
                 legend=unique(paste0(nodes$Term[1:total.terms]," (",nodes$Freq[1:total.terms],")")),
                 cex=.6,bg="transparent", #<---------------------
                 border = NA,pt.bg=many.colors,bty = "n",
                 x.intersp=.7,y.intersp=.9, #<---------------------
                 col=NA,text.font = 1,pch=22,
                 pt.cex=1)
          dev.off()}
      }
    } else if (total.terms >= 16) {
      print("grafico 2")
      for (i in -1:0) {
        #
        #
        #chord 2
        #
        go.term = drop_na(nodes[,c("Entry","Term")])
        colnames(go.term)[1] <- "GO"
        
        go.entry.term = merge(links,go.term,by="GO")#[c(3,2)]
        
        ##
        
        go.entry = data.frame("GO"=go.entry.term$GO,"Entry"=go.entry.term$Entry)
        go.entry=as.data.frame.matrix(table(go.entry)) 
        go.entry.mat = as.matrix(go.entry) 
        
        #########
        circos.clear()
        row_sum = sum(rowSums(abs(go.entry.mat)))
        col_sum = sum(colSums(abs(go.entry.mat)))
        small_gap = 1
        big_gap = 5
        nr = nrow(go.entry.mat)
        nc = ncol(go.entry.mat)
        n_sector = nr + nc
        row_sector_degree = (360 - small_gap*(n_sector - 2) - big_gap*2) * (row_sum/(row_sum + col_sum)) +
          small_gap*(nr-1)
        ###
        start_degree = 90 - (180 - row_sector_degree)/2
        #
        # gap entre cada cuerda
        circos.par(gap.after = c(rep(.5, nrow(go.entry.mat)-1), 5, rep(.5, ncol(go.entry.mat)-1),5),
                   #circos.par(gap.after = c(rep(0.5, nrow(a)-1), 20, rep(0.5, ncol(a)-1),20),
                   start.degree = start_degree)
        #circos.par(start.degree = 80, clock.wise = T)
        circos.par(track.margin = c(0,0))
        chord.plot.2 = function () {
          chordDiagram(go.entry.mat, order = nodes$Entry, # c(nodes$Term[1:total.terms],only.entrys)
                       grid.col = colors.fold.change,
                       transparency = 0,
                       directional = i,
                       annotationTrack = "grid",
                       #    c(% del datio del circulo, alto del track)
                       annotationTrackHeight = c(0.08, .01),
                       preAllocateTracks = 2)
          circos.track(track.index = 2,ylim = c(0, 1), track.height = uh(5, "cm"),
                       track.margin = c(5, uh(5, "cm")),
                       panel.fun = function(x, y) {
                         xcenter = get.cell.meta.data("xcenter")
                         circos.lines(c(xcenter, xcenter), c(0, uy(.2, "cm")),
                                      col = 1:4)
                       },bg.border =  NA)
          circos.track(track.index = 2, panel.fun = for(si in get.all.sector.index()[1:total.terms]) {
            xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 2)
            ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 2)
            circos.text(mean(xlim), mean(ylim), si, sector.index = si, track.index = 2,
                        facing = "clockwise",
                        niceFacing = TRUE,
                        adj = c(0.15,0.5),
                        cex = .5, #<---------------------
                        font = 2,
                        col = "black")
            
          }, bg.border =  NA)
          circos.track(track.index = 2, panel.fun = for(si in get.all.sector.index()[start.entry:total.nodes]) {
            xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 2)
            ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 2)
            circos.text(mean(xlim), mean(ylim), si, sector.index = si, track.index = 2,
                        facing = "clockwise",
                        niceFacing = TRUE,
                        adj = c(0.15,0.5),
                        cex = -0.003039*length(only.entrys) + 0.4913, #<---------------------
                        font = 2,
                        col = "black")
          }, bg.border = NA)
        }
        #
        {
          png(file=paste0("zxcvbnmNeVOmics_Plot_Chord_2_",i,".png"),
              width = 30,
              height = 10,
              units = "cm",
              res=900,bg="white")
          lgd_points = Legend(at = nodes$Term[1:total.terms], type = "grid",
                              size = unit(0, "mm"),
                              labels_gp = gpar(col = "white"),
                              legend_gp = gpar(col = "white"), 
                              title_position = "topleft",
                              title = "")
          lgd_list_vertical = packLegend(lgd_points)
          plot.new()
          circle_size = unit(1, "snpc") # snpc unit gives you a square region
          pushViewport(viewport(x = 0, y = 0.5, width = circle_size, height = circle_size,
                                just = c("left", "center")))
          par(omi = gridOMI(), new = TRUE)
          chord.plot.2()
          upViewport()
          pushViewport(viewport(x = circle_size, y = 0.5, width = grobWidth(lgd_list_vertical),
                                height = grobHeight(lgd_list_vertical), just = c("left", "center")))
          grid.draw(lgd_list_vertical)
          upViewport()
          rasterImage(legend_image, -1, .9, -1.08,.63)
          text(x=-1.03, y = seq(0.9,.63,l=5),
               labels = round(seq(max.value.exp,min.value.exp,length.out = 5),digits=0),font = 1,
               pos = 4, cex=.5) #<---------------------
          text(x=-1,y= .93, "logFC", adj = c(0,0) ,font = 2, pos = 3,
               cex = .6) #<---------------------
          text(x=1, y=1, "ASPECT", adj = c(0,0), 
               cex = .9,font = 2) #<---------------------
          legend(x=1.1,y=.95, 
                 legend=unique(paste0(nodes$Term[1:total.terms]," (",nodes$Freq[1:total.terms],")")),
                 cex=.5,bg="transparent", #<---------------------
                 border = NA,pt.bg=many.colors,bty = "n",
                 x.intersp=.7,y.intersp=.9, #<---------------------
                 col=NA,text.font = 1,pch=22,
                 pt.cex=1)
          dev.off()}
      }
    }
  } else {
    print('no chord plot')
  }
  #
  #
  #___________________________________________________________________________________________ UpSet plot
  #
  #
  # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  #
  # |||||||||     Upset Plot 
  #
  # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  #
  #
  ups=as.data.frame.matrix(t(table(links)))
  #
  # size of = set names
  if ((total.terms >= 1) && (total.terms <= 7)) {
    set.name.size=2
  } else if ((total.terms >= 8) && (total.terms <= 15)) {
    set.name.size=1.5
  } else if ((total.terms >= 16) && (total.terms <= 25)) {
    set.name.size=1.25
  } else if (total.terms >= 25) {
    set.name.size=0.9
  }
  #
  # size of = point.size and line.size (point.size)
  if ((total.terms >= 1) && (total.terms <= 10)) {
    point.size=5
  } else if ((total.terms >= 11) && (total.terms <= 20)) {
    point.size=4.16
  } else if ((total.terms >= 21) && (total.terms <= 30)) {
    point.size=3.3
  } else if (total.terms >= 31) {
    point.size=2.5
  }
  #
  #
  {png(file="zxcvbnmNeVOmics_Plot_UpSetR.png",width = 10,height = 7,units = "in",res=900,bg="white")
    upset(ups,sets=nodes$Entry[1:total.terms],
          sets.bar.color = "blue",order.by ="freq",empty.intersections = NULL,
          point.size=point.size,mainbar.y.label="poiuytrewq",sets.x.label = "",
          main.bar.color="black",matrix.color="black",shade.color="red",
          line.size=point.size/5,show.numbers = "yes",group.by = "degree",
          matrix.dot.alpha = 0,mb.ratio = c(0.5, 0.5),
          #c(intersection size title, intersection size tick labels, set size title, set size tick labels, set names, numbers above bars)
          #            1                           2                      3                   4                5                6           
          text.scale = c(2,2,1,1.5,set.name.size,2), 
          scale.sets = "identity",scale.intersections = "identity",
          shade.alpha = 0.35,set_size.angles = 0)
    dev.off()}
  #>
  print(warnings())
  #
  #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  #++++++++++++            there are no expression values             ++++++++++++++
  #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
} else {
  print("there are no expression values")
  #
  prot.col=rep("blue",1000)#9D9D9D
  col.set.node=c(rep("gray40",total.terms),
                 prot.col[nodes$num[start.entry:total.nodes]])
  col.set.link=c(many.colors[nodes$num[1:total.terms]],
                 prot.col[nodes$num[start.entry:total.nodes]])
  #
  layouts=c("nicely","kk","dh")
  #
  if ((total.terms > 1) && (total.terms <= 20)) {
    width.png=20
    print(paste("el ancho de la figura es:",width.png))
  } else if ((total.terms >= 21) && (total.terms <= 40)) {
    width.png=25
    print(paste("el ancho de la figura es:",width.png))
  } else if (total.terms > 41) {
    width.png=30
  }
  #
  #
  for (i in layouts){
    graf1 = ggraph(link_tbl,layout = "igraph", algorithm = i) +
      geom_edge_link(aes(colour = GO),alpha = 1,width = correl.width(-0.01296*length(only.entrys)+1.798))+ #<---------------------
    geom_node_point(colour= col.set.node,shape=c(rep(16,total.terms),rep(20,length(only.entrys))),
                    size= correl.point(-0.04259*length(only.entrys)+5.980),show.legend = NA, #<---------------------
                    position = "identity")+
      scale_edge_colour_manual(values=col.set.link,name = "ASPECT")+
      geom_node_text(aes(label=Label),
                     position ="identity",
                     size = correl.text(-0.02778*length(only.entrys)+4.639), #<---------------------
                     repel = T,
                     #vjust=3,
                     #hjust=.3,
                     fontface = "bold")+
      theme_classic() +
      theme(legend.key.width=unit(.4,"cm"),
            legend.key.height=unit(.4,"cm"),
            legend.title=element_text(size=10,face="bold"),
            legend.text=element_text(size=8),
            axis.line=element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.title = element_blank())
    graf1 = graf1 + scale_color_gradient(low = c("green","black"), high = "red") +
      labs(caption = paste("NeVOmics:",sep=" ",format(Sys.time(),"%a(%d)/%b/%Y   %X")))
    {
      png(file=paste0("zxcvbnmNeVOmics_Plot_Network_Des_Layout_",i,".png"),
          width = width.png,  #<---------------------
          height = 15, 
          units = "cm",
          res=900,bg="white")
      plot(graf1)
      dev.off()}
  }
  #
  for (i in layouts){
    graf1 = ggraph(link_tbl,layout = "igraph", algorithm = i) +
      geom_edge_link(aes(colour = GO),alpha = 1,width = correl.width(-0.01296*length(only.entrys)+1.798))+ #<---------------------
    geom_node_point(colour= col.set.node,shape=c(rep(16,total.terms),rep(20,length(only.entrys))),
                    size= correl.point(-0.04259*length(only.entrys)+5.980),show.legend = NA, #<---------------------
                    position = "identity")+
      scale_edge_colour_manual(values=col.set.link,name = "ASPECT")+
      geom_node_text(aes(label=name),
                     position ="identity",
                     size = correl.text(-0.02778*length(only.entrys)+4.639), #<---------------------
                     repel = T,
                     #vjust=3,
                     #hjust=.3,
                     fontface = "bold")+
      theme_classic() +
      theme(legend.key.width=unit(.4,"cm"),
            legend.key.height=unit(.4,"cm"),
            legend.title=element_text(size=10,face="bold"),
            legend.text=element_text(size=8),
            axis.line=element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.title = element_blank())
    graf1 = graf1 + scale_color_gradient(low = c("green","black"), high = "red") +
      labs(caption = paste("NeVOmics:",sep=" ",format(Sys.time(),"%a(%d)/%b/%Y   %X")))
    {
      png(file=paste0("zxcvbnmNeVOmics_Plot_Network_Id_Layout_",i,".png"),
          width = width.png,  #<---------------------
          height = 15, 
          units = "cm",
          res=900,bg="white")
      plot(graf1)
      dev.off()}
  }
  #
  for (i in layouts){
    graf1 = ggraph(link_tbl,layout = "igraph", algorithm = i) +
      geom_edge_link(aes(colour = GO),alpha = 1,width = correl.width(-0.01296*length(only.entrys)+1.798))+ #<---------------------
    geom_node_point(colour= col.set.node,shape=c(rep(16,total.terms),rep(20,length(only.entrys))),
                    size= correl.point(-0.04259*length(only.entrys)+5.980),show.legend = NA, #<---------------------
                    position = "identity")+
      scale_edge_colour_manual(values=col.set.link,name = "ASPECT")+
      geom_node_text(aes(label=Entry),
                     position ="identity",
                     size = correl.text(-0.02778*length(only.entrys)+4.639), #<---------------------
                     repel = T,
                     #vjust=3,
                     #hjust=.3,
                     fontface = "bold")+
      theme_classic() +
      theme(legend.key.width=unit(.4,"cm"),
            legend.key.height=unit(.4,"cm"),
            legend.title=element_text(size=10,face="bold"),
            legend.text=element_text(size=8),
            axis.line=element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.title = element_blank())
    graf1 = graf1 + scale_color_gradient(low = c("green","black"), high = "red") +
      labs(caption = paste("NeVOmics:",sep=" ",format(Sys.time(),"%a(%d)/%b/%Y   %X")))
    {
      png(file=paste0("zxcvbnmNeVOmics_Plot_Network_Layout_",i,".png"),
          width = width.png,  #<---------------------
          height = 15, 
          units = "cm",
          res=900,bg="white")
      plot(graf1)
      dev.off()}
  }
  #
  #
  #\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
  # loop si son mas de 80 proteinas
  if (length(only.entrys) > 80) {
    ##
    modelos = c('nicely','kk','dh', 'drl')
    V(link_tbl)$Etiqueta = c(as.character(nodes$Term[1:total.terms]),
                             as.character(rep(NA,total.nodes - total.terms)))
    ##
    for (i in modelos){
      graf1 = ggraph(link_tbl,layout = "igraph", algorithm = i) +
        geom_edge_link(aes(colour = GO),alpha = 1,width = correl.width(-0.01296*length(only.entrys)+1.798))+ #<---------------------
      geom_node_point(colour= col.set.node,shape=c(rep(16,total.terms),rep(20,length(only.entrys))),
                      size= correl.point(-0.04259*length(only.entrys)+5.980),show.legend = NA, #<---------------------
                      position = "identity")+
        scale_edge_colour_manual(values=col.set.link,name = "Biological Process")+
        geom_node_text(aes(label=Etiqueta),
                       position ="identity",
                       size = 3, #<---------------------
                       repel = T,
                       #vjust= fun.vjust(-0.04259*length(only.entrys)+5.980), # 3
                       #hjust= fun.hjust(-0.04259*length(only.entrys)+5.980), # .3
                       fontface = "bold")+
        theme_classic() +
        theme(legend.key.width=unit(.4,"cm"),
              legend.key.height=unit(.4,"cm"),
              legend.title=element_text(size=10,face="bold"),
              legend.text=element_text(size=8),
              axis.line=element_blank(),
              axis.text = element_blank(),
              axis.ticks = element_blank(),
              axis.title = element_blank())
      graf1 = graf1 + scale_color_gradient(low = c("green","black"), high = "red") +
        labs(caption = paste("NeVOmics:",sep=" ",format(Sys.time(),"%a(%d)/%b/%Y   %X")))
      {
        png(file=paste0("zxcvbnmNeVOmics_Plot_Network_Layout_XL_",i,".png"),
            width = width.png,  #<---------------------
            height = 15, 
            units = "cm",
            res=900,bg="white")
        plot(graf1)
        dev.off()}
    }
    #
  }
  #\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
  #___________________________________________________________________________________________ Complex plots
  #
  if (total.terms <= 30) {
    # --------------------------------------------------------------------- one bar configuration
    #
    # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    #
    # |||||||||      Bar plot
    #
    # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    #
    {
      for.bar = data.frame(c(data.frame("Entry"=c(nodes$Entry[1:total.terms],rep(NA,31-total.terms))),
                             data.frame("Freq"=c(nodes$Freq[1:total.terms],rep(NA,31-total.terms))),
                             data.frame("Fdr"=c(nodes$FDR[1:total.terms],rep(NA,31-total.terms))),
                             data.frame("Term"=c(nodes$Term[1:total.terms],rep(NA,31-total.terms))),
                             data.frame("P"=c(nodes$P[1:total.terms],rep(NA,31-total.terms)))))
      for.bar["log.pval"] = -log10(for.bar$P)
      for.bar["log.fdr"] = -log10(for.bar$Fdr)
    }
    ###
    #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    {
      barplot.enrichment = ggplot(for.bar) +
        geom_bar(aes(x = Entry, weight = log.pval),fill=rev(many.colors[1:(total.terms+1)]),width=.7) +
        geom_line(aes(x = c(1:31), y = -log10(0.05),
                      color=-log10(0.05)),
                  color = c(rep("transparent",31-total.terms),
                            rep("black",total.terms)), size = .2 )+
        geom_point(aes(x = c(1:31), y = -log10(0.05),color=-log10(0.05)),
                   position = position_dodge(0.3),
                   color = c(rep("transparent",31-total.terms),
                             rep("black",total.terms)), size = 1)+
        geom_line(aes(x = c(1:31), y = rev(log.fdr),
                      color=log.fdr),
                  color = c(rep("transparent",31-total.terms),
                            rep("blue",total.terms)), size = .2 )+
        geom_point(aes(x = c(1:31), y = rev(log.fdr),color=log.fdr),
                   position = position_dodge(0.3),shape=17,
                   color = c(rep("transparent",31-total.terms),
                             rep("blue",total.terms)), size = 1)+
        #scale_y_continuous(position = "right")+
        coord_flip()+
        scale_x_discrete(limits=rev(for.bar$Entry))+
        geom_text(aes(x=Entry,y=log.pval+.35,label=Freq),size=2.5,fontface ="bold")+ #<---------------------
      #ggtitle("-log(P-value) and\n -log(Corrected P-value)")+
      theme(axis.title.x =element_blank(),# element_text(colour="black",size=27,face="bold",hjust = 0,vjust = 10),
            axis.title.y = element_blank(),
            #plot.title = element_text(vjust = 4,hjust=.08,size = 12),
            axis.text.y = element_text(colour=c("white",rep("black",30)),size=8,
                                       face="bold",hjust = 1.3,vjust = .5), #element_text(size = 8,color="black"),#tama?o y color de texto de eje x
            axis.text.x = element_text(colour="black",size=7,face="bold",vjust = 0),
            axis.ticks.length =unit(.1, "cm"), #reduce la longitud de las marcas en la linea de escala de los ejes x y
            axis.ticks = element_line(color = "black",size = 1),
            axis.ticks.y.left = element_line(colour = "black",size=.5,unit(0, "cm")),
            #axis.ticks.x.top = element_line(colour = "black",size=.5,unit(3, "cm")),
            axis.line.y = element_blank(), axis.title.x.top =  element_blank(), 
            plot.margin = margin(1.5,1,0.5,0, "cm"))+
        scale_y_continuous(position = "right",
                           limit = c(0,round(max(for.bar$log.pval[1:total.terms]+1.5))),
                           breaks=seq(0,round(max(for.bar$log.pval[1:total.terms]+1.5),digits = 0),1))+
        labs(title = expression(bold("ASPECT")),
             subtitle = expression(bold("-log10("*italic("P")*"-value & adj."*italic("P")*"-value)")))+
        theme(plot.title = element_text(vjust = -.7,hjust=1.5,size = 12),
              plot.subtitle = element_text(vjust = -.7,hjust=.6,size = 8))
    }
    ###
    #
    # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    #
    # |||||||||      Complex Plots
    #
    # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    #
    layouts=c("nicely","kk","dh")
    #
    ###
    for (i in layouts){
      graf2 = ggraph(link_tbl,layout = "igraph",algorithm = i) +
        geom_edge_link(colour = col.edges.cluster,alpha = 1,
                       width = correl.width(-0.01296*length(only.entrys)+1.798))+ #<---------------------
      geom_node_point(colour= col.set.node,shape=c(rep(16,total.terms),rep(20,length(only.entrys))),
                      size= correl.point(-0.04259*length(only.entrys)+5.980), #<---------------------
                      position = "identity")+
        geom_node_text(aes(label=Label),
                       position ="identity",
                       size = correl.text(-0.02778*length(only.entrys)+4.639), #<---------------------
                       repel = T,
                       #vjust=3,
                       #hjust=.3,
                       fontface = "bold")+
        theme_classic()
      graf2 = graf2 +  theme(legend.key.width=unit(.4,"cm"),
                             legend.key.height=unit(.4,"cm"),
                             legend.position = c(1.06, .9),
                             legend.title=element_text(size=10,face="bold"),
                             legend.text=element_text(size=8),
                             axis.line=element_blank(),
                             axis.text = element_blank(),
                             axis.ticks = element_blank(),
                             axis.title = element_blank())
      graf2 = graf2 + scale_color_gradient(low = c("green","black"), high = "red") +
        labs(caption = paste("NeVOmics:",sep=" ",format(Sys.time(),"%a(%d)/%b/%Y   %X")))
      {
        # combine plots
        plot.with.inset = ggdraw() +
          draw_plot(graf2,width = .725, height = 1) +
          #                         x , y , width , height
          draw_plot(barplot.enrichment, .71,.01, .25, .87)
        ## Save plot
        png(file=paste("zxcvbnmNeVOmics_Plot_Network_Bar_Des_Layout_",i,".png"),
            width = 25,
            height = 15,
            units = "cm",
            res=900,bg="white")
        plot(plot.with.inset)
        dev.off()}
    }
    #
    for (i in layouts){
      graf2 = ggraph(link_tbl,layout = "igraph",algorithm = i) +
        geom_edge_link(colour = col.edges.cluster,alpha = 1,
                       width = correl.width(-0.01296*length(only.entrys)+1.798))+ #<---------------------
      geom_node_point(colour= col.set.node,shape=c(rep(16,total.terms),rep(20,length(only.entrys))),
                      size= correl.point(-0.04259*length(only.entrys)+5.980), #<---------------------
                      position = "identity")+
        geom_node_text(aes(label=name),
                       position ="identity",
                       size = correl.text(-0.02778*length(only.entrys)+4.639), #<---------------------
                       repel = T,
                       #vjust=3,
                       #hjust=.3,
                       fontface = "bold")+
        theme_classic()
      graf2 = graf2 +  theme(legend.key.width=unit(.4,"cm"),
                             legend.key.height=unit(.4,"cm"),
                             legend.position = c(1.06, .9),
                             legend.title=element_text(size=10,face="bold"),
                             legend.text=element_text(size=8),
                             axis.line=element_blank(),
                             axis.text = element_blank(),
                             axis.ticks = element_blank(),
                             axis.title = element_blank())
      graf2 = graf2 + scale_color_gradient(low = c("green","black"), high = "red") +
        labs(caption = paste("NeVOmics:",sep=" ",format(Sys.time(),"%a(%d)/%b/%Y   %X")))
      {
        # combine plots
        plot.with.inset = ggdraw() +
          draw_plot(graf2,width = .725, height = 1) +
          #                         x , y , width , height
          draw_plot(barplot.enrichment, .71,.01, .25, .87)
        ## Save plot
        png(file=paste("zxcvbnmNeVOmics_Plot_Network_Bar_Id_Layout_",i,".png"),
            width = 25,
            height = 15,
            units = "cm",
            res=900,bg="white")
        plot(plot.with.inset)
        dev.off()}
    }
    #
    for (i in layouts){
      graf2 = ggraph(link_tbl,layout = "igraph",algorithm = i) +
        geom_edge_link(colour = col.edges.cluster,alpha = 1,
                       width = correl.width(-0.01296*length(only.entrys)+1.798))+ #<---------------------
      geom_node_point(colour= col.set.node,shape=c(rep(16,total.terms),rep(20,length(only.entrys))),
                      size= correl.point(-0.04259*length(only.entrys)+5.980), #<---------------------
                      position = "identity")+
        geom_node_text(aes(label=Entry),
                       position ="identity",
                       size = correl.text(-0.02778*length(only.entrys)+4.639), #<---------------------
                       repel = T,
                       #vjust=3,
                       #hjust=.3,
                       fontface = "bold")+
        theme_classic()
      graf2 = graf2 +  theme(legend.key.width=unit(.4,"cm"),
                             legend.key.height=unit(.4,"cm"),
                             legend.position = c(1.06, .9),
                             legend.title=element_text(size=10,face="bold"),
                             legend.text=element_text(size=8),
                             axis.line=element_blank(),
                             axis.text = element_blank(),
                             axis.ticks = element_blank(),
                             axis.title = element_blank())
      graf2 = graf2 + scale_color_gradient(low = c("green","black"), high = "red") +
        labs(caption = paste("NeVOmics:",sep=" ",format(Sys.time(),"%a(%d)/%b/%Y   %X")))
      {
        # combine plots
        plot.with.inset = ggdraw() +
          draw_plot(graf2,width = .725, height = 1) +
          #                         x , y , width , height
          draw_plot(barplot.enrichment, .71,.01, .25, .87)
        ## Save plot
        png(file=paste("zxcvbnmNeVOmics_Plot_Network_Bar_Layout_",i,".png"),
            width = 25,
            height = 15,
            units = "cm",
            res=900,bg="white")
        plot(plot.with.inset)
        dev.off()}
    }
    #
    #
    #
    #\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    if (length(only.entrys) > 80) {
      modelos = c('nicely','kk','dh', 'drl')
      V(link_tbl)$Etiqueta = c(as.character(nodes$Term[1:total.terms]),
                               as.character(rep(NA,total.nodes - total.terms)))
      
      for (i in modelos){
        graf2 = ggraph(link_tbl,layout = "igraph",algorithm = i) +
          geom_edge_link(colour = col.edges.cluster,alpha = 1,
                         width = correl.width(-0.01296*length(only.entrys)+1.798))+ #<---------------------
        geom_node_point(colour= col.set.node,shape=c(rep(16,total.terms),rep(20,length(only.entrys))),
                        size= correl.point(-0.04259*length(only.entrys)+5.980), #<---------------------
                        position = "identity")+
          geom_node_text(aes(label=Etiqueta),
                         position ="identity",
                         size = 3, #<---------------------
                         repel = T,
                         #vjust=3,
                         #hjust=.3,
                         fontface = "bold")+
          theme_classic()
        graf2 = graf2 +  theme(legend.key.width=unit(.4,"cm"),
                               legend.key.height=unit(.4,"cm"),
                               legend.position = c(1.06, .9),
                               legend.title=element_text(size=10,face="bold"),
                               legend.text=element_text(size=8),
                               axis.line=element_blank(),
                               axis.text = element_blank(),
                               axis.ticks = element_blank(),
                               axis.title = element_blank())
        graf2 = graf2 + scale_color_gradient(low = c("green","black"), high = "red") +
          labs(caption = paste("NeVOmics:",sep=" ",format(Sys.time(),"%a(%d)/%b/%Y   %X")))
        {
          # combine plots
          plot.with.inset = ggdraw() +
            draw_plot(graf2,width = .725, height = 1) +
            #                         x , y , width , height
            draw_plot(barplot.enrichment, .71,.01, .25, .87)
          ## Save plot
          png(file=paste("zxcvbnmNeVOmics_Plot_Network_Bar_Des_Layout_XL_",i,".png"),
              width = 25,
              height = 15,
              units = "cm",
              res=900,bg="white")
          plot(plot.with.inset)
          dev.off()}
      }
    }
    #\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    #
    #
    #
  } else if ((total.terms >= 31) && (total.terms <= 60)) {
    # ------------------------------------------------------------------------  two bars configuration
    #
    #
    # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    #
    # |||||||||      Bar plots
    #
    # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    #
    nodes.1_30 = nodes[1:30,]
    #
    {
      for.bar.1_30 = data.frame(c(data.frame("Entry"=c(nodes.1_30$Entry,rep(NA,31-nrow(nodes.1_30)))),
                                  data.frame("Freq"=c(nodes.1_30$Freq,rep(NA,31-nrow(nodes.1_30)))),
                                  data.frame("Fdr"=c(nodes.1_30$FDR,rep(NA,31-nrow(nodes.1_30)))),
                                  data.frame("Term"=c(nodes.1_30$Term,rep(NA,31-nrow(nodes.1_30)))),
                                  data.frame("P"=c(nodes.1_30$P,rep(NA,31-nrow(nodes.1_30))))))
      for.bar.1_30["log.pval"] = -log10(for.bar.1_30$P)
      for.bar.1_30["log.fdr"] = -log10(for.bar.1_30$Fdr)
    }
    #
    barplot.p.val.1_30 = ggplot(for.bar.1_30) +
      geom_bar(aes(x = Entry, weight = log.pval),fill=rev(many.colors[1:(nrow(nodes.1_30)+1)]),width=.7) +
      geom_line(aes(x = c(1:31), y = -log10(0.05),
                    color=-log10(0.05)),
                color = c(rep("transparent",31-nrow(nodes.1_30)),
                          rep("black",nrow(nodes.1_30))), size = .2 )+
      geom_point(aes(x = c(1:31), y = -log10(0.05),color=-log10(0.05)),
                 position = position_dodge(0.3),
                 color = c(rep("transparent",31-nrow(nodes.1_30)),
                           rep("black",nrow(nodes.1_30))), size = 1)+
      geom_line(aes(x = c(1:31), y = rev(log.fdr),
                    color=log.fdr),
                color = c(rep("transparent",31-nrow(nodes.1_30)),
                          rep("blue",nrow(nodes.1_30))), size = .2 )+
      geom_point(aes(x = c(1:31), y = rev(log.fdr),color=log.fdr),
                 position = position_dodge(0.3),shape=17,
                 color = c(rep("transparent",31-nrow(nodes.1_30)),
                           rep("blue",nrow(nodes.1_30))), size = 1)+
      #scale_y_continuous(position = "right")+
      coord_flip()+
      scale_x_discrete(limits=rev(for.bar.1_30$Entry))+
      geom_text(aes(x=Entry,y=log.pval+.35,label=Freq),size=2.5,fontface ="bold")+ #<---------------------
    #ggtitle("-log(P-value) and\n -log(Corrected P-value)")+
    theme(plot.title = element_blank(),
          axis.title.x = element_blank(),# element_text(colour="black",size=27,face="bold",hjust = 0,vjust = 10),
          axis.title.y = element_blank(),
          axis.text.y = element_text(colour=c("white",rep("black",30)),size=8,
                                     face="bold",hjust = 1.3,vjust = .5), #element_text(size = 8,color="black"),#tama?o y color de texto de eje x
          axis.text.x = element_text(colour="black",size=7,face="bold",vjust = 0),
          axis.ticks.length =unit(.1, "cm"),
          axis.ticks = element_line(color = "black",size = 1),
          axis.ticks.y.left = element_line(colour = "black",size=.5,unit(0, "cm")),
          #axis.ticks.x.top = element_line(colour = "black",size=.5,unit(3, "cm")),
          axis.line.y = element_blank(), axis.title.x.top =  element_blank(), 
          plot.margin = margin(1.5,1,0.5,0, "cm"))+
      scale_y_continuous(position = "right",
                         limit = c(0,round(max(for.bar.1_30$log.pval[1:nrow(nodes.1_30)]+1.5))),
                         breaks=seq(0,round(max(for.bar.1_30$log.pval[1:nrow(nodes.1_30)]+1.5),digits = 0),1))+
      labs(title = expression(bold("ASPECT")),
           subtitle = expression(bold("-log10("*italic("P")*"-value & adj."*italic("P")*"-value)")))+
      theme(plot.title = element_text(vjust = 2,hjust=-.3,size = 15),
            plot.subtitle = element_text(vjust = -.7,hjust=.6,size = 8))
    #
    # Second part of the data
    #
    nodes.31_final = nodes[31:total.terms,]
    #
    {
      for.bar.31_final = data.frame(c(data.frame("Entry"=c(nodes.31_final$Entry,rep(NA,31-nrow(nodes.31_final)))),
                                      data.frame("Freq"=c(nodes.31_final$Freq,rep(NA,31-nrow(nodes.31_final)))),
                                      data.frame("Fdr"=c(nodes.31_final$FDR,rep(NA,31-nrow(nodes.31_final)))),
                                      data.frame("Term"=c(nodes.31_final$Term,rep(NA,31-nrow(nodes.31_final)))),
                                      data.frame("P"=c(nodes.31_final$P,rep(NA,31-nrow(nodes.31_final))))))
      for.bar.31_final["log.pval"] = -log10(for.bar.31_final$P)
      for.bar.31_final["log.fdr"] = -log10(for.bar.31_final$Fdr)
    }
    
    barplot.p.val.31_final = ggplot(for.bar.31_final) +
      geom_bar(aes(x = Entry, weight = log.pval),fill=rev(many.colors[(nrow(nodes.1_30)+1):(total.terms+1)]),width=.7) +
      geom_line(aes(x = c(1:31), y = -log10(0.05),
                    color=-log10(0.05)),
                color = c(rep("transparent",31-nrow(nodes.31_final)),
                          rep("black",nrow(nodes.31_final))), size = .2 )+
      geom_point(aes(x = c(1:31), y = -log10(0.05),color=-log10(0.05)),
                 position = position_dodge(0.3),
                 color = c(rep("transparent",31-nrow(nodes.31_final)),
                           rep("black",nrow(nodes.31_final))), size = 1)+
      geom_line(aes(x = c(1:31), y = rev(log.fdr),
                    color=log.fdr),
                color = c(rep("transparent",31-nrow(nodes.31_final)),
                          rep("blue",nrow(nodes.31_final))), size = .2 )+
      geom_point(aes(x = c(1:31), y = rev(log.fdr),color=log.fdr),
                 position = position_dodge(0.3),shape=17,
                 color = c(rep("transparent",31-nrow(nodes.31_final)),
                           rep("blue",nrow(nodes.31_final))), size = 1)+
      #scale_y_continuous(position = "right")+
      coord_flip()+
      scale_x_discrete(limits=rev(for.bar.31_final$Entry))+
      geom_text(aes(x=Entry,y=log.pval+.35,label=Freq),size=2.5,fontface ="bold")+ #<---------------------
    #ggtitle("-log(P-value) and\n -log(Corrected P-value)")+
    theme(plot.title = element_text(colour="white"),
          axis.title.x =element_blank(),# element_text(colour="black",size=27,face="bold",hjust = 0,vjust = 10),
          axis.title.y = element_blank(),
          axis.text.y = element_text(colour=c("white",rep("black",30)),size=8,
                                     face="bold",hjust = 1.3,vjust = .5), #element_text(size = 8,color="black"),#tama?o y color de texto de eje x
          axis.text.x = element_text(colour="black",size=7,face="bold",vjust = 0),
          axis.ticks.length =unit(.1, "cm"),
          axis.ticks = element_line(color = "black",size = 1),
          axis.ticks.y.left = element_line(colour = "black",size=.5,unit(0, "cm")),
          #axis.ticks.x.top = element_line(colour = "black",size=3,unit(3, "cm")),
          axis.line.y = element_blank(), axis.title.x.top =  element_blank(), 
          plot.margin = margin(1.5,1,0.5,0, "cm"))+
      scale_y_continuous(position = "right",
                         limit = c(0,round(max(for.bar.1_30$log.pval[1:nrow(nodes.1_30)]+1.5))),
                         breaks=seq(0,round(max(for.bar.1_30$log.pval[1:nrow(nodes.1_30)]+1.5),digits = 0),1))+
      labs(title = expression(bold("xxxxxxxx")),
           subtitle = expression(bold("-log10("*italic("P")*"-value & adj."*italic("P")*"-value)")))+
      theme(plot.title = element_text(vjust = 0.5,hjust=1.5,size = 12),
            plot.subtitle = element_text(vjust = -.7,hjust=.6,size = 8))
    #
    #
    #
    # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    #
    # |||||||||      Complex Plots
    #
    # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    #
    #
    graf2 = ggraph(link_tbl,layout = "igraph",algorithm = "kk") +
      geom_edge_link(colour = col.edges.cluster,alpha = 1,
                     width = correl.width(-0.01296*length(only.entrys)+1.798))+ #<---------------------
    geom_node_point(colour= col.set.node,shape=c(rep(16,total.terms),rep(20,length(only.entrys))),
                    size= correl.point(-0.04259*length(only.entrys)+5.980), #<---------------------
                    position = "identity")+
      geom_node_text(aes(label=Label),
                     position ="identity",
                     size = correl.text(-0.02778*length(only.entrys)+4.639), #<---------------------
                     repel = T,
                     #vjust=3,
                     #hjust=.3,
                     fontface = "bold")+
      theme_classic()
    graf2 = graf2 +  theme(legend.key.width=unit(.4,"cm"),
                           legend.key.height=unit(.4,"cm"),
                           legend.position = c(1.06, .9),
                           legend.title=element_text(size=10,face="bold"),
                           legend.text=element_text(size=8),
                           axis.line=element_blank(),
                           axis.text = element_blank(),
                           axis.ticks = element_blank(),
                           axis.title = element_blank())
    graf2 = graf2 + scale_color_gradient(low = c("green","black"), high = "red") +
      labs(caption = paste("NeVOmics:",sep=" ",format(Sys.time(),"%a(%d)/%b/%Y   %X")))
    #
    ###
    for (i in layouts){
      graf2 = ggraph(link_tbl,layout = "igraph",algorithm = i) +
        geom_edge_link(colour = col.edges.cluster,alpha = 1,
                       width = correl.width(-0.01296*length(only.entrys)+1.798))+ #<---------------------
      geom_node_point(colour= col.set.node,shape=c(rep(16,total.terms),rep(20,length(only.entrys))),
                      size= correl.point(-0.04259*length(only.entrys)+5.980), #<---------------------
                      position = "identity")+
        geom_node_text(aes(label=Label),
                       position ="identity",
                       size = correl.text(-0.02778*length(only.entrys)+4.639), #<---------------------
                       repel = T,
                       #vjust=3,
                       #hjust=.3,
                       fontface = "bold")+
        theme_classic()
      graf2 = graf2 +  theme(legend.key.width=unit(.4,"cm"),
                             legend.key.height=unit(.4,"cm"),
                             legend.position = c(1.06, .9),
                             legend.title=element_text(size=10,face="bold"),
                             legend.text=element_text(size=8),
                             axis.line=element_blank(),
                             axis.text = element_blank(),
                             axis.ticks = element_blank(),
                             axis.title = element_blank())
      graf2 = graf2 + scale_color_gradient(low = c("green","black"), high = "red") +
        labs(caption = paste("NeVOmics:",sep=" ",format(Sys.time(),"%a(%d)/%b/%Y   %X")))
      {
        plots_combined  = ggdraw() +
          draw_plot(barplot.p.val.1_30, x = 0, y = 0, width = .5, height = 1) +
          draw_plot(barplot.p.val.31_final, x = .45, y = 0, width = .5, height = 1)
        # combine plots
        plot.with.inset = ggdraw() +
          draw_plot(graf2,  width = .54, height = 1) +
          #                         x , y , width , height
          draw_plot(plots_combined, .53,.01, .45, .85) # 2 barplots
        
        ## Save plot
        png(file=paste("zxcvbnmNeVOmics_Plot_Network_Bar_Des_Layout_",i,".png"),
            width = 35,
            height = 15,
            units = "cm",
            res=900,bg="white")
        plot(plot.with.inset)
        dev.off()}
    }
    #
    for (i in layouts){
      graf2 = ggraph(link_tbl,layout = "igraph",algorithm = i) +
        geom_edge_link(colour = col.edges.cluster,alpha = 1,
                       width = correl.width(-0.01296*length(only.entrys)+1.798))+ #<---------------------
      geom_node_point(colour = col.set.node,shape=c(rep(16,total.terms),rep(20,length(only.entrys))),
                      size= correl.point(-0.04259*length(only.entrys)+5.980), #<---------------------
                      position = "identity")+
        geom_node_text(aes(label=name),
                       position ="identity",
                       size = correl.text(-0.02778*length(only.entrys)+4.639), #<---------------------
                       repel = T,
                       #vjust=3,
                       #hjust=.3,
                       fontface = "bold")+
        theme_classic()
      graf2 = graf2 +  theme(legend.key.width=unit(.4,"cm"),
                             legend.key.height=unit(.4,"cm"),
                             legend.position = c(1.06, .9),
                             legend.title=element_text(size=10,face="bold"),
                             legend.text=element_text(size=8),
                             axis.line=element_blank(),
                             axis.text = element_blank(),
                             axis.ticks = element_blank(),
                             axis.title = element_blank())
      graf2 = graf2 + scale_color_gradient(low = c("green","black"), high = "red") +
        labs(caption = paste("NeVOmics:",sep=" ",format(Sys.time(),"%a(%d)/%b/%Y   %X")))
      {
        plots_combined  = ggdraw() +
          draw_plot(barplot.p.val.1_30, x = 0, y = 0, width = .5, height = 1) +
          draw_plot(barplot.p.val.31_final, x = .45, y = 0, width = .5, height = 1)
        # combine plots
        plot.with.inset = ggdraw() +
          draw_plot(graf2,  width = .54, height = 1) +
          #                         x , y , width , height
          draw_plot(plots_combined, .53,.01, .45, .85) # 2 barplots
        
        ## Save plot
        png(file=paste("zxcvbnmNeVOmics_Plot_Network_Bar_Id_Layout_",i,".png"),
            width = 35,
            height = 15,
            units = "cm",
            res=900,bg="white")
        plot(plot.with.inset)
        dev.off()}
    }
    #
    for (i in layouts){
      graf2 = ggraph(link_tbl,layout = "igraph",algorithm = i) +
        geom_edge_link(colour = col.edges.cluster,alpha = 1,
                       width = correl.width(-0.01296*length(only.entrys)+1.798))+ #<---------------------
      geom_node_point(aes(colour= logFC),
                      size= correl.point(-0.04259*length(only.entrys)+5.980), #<---------------------
                      position = "identity")+
        geom_node_text(aes(label=Entry),
                       position ="identity",
                       size = correl.text(-0.02778*length(only.entrys)+4.639), #<---------------------
                       repel = T,
                       #vjust=3,
                       #hjust=.3,
                       fontface = "bold")+
        theme_classic()
      graf2 = graf2 +  theme(legend.key.width=unit(.4,"cm"),
                             legend.key.height=unit(.4,"cm"),
                             legend.position = c(1.06, .9),
                             legend.title=element_text(size=10,face="bold"),
                             legend.text=element_text(size=8),
                             axis.line=element_blank(),
                             axis.text = element_blank(),
                             axis.ticks = element_blank(),
                             axis.title = element_blank())
      graf2 = graf2 + scale_color_gradient(low = c("green","black"), high = "red") +
        labs(caption = paste("NeVOmics:",sep=" ",format(Sys.time(),"%a(%d)/%b/%Y   %X")))
      {
        plots_combined  = ggdraw() +
          draw_plot(barplot.p.val.1_30, x = 0, y = 0, width = .5, height = 1) +
          draw_plot(barplot.p.val.31_final, x = .45, y = 0, width = .5, height = 1)
        # combine plots
        plot.with.inset = ggdraw() +
          draw_plot(graf2,  width = .54, height = 1) +
          #                         x , y , width , height
          draw_plot(plots_combined, .53,.01, .45, .85) # 2 barplots
        
        ## Save plot
        png(file=paste("zxcvbnmNeVOmics_Plot_Network_Bar_Layout_",i,".png"),
            width = 35,
            height = 15,
            units = "cm",
            res=900,bg="white")
        plot(plot.with.inset)
        dev.off()}
    }
    #>
    #
    #
    #\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    if (length(only.entrys) > 80) {
      modelos = c('nicely','kk','dh', 'drl')
      V(link_tbl)$Etiqueta = c(as.character(nodes$Term[1:total.terms]),
                               as.character(rep(NA,total.nodes - total.terms)))
      for (i in modelos){
        graf2 = ggraph(link_tbl,layout = "igraph",algorithm = i) +
          geom_edge_link(colour = col.edges.cluster,alpha = 1,
                         width = correl.width(-0.01296*length(only.entrys)+1.798))+ #<---------------------
        geom_node_point(colour= col.set.node,shape=c(rep(16,total.terms),rep(20,length(only.entrys))),
                        size= correl.point(-0.04259*length(only.entrys)+5.980), #<---------------------
                        position = "identity")+
          geom_node_text(aes(label=Etiqueta),
                         position ="identity",
                         size = 3, #<---------------------
                         repel = T,
                         #vjust=3,
                         #hjust=.3,
                         fontface = "bold")+
          theme_classic()
        graf2 = graf2 +  theme(legend.key.width=unit(.4,"cm"),
                               legend.key.height=unit(.4,"cm"),
                               legend.position = c(1.06, .9),
                               legend.title=element_text(size=10,face="bold"),
                               legend.text=element_text(size=8),
                               axis.line=element_blank(),
                               axis.text = element_blank(),
                               axis.ticks = element_blank(),
                               axis.title = element_blank())
        graf2 = graf2 + scale_color_gradient(low = c("green","black"), high = "red") +
          labs(caption = paste("NeVOmics:",sep=" ",format(Sys.time(),"%a(%d)/%b/%Y   %X")))
        {
          plots_combined  = ggdraw() +
            draw_plot(barplot.p.val.1_30, x = 0, y = 0, width = .5, height = 1) +
            draw_plot(barplot.p.val.31_final, x = .45, y = 0, width = .5, height = 1)
          # combine plots
          plot.with.inset = ggdraw() +
            draw_plot(graf2,  width = .54, height = 1) +
            #                         x , y , width , height
            draw_plot(plots_combined, .53,.01, .45, .85) # 2 barplots
          
          ## Save plot
          png(file=paste("zxcvbnmNeVOmics_Plot_Network_Bar_Des_Layout_XL_",i,".png"),
              width = 35,
              height = 15,
              units = "cm",
              res=900,bg="white")
          plot(plot.with.inset)
          dev.off()}
      }
    }
    #\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    #
    #
    #
    #
  }
  #
  #___________________________________________________________________________________________ Chrod plots
  #
  #
  # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  #
  # |||||||||     Chord plots (with expression values)
  #
  # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  #
  # ---------------------------------------------------------- chord 1 
  
  if (length(only.entrys) <= 70) {
    if ((total.terms > 1) && (total.terms <= 15)) {
      print("grafico 1 y 2")
      for (i in -1:0) {
        #
        #chord 1
        #
        go.term = drop_na(nodes[,c("Entry","Term")])
        colnames(go.term)[1] <- "GO"
        #
        go.entry.term = merge(links,go.term,by="GO")#[c(3,2)]
        #
        go.paste.term=paste0(go.entry.term$GO," ~ ",go.entry.term$Term)
        #
        term.entry = data.frame("Term"=go.paste.term,"Entry"=go.entry.term$Entry)
        term.entry=as.data.frame.matrix(table(term.entry)) 
        term.entry.mat = as.matrix(term.entry) 
        #########
        row_sum = sum(rowSums(abs(term.entry.mat)))
        col_sum = sum(colSums(abs(term.entry.mat)))
        small_gap = 1
        big_gap = 15
        nr = nrow(term.entry.mat)
        nc = ncol(term.entry.mat)
        n_sector = nr + nc
        row_sector_degree = (360 - small_gap*(n_sector - 2) - big_gap*2) * (row_sum/(row_sum + col_sum)) +
          small_gap*(nr-1)
        ###
        start_degree = 90 - (180 - row_sector_degree)/2
        #
        circos.clear()
        # gap entre cada cuerda
        circos.par(gap.after = c(rep(1, nrow(term.entry.mat)-1), 15, rep(1, ncol(term.entry.mat)-1),15),
                   #circos.par(gap.after = c(rep(0.5, nrow(a)-1), 20, rep(0.5, ncol(a)-1),20),
                   start.degree = start_degree)
        #circos.par(start.degree = 80, clock.wise = T)
        circos.par(track.margin = c(0,0))
        chord.plot.1 = function () {
          chordDiagram(term.entry.mat, order = c(paste0(go.term$GO," ~ ",go.term$Term),only.entrys),
                       grid.col = col.set.link,
                       transparency = 0,
                       directional = i,
                       annotationTrack = "grid",
                       #    c(% del datio del circulo, alto del track)
                       annotationTrackHeight = c(0.08, .01),
                       preAllocateTracks = 2)
          circos.track(track.index = 2,ylim = c(0, 1), track.height = uh(5, "cm"),
                       track.margin = c(5, uh(5, "cm")),
                       panel.fun = function(x, y) {
                         xcenter = get.cell.meta.data("xcenter")
                         circos.lines(c(xcenter, xcenter), c(0, uy(.3, "cm")),
                                      col = 1:4)
                       },bg.border =  NA)
          circos.track(track.index = 2, panel.fun = for(si in get.all.sector.index()[1:total.terms]) {
            xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 2)
            ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 2)
            circos.text(mean(xlim), mean(ylim), si, sector.index = si, track.index = 2,
                        facing = "downward",
                        niceFacing = TRUE,
                        adj = c(0,0.5),
                        cex = .5, #<---------------------
                        font = 2,
                        col = NA)
            
          }, bg.border =  NA)
          circos.track(track.index = 2, panel.fun = for(si in get.all.sector.index()[start.entry:total.nodes]) {
            xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 2)
            ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 2)
            circos.text(mean(xlim), mean(ylim), si, sector.index = si, track.index = 2,
                        facing = "clockwise",
                        niceFacing = TRUE,
                        adj = c(0,0.5),
                        cex = -0.003039*length(only.entrys) + 0.4913, #<---------------------
                        font = 2,
                        col = "black")
          }, bg.border = NA)
        }
        #
        #
        {
          png(file=paste0("zxcvbnmNeVOmics_Plot_Chord_1_",i,".png"),
              width = 30,
              height = 10,
              units = "cm",
              res=1024,bg="white")
          chord.plot.1()
          text(x=.4, y=.9, "ASPECT", adj = c(0,0), 
               cex = .9,font = 2) #<---------------------
          dev.off()}
        #
        # ------------------------------------------------------------ chord 2
        #
        go.term = drop_na(nodes[,c("Entry","Term")])
        colnames(go.term)[1] <- "GO"
        
        go.entry.term = merge(links,go.term,by="GO")#[c(3,2)]
        
        ##
        
        go.entry = data.frame("GO"=go.entry.term$GO,"Entry"=go.entry.term$Entry)
        go.entry=as.data.frame.matrix(table(go.entry)) 
        go.entry.mat = as.matrix(go.entry) 
        
        #########
        circos.clear()
        row_sum = sum(rowSums(abs(go.entry.mat)))
        col_sum = sum(colSums(abs(go.entry.mat)))
        small_gap = 1
        big_gap = 5
        nr = nrow(go.entry.mat)
        nc = ncol(go.entry.mat)
        n_sector = nr + nc
        row_sector_degree = (360 - small_gap*(n_sector - 2) - big_gap*2) * (row_sum/(row_sum + col_sum)) +
          small_gap*(nr-1)
        ###
        start_degree = 90 - (180 - row_sector_degree)/2
        #
        # gap entre cada cuerda
        circos.par(gap.after = c(rep(.5, nrow(go.entry.mat)-1), 5, rep(.5, ncol(go.entry.mat)-1),5),
                   #circos.par(gap.after = c(rep(0.5, nrow(a)-1), 20, rep(0.5, ncol(a)-1),20),
                   start.degree = start_degree)
        #circos.par(start.degree = 80, clock.wise = T)
        circos.par(track.margin = c(0,0))
        chord.plot.2 = function () {
          chordDiagram(go.entry.mat, order = nodes$Entry, # c(nodes$Term[1:total.terms],only.entrys)
                       grid.col = col.set.link,
                       transparency = 0,
                       directional = i,
                       annotationTrack = "grid",
                       #    c(% del datio del circulo, alto del track)
                       annotationTrackHeight = c(0.08, .01),
                       preAllocateTracks = 2)
          circos.track(track.index = 2,ylim = c(0, 1), track.height = uh(5, "cm"),
                       track.margin = c(5, uh(5, "cm")),
                       panel.fun = function(x, y) {
                         xcenter = get.cell.meta.data("xcenter")
                         circos.lines(c(xcenter, xcenter), c(0, uy(.2, "cm")),
                                      col = 1:4)
                       },bg.border =  NA)
          circos.track(track.index = 2, panel.fun = for(si in get.all.sector.index()[1:total.terms]) {
            xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 2)
            ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 2)
            circos.text(mean(xlim), mean(ylim), si, sector.index = si, track.index = 2,
                        facing = "clockwise",
                        niceFacing = TRUE,
                        adj = c(0.15,0.5),
                        cex = .5, #<---------------------
                        font = 2,
                        col = "black")
            
          }, bg.border =  NA)
          circos.track(track.index = 2, panel.fun = for(si in get.all.sector.index()[start.entry:total.nodes]) {
            xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 2)
            ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 2)
            circos.text(mean(xlim), mean(ylim), si, sector.index = si, track.index = 2,
                        facing = "clockwise",
                        niceFacing = TRUE,
                        adj = c(0.15,0.5),
                        cex = -0.003039*length(only.entrys) + 0.4913, #<---------------------
                        font = 2,
                        col = "black")
          }, bg.border = NA)
        }
        #
        {
          png(file=paste0("zxcvbnmNeVOmics_Plot_Chord_2_",i,".png"),
              width = 30,
              height = 10,
              units = "cm",
              res=900,bg="white")
          lgd_points = Legend(at = nodes$Term[1:total.terms], type = "grid",
                              size = unit(0, "mm"),
                              labels_gp = gpar(col = "white"),
                              legend_gp = gpar(col = "white"), 
                              title_position = "topleft",
                              title = "")
          lgd_list_vertical = packLegend(lgd_points)
          plot.new()
          circle_size = unit(1, "snpc") # snpc unit gives you a square region
          pushViewport(viewport(x = 0, y = 0.5, width = circle_size, height = circle_size,
                                just = c("left", "center")))
          par(omi = gridOMI(), new = TRUE)
          chord.plot.2()
          upViewport()
          pushViewport(viewport(x = circle_size, y = 0.5, width = grobWidth(lgd_list_vertical),
                                height = grobHeight(lgd_list_vertical), just = c("left", "center")))
          grid.draw(lgd_list_vertical)
          upViewport()
          text(x=1, y=1, "ASPECT", adj = c(0,0), 
               cex = .9,font = 2) #<---------------------
          legend(x=1.1,y=.95, 
                 legend=unique(paste0(nodes$Term[1:total.terms]," (",nodes$Freq[1:total.terms],")")),
                 cex=.6,bg="transparent", #<---------------------
                 border = NA,pt.bg=many.colors,bty = "n",
                 x.intersp=.7,y.intersp=.9, #<---------------------
                 col=NA,text.font = 1,pch=22,
                 pt.cex=1)
          dev.off()}
      }
    } else if (total.terms >= 16) {
      print("grafico 2")
      for (i in -1:0) {
        #
        #
        # --------------------------------------------------- chord 2
        #
        go.term = drop_na(nodes[,c("Entry","Term")])
        colnames(go.term)[1] <- "GO"
        
        go.entry.term = merge(links,go.term,by="GO")#[c(3,2)]
        
        ##
        
        go.entry = data.frame("GO"=go.entry.term$GO,"Entry"=go.entry.term$Entry)
        go.entry=as.data.frame.matrix(table(go.entry)) 
        go.entry.mat = as.matrix(go.entry) 
        
        #########
        circos.clear()
        row_sum = sum(rowSums(abs(go.entry.mat)))
        col_sum = sum(colSums(abs(go.entry.mat)))
        small_gap = 1
        big_gap = 5
        nr = nrow(go.entry.mat)
        nc = ncol(go.entry.mat)
        n_sector = nr + nc
        row_sector_degree = (360 - small_gap*(n_sector - 2) - big_gap*2) * (row_sum/(row_sum + col_sum)) +
          small_gap*(nr-1)
        ###
        start_degree = 90 - (180 - row_sector_degree)/2
        #
        # gap entre cada cuerda
        circos.par(gap.after = c(rep(.5, nrow(go.entry.mat)-1), 5, rep(.5, ncol(go.entry.mat)-1),5),
                   #circos.par(gap.after = c(rep(0.5, nrow(a)-1), 20, rep(0.5, ncol(a)-1),20),
                   start.degree = start_degree)
        #circos.par(start.degree = 80, clock.wise = T)
        circos.par(track.margin = c(0,0))
        chord.plot.2 = function () {
          chordDiagram(go.entry.mat, order = nodes$Entry, # c(nodes$Term[1:total.terms],only.entrys)
                       grid.col = col.set.link,
                       transparency = 0,
                       directional = i,
                       annotationTrack = "grid",
                       #    c(% del datio del circulo, alto del track)
                       annotationTrackHeight = c(0.08, .01),
                       preAllocateTracks = 2)
          circos.track(track.index = 2,ylim = c(0, 1), track.height = uh(5, "cm"),
                       track.margin = c(5, uh(5, "cm")),
                       panel.fun = function(x, y) {
                         xcenter = get.cell.meta.data("xcenter")
                         circos.lines(c(xcenter, xcenter), c(0, uy(.2, "cm")),
                                      col = 1:4)
                       },bg.border =  NA)
          circos.track(track.index = 2, panel.fun = for(si in get.all.sector.index()[1:total.terms]) {
            xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 2)
            ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 2)
            circos.text(mean(xlim), mean(ylim), si, sector.index = si, track.index = 2,
                        facing = "clockwise",
                        niceFacing = TRUE,
                        adj = c(0.15,0.5),
                        cex = .5, #<---------------------
                        font = 2,
                        col = "black")
            
          }, bg.border =  NA)
          circos.track(track.index = 2, panel.fun = for(si in get.all.sector.index()[start.entry:total.nodes]) {
            xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 2)
            ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 2)
            circos.text(mean(xlim), mean(ylim), si, sector.index = si, track.index = 2,
                        facing = "clockwise",
                        niceFacing = TRUE,
                        adj = c(0.15,0.5),
                        cex =  -0.003039*length(only.entrys) + 0.4913, #<---------------------
                        font = 2,
                        col = "black")
          }, bg.border = NA)
        }
        #
        {
          png(file=paste0("zxcvbnmNeVOmics_Plot_Chord_2_",i,".png"),
              width = 30,
              height = 10,
              units = "cm",
              res=900,bg="white")
          lgd_points = Legend(at = nodes$Term[1:total.terms], type = "grid",
                              size = unit(0, "mm"),
                              labels_gp = gpar(col = "white"),
                              legend_gp = gpar(col = "white"), 
                              title_position = "topleft",
                              title = "")
          lgd_list_vertical = packLegend(lgd_points)
          plot.new()
          circle_size = unit(1, "snpc") # snpc unit gives you a square region
          pushViewport(viewport(x = 0, y = 0.5, width = circle_size, height = circle_size,
                                just = c("left", "center")))
          par(omi = gridOMI(), new = TRUE)
          chord.plot.2()
          upViewport()
          pushViewport(viewport(x = circle_size, y = 0.5, width = grobWidth(lgd_list_vertical),
                                height = grobHeight(lgd_list_vertical), just = c("left", "center")))
          grid.draw(lgd_list_vertical)
          upViewport()
          text(x=1, y=1, "ASPECT", adj = c(0,0), 
               cex = .9,font = 2) #<---------------------
          legend(x=1.1,y=.95, 
                 legend=unique(paste0(nodes$Term[1:total.terms]," (",nodes$Freq[1:total.terms],")")),
                 cex=.5,bg="transparent", #<---------------------
                 border = NA,pt.bg=many.colors,bty = "n",
                 x.intersp=.7,y.intersp=.9, #<---------------------
                 col=NA,text.font = 1,pch=22,
                 pt.cex=1)
          dev.off()}
      }
    }
  } else {
    print('no chord plot')
  }
  
  #
  #___________________________________________________________________________________________ UpSet plot
  #
  #
  # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  #
  # |||||||||     Upset Plot 
  #
  # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  #
  #
  ups=as.data.frame.matrix(t(table(links)))
  #
  # size of = set names
  if ((total.terms >= 1) && (total.terms <= 7)) {
    set.name.size=2
  } else if ((total.terms >= 8) && (total.terms <= 15)) {
    set.name.size=1.5
  } else if ((total.terms >= 16) && (total.terms <= 25)) {
    set.name.size=1.25
  } else if (total.terms >= 25) {
    set.name.size=0.9
  }
  #
  # size of = point.size and line.size (point.size)
  if ((total.terms >= 1) && (total.terms <= 10)) {
    point.size=5
  } else if ((total.terms >= 11) && (total.terms <= 20)) {
    point.size=4.16
  } else if ((total.terms >= 21) && (total.terms <= 30)) {
    point.size=3.3
  } else if (total.terms >= 31) {
    point.size=2.5
  }
  #
  #
  {png(file="zxcvbnmNeVOmics_Plot_UpSetR.png",width = 10,height = 7,units = "in",res=900,bg="white")
    upset(ups,sets=nodes$Entry[1:total.terms],
          sets.bar.color = "blue",order.by ="freq",empty.intersections = NULL,
          point.size=point.size,mainbar.y.label="poiuytrewq",sets.x.label = "",
          main.bar.color="black",matrix.color="black",shade.color="red",
          line.size=point.size/5,show.numbers = "yes",group.by = "degree",
          matrix.dot.alpha = 0,mb.ratio = c(0.5, 0.5),
          #c(intersection size title, intersection size tick labels, set size title, set size tick labels, set names, numbers above bars)
          #            1                           2                      3                   4                5                6           
          text.scale = c(2,2,1,1.5,set.name.size,2), 
          scale.sets = "identity",scale.intersections = "identity",
          shade.alpha = 0.35,set_size.angles = 0)
    dev.off()}
  #>
  print(warnings())
} 
