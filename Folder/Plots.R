

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
}
#
#
# width of png file for simple plot
#
if ((total.terms > 1) && (total.terms <= 20)) {
  width.png=20
  print(paste("el ancho de la figura es:",width.png))
} else if ((total.terms >= 21) && (total.terms <= 30)) {
  width.png=22.5
  print(paste("el ancho de la figura es:",width.png))
} else if ((total.terms >= 31) && (total.terms <= 54)) {
  width.png=25
  print(paste("el ancho de la figura es:",width.png))
} else if (total.terms > 55) {
  width.png=30
  print(paste("el ancho de la figura es:",width.png))
}
#
#
# Control of plots
#
#
if (ncol(nodes) == 7) {
  print("si hay valores de expresi?n")
  V(link_tbl)$Value = nodes$Exp
  #
  rbPal <- colorRampPalette(c('green','black','red'))
  legend_image <- as.raster(matrix(rbPal(100), ncol=1))
  xcol <- rbPal(100)[as.numeric(cut(nodes$Exp,breaks = 100))]
  xcol=xcol[!is.na(xcol)]
  colors.fold.change=c(many.colors[nodes$num[1:total.terms]],xcol)
  #
  min.value.exp=round(min(nodes$Exp[!is.na(nodes$Exp)]),digits=0)
  max.value.exp=round(max(nodes$Exp[!is.na(nodes$Exp)]),digits=0)
  #
  #
  if (total.terms <= 60) {
    print("si hay 60 o menos de 60, se puede hacer complex y simple plots")
    if (total.terms <= 30) {
      print("Son 30 o menos de 30 categor?as, usar configuraci?n de una barra")
      print(paste("el ancho de la figura es:",width.png))
      #
      # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
      # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
      #
      # |||||||||       Bar plot for Network and user
      #
      # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
      # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
      # linear model on geom_text: y = mx + b (y = 0.03x - 0.17)
      # Preparation of data frame (GO terms < 30)
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
      {
        barplot.enrichment = ggplot(for.bar) +
          geom_bar(aes(x = Entry, weight = log.pval),fill=rev(many.colors[1:(total.terms+1)]),width=.7) +
          geom_line(aes(x = c(1:31), y = -log10(0.05),
                        color=-log10(0.05)),
                    color = c(rep("transparent",31-total.terms),
                              rep("black",total.terms)), size = .7 )+
          geom_point(aes(x = c(1:31), y = -log10(0.05),color=-log10(0.05)),
                     position = position_dodge(0.3),
                     color = c(rep("transparent",31-total.terms),
                               rep("black",total.terms)), size = 4)+
          geom_line(aes(x = c(1:31), y = rev(log.fdr),
                        color=log.fdr),
                    color = c(rep("transparent",31-total.terms),
                              rep("blue",total.terms)), size = .7 )+
          geom_point(aes(x = c(1:31), y = rev(log.fdr),color=log.fdr),
                     position = position_dodge(0.3),shape=17,
                     color = c(rep("transparent",31-total.terms),
                               rep("blue",total.terms)), size = 4)+
          #scale_y_continuous(position = "right")+
          coord_flip()+
          scale_x_discrete(limits=rev(for.bar$Entry))+
          geom_text(aes(x=Entry,y=log.pval+.5,label=Freq),size=6.7,fontface ="bold")+ #<---------------------
        #ggtitle("-log(P-value) and\n -log(Corrected P-value)")+
          theme(plot.title = element_text(vjust = 4,hjust=.08,size = 20), # element_text(hjust = .1,size = 12,face="bold"),#posocion y tama?o del titulo sobre el grafico
                axis.title.x =element_blank(),# element_text(colour="black",size=27,face="bold",hjust = 0,vjust = 10),
                axis.title.y = element_blank(),
                axis.text.y = element_text(colour=c("white",rep("black",30)),size=19,
                                           face="bold",hjust = 1.3,vjust = .5), #element_text(size = 8,color="black"),#tama?o y color de texto de eje x
                axis.text.x = element_text(colour="black",size=19,face="bold",vjust = 0),
                axis.ticks.length =unit(.2, "cm"), #reduce la longitud de las marcas en la linea de escala de los ejes x y
                axis.ticks.y.left = element_line(colour = "black",size=3,unit(0, "cm")),
                axis.ticks.x.top = element_line(colour = "black",size=3,unit(3, "cm")),
                axis.line.y = element_blank(), axis.title.x.top =  element_blank(), 
                axis.ticks = element_line(color = "black",size = 1),
                plot.margin = margin(1.5,1,0,1, "cm"))+
          scale_y_continuous(position = "right",
                             limit = c(0,round(max(for.bar$log.pval[1:total.terms]+1.5))),
                             breaks=seq(0,round(max(for.bar$log.pval[1:total.terms]+1.5),digits = 0),1))+
          labs(title = expression(bold("ASPECT")),
               subtitle = expression(bold(" -log10("*italic("P")*"-value & adj."*italic("P")*"-value)")))+
          theme(plot.title = element_text(vjust = 7,hjust=-.3,size = 30),
                plot.subtitle = element_text(vjust = 3,hjust=.08,size = 20))
        {
          barplot.enrichment.user = barplot.enrichment + 
            labs(caption = paste("GOmics:",sep=" ",format(Sys.time(),"%a(%d)/%b/%Y   %X")))
          png(file="zxcvbnmPlot_Bar.png",
              width = 15,
              height = 30,
              units = "cm",
              res=900,bg="white")
          plot(barplot.enrichment.user)
          dev.off()
          }
      }
      #
      # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
      # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
      #
      # |||||||||         Size control
      #
      # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
      # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
      #
      #
      # geom_edge_link function: width
      # linear model: y = mx + b (y = -0.03x + 8.08) #<---------------------
      #
      # geom_node_point function: size
      # linear model: y = mx + b (y = -.012x + 24.32) #<---------------------
      #
      # geom_node_text function: size
      # linear model: y = mx + b (y = -0.01x + 8.36) #<---------------------
      #
      # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
      # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
      #
      # |||||||||     1)    4 simple Plots with Terms as label
      #
      # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
      # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
      #
      layouts=c("nicely","kk","dh")
      #
      for (i in layouts){
        graf1 = ggraph(link_tbl,layout = "igraph", algorithm = i) +
          geom_edge_link(aes(colour = GO),alpha = .7,width = -0.03*vcount(link_tbl)+8.08)+ #<---------------------
        geom_node_point(aes(colour= Value),
                        size=-.12*vcount(link_tbl)+24.32,show.legend = NA, #<---------------------
                        position = "identity")+
          scale_edge_colour_manual(values=colors.fold.change,name = "ASPECT")+
          geom_node_text(aes(label=Label),
                         position ="identity",
                         size = -0.01*vcount(link_tbl)+8.36, #<---------------------
                         repel = T,
                         vjust=3,
                         hjust=.3,
                         fontface = "bold")+
          theme_classic() +
          theme(legend.key.width=unit(1.2,"cm"),
                legend.key.height=unit(1.2,"cm"),
                legend.title=element_text(size=27,face="bold"),
                legend.text=element_text(size=25,face="bold"),
                axis.line=element_blank(),
                axis.text = element_blank(),
                axis.ticks = element_blank(),
                axis.title = element_blank())
        graf1 = graf1 + scale_color_gradient(low = c("green","black"), high = "red") +
          labs(caption = paste("GOmics:",sep=" ",format(Sys.time(),"%a(%d)/%b/%Y   %X")))
        {
          png(file=paste("zxcvbnmPlot_Network_Ter_Layout_",i,".png"),
              width = width.png,  #<---------------------
              height = 15, 
              units = "in",
              res=600,bg="white")
          plot(graf1)
          dev.off()
          }
      }
      #
      # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
      # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
      #
      # |||||||||    2)     4 simple Plots with GO as label
      #
      # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
      # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
      #
      for (i in layouts){
        graf1 = ggraph(link_tbl,layout = "igraph", algorithm = i) +
          geom_edge_link(aes(colour = GO),alpha = .7,width = -0.03*vcount(link_tbl)+8.08)+ #<---------------------
        geom_node_point(aes(colour= Value),
                        size=-.12*vcount(link_tbl)+24.32,show.legend = NA, #<---------------------
                        position = "identity")+
          scale_edge_colour_manual(values=colors.fold.change,name = "ASPECT")+
          geom_node_text(aes(label=name),
                         position ="identity",
                         size = -0.01*vcount(link_tbl)+8.36, #<---------------------
                         repel = T,
                         vjust=3,
                         hjust=.3,
                         fontface = "bold")+
          theme_classic() +
          theme(legend.key.width=unit(1.2,"cm"),
                legend.key.height=unit(1.2,"cm"),
                legend.title=element_text(size=27,face="bold"),
                legend.text=element_text(size=25,face="bold"),
                axis.line=element_blank(),
                axis.text = element_blank(),
                axis.ticks = element_blank(),
                axis.title = element_blank())
        graf1 = graf1 + scale_color_gradient(low = c("green","black"), high = "red") +
          labs(caption = paste("GOmics:",sep=" ",format(Sys.time(),"%a(%d)/%b/%Y   %X")))
        {
          png(file=paste("zxcvbnmPlot_Network_GO_Layout_",i,".png"),
              width = width.png,  #<--------------------- 
              height = 15, 
              units = "in",
              res=600,bg="white")
          plot(graf1)
          dev.off()
          }
      }
      #
      # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
      # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
      #
      # |||||||||    1)     4 Complex Plots with Term as label
      #
      # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
      # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
      #
      ## Construction of graphics with different designs
      #
      for (i in layouts){
        graf2 = ggraph(link_tbl,layout = "igraph",algorithm = i) +
          geom_edge_link(colour = col.edges.cluster,alpha = .7,
                         width = -0.03*vcount(link_tbl)+8.08)+ #<---------------------
        geom_node_point(aes(colour= Value),
                        size=-.12*vcount(link_tbl)+24.32, #<---------------------
                        position = "identity")+
          geom_node_text(aes(label=Label),
                         position ="identity",
                         size = -0.01*vcount(link_tbl)+8.36, #<---------------------
                         repel = T,
                         vjust=3,
                         hjust=.3,
                         fontface = "bold")+
          theme_classic()
        graf2 = graf2 +  theme(legend.key.width=unit(1.2,"cm"),
                               legend.position = c(1.05, .91),
                               legend.key.height=unit(1,"cm"),
                               legend.title=element_text(size=24,face="bold"),
                               legend.text=element_text(size=22),
                               axis.line=element_blank(),
                               axis.text = element_blank(),
                               axis.ticks = element_blank(),
                               axis.title = element_blank())
        graf2 = graf2 + scale_color_gradient(low = c("green","black"), high = "red") +
          labs(caption = paste("GOmics:",sep=" ",format(Sys.time(),"%a(%d)/%b/%Y   %X")))
        {
          # combine plots
          plot.with.inset = ggdraw() +
            draw_plot(graf2,  width = .72, height = 1) +
            #                         x , y , width , height
            draw_plot(barplot.enrichment, .71,.01, .25, .825)
          ## Save plot
          png(file=paste("zxcvbnmPlot_Network_Bar_Term_Layout_",i,".png"),
              width = 25,
              height = 15,
              units = "in",
              res=600,bg="white")
          plot(plot.with.inset)
          dev.off()
          }
      }
      #
      # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
      # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
      #
      # |||||||||     2)    4 Complex Plots with GO as label
      #
      # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
      # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
      #
      ## Construction of graphics with different designs
      #
      for (i in layouts){
        graf2 = ggraph(link_tbl,layout = "igraph",algorithm = i) +
          geom_edge_link(colour = col.edges.cluster,alpha = .7,
                         width = -0.03*vcount(link_tbl)+8.08)+ #<---------------------
        geom_node_point(aes(colour= Value),
                        size=-.12*vcount(link_tbl)+24.32, #<---------------------
                        position = "identity")+
          geom_node_text(aes(label= name),
                         position ="identity",
                         size = -0.01*vcount(link_tbl)+8.36, #<---------------------
                         repel = T,
                         vjust=3,
                         hjust=.3,
                         fontface = "bold")+
          theme_classic()
        graf2 = graf2 +  theme(legend.key.width=unit(1.2,"cm"),
                               legend.position = c(1.05, .91),
                               legend.key.height=unit(1,"cm"),
                               legend.title=element_text(size=24,face="bold"),
                               legend.text=element_text(size=22),
                               axis.line=element_blank(),
                               axis.text = element_blank(),
                               axis.ticks = element_blank(),
                               axis.title = element_blank())
        graf2 = graf2 + scale_color_gradient(low = c("green","black"), high = "red") +
          labs(caption = paste("GOmics:",sep=" ",format(Sys.time(),"%a(%d)/%b/%Y   %X")))
        {
          # combine plots
          plot.with.inset = ggdraw() +
            draw_plot(graf2,  width = .72, height = 1) +
            #                         x , y , width , height
            draw_plot(barplot.enrichment, .71,.01, .25, .825)
          ## Save plot
          png(file=paste("zxcvbnmPlot_Network_Bar_GO_Layout_",i,".png"),
              width = 25,
              height = 15,
              units = "in",
              res=600,bg="white")
          plot(plot.with.inset)
          dev.off()
          }
      }
      # >
      #
      # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
      # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
      #
      # |||||||||     Chord plots (with expression values)
      #
      # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
      # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
      #
      if ((vcount(link_tbl) > 1) && (vcount(link_tbl) <= 45)) {
        print("grafico 1 y 2")
        for (i in -1:0) {
          #
          # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
          # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
          #
          # |||||||||     1) Chord plots (with expression values)
          #
          # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
          # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
          #
          go.term = drop_na(nodes[,c("Entry","Term")])
          colnames(go.term)[1] <- "GO"
          
          go.entry.term = merge(links,go.term,by="GO")#[c(3,2)]
          
          term.entry = data.frame("Term"=go.entry.term$Term,"Entry"=go.entry.term$Entry)
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
            chordDiagram(term.entry.mat, order = c(nodes$Term[1:total.terms],only.entrys),
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
                           circos.lines(c(xcenter, xcenter), c(0, uy(.7, "cm")),
                                        col = 1:4)
                         },bg.border =  NA)
            circos.track(track.index = 2, panel.fun = for(si in get.all.sector.index()[1:total.terms]) {
              xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 2)
              ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 2)
              circos.text(mean(xlim), mean(ylim), si, sector.index = si, track.index = 2,
                          facing = "downward",
                          niceFacing = TRUE,
                          adj = c(0,0.5),
                          cex = -0.004*vcount(link_tbl)+1.48, #<---------------------
                          font = 2,
                          col = NA)
              
            }, bg.border =  NA)
            circos.track(track.index = 2, panel.fun = for(si in get.all.sector.index()[start.entry:vcount(link_tbl)]) {
              xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 2)
              ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 2)
              circos.text(mean(xlim), mean(ylim), si, sector.index = si, track.index = 2,
                          facing = "clockwise",
                          niceFacing = TRUE,
                          adj = c(0.15,0.5),
                          cex = -0.005*vcount(link_tbl)+1.35, #<---------------------
                          font = 2,
                          col = "black")
            }, bg.border = NA)
          }
          #
          {
            png(file=paste0("zxcvbnmPlot_Chord_1_",i,".png"),
                width = 30,
                height = 10,
                units = "in",
                res=900,bg="white")
            chord.plot.1()
            rasterImage(legend_image, -1.01, .7, -1.09,.43)
            text(x=-1, y = seq(0.7,.43,l=5),
                 labels = round(seq(max.value.exp,min.value.exp,length.out = 5),digits=0),font = 1,
                 pos = 4, cex=-0.005*vcount(link_tbl)+1.35) #<---------------------
            text(x=-.93,y= .75, "Value", adj = c(0,0) ,font = 2, pos = 3,
                 cex = -0.004*vcount(link_tbl)+1.48+.2) #<---------------------
            text(x=.4, y=.9, "ASPECT", adj = c(0,0), 
                 cex = -0.004*vcount(link_tbl)+1.48+.2,font = 2) #<---------------------
            dev.off()
          }
          #
          # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
          # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
          #
          # |||||||||     2) Chord plots (with expression values)
          #
          # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
          # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
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
                           circos.lines(c(xcenter, xcenter), c(0, uy(.5, "cm")),
                                        col = 1:4)
                         },bg.border =  NA)
            circos.track(track.index = 2, panel.fun = for(si in get.all.sector.index()[1:total.terms]) {
              xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 2)
              ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 2)
              circos.text(mean(xlim), mean(ylim), si, sector.index = si, track.index = 2,
                          facing = "clockwise",
                          niceFacing = TRUE,
                          adj = c(0.15,0.5),
                          cex = -0.004*vcount(link_tbl)+1.48, #<---------------------
                          font = 2,
                          col = "white")
              
            }, bg.border =  NA)
            circos.track(track.index = 2, panel.fun = for(si in get.all.sector.index()[start.entry:vcount(link_tbl)]) {
              xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 2)
              ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 2)
              circos.text(mean(xlim), mean(ylim), si, sector.index = si, track.index = 2,
                          facing = "clockwise",
                          niceFacing = TRUE,
                          adj = c(0.15,0.5),
                          cex = -0.005*vcount(link_tbl)+1.35, #<---------------------
                          font = 2,
                          col = "black")
            }, bg.border = NA)
          }
          #
          {
            png(file=paste0("zxcvbnmPlot_Chord_2_",i,".png"),
                width = 22.5,
                height = 10,
                units = "in",
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
            rasterImage(legend_image, -.93, .7, -1.01,.43)
            text(x=-.92, y = seq(0.7,.43,l=5),
                 labels = round(seq(max.value.exp,min.value.exp,length.out = 5),digits=0),font = 1,
                 pos = 4, cex=-0.005*vcount(link_tbl)+1.35) #<---------------------
            text(x=-.93,y= .75, "Value", adj = c(0,0) ,font = 2, pos = 3,
                 cex = -0.004*vcount(link_tbl)+1.48+.2) #<---------------------
            text(x=.65, y=.75, "ASPECT", adj = c(0,0), 
                 cex = -0.004*vcount(link_tbl)+1.48+.2,font = 2) #<---------------------
            legend(x=.65,y=.7, 
                   legend=unique(paste0(nodes$Term[1:total.terms]," (",nodes$Freq[1:total.terms],")")),
                   cex=-0.004*vcount(link_tbl)+1.48,bg="transparent", #<---------------------
                   border = NA,pt.bg=many.colors,bty = "n",
                   x.intersp=1.3,y.intersp=-0.003*vcount(link_tbl)+1.21, #<---------------------
                   col=NA,text.font = 2,pch=22,
                   pt.cex=-0.005*vcount(link_tbl)+2.35)
            dev.off()
          }
          # >
        }
      } else if ((vcount(link_tbl) >= 46) && (vcount(link_tbl) <= 60)) {
        print("grafico 2")
        for (i in -1:0) {
          #
          # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
          # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
          #
          # |||||||||     2) Chord plots (with expression values)
          #
          # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
          # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
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
                         #    c(% del datio del c?rculo, alto del track)
                         annotationTrackHeight = c(0.08, .01),
                         preAllocateTracks = 2)
            circos.track(track.index = 2,ylim = c(0, 1), track.height = uh(5, "cm"),
                         track.margin = c(5, uh(5, "cm")),
                         panel.fun = function(x, y) {
                           xcenter = get.cell.meta.data("xcenter")
                           circos.lines(c(xcenter, xcenter), c(0, uy(.5, "cm")),
                                        col = 1:4)
                         },bg.border =  NA)
            circos.track(track.index = 2, panel.fun = for(si in get.all.sector.index()[1:total.terms]) {
              xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 2)
              ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 2)
              circos.text(mean(xlim), mean(ylim), si, sector.index = si, track.index = 2,
                          facing = "clockwise",
                          niceFacing = TRUE,
                          adj = c(0.15,0.5),
                          cex = -0.004*vcount(link_tbl)+1.48, #<---------------------
                          font = 2,
                          col = "white")
              
            }, bg.border =  NA)
            circos.track(track.index = 2, panel.fun = for(si in get.all.sector.index()[start.entry:vcount(link_tbl)]) {
              xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 2)
              ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 2)
              circos.text(mean(xlim), mean(ylim), si, sector.index = si, track.index = 2,
                          facing = "clockwise",
                          niceFacing = TRUE,
                          adj = c(0.15,0.5),
                          cex = -0.005*vcount(link_tbl)+1.35, #<---------------------
                          font = 2,
                          col = "black")
            }, bg.border = NA)
          }
          #
          {
            png(file=paste0("zxcvbnmPlot_Chord_2_",i,".png"),
                width = 22.5,
                height = 10,
                units = "in",
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
            rasterImage(legend_image, -.93, .7, -1.01,.43)
            text(x=-.92, y = seq(0.7,.43,l=5),
                 labels = round(seq(max.value.exp,min.value.exp,length.out = 5),digits=0),font = 1,
                 pos = 4, cex=-0.005*vcount(link_tbl)+1.35) #<---------------------
            text(x=-.93,y= .75, "Value", adj = c(0,0) ,font = 2, pos = 3,
                 cex = -0.004*vcount(link_tbl)+1.48+.2) #<---------------------
            text(x=.65, y=.75, "ASPECT", adj = c(0,0), 
                 cex = -0.004*vcount(link_tbl)+1.48+.2,font = 2) #<---------------------
            legend(x=.65,y=.7, 
                   legend=unique(paste0(nodes$Term[1:total.terms]," (",nodes$Freq[1:total.terms],")")),
                   cex=-0.004*vcount(link_tbl)+1.48,bg="transparent", #<---------------------
                   border = NA,pt.bg=many.colors,bty = "n",
                   x.intersp=1.3,y.intersp=-0.003*vcount(link_tbl)+1.21, #<---------------------
                   col=NA,text.font = 2,pch=22,
                   pt.cex=-0.005*vcount(link_tbl)+2.35)
            dev.off()
          }
          # >
        }
      } else if (vcount(link_tbl) >= 61) {
        print("grafico 3")
        for (i in -1:0) {
          #
          # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
          # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
          #
          # |||||||||     3) Chord plots (with expression values)
          #
          # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
          # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
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
          circos.par(start.degree = 88, clock.wise = T)
          circos.par(track.margin = c(0, 0)) # 1
          #
          chord.plot.3 = function () {
            chordDiagram(go.entry.mat, order = nodes$Entry, # c(nodes$Term[1:total.terms],only.entrys)
                         grid.col = colors.fold.change,
                         transparency = i,
                         directional = -1,
                         annotationTrack = "grid",
                         #    c(% del datio del c?rculo, alto del track)
                         annotationTrackHeight = c(0.08, .01),
                         preAllocateTracks = 2)
            circos.track(track.index = 2,ylim = c(0, 1), track.height = uh(5, "cm"),
                         track.margin = c(5, uh(5, "cm")),
                         panel.fun = function(x, y) {
                           xcenter = get.cell.meta.data("xcenter")
                           circos.lines(c(xcenter, xcenter), c(0, uy(.5, "cm")),
                                        col = 1:4)
                         },bg.border =  NA)
            circos.track(track.index = 2, panel.fun = for(si in get.all.sector.index()[1:total.terms]) {
              xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 2)
              ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 2)
              circos.text(mean(xlim), mean(ylim), si, sector.index = si, track.index = 2,
                          facing = "clockwise",
                          niceFacing = TRUE,
                          adj = c(0.15,0.5),
                          cex = -0.004*vcount(link_tbl)+1.48, #<---------------------
                          font = 2,
                          col = "white")
              
            }, bg.border =  NA)
            circos.track(track.index = 2, panel.fun = for(si in get.all.sector.index()[start.entry:vcount(link_tbl)]) {
              xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 2)
              ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 2)
              circos.text(mean(xlim), mean(ylim), si, sector.index = si, track.index = 2,
                          facing = "clockwise",
                          niceFacing = TRUE,
                          adj = c(0.15,0.5),
                          cex = -0.005*vcount(link_tbl)+1.35, #<---------------------
                          font = 2,
                          col = "black")
            }, bg.border = NA)
          }
          #
          {
            png(file=paste0("zxcvbnmPlot_Chord_3_",i,".png"),
                width = 22.5,
                height = 10,
                units = "in",
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
            chord.plot.3()
            upViewport()
            pushViewport(viewport(x = circle_size, y = 0.5, width = grobWidth(lgd_list_vertical),
                                  height = grobHeight(lgd_list_vertical), just = c("left", "center")))
            grid.draw(lgd_list_vertical)
            upViewport()
            rasterImage(legend_image, -.93, .7, -1.01,.43)
            text(x=-.92, y = seq(0.7,.43,l=5),
                 labels = round(seq(max.value.exp,min.value.exp,length.out = 5),digits=0),font = 1,
                 pos = 4, cex=-0.005*vcount(link_tbl)+1.35) #<---------------------
            text(x=-.93,y= .75, "Value", adj = c(0,0) ,font = 2, pos = 3,
                 cex = -0.004*vcount(link_tbl)+1.48+.2) #<---------------------
            text(x=.65, y=.75, "ASPECT", adj = c(0,0), 
                 cex = -0.004*vcount(link_tbl)+1.48+.2,font = 2) #<---------------------
            legend(x=.65,y=.7, 
                   legend=unique(paste0(nodes$Term[1:total.terms]," (",nodes$Freq[1:total.terms],")")),
                   cex=-0.004*vcount(link_tbl)+1.48,bg="transparent", #<---------------------
                   border = NA,pt.bg=many.colors,bty = "n",
                   x.intersp=1.3,y.intersp=-0.003*vcount(link_tbl)+1.21, #<---------------------
                   col=NA,text.font = 2,pch=22,
                   pt.cex=-0.005*vcount(link_tbl)+2.35)
            dev.off()
          }
          # >
        }
      } else {
        print("no tiene que imprimirse este mensaje")
      }
      # >
      #
      # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
      # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
      #
      # |||||||||     Upset Plot 
      #
      # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
      # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
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
          set.name.size=1.35
          #print(gnp.node.size)
        } else if ((length(unique(links$Entry)) >= 9.1) && (length(unique(links$Entry)) <= 13)) {
          set.name.size=1.2
          #print(gnp.node.size)
        } else if (length(unique(links$Entry)) >= 13.1) {
          set.name.size=0.9
          #print(gnp.node.size)
        }
        #
        {png(file="zxcvbnmPlot_UpSetR.png",width = 10,height = 7,units = "in",res=700,bg="white")
          upset(ups,sets=nodes$Entry[1:total.terms],
                sets.bar.color = "darkgoldenrod2",order.by ="freq",empty.intersections = NULL,
                point.size=2,mainbar.y.label="poiuytrewq",sets.x.label = "",
                main.bar.color="black",matrix.color="black",shade.color="wheat3",
                line.size=0.5,show.numbers = "yes",group.by = "degree",
                matrix.dot.alpha = 0.5,mb.ratio = c(0.7, 0.3),
                #c(intersection size title, intersection size tick labels, set size title, set size tick labels, set names, numbers above bars)
                text.scale = c(1.5,1.5,0,1.5,set.name.size,1.5), 
                scale.sets = "identity",scale.intersections = "identity",
                shade.alpha = 0.35,set_size.angles = 0,
                attribute.plots = list(gridrows = 45, 
                                       plots = list(list(plot = myplot2,
                                                         x="term",y = "log.Pval", 
                                                         queries = F)),ncols=2))
          dev.off()}
      }
      # >
    } else {
      print("son 31 o m?s de 31 y menos de 60, usar la configuraci?n de dos barras")
      print(paste("el ancho de la figura es:",width.png))
      # aqu? usar la nueva configuraci?n para agregar dos barras
      # controlar el ancho de la figura para esta figura don dos barras
      # First part of the data
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
                            rep("black",nrow(nodes.1_30))), size = .7 )+
        geom_point(aes(x = c(1:31), y = -log10(0.05),color=-log10(0.05)),
                   position = position_dodge(0.3),
                   color = c(rep("transparent",31-nrow(nodes.1_30)),
                             rep("black",nrow(nodes.1_30))), size = 4)+
        geom_line(aes(x = c(1:31), y = rev(log.fdr),
                      color=log.fdr),
                  color = c(rep("transparent",31-nrow(nodes.1_30)),
                            rep("blue",nrow(nodes.1_30))), size = .7 )+
        geom_point(aes(x = c(1:31), y = rev(log.fdr),color=log.fdr),
                   position = position_dodge(0.3),shape=17,
                   color = c(rep("transparent",31-nrow(nodes.1_30)),
                             rep("blue",nrow(nodes.1_30))), size = 4)+
        #scale_y_continuous(position = "right")+
        coord_flip()+
        scale_x_discrete(limits=rev(for.bar.1_30$Entry))+
        geom_text(aes(x=Entry,y=log.pval+.5,label=Freq),size=6.7,fontface ="bold")+ #<---------------------
      #ggtitle("-log(P-value) and\n -log(Corrected P-value)")+
        theme(plot.title = element_text(vjust = 4,hjust=.08,size = 20), # element_text(hjust = .1,size = 12,face="bold"),#posocion y tama?o del titulo sobre el grafico
              axis.title.x =element_blank(),# element_text(colour="black",size=27,face="bold",hjust = 0,vjust = 10),
              axis.title.y = element_blank(),
              axis.text.y = element_text(colour=c("white",rep("black",30)),size=19,
                                         face="bold",hjust = 1.3,vjust = .5), #element_text(size = 8,color="black"),#tama?o y color de texto de eje x
              axis.text.x = element_text(colour="black",size=19,face="bold",vjust = 0),
              axis.ticks.length =unit(.2, "cm"), #reduce la longitud de las marcas en la linea de escala de los ejes x y
              axis.ticks.y.left = element_line(colour = "black",size=3,unit(0, "cm")),
              axis.ticks.x.top = element_line(colour = "black",size=3,unit(3, "cm")),
              axis.line.y = element_blank(), axis.title.x.top =  element_blank(), 
              axis.ticks = element_line(color = "black",size = 1),
              plot.margin = margin(1.5,1,0,1, "cm"))+
        scale_y_continuous(position = "right",
                           limit = c(0,round(max(for.bar.1_30$log.pval[1:nrow(nodes.1_30)]+1.5))),
                           breaks=seq(0,round(max(for.bar.1_30$log.pval[1:nrow(nodes.1_30)]+1.5),digits = 0),1))+
        labs(title = expression(bold("ASPECT")),
             subtitle = expression(bold(" -log10("*italic("P")*"-value & adj."*italic("P")*"-value)")))+
        theme(plot.title = element_text(vjust = 7,hjust=-.3,size = 30),
              plot.subtitle = element_text(vjust = 3,hjust=.08,size = 20))
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
                            rep("black",nrow(nodes.31_final))), size = .7 )+
        geom_point(aes(x = c(1:31), y = -log10(0.05),color=-log10(0.05)),
                   position = position_dodge(0.3),
                   color = c(rep("transparent",31-nrow(nodes.31_final)),
                             rep("black",nrow(nodes.31_final))), size = 4)+
        geom_line(aes(x = c(1:31), y = rev(log.fdr),
                      color=log.fdr),
                  color = c(rep("transparent",31-nrow(nodes.31_final)),
                            rep("blue",nrow(nodes.31_final))), size = .7 )+
        geom_point(aes(x = c(1:31), y = rev(log.fdr),color=log.fdr),
                   position = position_dodge(0.3),shape=17,
                   color = c(rep("transparent",31-nrow(nodes.31_final)),
                             rep("blue",nrow(nodes.31_final))), size = 4)+
        #scale_y_continuous(position = "right")+
        coord_flip()+
        scale_x_discrete(limits=rev(for.bar.31_final$Entry))+
        geom_text(aes(x=Entry,y=log.pval+.5,label=Freq),size=6.7,fontface ="bold")+ #<---------------------
      #ggtitle("-log(P-value) and\n -log(Corrected P-value)")+
        theme(plot.title = element_text(vjust = 4,hjust=.08,size = 20), # element_text(hjust = .1,size = 12,face="bold"),#posocion y tama?o del titulo sobre el grafico
              axis.title.x =element_blank(),# element_text(colour="black",size=27,face="bold",hjust = 0,vjust = 10),
              axis.title.y = element_blank(),
              axis.text.y = element_text(colour=c("white",rep("black",30)),size=19,
                                         face="bold",hjust = 1.3,vjust = .5), #element_text(size = 8,color="black"),#tama?o y color de texto de eje x
              axis.text.x = element_text(colour="black",size=19,face="bold",vjust = 0),
              axis.ticks.length =unit(.2, "cm"), #reduce la longitud de las marcas en la linea de escala de los ejes x y
              axis.ticks.y.left = element_line(colour = "black",size=3,unit(0, "cm")),
              axis.ticks.x.top = element_line(colour = "black",size=3,unit(3, "cm")),
              axis.line.y = element_blank(), axis.title.x.top =  element_blank(), 
              axis.ticks = element_line(color = "black",size = 1),
              plot.margin = margin(1.5,1,0,1, "cm"))+
        scale_y_continuous(position = "right",
                           limit = c(0,round(max(for.bar.1_30$log.pval[1:nrow(nodes.1_30)]+1.5))),
                           breaks=seq(0,round(max(for.bar.1_30$log.pval[1:nrow(nodes.1_30)]+1.5),digits = 0),1))+
        labs(title = expression(bold("")),
             subtitle = expression(bold(" -log10("*italic("P")*"-value & adj."*italic("P")*"-value)")))+
        theme(plot.title = element_text(vjust = 7,hjust=-.3,size = 30),
              plot.subtitle = element_text(vjust = 3,hjust=.08,size = 20))
      #
      { # plots_combined used for network
        barplot.p.val.1_30.user = barplot.p.val.1_30 +
          labs(caption = paste("GOmics:",sep=" ",format(Sys.time(),"%a(%d)/%b/%Y   %X")))
        plots_combined  = ggdraw() +
          draw_plot(barplot.p.val.1_30.user, x = 0, y = 0, width = .5, height = 1) +
          draw_plot(barplot.p.val.31_final, x = .45, y = 0, width = .5, height = 1)
        # Barplot for user
        png(file="zxcvbnmPlot_Bar.png",
            width = 30,
            height = 30,
            units = "cm",
            res=900,bg="white")
        plot(plots_combined)
        dev.off()
        }
      #
    #
    #
    # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    #
    # |||||||||     1)    4 simple Plots with Terms as label
    #
    # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    #
    layouts=c("nicely","kk","dh")
    #
    for (i in layouts){
      graf1 = ggraph(link_tbl,layout = "igraph", algorithm = i) +
        geom_edge_link(aes(colour = GO),alpha = .7,width = -0.03*vcount(link_tbl)+8.08)+ #<---------------------
      geom_node_point(aes(colour= Value),
                      size=-.12*vcount(link_tbl)+24.32,show.legend = NA, #<---------------------
                      position = "identity")+
        scale_edge_colour_manual(values=colors.fold.change,name = "ASPECT")+
        geom_node_text(aes(label=Label),
                       position ="identity",
                       size = -0.01*vcount(link_tbl)+8.36, #<---------------------
                       repel = T,
                       vjust=3,
                       hjust=.3,
                       fontface = "bold")+
        theme_classic() +
        theme(legend.key.width=unit(1.2,"cm"),
              legend.key.height=unit(1.2,"cm"),
              legend.title=element_text(size=27,face="bold"),
              legend.text=element_text(size=25,face="bold"),
              axis.line=element_blank(),
              axis.text = element_blank(),
              axis.ticks = element_blank(),
              axis.title = element_blank())
      graf1 = graf1 + scale_color_gradient(low = c("green","black"), high = "red") +
        labs(caption = paste("GOmics:",sep=" ",format(Sys.time(),"%a(%d)/%b/%Y   %X")))
      {
        ## Save plot
        png(file=paste("zxcvbnmPlot_Network_Term_Layout_",i,".png"),
            width = width.png,
            height = 15,
            units = "in",
            res=600,bg="white")
        plot(graf1)
        dev.off()
        }
    }
    #
    # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    #
    # |||||||||    2)     4 simple Plots with GO as label
    #
    # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    #
    for (i in layouts){
      graf1 = ggraph(link_tbl,layout = "igraph", algorithm = i) +
        geom_edge_link(aes(colour = GO),alpha = .7,width = -0.03*vcount(link_tbl)+8.08)+ #<---------------------
      geom_node_point(aes(colour= Value),
                      size=-.12*vcount(link_tbl)+24.32,show.legend = NA, #<---------------------
                      position = "identity")+
        scale_edge_colour_manual(values=colors.fold.change,name = "ASPECT")+
        geom_node_text(aes(label=name),
                       position ="identity",
                       size = -0.01*vcount(link_tbl)+8.36, #<---------------------
                       repel = T,
                       vjust=3,
                       hjust=.3,
                       fontface = "bold")+
        theme_classic() +
        theme(legend.key.width=unit(1.2,"cm"),
              legend.key.height=unit(1.2,"cm"),
              legend.title=element_text(size=27,face="bold"),
              legend.text=element_text(size=25,face="bold"),
              axis.line=element_blank(),
              axis.text = element_blank(),
              axis.ticks = element_blank(),
              axis.title = element_blank())
      graf1 = graf1 + scale_color_gradient(low = c("green","black"), high = "red") +
        labs(caption = paste("GOmics:",sep=" ",format(Sys.time(),"%a(%d)/%b/%Y   %X")))
      {
        ## Save plot
        png(file=paste("zxcvbnmPlot_Network_GO_Layout_",i,".png"),
            width = width.png,
            height = 15,
            units = "in",
            res=600,bg="white")
        plot(graf1)
        dev.off()
        }
    }
    #
    # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    #
    # |||||||||    1)     4 Complex Plots with Term as label
    #
    # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    #
    ## Construction of graphics with different designs
    #
    for (i in layouts){
      graf2 = ggraph(link_tbl,layout = "igraph",algorithm = i) +
        geom_edge_link(colour = col.edges.cluster,alpha = .7,
                       width = -0.03*vcount(link_tbl)+8.08)+ #<---------------------
      geom_node_point(aes(colour= Value),
                      size=-.12*vcount(link_tbl)+24.32, #<---------------------
                      position = "identity")+
        geom_node_text(aes(label=Label),
                       position ="identity",
                       size = -0.01*vcount(link_tbl)+8.36, #<---------------------
                       repel = T,
                       vjust=3,
                       hjust=.3,
                       fontface = "bold")+
        theme_classic()
      graf2 = graf2 +  theme(legend.key.width=unit(1.2,"cm"),
                             legend.position = c(1.05, .91),
                             legend.key.height=unit(1,"cm"),
                             legend.title=element_text(size=24,face="bold"),
                             legend.text=element_text(size=22),
                             axis.line=element_blank(),
                             axis.text = element_blank(),
                             axis.ticks = element_blank(),
                             axis.title = element_blank())
      graf2 = graf2 + scale_color_gradient(low = c("green","black"), high = "red") +
        labs(caption = paste("GOmics:",sep=" ",format(Sys.time(),"%a(%d)/%b/%Y   %X")))
      {
        plots_combined  = ggdraw() +
          draw_plot(barplot.p.val.1_30, x = 0, y = 0, width = .5, height = 1) +
          draw_plot(barplot.p.val.31_final, x = .45, y = 0, width = .5, height = 1)
        # combine plots
        plot.with.inset = ggdraw() +
          draw_plot(graf2,  width = .54, height = 1) +
          #                         x , y , width , height
          draw_plot(plots_combined, .52,.01, .45, .825) # 2 barplots
        
        ## Save plot
        png(file=paste("zxcvbnmPlot_Network_Bar_Term_Layout_",i,".png"),
            width = 35,
            height = 15,
            units = "in",
            res=600,bg="white")
        plot(plot.with.inset)
        dev.off()
        }
    }
    #
    # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    #
    # |||||||||     2)    4 Complex Plots with GO as label
    #
    # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    #
    ## Construction of graphics with different designs
    #
    for (i in layouts){
      graf2 = ggraph(link_tbl,layout = "igraph",algorithm = i) +
        geom_edge_link(colour = col.edges.cluster,alpha = .7,
                       width = -0.03*vcount(link_tbl)+8.08)+ #<---------------------
      geom_node_point(aes(colour= Value),
                      size=-.12*vcount(link_tbl)+24.32, #<---------------------
                      position = "identity")+
        geom_node_text(aes(label= name),
                       position ="identity",
                       size = -0.01*vcount(link_tbl)+8.36, #<---------------------
                       repel = T,
                       vjust=3,
                       hjust=.3,
                       fontface = "bold")+
        theme_classic()
      graf2 = graf2 +  theme(legend.key.width=unit(1.2,"cm"),
                             legend.position = c(1.05, .91),
                             legend.key.height=unit(1,"cm"),
                             legend.title=element_text(size=24,face="bold"),
                             legend.text=element_text(size=22),
                             axis.line=element_blank(),
                             axis.text = element_blank(),
                             axis.ticks = element_blank(),
                             axis.title = element_blank())
      graf2 = graf2 + scale_color_gradient(low = c("green","black"), high = "red") +
        labs(caption = paste("GOmics:",sep=" ",format(Sys.time(),"%a(%d)/%b/%Y   %X")))
      {
        plots_combined  = ggdraw() +
          draw_plot(barplot.p.val.1_30, x = 0, y = 0, width = .5, height = 1) +
          draw_plot(barplot.p.val.31_final, x = .45, y = 0, width = .5, height = 1)
        # combine plots
        plot.with.inset = ggdraw() +
          draw_plot(graf2,  width = .54, height = 1) +
          #                         x , y , width , height
          draw_plot(plots_combined, .52,.01, .45, .825) # 2 barplots
        
        ## Save plot
        png(file=paste("zxcvbnmPlot_Network_Bar_GO_Layout_",i,".png"),
            width = 35,
            height = 15,
            units = "in",
            res=600,bg="white")
        plot(plot.with.inset)
        dev.off()
        }
    }
    # >
    #
    # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    #
    # |||||||||     Chord plots (with expression values)
    #
    # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    #
    if ((vcount(link_tbl) > 1) && (vcount(link_tbl) <= 45)) {
      print("grafico 1 y 2")
      for (i in -1:0) {
        #
        # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
        # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
        #
        # |||||||||     1) Chord plots (with expression values)
        #
        # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
        # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
        #
        go.term = drop_na(nodes[,c("Entry","Term")])
        colnames(go.term)[1] <- "GO"
        
        go.entry.term = merge(links,go.term,by="GO")#[c(3,2)]
        
        term.entry = data.frame("Term"=go.entry.term$Term,"Entry"=go.entry.term$Entry)
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
          chordDiagram(term.entry.mat, order = c(nodes$Term[1:total.terms],only.entrys),
                       grid.col = colors.fold.change,
                       transparency = 0,
                       directional = i,
                       annotationTrack = "grid",
                       #    c(% del datio del c?rculo, alto del track)
                       annotationTrackHeight = c(0.08, .01),
                       preAllocateTracks = 2)
          circos.track(track.index = 2,ylim = c(0, 1), track.height = uh(5, "cm"),
                       track.margin = c(5, uh(5, "cm")),
                       panel.fun = function(x, y) {
                         xcenter = get.cell.meta.data("xcenter")
                         circos.lines(c(xcenter, xcenter), c(0, uy(.7, "cm")),
                                      col = 1:4)
                       },bg.border =  NA)
          circos.track(track.index = 2, panel.fun = for(si in get.all.sector.index()[1:total.terms]) {
            xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 2)
            ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 2)
            circos.text(mean(xlim), mean(ylim), si, sector.index = si, track.index = 2,
                        facing = "downward",
                        niceFacing = TRUE,
                        adj = c(0,0.5),
                        cex = -0.004*vcount(link_tbl)+1.48, #<---------------------
                        font = 2,
                        col = NA)
            
          }, bg.border =  NA)
          circos.track(track.index = 2, panel.fun = for(si in get.all.sector.index()[start.entry:vcount(link_tbl)]) {
            xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 2)
            ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 2)
            circos.text(mean(xlim), mean(ylim), si, sector.index = si, track.index = 2,
                        facing = "clockwise",
                        niceFacing = TRUE,
                        adj = c(0.15,0.5),
                        cex = -0.005*vcount(link_tbl)+1.35, #<---------------------
                        font = 2,
                        col = "black")
          }, bg.border = NA)
        }
        #
        {
          png(file=paste0("zxcvbnmPlot_Chord_1_",i,".png"),
              width = 30,
              height = 10,
              units = "in",
              res=900,bg="white")
          chord.plot.1()
          rasterImage(legend_image, -1.01, .7, -1.09,.43)
          text(x=-1, y = seq(0.7,.43,l=5),
               labels = round(seq(max.value.exp,min.value.exp,length.out = 5),digits=0),font = 1,
               pos = 4, cex=-0.005*vcount(link_tbl)+1.35) #<---------------------
          text(x=-.93,y= .75, "Value", adj = c(0,0) ,font = 2, pos = 3,
               cex = -0.004*vcount(link_tbl)+1.48+.2) #<---------------------
          text(x=.4, y=.9, "ASPECT", adj = c(0,0), 
               cex = -0.004*vcount(link_tbl)+1.48+.2,font = 2) #<---------------------
          dev.off()
        }
        #
        # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
        # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
        #
        # |||||||||     2) Chord plots (with expression values)
        #
        # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
        # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
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
                       #    c(% del datio del c?rculo, alto del track)
                       annotationTrackHeight = c(0.08, .01),
                       preAllocateTracks = 2)
          circos.track(track.index = 2,ylim = c(0, 1), track.height = uh(5, "cm"),
                       track.margin = c(5, uh(5, "cm")),
                       panel.fun = function(x, y) {
                         xcenter = get.cell.meta.data("xcenter")
                         circos.lines(c(xcenter, xcenter), c(0, uy(.5, "cm")),
                                      col = 1:4)
                       },bg.border =  NA)
          circos.track(track.index = 2, panel.fun = for(si in get.all.sector.index()[1:total.terms]) {
            xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 2)
            ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 2)
            circos.text(mean(xlim), mean(ylim), si, sector.index = si, track.index = 2,
                        facing = "clockwise",
                        niceFacing = TRUE,
                        adj = c(0.15,0.5),
                        cex = -0.004*vcount(link_tbl)+1.48, #<---------------------
                        font = 2,
                        col = "white")
            
          }, bg.border =  NA)
          circos.track(track.index = 2, panel.fun = for(si in get.all.sector.index()[start.entry:vcount(link_tbl)]) {
            xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 2)
            ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 2)
            circos.text(mean(xlim), mean(ylim), si, sector.index = si, track.index = 2,
                        facing = "clockwise",
                        niceFacing = TRUE,
                        adj = c(0.15,0.5),
                        cex = -0.005*vcount(link_tbl)+1.35, #<---------------------
                        font = 2,
                        col = "black")
          }, bg.border = NA)
        }
        #
        {
          png(file=paste0("zxcvbnmPlot_Chord_2_",i,".png"),
              width = 22.5,
              height = 10,
              units = "in",
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
          rasterImage(legend_image, -.93, .7, -1.01,.43)
          text(x=-.92, y = seq(0.7,.43,l=5),
               labels = round(seq(max.value.exp,min.value.exp,length.out = 5),digits=0),font = 1,
               pos = 4, cex=-0.005*vcount(link_tbl)+1.35) #<---------------------
          text(x=-.93,y= .75, "Value", adj = c(0,0) ,font = 2, pos = 3,
               cex = -0.004*vcount(link_tbl)+1.48+.2) #<---------------------
          text(x=.65, y=.75, "ASPECT", adj = c(0,0), 
               cex = -0.004*vcount(link_tbl)+1.48+.2,font = 2) #<---------------------
          legend(x=.65,y=.7, 
                 legend=unique(paste0(nodes$Term[1:total.terms]," (",nodes$Freq[1:total.terms],")")),
                 cex=-0.004*vcount(link_tbl)+1.48,bg="transparent", #<---------------------
                 border = NA,pt.bg=many.colors,bty = "n",
                 x.intersp=1.3,y.intersp=-0.003*vcount(link_tbl)+1.21, #<---------------------
                 col=NA,text.font = 2,pch=22,
                 pt.cex=-0.005*vcount(link_tbl)+2.35)
          dev.off()
        }
        # >
      }
    } else if ((vcount(link_tbl) >= 46) && (vcount(link_tbl) <= 60)) {
      print("grafico 2")
      for (i in -1:0) {
        #
        # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
        # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
        #
        # |||||||||     2) Chord plots (with expression values)
        #
        # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
        # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
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
                       #    c(% del datio del c?rculo, alto del track)
                       annotationTrackHeight = c(0.08, .01),
                       preAllocateTracks = 2)
          circos.track(track.index = 2,ylim = c(0, 1), track.height = uh(5, "cm"),
                       track.margin = c(5, uh(5, "cm")),
                       panel.fun = function(x, y) {
                         xcenter = get.cell.meta.data("xcenter")
                         circos.lines(c(xcenter, xcenter), c(0, uy(.5, "cm")),
                                      col = 1:4)
                       },bg.border =  NA)
          circos.track(track.index = 2, panel.fun = for(si in get.all.sector.index()[1:total.terms]) {
            xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 2)
            ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 2)
            circos.text(mean(xlim), mean(ylim), si, sector.index = si, track.index = 2,
                        facing = "clockwise",
                        niceFacing = TRUE,
                        adj = c(0.15,0.5),
                        cex = -0.004*vcount(link_tbl)+1.48, #<---------------------
                        font = 2,
                        col = "white")
            
          }, bg.border =  NA)
          circos.track(track.index = 2, panel.fun = for(si in get.all.sector.index()[start.entry:vcount(link_tbl)]) {
            xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 2)
            ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 2)
            circos.text(mean(xlim), mean(ylim), si, sector.index = si, track.index = 2,
                        facing = "clockwise",
                        niceFacing = TRUE,
                        adj = c(0.15,0.5),
                        cex = -0.005*vcount(link_tbl)+1.35, #<---------------------
                        font = 2,
                        col = "black")
          }, bg.border = NA)
        }
        #
        {
          png(file=paste0("zxcvbnmPlot_Chord_2_",i,".png"),
              width = 22.5,
              height = 10,
              units = "in",
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
          rasterImage(legend_image, -.93, .7, -1.01,.43)
          text(x=-.92, y = seq(0.7,.43,l=5),
               labels = round(seq(max.value.exp,min.value.exp,length.out = 5),digits=0),font = 1,
               pos = 4, cex=-0.005*vcount(link_tbl)+1.35) #<---------------------
          text(x=-.93,y= .75, "Value", adj = c(0,0) ,font = 2, pos = 3,
               cex = -0.004*vcount(link_tbl)+1.48+.2) #<---------------------
          text(x=.65, y=.75, "ASPECT", adj = c(0,0), 
               cex = -0.004*vcount(link_tbl)+1.48+.2,font = 2) #<---------------------
          legend(x=.65,y=.7, 
                 legend=unique(paste0(nodes$Term[1:total.terms]," (",nodes$Freq[1:total.terms],")")),
                 cex=-0.004*vcount(link_tbl)+1.48,bg="transparent", #<---------------------
                 border = NA,pt.bg=many.colors,bty = "n",
                 x.intersp=1.3,y.intersp=-0.003*vcount(link_tbl)+1.21, #<---------------------
                 col=NA,text.font = 2,pch=22,
                 pt.cex=-0.005*vcount(link_tbl)+2.35)
          dev.off()
        }
        # >
      }
    } else if (vcount(link_tbl) >= 61) {
      print("grafico 3")
      for (i in -1:0) {
        #
        # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
        # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
        #
        # |||||||||     3) Chord plots (with expression values)
        #
        # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
        # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
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
        circos.par(start.degree = 88, clock.wise = T)
        circos.par(track.margin = c(0, 0)) # 1
        #
        chord.plot.3 = function () {
          chordDiagram(go.entry.mat, order = nodes$Entry, # c(nodes$Term[1:total.terms],only.entrys)
                       grid.col = colors.fold.change,
                       transparency = i,
                       directional = -1,
                       annotationTrack = "grid",
                       #    c(% del datio del c?rculo, alto del track)
                       annotationTrackHeight = c(0.08, .01),
                       preAllocateTracks = 2)
          circos.track(track.index = 2,ylim = c(0, 1), track.height = uh(5, "cm"),
                       track.margin = c(5, uh(5, "cm")),
                       panel.fun = function(x, y) {
                         xcenter = get.cell.meta.data("xcenter")
                         circos.lines(c(xcenter, xcenter), c(0, uy(.5, "cm")),
                                      col = 1:4)
                       },bg.border =  NA)
          circos.track(track.index = 2, panel.fun = for(si in get.all.sector.index()[1:total.terms]) {
            xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 2)
            ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 2)
            circos.text(mean(xlim), mean(ylim), si, sector.index = si, track.index = 2,
                        facing = "clockwise",
                        niceFacing = TRUE,
                        adj = c(0.15,0.5),
                        cex = -0.004*vcount(link_tbl)+1.48, #<---------------------
                        font = 2,
                        col = "white")
            
          }, bg.border =  NA)
          circos.track(track.index = 2, panel.fun = for(si in get.all.sector.index()[start.entry:vcount(link_tbl)]) {
            xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 2)
            ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 2)
            circos.text(mean(xlim), mean(ylim), si, sector.index = si, track.index = 2,
                        facing = "clockwise",
                        niceFacing = TRUE,
                        adj = c(0.15,0.5),
                        cex = -0.005*vcount(link_tbl)+1.35, #<---------------------
                        font = 2,
                        col = "black")
          }, bg.border = NA)
        }
        #
        {
          png(file=paste0("zxcvbnmPlot_Chord_3_",i,".png"),
              width = 22.5,
              height = 10,
              units = "in",
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
          chord.plot.3()
          upViewport()
          pushViewport(viewport(x = circle_size, y = 0.5, width = grobWidth(lgd_list_vertical),
                                height = grobHeight(lgd_list_vertical), just = c("left", "center")))
          grid.draw(lgd_list_vertical)
          upViewport()
          rasterImage(legend_image, -.93, .7, -1.01,.43)
          text(x=-.92, y = seq(0.7,.43,l=5),
               labels = round(seq(max.value.exp,min.value.exp,length.out = 5),digits=0),font = 1,
               pos = 4, cex=-0.005*vcount(link_tbl)+1.35) #<---------------------
          text(x=-.93,y= .75, "Value", adj = c(0,0) ,font = 2, pos = 3,
               cex = -0.004*vcount(link_tbl)+1.48+.2) #<---------------------
          text(x=.65, y=.75, "ASPECT", adj = c(0,0), 
               cex = -0.004*vcount(link_tbl)+1.48+.2,font = 2) #<---------------------
          legend(x=.65,y=.7, 
                 legend=unique(paste0(nodes$Term[1:total.terms]," (",nodes$Freq[1:total.terms],")")),
                 cex=-0.004*vcount(link_tbl)+1.48,bg="transparent", #<---------------------
                 border = NA,pt.bg=many.colors,bty = "n",
                 x.intersp=1.3,y.intersp=-0.003*vcount(link_tbl)+1.21, #<---------------------
                 col=NA,text.font = 2,pch=22,
                 pt.cex=-0.005*vcount(link_tbl)+2.35)
          dev.off()
        }
        # >
      }
    } else {
      print("no tiene que imprimirse este mensaje")
    }
    # >
    #
    # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    #
    # |||||||||     Upset Plot (with expression values)
    #
    # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
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
        set.name.size=1.35
        #print(gnp.node.size)
      } else if ((length(unique(links$Entry)) >= 9.1) && (length(unique(links$Entry)) <= 13)) {
        set.name.size=1.2
        #print(gnp.node.size)
      } else if (length(unique(links$Entry)) >= 13.1) {
        set.name.size=0.9
        #print(gnp.node.size)
      }
      #
      {png(file="zxcvbnmPlot_UpSetR.png",width = 10,height = 7,units = "in",res=700,bg="white")
        upset(ups,sets=nodes$Entry[1:total.terms],
              sets.bar.color = "darkgoldenrod2",order.by ="freq",empty.intersections = NULL,
              point.size=2,mainbar.y.label="poiuytrewq",sets.x.label = "",
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
    # >
    }
  } else {
    print("hay 61 o m?s de 61, se puede hacer ?nicamente plots simples")
    print(paste("el ancho de la figura es:",width.png))
    #
    # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    #
    # |||||||||     1)    4 simple Plots with Terms as label
    #
    # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    #
    layouts=c("nicely","kk","dh")
    #
    for (i in layouts){
      graf1 = ggraph(link_tbl,layout = "igraph", algorithm = i) +
        geom_edge_link(aes(colour = GO),alpha = .7,width = -0.03*vcount(link_tbl)+8.08)+ #<---------------------
      geom_node_point(aes(colour= Value),
                      size=-.12*vcount(link_tbl)+24.32,show.legend = NA, #<---------------------
                      position = "identity")+
        scale_edge_colour_manual(values=colors.fold.change,name = "ASPECT")+
        geom_node_text(aes(label=Label),
                       position ="identity",
                       size = -0.01*vcount(link_tbl)+8.36, #<---------------------
                       repel = T,
                       vjust=3,
                       hjust=.3,
                       fontface = "bold")+
        theme_classic() +
        theme(legend.key.width=unit(1.2,"cm"),
              legend.key.height=unit(1.2,"cm"),
              legend.title=element_text(size=27,face="bold"),
              legend.text=element_text(size=25,face="bold"),
              axis.line=element_blank(),
              axis.text = element_blank(),
              axis.ticks = element_blank(),
              axis.title = element_blank())
      graf1 = graf1 + scale_color_gradient(low = c("green","black"), high = "red") +
        labs(caption = paste("GOmics:",sep=" ",format(Sys.time(),"%a(%d)/%b/%Y   %X")))
      {
        png(file=paste("zxcvbnmPlot_Network_Ter_Layout_",i,".png"),
            width = width.png,  #<---------------------
            height = 15, 
            units = "in",
            res=600,bg="white")
        plot(graf1)
        dev.off()
        }
    }
    #
    # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    #
    # |||||||||    2)     4 simple Plots with GO as label
    #
    # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    #
    for (i in layouts){
      graf1 = ggraph(link_tbl,layout = "igraph", algorithm = i) +
        geom_edge_link(aes(colour = GO),alpha = .7,width = -0.03*vcount(link_tbl)+8.08)+ #<---------------------
      geom_node_point(aes(colour= Value),
                      size=-.12*vcount(link_tbl)+24.32,show.legend = NA, #<---------------------
                      position = "identity")+
        scale_edge_colour_manual(values=colors.fold.change,name = "ASPECT")+
        geom_node_text(aes(label=name),
                       position ="identity",
                       size = -0.01*vcount(link_tbl)+8.36, #<---------------------
                       repel = T,
                       vjust=3,
                       hjust=.3,
                       fontface = "bold")+
        theme_classic() +
        theme(legend.key.width=unit(1.2,"cm"),
              legend.key.height=unit(1.2,"cm"),
              legend.title=element_text(size=27,face="bold"),
              legend.text=element_text(size=25,face="bold"),
              axis.line=element_blank(),
              axis.text = element_blank(),
              axis.ticks = element_blank(),
              axis.title = element_blank())
      graf1 = graf1 + scale_color_gradient(low = c("green","black"), high = "red") +
        labs(caption = paste("GOmics:",sep=" ",format(Sys.time(),"%a(%d)/%b/%Y   %X")))
      {
        png(file=paste("zxcvbnmPlot_Network_GO_Layout_",i,".png"),
            width = width.png,  #<--------------------- 
            height = 15, 
            units = "in",
            res=600,bg="white")
        plot(graf1)
        dev.off()
        }
    }
    # >
    #
    # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    #
    # |||||||||     Chord plots (with expression values)
    #
    # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    #
    if ((vcount(link_tbl) > 1) && (vcount(link_tbl) <= 45)) {
      print("grafico 1 y 2")
      for (i in -1:0) {
        #
        # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
        # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
        #
        # |||||||||     1) Chord plots (with expression values)
        #
        # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
        # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
        #
        go.term = drop_na(nodes[,c("Entry","Term")])
        colnames(go.term)[1] <- "GO"
        
        go.entry.term = merge(links,go.term,by="GO")#[c(3,2)]
        
        term.entry = data.frame("Term"=go.entry.term$Term,"Entry"=go.entry.term$Entry)
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
          chordDiagram(term.entry.mat, order = c(nodes$Term[1:total.terms],only.entrys),
                       grid.col = colors.fold.change,
                       transparency = 0,
                       directional = i,
                       annotationTrack = "grid",
                       #    c(% del datio del c?rculo, alto del track)
                       annotationTrackHeight = c(0.08, .01),
                       preAllocateTracks = 2)
          circos.track(track.index = 2,ylim = c(0, 1), track.height = uh(5, "cm"),
                       track.margin = c(5, uh(5, "cm")),
                       panel.fun = function(x, y) {
                         xcenter = get.cell.meta.data("xcenter")
                         circos.lines(c(xcenter, xcenter), c(0, uy(.7, "cm")),
                                      col = 1:4)
                       },bg.border =  NA)
          circos.track(track.index = 2, panel.fun = for(si in get.all.sector.index()[1:total.terms]) {
            xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 2)
            ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 2)
            circos.text(mean(xlim), mean(ylim), si, sector.index = si, track.index = 2,
                        facing = "downward",
                        niceFacing = TRUE,
                        adj = c(0,0.5),
                        cex = -0.004*vcount(link_tbl)+1.48, #<---------------------
                        font = 2,
                        col = NA)
            
          }, bg.border =  NA)
          circos.track(track.index = 2, panel.fun = for(si in get.all.sector.index()[start.entry:vcount(link_tbl)]) {
            xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 2)
            ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 2)
            circos.text(mean(xlim), mean(ylim), si, sector.index = si, track.index = 2,
                        facing = "clockwise",
                        niceFacing = TRUE,
                        adj = c(0.15,0.5),
                        cex = -0.005*vcount(link_tbl)+1.35, #<---------------------
                        font = 2,
                        col = "black")
          }, bg.border = NA)
        }
        #
        {
          png(file=paste0("zxcvbnmPlot_Chord_1_",i,".png"),
              width = 30,
              height = 10,
              units = "in",
              res=900,bg="white")
          chord.plot.1()
          rasterImage(legend_image, -1.01, .7, -1.09,.43)
          text(x=-1, y = seq(0.7,.43,l=5),
               labels = round(seq(max.value.exp,min.value.exp,length.out = 5),digits=0),font = 1,
               pos = 4, cex=-0.005*vcount(link_tbl)+1.35) #<---------------------
          text(x=-.93,y= .75, "Value", adj = c(0,0) ,font = 2, pos = 3,
               cex = -0.004*vcount(link_tbl)+1.48+.2) #<---------------------
          text(x=.4, y=.9, "ASPECT", adj = c(0,0), 
               cex = -0.004*vcount(link_tbl)+1.48+.2,font = 2) #<---------------------
          dev.off()
        }
        #
        # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
        # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
        #
        # |||||||||     2) Chord plots (with expression values)
        #
        # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
        # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
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
                       #    c(% del datio del c?rculo, alto del track)
                       annotationTrackHeight = c(0.08, .01),
                       preAllocateTracks = 2)
          circos.track(track.index = 2,ylim = c(0, 1), track.height = uh(5, "cm"),
                       track.margin = c(5, uh(5, "cm")),
                       panel.fun = function(x, y) {
                         xcenter = get.cell.meta.data("xcenter")
                         circos.lines(c(xcenter, xcenter), c(0, uy(.5, "cm")),
                                      col = 1:4)
                       },bg.border =  NA)
          circos.track(track.index = 2, panel.fun = for(si in get.all.sector.index()[1:total.terms]) {
            xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 2)
            ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 2)
            circos.text(mean(xlim), mean(ylim), si, sector.index = si, track.index = 2,
                        facing = "clockwise",
                        niceFacing = TRUE,
                        adj = c(0.15,0.5),
                        cex = -0.004*vcount(link_tbl)+1.48, #<---------------------
                        font = 2,
                        col = "white")
            
          }, bg.border =  NA)
          circos.track(track.index = 2, panel.fun = for(si in get.all.sector.index()[start.entry:vcount(link_tbl)]) {
            xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 2)
            ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 2)
            circos.text(mean(xlim), mean(ylim), si, sector.index = si, track.index = 2,
                        facing = "clockwise",
                        niceFacing = TRUE,
                        adj = c(0.15,0.5),
                        cex = -0.005*vcount(link_tbl)+1.35, #<---------------------
                        font = 2,
                        col = "black")
          }, bg.border = NA)
        }
        #
        {
          png(file=paste0("zxcvbnmPlot_Chord_2_",i,".png"),
              width = 22.5,
              height = 10,
              units = "in",
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
          rasterImage(legend_image, -.93, .7, -1.01,.43)
          text(x=-.92, y = seq(0.7,.43,l=5),
               labels = round(seq(max.value.exp,min.value.exp,length.out = 5),digits=0),font = 1,
               pos = 4, cex=-0.005*vcount(link_tbl)+1.35) #<---------------------
          text(x=-.93,y= .75, "Value", adj = c(0,0) ,font = 2, pos = 3,
               cex = -0.004*vcount(link_tbl)+1.48+.2) #<---------------------
          text(x=.65, y=.75, "ASPECT", adj = c(0,0), 
               cex = -0.004*vcount(link_tbl)+1.48+.2,font = 2) #<---------------------
          legend(x=.65,y=.7, 
                 legend=unique(paste0(nodes$Term[1:total.terms]," (",nodes$Freq[1:total.terms],")")),
                 cex=-0.004*vcount(link_tbl)+1.48,bg="transparent", #<---------------------
                 border = NA,pt.bg=many.colors,bty = "n",
                 x.intersp=1.3,y.intersp=-0.003*vcount(link_tbl)+1.21, #<---------------------
                 col=NA,text.font = 2,pch=22,
                 pt.cex=-0.005*vcount(link_tbl)+2.35)
          dev.off()
        }
        # >
      }
    } else if ((vcount(link_tbl) >= 46) && (vcount(link_tbl) <= 60)) {
      print("grafico 2")
      for (i in -1:0) {
        #
        # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
        # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
        #
        # |||||||||     2) Chord plots (with expression values)
        #
        # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
        # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
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
                       #    c(% del datio del c?rculo, alto del track)
                       annotationTrackHeight = c(0.08, .01),
                       preAllocateTracks = 2)
          circos.track(track.index = 2,ylim = c(0, 1), track.height = uh(5, "cm"),
                       track.margin = c(5, uh(5, "cm")),
                       panel.fun = function(x, y) {
                         xcenter = get.cell.meta.data("xcenter")
                         circos.lines(c(xcenter, xcenter), c(0, uy(.5, "cm")),
                                      col = 1:4)
                       },bg.border =  NA)
          circos.track(track.index = 2, panel.fun = for(si in get.all.sector.index()[1:total.terms]) {
            xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 2)
            ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 2)
            circos.text(mean(xlim), mean(ylim), si, sector.index = si, track.index = 2,
                        facing = "clockwise",
                        niceFacing = TRUE,
                        adj = c(0.15,0.5),
                        cex = -0.004*vcount(link_tbl)+1.48, #<---------------------
                        font = 2,
                        col = "white")
            
          }, bg.border =  NA)
          circos.track(track.index = 2, panel.fun = for(si in get.all.sector.index()[start.entry:vcount(link_tbl)]) {
            xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 2)
            ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 2)
            circos.text(mean(xlim), mean(ylim), si, sector.index = si, track.index = 2,
                        facing = "clockwise",
                        niceFacing = TRUE,
                        adj = c(0.15,0.5),
                        cex = -0.005*vcount(link_tbl)+1.35, #<---------------------
                        font = 2,
                        col = "black")
          }, bg.border = NA)
        }
        #
        {
          png(file=paste0("zxcvbnmPlot_Chord_2_",i,".png"),
              width = 22.5,
              height = 10,
              units = "in",
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
          rasterImage(legend_image, -.93, .7, -1.01,.43)
          text(x=-.92, y = seq(0.7,.43,l=5),
               labels = round(seq(max.value.exp,min.value.exp,length.out = 5),digits=0),font = 1,
               pos = 4, cex=-0.005*vcount(link_tbl)+1.35) #<---------------------
          text(x=-.93,y= .75, "Value", adj = c(0,0) ,font = 2, pos = 3,
               cex = -0.004*vcount(link_tbl)+1.48+.2) #<---------------------
          text(x=.65, y=.75, "ASPECT", adj = c(0,0), 
               cex = -0.004*vcount(link_tbl)+1.48+.2,font = 2) #<---------------------
          legend(x=.65,y=.7, 
                 legend=unique(paste0(nodes$Term[1:total.terms]," (",nodes$Freq[1:total.terms],")")),
                 cex=-0.004*vcount(link_tbl)+1.48,bg="transparent", #<---------------------
                 border = NA,pt.bg=many.colors,bty = "n",
                 x.intersp=1.3,y.intersp=-0.003*vcount(link_tbl)+1.21, #<---------------------
                 col=NA,text.font = 2,pch=22,
                 pt.cex=-0.005*vcount(link_tbl)+2.35)
          dev.off()
        }
        # >
      }
    } else if (vcount(link_tbl) >= 61) {
      print("grafico 3")
      for (i in -1:0) {
        #
        # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
        # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
        #
        # |||||||||     3) Chord plots (with expression values)
        #
        # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
        # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
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
        circos.par(start.degree = 88, clock.wise = T)
        circos.par(track.margin = c(0, 0)) # 1
        #
        chord.plot.3 = function () {
          chordDiagram(go.entry.mat, order = nodes$Entry, # c(nodes$Term[1:total.terms],only.entrys)
                       grid.col = colors.fold.change,
                       transparency = i,
                       directional = -1,
                       annotationTrack = "grid",
                       #    c(% del datio del c?rculo, alto del track)
                       annotationTrackHeight = c(0.08, .01),
                       preAllocateTracks = 2)
          circos.track(track.index = 2,ylim = c(0, 1), track.height = uh(5, "cm"),
                       track.margin = c(5, uh(5, "cm")),
                       panel.fun = function(x, y) {
                         xcenter = get.cell.meta.data("xcenter")
                         circos.lines(c(xcenter, xcenter), c(0, uy(.5, "cm")),
                                      col = 1:4)
                       },bg.border =  NA)
          circos.track(track.index = 2, panel.fun = for(si in get.all.sector.index()[1:total.terms]) {
            xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 2)
            ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 2)
            circos.text(mean(xlim), mean(ylim), si, sector.index = si, track.index = 2,
                        facing = "clockwise",
                        niceFacing = TRUE,
                        adj = c(0.15,0.5),
                        cex = -0.004*vcount(link_tbl)+1.48, #<---------------------
                        font = 2,
                        col = "white")
            
          }, bg.border =  NA)
          circos.track(track.index = 2, panel.fun = for(si in get.all.sector.index()[start.entry:vcount(link_tbl)]) {
            xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 2)
            ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 2)
            circos.text(mean(xlim), mean(ylim), si, sector.index = si, track.index = 2,
                        facing = "clockwise",
                        niceFacing = TRUE,
                        adj = c(0.15,0.5),
                        cex = -0.005*vcount(link_tbl)+1.35, #<---------------------
                        font = 2,
                        col = "black")
          }, bg.border = NA)
        }
        #
        {
          png(file=paste0("zxcvbnmPlot_Chord_3_",i,".png"),
              width = 22.5,
              height = 10,
              units = "in",
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
          chord.plot.3()
          upViewport()
          pushViewport(viewport(x = circle_size, y = 0.5, width = grobWidth(lgd_list_vertical),
                                height = grobHeight(lgd_list_vertical), just = c("left", "center")))
          grid.draw(lgd_list_vertical)
          upViewport()
          rasterImage(legend_image, -.93, .7, -1.01,.43)
          text(x=-.92, y = seq(0.7,.43,l=5),
               labels = round(seq(max.value.exp,min.value.exp,length.out = 5),digits=0),font = 1,
               pos = 4, cex=-0.005*vcount(link_tbl)+1.35) #<---------------------
          text(x=-.93,y= .75, "Value", adj = c(0,0) ,font = 2, pos = 3,
               cex = -0.004*vcount(link_tbl)+1.48+.2) #<---------------------
          text(x=.65, y=.75, "ASPECT", adj = c(0,0), 
               cex = -0.004*vcount(link_tbl)+1.48+.2,font = 2) #<---------------------
          legend(x=.65,y=.7, 
                 legend=unique(paste0(nodes$Term[1:total.terms]," (",nodes$Freq[1:total.terms],")")),
                 cex=-0.004*vcount(link_tbl)+1.48,bg="transparent", #<---------------------
                 border = NA,pt.bg=many.colors,bty = "n",
                 x.intersp=1.3,y.intersp=-0.003*vcount(link_tbl)+1.21, #<---------------------
                 col=NA,text.font = 2,pch=22,
                 pt.cex=-0.005*vcount(link_tbl)+2.35)
          dev.off()
        }
        # >
      }
    } else {
      print("no tiene que imprimirse este mensaje")
    }
    # >
    #
    # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    #
    # |||||||||     Upset Plot (with expression values)
    #
    # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
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
        set.name.size=1.35
        #print(gnp.node.size)
      } else if ((length(unique(links$Entry)) >= 9.1) && (length(unique(links$Entry)) <= 13)) {
        set.name.size=1.2
        #print(gnp.node.size)
      } else if (length(unique(links$Entry)) >= 13.1) {
        set.name.size=0.9
        #print(gnp.node.size)
      }
      #
      {png(file="zxcvbnmPlot_UpSetR.png",width = 10,height = 7,units = "in",res=700,bg="white")
        upset(ups,sets=nodes$Entry[1:total.terms],
              sets.bar.color = "darkgoldenrod2",order.by ="freq",empty.intersections = NULL,
              point.size=2,mainbar.y.label="poiuytrewq",sets.x.label = "",
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
    # >
    # >
  }
# =================================================================================================
# =================================================================================================
# =================================================================================================
# =================================================================================================
# =================================================================================================
# =================================================================================================
} else {
  print("no hay valores de expresi?n")
  prot.col=rep("blue",1000)#9D9D9D
  col.set=c(many.colors[nodes$num[1:total.terms]],
            prot.col[nodes$num[start.entry:total.nodes]])
  #
  #
  if (total.terms <= 60) {
    print("si hay 60 o menos de 60, se puede hacer complex y simple plots")
    if (total.terms <= 30) {
      print("Son 30 o menos de 30 categor?as, usar configuraci?n de una barra")
      print(paste("el ancho de la figura es:",width.png))
      #
      # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
      # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
      #
      # |||||||||       Bar plot for Network and user
      #
      # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
      # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
      # linear model on geom_text: y = mx + b (y = 0.03x - 0.17)
      # Preparation of data frame (GO terms < 30)
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
      {
        barplot.enrichment = ggplot(for.bar) +
          geom_bar(aes(x = Entry, weight = log.pval),fill=rev(many.colors[1:(total.terms+1)]),width=.7) +
          geom_line(aes(x = c(1:31), y = -log10(0.05),
                        color=-log10(0.05)),
                    color = c(rep("transparent",31-total.terms),
                              rep("black",total.terms)), size = .7 )+
          geom_point(aes(x = c(1:31), y = -log10(0.05),color=-log10(0.05)),
                     position = position_dodge(0.3),
                     color = c(rep("transparent",31-total.terms),
                               rep("black",total.terms)), size = 4)+
          geom_line(aes(x = c(1:31), y = rev(log.fdr),
                        color=log.fdr),
                    color = c(rep("transparent",31-total.terms),
                              rep("blue",total.terms)), size = .7 )+
          geom_point(aes(x = c(1:31), y = rev(log.fdr),color=log.fdr),
                     position = position_dodge(0.3),shape=17,
                     color = c(rep("transparent",31-total.terms),
                               rep("blue",total.terms)), size = 4)+
          #scale_y_continuous(position = "right")+
          coord_flip()+
          scale_x_discrete(limits=rev(for.bar$Entry))+
          geom_text(aes(x=Entry,y=log.pval+.5,label=Freq),size=6.7,fontface ="bold")+ #<---------------------
        #ggtitle("-log(P-value) and\n -log(Corrected P-value)")+
          theme(plot.title = element_text(vjust = 4,hjust=.08,size = 20), # element_text(hjust = .1,size = 12,face="bold"),#posocion y tama?o del titulo sobre el grafico
                axis.title.x =element_blank(),# element_text(colour="black",size=27,face="bold",hjust = 0,vjust = 10),
                axis.title.y = element_blank(),
                axis.text.y = element_text(colour=c("white",rep("black",30)),size=19,
                                           face="bold",hjust = 1.3,vjust = .5), #element_text(size = 8,color="black"),#tama?o y color de texto de eje x
                axis.text.x = element_text(colour="black",size=19,face="bold",vjust = 0),
                axis.ticks.length =unit(.2, "cm"), #reduce la longitud de las marcas en la linea de escala de los ejes x y
                axis.ticks.y.left = element_line(colour = "black",size=3,unit(0, "cm")),
                axis.ticks.x.top = element_line(colour = "black",size=3,unit(3, "cm")),
                axis.line.y = element_blank(), axis.title.x.top =  element_blank(), 
                axis.ticks = element_line(color = "black",size = 1),
                plot.margin = margin(1.5,1,0,1, "cm"))+
          scale_y_continuous(position = "right",
                             limit = c(0,round(max(for.bar$log.pval[1:total.terms]+1.5))),
                             breaks=seq(0,round(max(for.bar$log.pval[1:total.terms]+1.5),digits = 0),1))+
          labs(title = expression(bold("ASPECT")),
               subtitle = expression(bold(" -log10("*italic("P")*"-value & adj."*italic("P")*"-value)")))+
          theme(plot.title = element_text(vjust = 7,hjust=-.3,size = 30),
                plot.subtitle = element_text(vjust = 3,hjust=.08,size = 20))
        {
          barplot.enrichment.user = barplot.enrichment + 
            labs(caption = paste("GOmics:",sep=" ",format(Sys.time(),"%a(%d)/%b/%Y   %X")))
          png(file="zxcvbnmPlot_Bar.png",
              width = 15,
              height = 30,
              units = "cm",
              res=900,bg="white")
          plot(barplot.enrichment.user)
          dev.off()
          }
      }
      #
      # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
      # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
      #
      # |||||||||         Size control
      #
      # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
      # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
      #
      #
      # geom_edge_link function: width
      # linear model: y = mx + b (y = -0.03x + 8.08) #<---------------------
      #
      # geom_node_point function: size
      # linear model: y = mx + b (y = -.012x + 24.32) #<---------------------
      #
      # geom_node_text function: size
      # linear model: y = mx + b (y = -0.01x + 8.36) #<---------------------
      #
      # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
      # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
      #
      # |||||||||     1)    4 simple Plots with Terms as label
      #
      # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
      # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
      #
      layouts=c("nicely","kk","dh")
      #
      for (i in layouts){
        graf1 = ggraph(link_tbl,layout = "igraph", algorithm = i) +
          geom_edge_link(aes(colour = GO),alpha = .7,width = -0.03*vcount(link_tbl)+8.08)+ #<---------------------
        geom_node_point(colour= col.set,shape=c(rep(16,total.terms),rep(20,length(only.entrys))),
                        size=-.12*vcount(link_tbl)+24.32,show.legend = NA, #<---------------------
                        position = "identity")+
          scale_edge_colour_manual(values=col.set,name = "ASPECT")+
          #scale_edge_colour_manual(values=colors.fold.change)+
          geom_node_text(aes(label=Label),
                         position ="identity",
                         size = -0.01*vcount(link_tbl)+8.36, #<---------------------
                         repel = T,
                         vjust=3,
                         hjust=.3,
                         fontface = "bold")+
          theme_classic() +
          theme(legend.key.width=unit(1.2,"cm"),
                legend.key.height=unit(1.2,"cm"),
                legend.title=element_text(size=27,face="bold"),
                legend.text=element_text(size=25,face="bold"),
                axis.line=element_blank(),
                axis.text = element_blank(),
                axis.ticks = element_blank(),
                axis.title = element_blank())
        graf1 = graf1 + 
          labs(caption = paste("GOmics:",sep=" ",format(Sys.time(),"%a(%d)/%b/%Y   %X")))
        {
          png(file=paste("zxcvbnmPlot_Network_Ter_Layout_",i,".png"),
              width = width.png,  #<---------------------
              height = 15, 
              units = "in",
              res=600,bg="white")
          plot(graf1)
          dev.off()
          }
      }
      #
      # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
      # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
      #
      # |||||||||    2)     4 simple Plots with GO as label
      #
      # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
      # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
      #
      for (i in layouts){
        graf1 = ggraph(link_tbl,layout = "igraph", algorithm = i) +
          geom_edge_link(aes(colour = GO),alpha = .7,width = -0.03*vcount(link_tbl)+8.08)+ #<---------------------
        geom_node_point(colour= col.set,shape=c(rep(16,total.terms),rep(20,length(only.entrys))),
                        size=-.12*vcount(link_tbl)+24.32,show.legend = NA, #<---------------------
                        position = "identity")+
          scale_edge_colour_manual(values=col.set,name = "ASPECT")+
          #scale_edge_colour_manual(values=colors.fold.change)+
          geom_node_text(aes(label=name),
                         position ="identity",
                         size = -0.01*vcount(link_tbl)+8.36, #<---------------------
                         repel = T,
                         vjust=3,
                         hjust=.3,
                         fontface = "bold")+
          theme_classic() +
          theme(legend.key.width=unit(1.2,"cm"),
                legend.key.height=unit(1.2,"cm"),
                legend.title=element_text(size=27,face="bold"),
                legend.text=element_text(size=25,face="bold"),
                axis.line=element_blank(),
                axis.text = element_blank(),
                axis.ticks = element_blank(),
                axis.title = element_blank())
        graf1 = graf1 +
          labs(caption = paste("GOmics:",sep=" ",format(Sys.time(),"%a(%d)/%b/%Y   %X")))
        {
          png(file=paste("zxcvbnmPlot_Network_GO_Layout_",i,".png"),
              width = width.png,  #<--------------------- 
              height = 15, 
              units = "in",
              res=600,bg="white")
          plot(graf1)
          dev.off()
          }
      }
      #
      # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
      # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
      #
      # |||||||||    1)     Complex Plots with Term as label
      #
      # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
      # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
      #
      ## Construction of graphics with different designs
      #
      for (i in layouts){
        graf2 = ggraph(link_tbl,layout = "igraph",algorithm = i) +
          geom_edge_link(colour = col.edges.cluster,alpha = .7,
                         width = -0.03*vcount(link_tbl)+8.08)+ #<---------------------
        geom_node_point(colour= col.set,shape=c(rep(16,total.terms),rep(20,length(only.entrys))),
                        size=-.12*vcount(link_tbl)+24.32, #<---------------------
                        position = "identity")+
          geom_node_text(aes(label=Label),
                         position ="identity",
                         size = -0.01*vcount(link_tbl)+8.36, #<---------------------
                         repel = T,
                         vjust=3,
                         hjust=.3,
                         fontface = "bold")+
          theme_classic()
        graf2 = graf2 +  theme(legend.key.width=unit(1.2,"cm"),
                               legend.position = c(1.05, .91),
                               legend.key.height=unit(1,"cm"),
                               legend.title=element_text(size=24,face="bold"),
                               legend.text=element_text(size=22),
                               axis.line=element_blank(),
                               axis.text = element_blank(),
                               axis.ticks = element_blank(),
                               axis.title = element_blank())
        graf2 = graf2 + scale_color_gradient(low = c("green","black"), high = "red") +
          labs(caption = paste("GOmics:",sep=" ",format(Sys.time(),"%a(%d)/%b/%Y   %X")))
        {
          # combine plots
          plot.with.inset = ggdraw() +
            draw_plot(graf2,  width = .72, height = 1) +
            #                         x , y , width , height
            draw_plot(barplot.enrichment, .71,.01, .25, .825)
          ## Save plot
          png(file=paste("zxcvbnmPlot_Network_Bar_Term_Layout_",i,".png"),
              width = 25,
              height = 15,
              units = "in",
              res=600,bg="white")
          plot(plot.with.inset)
          dev.off()
          }
      }
      #
      # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
      # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
      #
      # |||||||||     2)    4 Complex Plots with GO as label
      #
      # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
      # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
      #
      ## Construction of graphics with different designs
      #
      for (i in layouts){
        graf2 = ggraph(link_tbl,layout = "igraph",algorithm = i) +
          geom_edge_link(colour = col.edges.cluster,alpha = .7,
                         width = -0.03*vcount(link_tbl)+8.08)+ #<---------------------
        geom_node_point(colour = col.set,shape=c(rep(16,total.terms),rep(20,length(only.entrys))),
                        size=-.12*vcount(link_tbl)+24.32, #<---------------------
                        position = "identity")+
          geom_node_text(aes(label= name),
                         position ="identity",
                         size = -0.01*vcount(link_tbl)+8.36, #<---------------------
                         repel = T,
                         vjust=3,
                         hjust=.3,
                         fontface = "bold")+
          theme_classic()
        graf2 = graf2 +  theme(legend.key.width = unit(1.2,"cm"),
                               legend.position = c(1.05, .91),
                               legend.key.height = unit(1,"cm"),
                               legend.title = element_text(size=24,face="bold"),
                               legend.text = element_text(size=22),
                               axis.line = element_blank(),
                               axis.text = element_blank(),
                               axis.ticks = element_blank(),
                               axis.title = element_blank())
        graf2 = graf2 + scale_color_gradient(low = c("green","black"), high = "red") +
          labs(caption = paste("GOmics:",sep =" ",format(Sys.time(),"%a(%d)/%b/%Y   %X")))
        {
          # combine plots
          plot.with.inset = ggdraw() +
            draw_plot(graf2,  width = .72, height = 1) +
            #                         x , y , width , height
            draw_plot(barplot.enrichment, .71,.01, .25, .825)
          ## Save plot
          png(file=paste("zxcvbnmPlot_Network_Bar_GO_Layout_",i,".png"),
              width = 25,
              height = 15,
              units = "in",
              res=600,bg="white")
          plot(plot.with.inset)
          dev.off()
          }
      }
      # >
      #
      # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
      # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
      #
      # |||||||||     1) Chord plots (without expression values)
      #
      # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
      # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
      #
      if ((vcount(link_tbl) > 1) && (vcount(link_tbl) <= 45)) {
        print("grafico 1 y 2")
        for (i in -1:0) {
          #
          # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
          # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
          #
          # |||||||||     1) Chord plots (without expression values)
          #
          # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
          # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
          #
          prot.col=rep("blue",1000)#9D9D9D
          col.set=c(many.colors[nodes$num[1:total.terms]],
                    prot.col[nodes$num[start.entry:total.nodes]])
          #
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
            chordDiagram(term.entry.mat, order = c(nodes$Term[1:total.terms],only.entrys),
                         grid.col = col.set,
                         transparency = 0,
                         directional = i,
                         annotationTrack = "grid",
                         #    c(% del datio del c?rculo, alto del track)
                         annotationTrackHeight = c(0.08, .01),
                         preAllocateTracks = 2)
            circos.track(track.index = 2,ylim = c(0, 1), track.height = uh(5, "cm"),
                         track.margin = c(5, uh(5, "cm")),
                         panel.fun = function(x, y) {
                           xcenter = get.cell.meta.data("xcenter")
                           circos.lines(c(xcenter, xcenter), c(0, uy(.7, "cm")),
                                        col = 1:4)
                         },bg.border =  NA)
            circos.track(track.index = 2, panel.fun = for(si in get.all.sector.index()[1:total.terms]) {
              xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 2)
              ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 2)
              circos.text(mean(xlim), mean(ylim), si, sector.index = si, track.index = 2,
                          facing = "downward",
                          niceFacing = TRUE,
                          adj = c(0,0.5),
                          cex = -0.004*vcount(link_tbl)+1.48, #<---------------------
                          font = 2,
                          col = NA)
              
            }, bg.border =  NA)
            circos.track(track.index = 2, panel.fun = for(si in get.all.sector.index()[start.entry:vcount(link_tbl)]) {
              xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 2)
              ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 2)
              circos.text(mean(xlim), mean(ylim), si, sector.index = si, track.index = 2,
                          facing = "clockwise",
                          niceFacing = TRUE,
                          adj = c(0.15,0.5),
                          cex = -0.005*vcount(link_tbl)+1.35, #<---------------------
                          font = 2,
                          col = "black")
            }, bg.border = NA)
          }
          #
          {
            png(file=paste0("Plot_Chord_1_",i,".png"),
                width = 30,
                height = 10,
                units = "in",
                res=900,bg="white")
            chord.plot.1()
            text(x=.4, y=.9, "ASPECT", adj = c(0,0), 
                 cex = -0.004*vcount(link_tbl)+1.48+.2,font = 2) #<---------------------
            dev.off()
          }
          #
          #
          # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
          # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
          #
          # |||||||||     2) Chord plots (with expression values)
          #
          # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
          # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
          #
          go.term = drop_na(nodes[,c("Entry","Term")])
          colnames(go.term)[1] <- "GO"
          #
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
          chord.plot.1 = function () {
            chordDiagram(go.entry.mat, order = nodes$Entry, # c(nodes$Term[1:total.terms],only.entrys)
                         grid.col = col.set,
                         transparency = 0,
                         directional = i,
                         annotationTrack = "grid",
                         #    c(% del datio del c?rculo, alto del track)
                         annotationTrackHeight = c(0.08, .01),
                         preAllocateTracks = 2)
            circos.track(track.index = 2,ylim = c(0, 1), track.height = uh(5, "cm"),
                         track.margin = c(5, uh(5, "cm")),
                         panel.fun = function(x, y) {
                           xcenter = get.cell.meta.data("xcenter")
                           circos.lines(c(xcenter, xcenter), c(0, uy(.5, "cm")),
                                        col = 1:4)
                         },bg.border =  NA)
            circos.track(track.index = 2, panel.fun = for(si in get.all.sector.index()[1:total.terms]) {
              xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 2)
              ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 2)
              circos.text(mean(xlim), mean(ylim), si, sector.index = si, track.index = 2,
                          facing = "clockwise",
                          niceFacing = TRUE,
                          adj = c(0.15,0.5),
                          cex = -0.004*vcount(link_tbl)+1.48, #<---------------------
                          font = 2,
                          col = "white")
            }, bg.border =  NA)
            circos.track(track.index = 2, panel.fun = for(si in get.all.sector.index()[start.entry:vcount(link_tbl)]) {
              xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 2)
              ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 2)
              circos.text(mean(xlim), mean(ylim), si, sector.index = si, track.index = 2,
                          facing = "clockwise",
                          niceFacing = TRUE,
                          adj = c(0.15,0.5),
                          cex = -0.005*vcount(link_tbl)+1.35, #<---------------------
                          font = 2,
                          col = "black")
            }, bg.border = NA)
          }
          #
          {
            png(file=paste0("Plot_Chord_2_",i,".png"),
                width = 22.5,
                height = 10,
                units = "in",
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
            chord.plot.1()
            upViewport()
            pushViewport(viewport(x = circle_size, y = 0.5, width = grobWidth(lgd_list_vertical),
                                  height = grobHeight(lgd_list_vertical), just = c("left", "center")))
            grid.draw(lgd_list_vertical)
            upViewport()
            text(x=.65, y=.75, "ASPECT", adj = c(0,0), 
                 cex = -0.004*vcount(link_tbl)+1.48+.2,font = 2) #<---------------------
            legend(x=.65,y=.7, 
                   legend=unique(paste0(nodes$Term[1:total.terms]," (",nodes$Freq[1:total.terms],")")),
                   cex=-0.004*vcount(link_tbl)+1.48,bg="transparent", #<---------------------
                   border = NA,pt.bg=many.colors,bty = "n",
                   x.intersp=1.3,y.intersp=-0.003*vcount(link_tbl)+1.21, #<---------------------
                   col=NA,text.font = 2,pch=22,
                   pt.cex=-0.005*vcount(link_tbl)+2.35)
            dev.off()
          }
          # >
        }
      } else if ((vcount(link_tbl) >= 46) && (vcount(link_tbl) <= 60)) {
        print("grafico 2")
        for (i in -1:0) {
          #
          # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
          # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
          #
          # |||||||||     2) Chord plots (with expression values)
          #
          # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
          # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
          #
          go.term = drop_na(nodes[,c("Entry","Term")])
          colnames(go.term)[1] <- "GO"
          #
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
          chord.plot.1 = function () {
            chordDiagram(go.entry.mat, order = nodes$Entry, # c(nodes$Term[1:total.terms],only.entrys)
                         grid.col = col.set,
                         transparency = 0,
                         directional = i,
                         annotationTrack = "grid",
                         #    c(% del datio del c?rculo, alto del track)
                         annotationTrackHeight = c(0.08, .01),
                         preAllocateTracks = 2)
            circos.track(track.index = 2,ylim = c(0, 1), track.height = uh(5, "cm"),
                         track.margin = c(5, uh(5, "cm")),
                         panel.fun = function(x, y) {
                           xcenter = get.cell.meta.data("xcenter")
                           circos.lines(c(xcenter, xcenter), c(0, uy(.5, "cm")),
                                        col = 1:4)
                         },bg.border =  NA)
            circos.track(track.index = 2, panel.fun = for(si in get.all.sector.index()[1:total.terms]) {
              xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 2)
              ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 2)
              circos.text(mean(xlim), mean(ylim), si, sector.index = si, track.index = 2,
                          facing = "clockwise",
                          niceFacing = TRUE,
                          adj = c(0.15,0.5),
                          cex = -0.004*vcount(link_tbl)+1.48, #<---------------------
                          font = 2,
                          col = "white")
            }, bg.border =  NA)
            circos.track(track.index = 2, panel.fun = for(si in get.all.sector.index()[start.entry:vcount(link_tbl)]) {
              xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 2)
              ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 2)
              circos.text(mean(xlim), mean(ylim), si, sector.index = si, track.index = 2,
                          facing = "clockwise",
                          niceFacing = TRUE,
                          adj = c(0.15,0.5),
                          cex = -0.005*vcount(link_tbl)+1.35, #<---------------------
                          font = 2,
                          col = "black")
            }, bg.border = NA)
          }
          #
          {
            png(file=paste0("Plot_Chord_2_",i,".png"),
                width = 22.5,
                height = 10,
                units = "in",
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
            chord.plot.1()
            upViewport()
            pushViewport(viewport(x = circle_size, y = 0.5, width = grobWidth(lgd_list_vertical),
                                  height = grobHeight(lgd_list_vertical), just = c("left", "center")))
            grid.draw(lgd_list_vertical)
            upViewport()
            text(x=.65, y=.75, "ASPECT", adj = c(0,0), 
                 cex = -0.004*vcount(link_tbl)+1.48+.2,font = 2) #<---------------------
            legend(x=.65,y=.7, 
                   legend=unique(paste0(nodes$Term[1:total.terms]," (",nodes$Freq[1:total.terms],")")),
                   cex=-0.004*vcount(link_tbl)+1.48,bg="transparent", #<---------------------
                   border = NA,pt.bg=many.colors,bty = "n",
                   x.intersp=1.3,y.intersp=-0.003*vcount(link_tbl)+1.21, #<---------------------
                   col=NA,text.font = 2,pch=22,
                   pt.cex=-0.005*vcount(link_tbl)+2.35)
            dev.off()
          }
          # >
        }
      } else if (vcount(link_tbl) >= 61) {
        print("grafico 3")
        for (i in -1:0) {
          #
          # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
          # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
          #
          # |||||||||     3) Chord plots (without expression values)
          #
          # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
          # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
          #
          go.term = drop_na(nodes[,c("Entry","Term")])
          colnames(go.term)[1] <- "GO"
          #
          go.entry.term = merge(links,go.term,by="GO")#[c(3,2)]
          ##
          go.entry = data.frame("GO"=go.entry.term$GO,"Entry"=go.entry.term$Entry)
          go.entry=as.data.frame.matrix(table(go.entry)) 
          go.entry.mat = as.matrix(go.entry) 
          #########
          circos.clear()
          circos.par(start.degree = 88, clock.wise = T)
          circos.par(track.margin = c(0, 0)) # 1
          #
          chord.plot.1 = function () {
            chordDiagram(go.entry.mat, order = nodes$Entry, # c(nodes$Term[1:total.terms],only.entrys)
                         grid.col = col.set,
                         transparency = 0,
                         directional = i,
                         annotationTrack = "grid",
                         #    c(% del datio del c?rculo, alto del track)
                         annotationTrackHeight = c(0.08, .01),
                         preAllocateTracks = 2)
            circos.track(track.index = 2,ylim = c(0, 1), track.height = uh(5, "cm"),
                         track.margin = c(5, uh(5, "cm")),
                         panel.fun = function(x, y) {
                           xcenter = get.cell.meta.data("xcenter")
                           circos.lines(c(xcenter, xcenter), c(0, uy(.5, "cm")),
                                        col = 1:4)
                         },bg.border =  NA)
            circos.track(track.index = 2, panel.fun = for(si in get.all.sector.index()[1:total.terms]) {
              xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 2)
              ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 2)
              circos.text(mean(xlim), mean(ylim), si, sector.index = si, track.index = 2,
                          facing = "clockwise",
                          niceFacing = TRUE,
                          adj = c(0.15,0.5),
                          cex = -0.004*vcount(link_tbl)+1.48, #<---------------------
                          font = 2,
                          col = "white")
            }, bg.border =  NA)
            circos.track(track.index = 2, panel.fun = for(si in get.all.sector.index()[start.entry:vcount(link_tbl)]) {
              xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 2)
              ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 2)
              circos.text(mean(xlim), mean(ylim), si, sector.index = si, track.index = 2,
                          facing = "clockwise",
                          niceFacing = TRUE,
                          adj = c(0.15,0.5),
                          cex = -0.005*vcount(link_tbl)+1.35, #<---------------------
                          font = 2,
                          col = "black")
            }, bg.border = NA)
          }
          #
          {
            png(file=paste0("Plot_Chord_3_",i,".png"),
                width = 22.5,
                height = 10,
                units = "in",
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
            chord.plot.1()
            upViewport()
            pushViewport(viewport(x = circle_size, y = 0.5, width = grobWidth(lgd_list_vertical),
                                  height = grobHeight(lgd_list_vertical), just = c("left", "center")))
            grid.draw(lgd_list_vertical)
            upViewport()
            text(x=.65, y=.75, "ASPECT", adj = c(0,0), 
                 cex = -0.004*vcount(link_tbl)+1.48+.2,font = 2) #<---------------------
            legend(x=.65,y=.7, 
                   legend=unique(paste0(nodes$Term[1:total.terms]," (",nodes$Freq[1:total.terms],")")),
                   cex=-0.004*vcount(link_tbl)+1.48,bg="transparent", #<---------------------
                   border = NA,pt.bg=many.colors,bty = "n",
                   x.intersp=1.3,y.intersp=-0.003*vcount(link_tbl)+1.21, #<---------------------
                   col=NA,text.font = 2,pch=22,
                   pt.cex=-0.005*vcount(link_tbl)+2.35)
            dev.off()
          }
          # >
        }
      } else {
        print("este mensaje no tiene que imprimirse")
      }
      # >
      # >
      #
      # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
      # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
      #
      # |||||||||     Upset Plot (without expression values)
      #
      # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
      # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
      #
      {
        val=length(only.entrys) # nrow number (total genes)
        #
        ups=as.data.frame.matrix(t(table(links)))
        ups[["log.Pval"]]=c(-log10(nodes$P[1:total.terms]),rep(NA,val-total.terms))
        ups[["term"]]=c(nodes$Entry[1:total.terms],rep(NA,(val-total.terms)))
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
          set.name.size=1.35
          #print(gnp.node.size)
        } else if ((length(unique(links$Entry)) >= 9.1) && (length(unique(links$Entry)) <= 13)) {
          set.name.size=1.2
          #print(gnp.node.size)
        } else if (length(unique(links$Entry)) >= 13.1) {
          set.name.size=0.9
          #print(gnp.node.size)
        }
        #
        png(file="zxcvbnmPlot_UpSetR.png",width = 10,height = 7,units = "in",res=700,bg="white")
        upset(ups,sets=nodes$Entry[1:total.terms],
              sets.bar.color = "darkgoldenrod2",order.by ="freq",empty.intersections = NULL,
              point.size=2,mainbar.y.label="poiuytrewq",sets.x.label = "",
              main.bar.color="black",matrix.color="black",shade.color="wheat3",
              line.size=0.5,show.numbers = "yes",group.by = "degree",
              matrix.dot.alpha = 0.5,mb.ratio = c(0.7, 0.3),
              #c(intersection size title, intersection size tick labels, set size title, set size tick labels, set names, numbers above bars)
              text.scale = c(1.5,1.5,0,1.5,set.name.size,1.5), 
              scale.sets = "identity",scale.intersections = "identity",
              shade.alpha = 0.35,set_size.angles = 0,
              attribute.plots = list(gridrows = 45, plots = list(list(plot = myplot2,x="term",y = "log.Pval", queries = F)),ncols=2))
        dev.off() 
      }
      #>
    } else {
      print("son 31 o m?s de 31, usar la configuraci?n de dos barras")
      print(paste("el ancho de la figura es:",width.png))
      # aqu? usar la nueva configuraci?n para agregar dos barras
      # controlar el ancho de la figura para esta figura don dos barras
      # >>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      # aqu? usar la nueva configuraci?n para agregar dos barras
      # controlar el ancho de la figura para esta figura don dos barras
      # First part of the data
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
                            rep("black",nrow(nodes.1_30))), size = .7 )+
        geom_point(aes(x = c(1:31), y = -log10(0.05),color=-log10(0.05)),
                   position = position_dodge(0.3),
                   color = c(rep("transparent",31-nrow(nodes.1_30)),
                             rep("black",nrow(nodes.1_30))), size = 4)+
        geom_line(aes(x = c(1:31), y = rev(log.fdr),
                      color=log.fdr),
                  color = c(rep("transparent",31-nrow(nodes.1_30)),
                            rep("blue",nrow(nodes.1_30))), size = .7 )+
        geom_point(aes(x = c(1:31), y = rev(log.fdr),color=log.fdr),
                   position = position_dodge(0.3),shape=17,
                   color = c(rep("transparent",31-nrow(nodes.1_30)),
                             rep("blue",nrow(nodes.1_30))), size = 4)+
        #scale_y_continuous(position = "right")+
        coord_flip()+
        scale_x_discrete(limits=rev(for.bar.1_30$Entry))+
        geom_text(aes(x=Entry,y=log.pval+.5,label=Freq),size=6.7,fontface ="bold")+ #<---------------------
      ggtitle("-log(P-value) and\n -log(Corrected P-value)")+
        theme(plot.title = element_text(vjust = 4,hjust=.08,size = 20), # element_text(hjust = .1,size = 12,face="bold"),#posocion y tama?o del titulo sobre el grafico
              axis.title.x =element_blank(),# element_text(colour="black",size=27,face="bold",hjust = 0,vjust = 10),
              axis.title.y = element_blank(),
              axis.text.y = element_text(colour=c("white",rep("black",30)),size=19,
                                         face="bold",hjust = 1.3,vjust = .5), #element_text(size = 8,color="black"),#tama?o y color de texto de eje x
              axis.text.x = element_text(colour="black",size=19,face="bold",vjust = 0),
              axis.ticks.length =unit(.2, "cm"), #reduce la longitud de las marcas en la linea de escala de los ejes x y
              axis.ticks.y.left = element_line(colour = "black",size=3,unit(0, "cm")),
              axis.ticks.x.top = element_line(colour = "black",size=3,unit(3, "cm")),
              axis.line.y = element_blank(), axis.title.x.top =  element_blank(), 
              axis.ticks = element_line(color = "black",size = 1),
              plot.margin = margin(1.5,1,0,1, "cm"))+
        scale_y_continuous(position = "right",
                           limit = c(0,round(max(for.bar.1_30$log.pval[1:nrow(nodes.1_30)]+1.5))),
                           breaks=seq(0,round(max(for.bar.1_30$log.pval[1:nrow(nodes.1_30)]+1.5),digits = 0),1))+
        labs(title = expression(bold("ASPECT")),
             subtitle = expression(bold(" -log10("*italic("P")*"-value & adj."*italic("P")*"-value)")))+
        theme(plot.title = element_text(vjust = 7,hjust=-.3,size = 30),
              plot.subtitle = element_text(vjust = 3,hjust=.08,size = 20))
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
                            rep("black",nrow(nodes.31_final))), size = .7 )+
        geom_point(aes(x = c(1:31), y = -log10(0.05),color=-log10(0.05)),
                   position = position_dodge(0.3),
                   color = c(rep("transparent",31-nrow(nodes.31_final)),
                             rep("black",nrow(nodes.31_final))), size = 4)+
        geom_line(aes(x = c(1:31), y = rev(log.fdr),
                      color=log.fdr),
                  color = c(rep("transparent",31-nrow(nodes.31_final)),
                            rep("blue",nrow(nodes.31_final))), size = .7 )+
        geom_point(aes(x = c(1:31), y = rev(log.fdr),color=log.fdr),
                   position = position_dodge(0.3),shape=17,
                   color = c(rep("transparent",31-nrow(nodes.31_final)),
                             rep("blue",nrow(nodes.31_final))), size = 4)+
        #scale_y_continuous(position = "right")+
        coord_flip()+
        scale_x_discrete(limits=rev(for.bar.31_final$Entry))+
        geom_text(aes(x=Entry,y=log.pval+.5,label=Freq),size=6.7,fontface ="bold")+ #<---------------------
      ggtitle("-log(P-value) and\n -log(Corrected P-value)")+
        theme(plot.title = element_text(vjust = 4,hjust=.08,size = 20), # element_text(hjust = .1,size = 12,face="bold"),#posocion y tama?o del titulo sobre el grafico
              axis.title.x =element_blank(),# element_text(colour="black",size=27,face="bold",hjust = 0,vjust = 10),
              axis.title.y = element_blank(),
              axis.text.y = element_text(colour=c("white",rep("black",30)),size=19,
                                         face="bold",hjust = 1.3,vjust = .5), #element_text(size = 8,color="black"),#tama?o y color de texto de eje x
              axis.text.x = element_text(colour="black",size=19,face="bold",vjust = 0),
              axis.ticks.length =unit(.2, "cm"), #reduce la longitud de las marcas en la linea de escala de los ejes x y
              axis.ticks.y.left = element_line(colour = "black",size=3,unit(0, "cm")),
              axis.ticks.x.top = element_line(colour = "black",size=3,unit(3, "cm")),
              axis.line.y = element_blank(), axis.title.x.top =  element_blank(), 
              axis.ticks = element_line(color = "black",size = 1),
              plot.margin = margin(1.5,1,0,1, "cm"))+
        scale_y_continuous(position = "right",
                           limit = c(0,round(max(for.bar.1_30$log.pval[1:nrow(nodes.1_30)]+1.5))),
                           breaks=seq(0,round(max(for.bar.1_30$log.pval[1:nrow(nodes.1_30)]+1.5),digits = 0),1))+
        labs(title = expression(bold("")),
             subtitle = expression(bold(" -log10("*italic("P")*"-value & adj."*italic("P")*"-value)")))+
        theme(plot.title = element_text(vjust = 7,hjust=-.3,size = 30),
              plot.subtitle = element_text(vjust = 3,hjust=.08,size = 20))
      #
      { # plots_combined used for network
        barplot.p.val.1_30.user = barplot.p.val.1_30 +
          labs(caption = paste("GOmics:",sep=" ",format(Sys.time(),"%a(%d)/%b/%Y   %X")))
        plots_combined  = ggdraw() +
          draw_plot(barplot.p.val.1_30.user, x = 0, y = 0, width = .5, height = 1) +
          draw_plot(barplot.p.val.31_final, x = .45, y = 0, width = .5, height = 1)
        # Barplot for user
        png(file="zxcvbnmPlot_Bar.png",
            width = 30,
            height = 30,
            units = "cm",
            res=900,bg="white")
        plot(plots_combined)
        dev.off()
        }
      #
    #
    # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    #
    # |||||||||     1)    4 simple Plots with Terms as label
    #
    # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    #
    layouts=c("nicely","kk","dh")
    #
    for (i in layouts){
      graf1 = ggraph(link_tbl,layout = "igraph", algorithm = i) +
        geom_edge_link(aes(colour = GO),alpha = .7,width = -0.03*vcount(link_tbl)+8.08)+ #<---------------------
      geom_node_point(colour = col.set,shape=c(rep(16,total.terms),rep(20,length(only.entrys))),
                      size=-.12*vcount(link_tbl)+24.32,show.legend = NA, #<---------------------
                      position = "identity")+
        scale_edge_colour_manual(values=col.set,name = "ASPECT")+
        #scale_edge_colour_manual(values=colors.fold.change)+
        geom_node_text(aes(label=Label),
                       position ="identity",
                       size = -0.01*vcount(link_tbl)+8.36, #<---------------------
                       repel = T,
                       vjust=3,
                       hjust=.3,
                       fontface = "bold")+
        theme_classic() +
        theme(legend.key.width=unit(1.2,"cm"),
              legend.key.height=unit(1.2,"cm"),
              legend.title=element_text(size=27,face="bold"),
              legend.text=element_text(size=25,face="bold"),
              axis.line=element_blank(),
              axis.text = element_blank(),
              axis.ticks = element_blank(),
              axis.title = element_blank())
      graf1 = graf1 +
        labs(caption = paste("GOmics:",sep=" ",format(Sys.time(),"%a(%d)/%b/%Y   %X")))
      {
        ## Save plot
        png(file=paste("zxcvbnmPlot_Network_Term_Layout_",i,".png"),
            width = width.png,
            height = 15,
            units = "in",
            res=600,bg="white")
        plot(graf1)
        dev.off()
        }
    }
    #
    # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    #
    # |||||||||    2)     4 simple Plots with GO as label
    #
    # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    #
    for (i in layouts){
      graf1 = ggraph(link_tbl,layout = "igraph", algorithm = i) +
        geom_edge_link(aes(colour = GO),alpha = .7,width = -0.03*vcount(link_tbl)+8.08)+ #<---------------------
      geom_node_point(colour = col.set,shape=c(rep(16,total.terms),rep(20,length(only.entrys))),
                      size=-.12*vcount(link_tbl)+24.32,show.legend = NA, #<---------------------
                      position = "identity")+
        scale_edge_colour_manual(values=col.set,name = "ASPECT")+
        #scale_edge_colour_manual(values=colors.fold.change)+
        geom_node_text(aes(label=name),
                       position ="identity",
                       size = -0.01*vcount(link_tbl)+8.36, #<---------------------
                       repel = T,
                       vjust=3,
                       hjust=.3,
                       fontface = "bold")+
        theme_classic() +
        theme(legend.key.width=unit(1.2,"cm"),
              legend.key.height=unit(1.2,"cm"),
              legend.title=element_text(size=27,face="bold"),
              legend.text=element_text(size=25,face="bold"),
              axis.line=element_blank(),
              axis.text = element_blank(),
              axis.ticks = element_blank(),
              axis.title = element_blank())
      graf1 = graf1 +
        labs(caption = paste("GOmics:",sep=" ",format(Sys.time(),"%a(%d)/%b/%Y   %X")))
      {
        ## Save plot
        png(file=paste("zxcvbnmPlot_Network_GO_Layout_",i,".png"),
            width = width.png,
            height = 15,
            units = "in",
            res=600,bg="white")
        plot(graf1)
        dev.off()
        }
    }
    #
    # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    #
    # |||||||||    1)     4 Complex Plots with Term as label
    #
    # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    #
    ## Construction of graphics with different designs
    #
    for (i in layouts){
      graf2 = ggraph(link_tbl,layout = "igraph",algorithm = i) +
        geom_edge_link(colour = col.edges.cluster,alpha = .7,
                       width = -0.03*vcount(link_tbl)+8.08)+ #<---------------------
      geom_node_point(colour = col.set,shape=c(rep(16,total.terms),rep(20,length(only.entrys))),
                      size=-.12*vcount(link_tbl)+24.32, #<---------------------
                      position = "identity")+
        geom_node_text(aes(label=Label),
                       position ="identity",
                       size = -0.01*vcount(link_tbl)+8.36, #<---------------------
                       repel = T,
                       vjust=3,
                       hjust=.3,
                       fontface = "bold")+
        theme_classic()
      graf2 = graf2 +  theme(legend.key.width=unit(1.2,"cm"),
                             legend.position = c(1.05, .91),
                             legend.key.height=unit(1,"cm"),
                             legend.title=element_text(size=24,face="bold"),
                             legend.text=element_text(size=22),
                             axis.line=element_blank(),
                             axis.text = element_blank(),
                             axis.ticks = element_blank(),
                             axis.title = element_blank())
      graf2 = graf2 + scale_color_gradient(low = c("green","black"), high = "red") +
        labs(caption = paste("GOmics:",sep=" ",format(Sys.time(),"%a(%d)/%b/%Y   %X")))
      {
        plots_combined  = ggdraw() +
          draw_plot(barplot.p.val.1_30, x = 0, y = 0, width = .5, height = 1) +
          draw_plot(barplot.p.val.31_final, x = .45, y = 0, width = .5, height = 1)
        # combine plots
        plot.with.inset = ggdraw() +
          draw_plot(graf2,  width = .54, height = 1) +
          #                         x , y , width , height
          draw_plot(plots_combined, .52,.01, .45, .825) # 2 barplots
        
        ## Save plot
        png(file=paste("zxcvbnmPlot_Network_Bar_Term_Layout_",i,".png"),
            width = 35,
            height = 15,
            units = "in",
            res=600,bg="white")
        plot(plot.with.inset)
        dev.off()
        }
    }
    #
    # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    #
    # |||||||||     2)    4 Complex Plots with GO as label
    #
    # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    #
    ## Construction of graphics with different designs
    #
    for (i in layouts){
      graf2 = ggraph(link_tbl,layout = "igraph",algorithm = i) +
        geom_edge_link(colour = col.edges.cluster,alpha = .7,
                       width = -0.03*vcount(link_tbl)+8.08)+ #<---------------------
      geom_node_point(colour = col.set,shape=c(rep(16,total.terms),rep(20,length(only.entrys))),
                      size=-.12*vcount(link_tbl)+24.32, #<---------------------
                      position = "identity")+
        geom_node_text(aes(label= name),
                       position ="identity",
                       size = -0.01*vcount(link_tbl)+8.36, #<---------------------
                       repel = T,
                       vjust=3,
                       hjust=.3,
                       fontface = "bold")+
        theme_classic()
      graf2 = graf2 +  theme(legend.key.width=unit(1.2,"cm"),
                             legend.position = c(1.05, .91),
                             legend.key.height=unit(1,"cm"),
                             legend.title=element_text(size=24,face="bold"),
                             legend.text=element_text(size=22),
                             axis.line=element_blank(),
                             axis.text = element_blank(),
                             axis.ticks = element_blank(),
                             axis.title = element_blank())
      graf2 = graf2 + scale_color_gradient(low = c("green","black"), high = "red") +
        labs(caption = paste("GOmics:",sep=" ",format(Sys.time(),"%a(%d)/%b/%Y   %X")))
      {
        plots_combined  = ggdraw() +
          draw_plot(barplot.p.val.1_30, x = 0, y = 0, width = .5, height = 1) +
          draw_plot(barplot.p.val.31_final, x = .45, y = 0, width = .5, height = 1)
        # combine plots
        plot.with.inset = ggdraw() +
          draw_plot(graf2,  width = .54, height = 1) +
          #                         x , y , width , height
          draw_plot(plots_combined, .52,.01, .45, .825) # 2 barplots
        
        ## Save plot
        png(file=paste("zxcvbnmPlot_Network_Bar_GO_Layout_",i,".png"),
            width = 35,
            height = 15,
            units = "in",
            res=600,bg="white")
        plot(plot.with.inset)
        dev.off()
        }
    }
    # >
    #
    # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    #
    # |||||||||     1) Chord plots (without expression values)
    #
    # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    #
    if ((vcount(link_tbl) > 1) && (vcount(link_tbl) <= 45)) {
      print("grafico 1 y 2")
      for (i in -1:0) {
        #
        # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
        # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
        #
        # |||||||||     1) Chord plots (without expression values)
        #
        # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
        # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
        #
        prot.col=rep("blue",1000)#9D9D9D
        col.set=c(many.colors[nodes$num[1:total.terms]],
                  prot.col[nodes$num[start.entry:total.nodes]])
        #
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
          chordDiagram(term.entry.mat, order = c(nodes$Term[1:total.terms],only.entrys),
                       grid.col = col.set,
                       transparency = 0,
                       directional = i,
                       annotationTrack = "grid",
                       #    c(% del datio del c?rculo, alto del track)
                       annotationTrackHeight = c(0.08, .01),
                       preAllocateTracks = 2)
          circos.track(track.index = 2,ylim = c(0, 1), track.height = uh(5, "cm"),
                       track.margin = c(5, uh(5, "cm")),
                       panel.fun = function(x, y) {
                         xcenter = get.cell.meta.data("xcenter")
                         circos.lines(c(xcenter, xcenter), c(0, uy(.7, "cm")),
                                      col = 1:4)
                       },bg.border =  NA)
          circos.track(track.index = 2, panel.fun = for(si in get.all.sector.index()[1:total.terms]) {
            xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 2)
            ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 2)
            circos.text(mean(xlim), mean(ylim), si, sector.index = si, track.index = 2,
                        facing = "downward",
                        niceFacing = TRUE,
                        adj = c(0,0.5),
                        cex = -0.004*vcount(link_tbl)+1.48, #<---------------------
                        font = 2,
                        col = NA)
            
          }, bg.border =  NA)
          circos.track(track.index = 2, panel.fun = for(si in get.all.sector.index()[start.entry:vcount(link_tbl)]) {
            xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 2)
            ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 2)
            circos.text(mean(xlim), mean(ylim), si, sector.index = si, track.index = 2,
                        facing = "clockwise",
                        niceFacing = TRUE,
                        adj = c(0.15,0.5),
                        cex = -0.005*vcount(link_tbl)+1.35, #<---------------------
                        font = 2,
                        col = "black")
          }, bg.border = NA)
        }
        #
        {
          png(file=paste0("Plot_Chord_1_",i,".png"),
              width = 30,
              height = 10,
              units = "in",
              res=900,bg="white")
          chord.plot.1()
          text(x=.4, y=.9, "ASPECT", adj = c(0,0), 
               cex = -0.004*vcount(link_tbl)+1.48+.2,font = 2) #<---------------------
          dev.off()
        }
        #
        #
        # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
        # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
        #
        # |||||||||     2) Chord plots (with expression values)
        #
        # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
        # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
        #
        go.term = drop_na(nodes[,c("Entry","Term")])
        colnames(go.term)[1] <- "GO"
        #
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
        chord.plot.1 = function () {
          chordDiagram(go.entry.mat, order = nodes$Entry, # c(nodes$Term[1:total.terms],only.entrys)
                       grid.col = col.set,
                       transparency = 0,
                       directional = i,
                       annotationTrack = "grid",
                       #    c(% del datio del c?rculo, alto del track)
                       annotationTrackHeight = c(0.08, .01),
                       preAllocateTracks = 2)
          circos.track(track.index = 2,ylim = c(0, 1), track.height = uh(5, "cm"),
                       track.margin = c(5, uh(5, "cm")),
                       panel.fun = function(x, y) {
                         xcenter = get.cell.meta.data("xcenter")
                         circos.lines(c(xcenter, xcenter), c(0, uy(.5, "cm")),
                                      col = 1:4)
                       },bg.border =  NA)
          circos.track(track.index = 2, panel.fun = for(si in get.all.sector.index()[1:total.terms]) {
            xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 2)
            ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 2)
            circos.text(mean(xlim), mean(ylim), si, sector.index = si, track.index = 2,
                        facing = "clockwise",
                        niceFacing = TRUE,
                        adj = c(0.15,0.5),
                        cex = -0.004*vcount(link_tbl)+1.48, #<---------------------
                        font = 2,
                        col = "white")
          }, bg.border =  NA)
          circos.track(track.index = 2, panel.fun = for(si in get.all.sector.index()[start.entry:vcount(link_tbl)]) {
            xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 2)
            ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 2)
            circos.text(mean(xlim), mean(ylim), si, sector.index = si, track.index = 2,
                        facing = "clockwise",
                        niceFacing = TRUE,
                        adj = c(0.15,0.5),
                        cex = -0.005*vcount(link_tbl)+1.35, #<---------------------
                        font = 2,
                        col = "black")
          }, bg.border = NA)
        }
        #
        {
          png(file=paste0("Plot_Chord_2_",i,".png"),
              width = 22.5,
              height = 10,
              units = "in",
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
          chord.plot.1()
          upViewport()
          pushViewport(viewport(x = circle_size, y = 0.5, width = grobWidth(lgd_list_vertical),
                                height = grobHeight(lgd_list_vertical), just = c("left", "center")))
          grid.draw(lgd_list_vertical)
          upViewport()
          text(x=.65, y=.75, "ASPECT", adj = c(0,0), 
               cex = -0.004*vcount(link_tbl)+1.48+.2,font = 2) #<---------------------
          legend(x=.65,y=.7, 
                 legend=unique(paste0(nodes$Term[1:total.terms]," (",nodes$Freq[1:total.terms],")")),
                 cex=-0.004*vcount(link_tbl)+1.48,bg="transparent", #<---------------------
                 border = NA,pt.bg=many.colors,bty = "n",
                 x.intersp=1.3,y.intersp=-0.003*vcount(link_tbl)+1.21, #<---------------------
                 col=NA,text.font = 2,pch=22,
                 pt.cex=-0.005*vcount(link_tbl)+2.35)
          dev.off()
        }
        # >
      }
    } else if ((vcount(link_tbl) >= 46) && (vcount(link_tbl) <= 60)) {
      print("grafico 2")
      for (i in -1:0) {
        #
        # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
        # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
        #
        # |||||||||     2) Chord plots (with expression values)
        #
        # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
        # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
        #
        go.term = drop_na(nodes[,c("Entry","Term")])
        colnames(go.term)[1] <- "GO"
        #
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
        chord.plot.1 = function () {
          chordDiagram(go.entry.mat, order = nodes$Entry, # c(nodes$Term[1:total.terms],only.entrys)
                       grid.col = col.set,
                       transparency = 0,
                       directional = i,
                       annotationTrack = "grid",
                       #    c(% del datio del c?rculo, alto del track)
                       annotationTrackHeight = c(0.08, .01),
                       preAllocateTracks = 2)
          circos.track(track.index = 2,ylim = c(0, 1), track.height = uh(5, "cm"),
                       track.margin = c(5, uh(5, "cm")),
                       panel.fun = function(x, y) {
                         xcenter = get.cell.meta.data("xcenter")
                         circos.lines(c(xcenter, xcenter), c(0, uy(.5, "cm")),
                                      col = 1:4)
                       },bg.border =  NA)
          circos.track(track.index = 2, panel.fun = for(si in get.all.sector.index()[1:total.terms]) {
            xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 2)
            ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 2)
            circos.text(mean(xlim), mean(ylim), si, sector.index = si, track.index = 2,
                        facing = "clockwise",
                        niceFacing = TRUE,
                        adj = c(0.15,0.5),
                        cex = -0.004*vcount(link_tbl)+1.48, #<---------------------
                        font = 2,
                        col = "white")
          }, bg.border =  NA)
          circos.track(track.index = 2, panel.fun = for(si in get.all.sector.index()[start.entry:vcount(link_tbl)]) {
            xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 2)
            ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 2)
            circos.text(mean(xlim), mean(ylim), si, sector.index = si, track.index = 2,
                        facing = "clockwise",
                        niceFacing = TRUE,
                        adj = c(0.15,0.5),
                        cex = -0.005*vcount(link_tbl)+1.35, #<---------------------
                        font = 2,
                        col = "black")
          }, bg.border = NA)
        }
        #
        {
          png(file=paste0("Plot_Chord_2_",i,".png"),
              width = 22.5,
              height = 10,
              units = "in",
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
          chord.plot.1()
          upViewport()
          pushViewport(viewport(x = circle_size, y = 0.5, width = grobWidth(lgd_list_vertical),
                                height = grobHeight(lgd_list_vertical), just = c("left", "center")))
          grid.draw(lgd_list_vertical)
          upViewport()
          text(x=.65, y=.75, "ASPECT", adj = c(0,0), 
               cex = -0.004*vcount(link_tbl)+1.48+.2,font = 2) #<---------------------
          legend(x=.65,y=.7, 
                 legend=unique(paste0(nodes$Term[1:total.terms]," (",nodes$Freq[1:total.terms],")")),
                 cex=-0.004*vcount(link_tbl)+1.48,bg="transparent", #<---------------------
                 border = NA,pt.bg=many.colors,bty = "n",
                 x.intersp=1.3,y.intersp=-0.003*vcount(link_tbl)+1.21, #<---------------------
                 col=NA,text.font = 2,pch=22,
                 pt.cex=-0.005*vcount(link_tbl)+2.35)
          dev.off()
        }
        # >
      }
    } else if (vcount(link_tbl) >= 61) {
      print("grafico 3")
      for (i in -1:0) {
        #
        # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
        # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
        #
        # |||||||||     3) Chord plots (without expression values)
        #
        # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
        # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
        #
        go.term = drop_na(nodes[,c("Entry","Term")])
        colnames(go.term)[1] <- "GO"
        #
        go.entry.term = merge(links,go.term,by="GO")#[c(3,2)]
        ##
        go.entry = data.frame("GO"=go.entry.term$GO,"Entry"=go.entry.term$Entry)
        go.entry=as.data.frame.matrix(table(go.entry)) 
        go.entry.mat = as.matrix(go.entry) 
        #########
        circos.clear()
        circos.par(start.degree = 88, clock.wise = T)
        circos.par(track.margin = c(0, 0)) # 1
        #
        chord.plot.1 = function () {
          chordDiagram(go.entry.mat, order = nodes$Entry, # c(nodes$Term[1:total.terms],only.entrys)
                       grid.col = col.set,
                       transparency = 0,
                       directional = i,
                       annotationTrack = "grid",
                       #    c(% del datio del c?rculo, alto del track)
                       annotationTrackHeight = c(0.08, .01),
                       preAllocateTracks = 2)
          circos.track(track.index = 2,ylim = c(0, 1), track.height = uh(5, "cm"),
                       track.margin = c(5, uh(5, "cm")),
                       panel.fun = function(x, y) {
                         xcenter = get.cell.meta.data("xcenter")
                         circos.lines(c(xcenter, xcenter), c(0, uy(.5, "cm")),
                                      col = 1:4)
                       },bg.border =  NA)
          circos.track(track.index = 2, panel.fun = for(si in get.all.sector.index()[1:total.terms]) {
            xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 2)
            ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 2)
            circos.text(mean(xlim), mean(ylim), si, sector.index = si, track.index = 2,
                        facing = "clockwise",
                        niceFacing = TRUE,
                        adj = c(0.15,0.5),
                        cex = -0.004*vcount(link_tbl)+1.48, #<---------------------
                        font = 2,
                        col = "white")
          }, bg.border =  NA)
          circos.track(track.index = 2, panel.fun = for(si in get.all.sector.index()[start.entry:vcount(link_tbl)]) {
            xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 2)
            ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 2)
            circos.text(mean(xlim), mean(ylim), si, sector.index = si, track.index = 2,
                        facing = "clockwise",
                        niceFacing = TRUE,
                        adj = c(0.15,0.5),
                        cex = -0.005*vcount(link_tbl)+1.35, #<---------------------
                        font = 2,
                        col = "black")
          }, bg.border = NA)
        }
        #
        {
          png(file=paste0("Plot_Chord_3_",i,".png"),
              width = 22.5,
              height = 10,
              units = "in",
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
          chord.plot.1()
          upViewport()
          pushViewport(viewport(x = circle_size, y = 0.5, width = grobWidth(lgd_list_vertical),
                                height = grobHeight(lgd_list_vertical), just = c("left", "center")))
          grid.draw(lgd_list_vertical)
          upViewport()
          text(x=.65, y=.75, "ASPECT", adj = c(0,0), 
               cex = -0.004*vcount(link_tbl)+1.48+.2,font = 2) #<---------------------
          legend(x=.65,y=.7, 
                 legend=unique(paste0(nodes$Term[1:total.terms]," (",nodes$Freq[1:total.terms],")")),
                 cex=-0.004*vcount(link_tbl)+1.48,bg="transparent", #<---------------------
                 border = NA,pt.bg=many.colors,bty = "n",
                 x.intersp=1.3,y.intersp=-0.003*vcount(link_tbl)+1.21, #<---------------------
                 col=NA,text.font = 2,pch=22,
                 pt.cex=-0.005*vcount(link_tbl)+2.35)
          dev.off()
        }
        # >
      }
    } else {
      print("este mensaje no tiene que imprimirse")
    }
    # >
    # >
    #
    # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    #
    # |||||||||     Upset Plot (without expression values)
    #
    # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    #
    {
      val=length(only.entrys) # nrow number (total genes)
      #
      ups=as.data.frame.matrix(t(table(links)))
      ups[["log.Pval"]]=c(-log10(nodes$P[1:total.terms]),rep(NA,val-total.terms))
      ups[["term"]]=c(nodes$Entry[1:total.terms],rep(NA,(val-total.terms)))
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
        set.name.size=1.35
        #print(gnp.node.size)
      } else if ((length(unique(links$Entry)) >= 9.1) && (length(unique(links$Entry)) <= 13)) {
        set.name.size=1.2
        #print(gnp.node.size)
      } else if (length(unique(links$Entry)) >= 13.1) {
        set.name.size=0.9
        #print(gnp.node.size)
      }
      #
      png(file="zxcvbnmPlot_UpSetR.png",width = 10,height = 7,units = "in",res=700,bg="white")
      upset(ups,sets=nodes$Entry[1:total.terms],
            sets.bar.color = "darkgoldenrod2",order.by ="freq",empty.intersections = NULL,
            point.size=2,mainbar.y.label="poiuytrewq",sets.x.label = "",
            main.bar.color="black",matrix.color="black",shade.color="wheat3",
            line.size=0.5,show.numbers = "yes",group.by = "degree",
            matrix.dot.alpha = 0.5,mb.ratio = c(0.7, 0.3),
            #c(intersection size title, intersection size tick labels, set size title, set size tick labels, set names, numbers above bars)
            text.scale = c(1.5,1.5,0,1.5,set.name.size,1.5), 
            scale.sets = "identity",scale.intersections = "identity",
            shade.alpha = 0.35,set_size.angles = 0,
            attribute.plots = list(gridrows = 45, plots = list(list(plot = myplot2,x="term",y = "log.Pval", queries = F)),ncols=2))
      dev.off() 
    }
    #>
    }
  } else {
    print("hay 61 o m?s de 61, se puede hacer ?nicamente plots simples")
    print(paste("el ancho de la figura es:",width.png))
    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>
    #
    # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    #
    # |||||||||     1)    4 simple Plots with Terms as label
    #
    # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    #
    layouts=c("nicely","kk","dh")
    #
    for (i in layouts){
      graf1 = ggraph(link_tbl,layout = "igraph", algorithm = i) +
        geom_edge_link(aes(colour = GO),alpha = .7,width = -0.03*vcount(link_tbl)+8.08)+ #<---------------------
      geom_node_point(colour = col.set,shape=c(rep(16,total.terms),rep(20,length(only.entrys))),
                      size=-.12*vcount(link_tbl)+24.32,show.legend = NA, #<---------------------
                      position = "identity")+
        scale_edge_colour_manual(values=col.set,name = "ASPECT")+
        #scale_edge_colour_manual(values=colors.fold.change)+
        geom_node_text(aes(label=Label),
                       position ="identity",
                       size = -0.01*vcount(link_tbl)+8.36, #<---------------------
                       repel = T,
                       vjust=3,
                       hjust=.3,
                       fontface = "bold")+
        theme_classic() +
        theme(legend.key.width=unit(1.2,"cm"),
              legend.key.height=unit(1.2,"cm"),
              legend.title=element_text(size=27,face="bold"),
              legend.text=element_text(size=25,face="bold"),
              axis.line=element_blank(),
              axis.text = element_blank(),
              axis.ticks = element_blank(),
              axis.title = element_blank())
      graf1 = graf1 +
        labs(caption = paste("GOmics:",sep=" ",format(Sys.time(),"%a(%d)/%b/%Y   %X")))
      {
        png(file=paste("zxcvbnmPlot_Network_Ter_Layout_",i,".png"),
            width = width.png,  #<---------------------
            height = 15, 
            units = "in",
            res=600,bg="white")
        plot(graf1)
        dev.off()
        }
    }
    #
    # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    #
    # |||||||||    2)     4 simple Plots with GO as label
    #
    # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    #
    for (i in layouts){
      graf1 = ggraph(link_tbl,layout = "igraph", algorithm = i) +
        geom_edge_link(aes(colour = GO),alpha = .7,width = -0.03*vcount(link_tbl)+8.08)+ #<---------------------
      geom_node_point(colour = col.set,shape=c(rep(16,total.terms),rep(20,length(only.entrys))),
                      size=-.12*vcount(link_tbl)+24.32,show.legend = NA, #<---------------------
                      position = "identity")+
        scale_edge_colour_manual(values=col.set,name = "ASPECT")+
        #scale_edge_colour_manual(values=colors.fold.change)+
        geom_node_text(aes(label=name),
                       position ="identity",
                       size = -0.01*vcount(link_tbl)+8.36, #<---------------------
                       repel = T,
                       vjust=3,
                       hjust=.3,
                       fontface = "bold")+
        theme_classic() +
        theme(legend.key.width=unit(1.2,"cm"),
              legend.key.height=unit(1.2,"cm"),
              legend.title=element_text(size=27,face="bold"),
              legend.text=element_text(size=25,face="bold"),
              axis.line=element_blank(),
              axis.text = element_blank(),
              axis.ticks = element_blank(),
              axis.title = element_blank())
      graf1 = graf1 +
        labs(caption = paste("GOmics:",sep=" ",format(Sys.time(),"%a(%d)/%b/%Y   %X")))
      {
        png(file=paste("zxcvbnmPlot_Network_GO_Layout_",i,".png"),
            width = width.png,  #<--------------------- 
            height = 15, 
            units = "in",
            res=600,bg="white")
        plot(graf1)
        dev.off()
        }
    }
    # >
    #
    # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    #
    # |||||||||     1) Chord plots (without expression values)
    #
    # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    #
    prot.col=rep("blue",1000)#9D9D9D
    col.set=c(many.colors[nodes$num[1:total.terms]],
              prot.col[nodes$num[start.entry:total.nodes]])
    #
    if ((vcount(link_tbl) > 1) && (vcount(link_tbl) <= 45)) {
      print("grafico 1 y 2")
      for (i in -1:0) {
        #
        # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
        # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
        #
        # |||||||||     1) Chord plots (without expression values)
        #
        # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
        # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
        #
        prot.col=rep("blue",1000)#9D9D9D
        col.set=c(many.colors[nodes$num[1:total.terms]],
                  prot.col[nodes$num[start.entry:total.nodes]])
        #
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
          chordDiagram(term.entry.mat, order = c(nodes$Term[1:total.terms],only.entrys),
                       grid.col = col.set,
                       transparency = 0,
                       directional = i,
                       annotationTrack = "grid",
                       #    c(% del datio del c?rculo, alto del track)
                       annotationTrackHeight = c(0.08, .01),
                       preAllocateTracks = 2)
          circos.track(track.index = 2,ylim = c(0, 1), track.height = uh(5, "cm"),
                       track.margin = c(5, uh(5, "cm")),
                       panel.fun = function(x, y) {
                         xcenter = get.cell.meta.data("xcenter")
                         circos.lines(c(xcenter, xcenter), c(0, uy(.7, "cm")),
                                      col = 1:4)
                       },bg.border =  NA)
          circos.track(track.index = 2, panel.fun = for(si in get.all.sector.index()[1:total.terms]) {
            xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 2)
            ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 2)
            circos.text(mean(xlim), mean(ylim), si, sector.index = si, track.index = 2,
                        facing = "downward",
                        niceFacing = TRUE,
                        adj = c(0,0.5),
                        cex = -0.004*vcount(link_tbl)+1.48, #<---------------------
                        font = 2,
                        col = NA)
            
          }, bg.border =  NA)
          circos.track(track.index = 2, panel.fun = for(si in get.all.sector.index()[start.entry:vcount(link_tbl)]) {
            xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 2)
            ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 2)
            circos.text(mean(xlim), mean(ylim), si, sector.index = si, track.index = 2,
                        facing = "clockwise",
                        niceFacing = TRUE,
                        adj = c(0.15,0.5),
                        cex = -0.005*vcount(link_tbl)+1.35, #<---------------------
                        font = 2,
                        col = "black")
          }, bg.border = NA)
        }
        #
        {
          png(file=paste0("Plot_Chord_1_",i,".png"),
              width = 30,
              height = 10,
              units = "in",
              res=900,bg="white")
          chord.plot.1()
          text(x=.4, y=.9, "ASPECT", adj = c(0,0), 
               cex = -0.004*vcount(link_tbl)+1.48+.2,font = 2) #<---------------------
          dev.off()
        }
        #
        #
        # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
        # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
        #
        # |||||||||     2) Chord plots (without expression values)
        #
        # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
        # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
        #
        go.term = drop_na(nodes[,c("Entry","Term")])
        colnames(go.term)[1] <- "GO"
        #
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
        chord.plot.1 = function () {
          chordDiagram(go.entry.mat, order = nodes$Entry, # c(nodes$Term[1:total.terms],only.entrys)
                       grid.col = col.set,
                       transparency = 0,
                       directional = i,
                       annotationTrack = "grid",
                       #    c(% del datio del c?rculo, alto del track)
                       annotationTrackHeight = c(0.08, .01),
                       preAllocateTracks = 2)
          circos.track(track.index = 2,ylim = c(0, 1), track.height = uh(5, "cm"),
                       track.margin = c(5, uh(5, "cm")),
                       panel.fun = function(x, y) {
                         xcenter = get.cell.meta.data("xcenter")
                         circos.lines(c(xcenter, xcenter), c(0, uy(.5, "cm")),
                                      col = 1:4)
                       },bg.border =  NA)
          circos.track(track.index = 2, panel.fun = for(si in get.all.sector.index()[1:total.terms]) {
            xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 2)
            ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 2)
            circos.text(mean(xlim), mean(ylim), si, sector.index = si, track.index = 2,
                        facing = "clockwise",
                        niceFacing = TRUE,
                        adj = c(0.15,0.5),
                        cex = -0.004*vcount(link_tbl)+1.48, #<---------------------
                        font = 2,
                        col = "white")
          }, bg.border =  NA)
          circos.track(track.index = 2, panel.fun = for(si in get.all.sector.index()[start.entry:vcount(link_tbl)]) {
            xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 2)
            ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 2)
            circos.text(mean(xlim), mean(ylim), si, sector.index = si, track.index = 2,
                        facing = "clockwise",
                        niceFacing = TRUE,
                        adj = c(0.15,0.5),
                        cex = -0.005*vcount(link_tbl)+1.35, #<---------------------
                        font = 2,
                        col = "black")
          }, bg.border = NA)
        }
        #
        {
          png(file=paste0("Plot_Chord_2_",i,".png"),
              width = 22.5,
              height = 10,
              units = "in",
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
          chord.plot.1()
          upViewport()
          pushViewport(viewport(x = circle_size, y = 0.5, width = grobWidth(lgd_list_vertical),
                                height = grobHeight(lgd_list_vertical), just = c("left", "center")))
          grid.draw(lgd_list_vertical)
          upViewport()
          text(x=.65, y=.75, "ASPECT", adj = c(0,0), 
               cex = -0.004*vcount(link_tbl)+1.48+.2,font = 2) #<---------------------
          legend(x=.65,y=.7, 
                 legend=unique(paste0(nodes$Term[1:total.terms]," (",nodes$Freq[1:total.terms],")")),
                 cex=-0.004*vcount(link_tbl)+1.48,bg="transparent", #<---------------------
                 border = NA,pt.bg=many.colors,bty = "n",
                 x.intersp=1.3,y.intersp=-0.003*vcount(link_tbl)+1.21, #<---------------------
                 col=NA,text.font = 2,pch=22,
                 pt.cex=-0.005*vcount(link_tbl)+2.35)
          dev.off()
        }
        # >
      }
    } else if ((vcount(link_tbl) >= 46) && (vcount(link_tbl) <= 60)) {
      print("grafico 2")
      for (i in -1:0) {
        #
        # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
        # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
        #
        # |||||||||     2) Chord plots (without expression values)
        #
        # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
        # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
        #
        go.term = drop_na(nodes[,c("Entry","Term")])
        colnames(go.term)[1] <- "GO"
        #
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
        chord.plot.1 = function () {
          chordDiagram(go.entry.mat, order = nodes$Entry, # c(nodes$Term[1:total.terms],only.entrys)
                       grid.col = col.set,
                       transparency = 0,
                       directional = i,
                       annotationTrack = "grid",
                       #    c(% del datio del c?rculo, alto del track)
                       annotationTrackHeight = c(0.08, .01),
                       preAllocateTracks = 2)
          circos.track(track.index = 2,ylim = c(0, 1), track.height = uh(5, "cm"),
                       track.margin = c(5, uh(5, "cm")),
                       panel.fun = function(x, y) {
                         xcenter = get.cell.meta.data("xcenter")
                         circos.lines(c(xcenter, xcenter), c(0, uy(.5, "cm")),
                                      col = 1:4)
                       },bg.border =  NA)
          circos.track(track.index = 2, panel.fun = for(si in get.all.sector.index()[1:total.terms]) {
            xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 2)
            ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 2)
            circos.text(mean(xlim), mean(ylim), si, sector.index = si, track.index = 2,
                        facing = "clockwise",
                        niceFacing = TRUE,
                        adj = c(0.15,0.5),
                        cex = -0.004*vcount(link_tbl)+1.48, #<---------------------
                        font = 2,
                        col = "white")
          }, bg.border =  NA)
          circos.track(track.index = 2, panel.fun = for(si in get.all.sector.index()[start.entry:vcount(link_tbl)]) {
            xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 2)
            ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 2)
            circos.text(mean(xlim), mean(ylim), si, sector.index = si, track.index = 2,
                        facing = "clockwise",
                        niceFacing = TRUE,
                        adj = c(0.15,0.5),
                        cex = -0.005*vcount(link_tbl)+1.35, #<---------------------
                        font = 2,
                        col = "black")
          }, bg.border = NA)
        }
        #
        {
          png(file=paste0("Plot_Chord_2_",i,".png"),
              width = 22.5,
              height = 10,
              units = "in",
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
          chord.plot.1()
          upViewport()
          pushViewport(viewport(x = circle_size, y = 0.5, width = grobWidth(lgd_list_vertical),
                                height = grobHeight(lgd_list_vertical), just = c("left", "center")))
          grid.draw(lgd_list_vertical)
          upViewport()
          text(x=.65, y=.75, "ASPECT", adj = c(0,0), 
               cex = -0.004*vcount(link_tbl)+1.48+.2,font = 2) #<---------------------
          legend(x=.65,y=.7, 
                 legend=unique(paste0(nodes$Term[1:total.terms]," (",nodes$Freq[1:total.terms],")")),
                 cex=-0.004*vcount(link_tbl)+1.48,bg="transparent", #<---------------------
                 border = NA,pt.bg=many.colors,bty = "n",
                 x.intersp=1.3,y.intersp=-0.003*vcount(link_tbl)+1.21, #<---------------------
                 col=NA,text.font = 2,pch=22,
                 pt.cex=-0.005*vcount(link_tbl)+2.35)
          dev.off()
        }
        # >
      }
    } else if (vcount(link_tbl) >= 61) {
      print("grafico 3")
      for (i in -1:0) {
        #
        # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
        # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
        #
        # |||||||||     3) Chord plots (without expression values)
        #
        # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
        # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
        #
        go.term = drop_na(nodes[,c("Entry","Term")])
        colnames(go.term)[1] <- "GO"
        #
        go.entry.term = merge(links,go.term,by="GO")#[c(3,2)]
        ##
        go.entry = data.frame("GO"=go.entry.term$GO,"Entry"=go.entry.term$Entry)
        go.entry=as.data.frame.matrix(table(go.entry)) 
        go.entry.mat = as.matrix(go.entry) 
        #########
        circos.clear()
        circos.par(start.degree = 88, clock.wise = T)
        circos.par(track.margin = c(0, 0)) # 1
        #
        chord.plot.1 = function () {
          chordDiagram(go.entry.mat, order = nodes$Entry, # c(nodes$Term[1:total.terms],only.entrys)
                       grid.col = col.set,
                       transparency = 0,
                       directional = i,
                       annotationTrack = "grid",
                       #    c(% del datio del c?rculo, alto del track)
                       annotationTrackHeight = c(0.08, .01),
                       preAllocateTracks = 2)
          circos.track(track.index = 2,ylim = c(0, 1), track.height = uh(5, "cm"),
                       track.margin = c(5, uh(5, "cm")),
                       panel.fun = function(x, y) {
                         xcenter = get.cell.meta.data("xcenter")
                         circos.lines(c(xcenter, xcenter), c(0, uy(.5, "cm")),
                                      col = 1:4)
                       },bg.border =  NA)
          circos.track(track.index = 2, panel.fun = for(si in get.all.sector.index()[1:total.terms]) {
            xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 2)
            ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 2)
            circos.text(mean(xlim), mean(ylim), si, sector.index = si, track.index = 2,
                        facing = "clockwise",
                        niceFacing = TRUE,
                        adj = c(0.15,0.5),
                        cex = -0.004*vcount(link_tbl)+1.48, #<---------------------
                        font = 2,
                        col = "white")
          }, bg.border =  NA)
          circos.track(track.index = 2, panel.fun = for(si in get.all.sector.index()[start.entry:vcount(link_tbl)]) {
            xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 2)
            ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 2)
            circos.text(mean(xlim), mean(ylim), si, sector.index = si, track.index = 2,
                        facing = "clockwise",
                        niceFacing = TRUE,
                        adj = c(0.15,0.5),
                        cex = -0.005*vcount(link_tbl)+1.35, #<---------------------
                        font = 2,
                        col = "black")
          }, bg.border = NA)
        }
        #
        {
          png(file=paste0("Plot_Chord_3_",i,".png"),
              width = 22.5,
              height = 10,
              units = "in",
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
          chord.plot.1()
          upViewport()
          pushViewport(viewport(x = circle_size, y = 0.5, width = grobWidth(lgd_list_vertical),
                                height = grobHeight(lgd_list_vertical), just = c("left", "center")))
          grid.draw(lgd_list_vertical)
          upViewport()
          text(x=.65, y=.75, "ASPECT", adj = c(0,0), 
               cex = -0.004*vcount(link_tbl)+1.48+.2,font = 2) #<---------------------
          legend(x=.65,y=.7, 
                 legend=unique(paste0(nodes$Term[1:total.terms]," (",nodes$Freq[1:total.terms],")")),
                 cex=-0.004*vcount(link_tbl)+1.48,bg="transparent", #<---------------------
                 border = NA,pt.bg=many.colors,bty = "n",
                 x.intersp=1.3,y.intersp=-0.003*vcount(link_tbl)+1.21, #<---------------------
                 col=NA,text.font = 2,pch=22,
                 pt.cex=-0.005*vcount(link_tbl)+2.35)
          dev.off()
        }
        # >
      }
    } else {
      print("este mensaje no tiene que imprimirse")
    }
    # >
    # >
    #
    # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    #
    # |||||||||     Upset Plot (without expression values)
    #
    # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    #
    {
      val=length(only.entrys) # nrow number (total genes)
      #
      ups=as.data.frame.matrix(t(table(links)))
      ups[["log.Pval"]]=c(-log10(nodes$P[1:total.terms]),rep(NA,val-total.terms))
      ups[["term"]]=c(nodes$Entry[1:total.terms],rep(NA,(val-total.terms)))
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
        set.name.size=1.35
        #print(gnp.node.size)
      } else if ((length(unique(links$Entry)) >= 9.1) && (length(unique(links$Entry)) <= 13)) {
        set.name.size=1.2
        #print(gnp.node.size)
      } else if (length(unique(links$Entry)) >= 13.1) {
        set.name.size=0.9
        #print(gnp.node.size)
      }
      #
      png(file="zxcvbnmPlot_UpSetR.png",width = 10,height = 7,units = "in",res=700,bg="white")
      upset(ups,sets=nodes$Entry[1:total.terms],
            sets.bar.color = "darkgoldenrod2",order.by ="freq",empty.intersections = NULL,
            point.size=2,mainbar.y.label="poiuytrewq",sets.x.label = "",
            main.bar.color="black",matrix.color="black",shade.color="wheat3",
            line.size=0.5,show.numbers = "yes",group.by = "degree",
            matrix.dot.alpha = 0.5,mb.ratio = c(0.7, 0.3),
            #c(intersection size title, intersection size tick labels, set size title, set size tick labels, set names, numbers above bars)
            text.scale = c(1.5,1.5,0,1.5,set.name.size,1.5), 
            scale.sets = "identity",scale.intersections = "identity",
            shade.alpha = 0.35,set_size.angles = 0,
            attribute.plots = list(gridrows = 45, plots = list(list(plot = myplot2,x="term",y = "log.Pval", queries = F)),ncols=2))
      dev.off() 
    }
    #>
  }
}
#
