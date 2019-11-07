
{
  library(tidyverse)
  library(circlize)
  library(RColorBrewer)
  library(grid)
  library(ComplexHeatmap)
  library(gridBase)
}


nodes = read_csv(list.files(pattern = "nodes"))

links = read_csv(list.files(pattern = "edges"))

# la configuracion de colores viene desde python
barcolortitle = links$bar_title[!is.na(links$bar_title)]

links$bar_title <- NULL # extraigo la informacion de la columna "bar_title" y la elimino
#----------------------------------------------

colors.for.bar.rasterimage = links$bar_color_R[!is.na(links$bar_color_R)]

colors.for.bar.rasterimage = rev(colors.for.bar.rasterimage)

links$bar_color_R <- NULL # extraigo la informacion de la columna "bar_color_R" y la elimino
#-----------------------------------------------
links = drop_na(links)

columnas.links = colnames(links)


if (columnas.links[1] == "Path"){
  title.legend = ""#"KEGG Pathways"
} else if (columnas.links[1] == "GObp"){
  title.legend = ""#"Biological Process"
} else if (columnas.links[1] == "GOmf"){
  title.legend = ""#"Molecular Function"
} else if (columnas.links[1] == "GOcc"){
  title.legend = ""#"Cellular Component"
}




{
  # cambio del nombre de las columnas
  names(links)[4] = "source"
  names(links)[2] = "target"
  
  total.terms = nrow(nodes)
  start.entry = total.terms + 1
  
  entry.colors.unicos = unique(links[c("target", "entry_colors")])
  
  only.entrys = entry.colors.unicos$target
  id = rep("",total.terms)
  nodes.label = c(id,only.entrys)
  
  # colores ya definidos desde python
  prot.col = entry.colors.unicos$entry_colors
  
  col.set.link=   c(nodes$term_colors, prot.col)
  
  
  ids.genes = links[,c("source", "target")]
  # crear una matriz, dentro estaun dataframe = 
  ids.genes.mat = as.matrix(as.data.frame.matrix(table(ids.genes)))
}

# get information for index
indice = rownames(ids.genes.mat)
indice2 = data.frame("Term"=indice)
indice3 = merge(indice2,nodes,by="Term")

# change index
term.total_counts = paste0(indice3$Term, " (", indice3$list_count, ")")
rownames(ids.genes.mat)  = term.total_counts


{
  gap_term_genes = 7
  gap_temrs = 2
  gap_genes = 0.5
  terminos = nrow(ids.genes.mat)
  genes = ncol(ids.genes.mat)
  suma_espacios_terms = (terminos - 1) * gap_temrs
  suma_espacios_genes = (genes - 1) * gap_genes
  total.puntos = (gap_term_genes * 2) + terminos + suma_espacios_terms + genes + suma_espacios_genes
  partes = 360 / total.puntos
  seccion.terms = partes * (terminos + suma_espacios_terms - 7) # los 7 son los de un lado, gap
  angulo = seccion.terms / 2
  
}


total.nodes = terminos + genes


if (total.terms <= 15) {
  cex.terms = 0.45
} else if ((total.terms >= 16) && (total.terms <= 25)) {
  cex.terms = 0.35
} else if (total.terms >= 26) {
  cex.terms = 0.25
}

if (total.terms <= 20) {
  cex.legend.box = 0.55
  punto.legend.box = 1
} else if (total.terms >= 21) {
  cex.legend.box = 0.4
  punto.legend.box = 0.7
}

if ((-0.003039*length(only.entrys) + 0.47) < 0) {
  cex.genes = 0.1
} else if ((-0.003039*length(only.entrys) + 0.47) > 0) {
  cex.genes = -0.003039*length(only.entrys) + 0.47
}



if (length(unique(links$entry_colors)) > 1) {
  # datos con valores asociados a los genes/proteinas
  print('hay mas de un color, crear graficos con rampa')
  
  rbPal <- colorRampPalette(colors.for.bar.rasterimage)
  legend_image <- as.raster(matrix(rbPal(length(only.entrys)), ncol=1))
  
  # aqui agregar el cofigo de R para los graficos
  
  circos.clear()
  # gap entre cada cuerda
  circos.par(gap.after = c(rep(2, nrow(ids.genes.mat)-1), gap_term_genes,
                           rep(0.5, ncol(ids.genes.mat)-1), gap_term_genes),
             start.degree = angulo, clock.wise = T)
  
  circos.par(track.margin = c(0,0))
  
  
  
  #------------------------------------------
  # posicion de la barra
  ## bar
  xleft = -1
  ybottom = 0.375
  xright = -1.075
  ytop = 0.825
  #------------------------------------------
  
  # grafico 1
  
  for (i in -1:0){
    chord.plot.1 = function () {
      chordDiagram(ids.genes.mat, order = c(paste0(nodes$Term, " (", nodes$list_count, ")"), only.entrys),
                   grid.col = col.set.link, scale = TRUE,
                   transparency = 0.25,
                   directional = i,
                   annotationTrack = NULL,
                   #    c(% del datio del circulo, alto del track)
                   annotationTrackHeight = c(0.1, .01),
                   preAllocateTracks = 2)
      circos.track(track.index = 2,ylim = c(0, 1), track.height = uh(5, "cm"),
                   track.margin = c(0, uh(0, "cm")),
                   panel.fun = function(x, y) {
                     xcenter = get.cell.meta.data("xcenter")
                     circos.lines(c(xcenter, xcenter), c(0, 0.3),
                                  lwd = .3, col = "grey", type = "-")
                   },bg.border =  NA)
      circos.track(track.index = 2, panel.fun = for(si in get.all.sector.index()[1:total.terms]) {
        xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 2)
        ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 2)
        circos.text(mean(xlim), mean(ylim), si, sector.index = si, track.index = 2,
                    facing = "clockwise",
                    niceFacing = TRUE,
                    adj = c(0,0.5),
                    cex = cex.terms, 
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
                    cex = cex.genes, #<---------------------
                    font = 2,
                    col = "black")
      }, bg.border = NA)
    }
    
    {
      png(file=paste0("job_KEGG_plots/NeVOmics_Plot_Chord_NULL_clockwis_",i,".png"),
          width = 30,
          height = 10,
          units = "cm",
          res=1024,bg="white")
      chord.plot.1()
      text(x=.4, y=.9, title.legend, adj = c(0,0), 
           cex = .8,font = 2)

      dev.off()}
  }
  
  # grafico 1.1
  
  for (i in -1:0){
    chord.plot.1 = function () {
      chordDiagram(ids.genes.mat, order = c(paste0(nodes$Term, " (", nodes$list_count, ")"), only.entrys),
                   grid.col = col.set.link, scale = TRUE,
                   transparency = 0.25,
                   directional = i,
                   annotationTrack = NULL,
                   #    c(% del datio del circulo, alto del track)
                   annotationTrackHeight = c(0.1, .01),
                   preAllocateTracks = 2)
      circos.track(track.index = 2,ylim = c(0, 1), track.height = uh(5, "cm"),
                   track.margin = c(0, uh(0, "cm")),
                   panel.fun = function(x, y) {
                     xcenter = get.cell.meta.data("xcenter")
                     circos.lines(c(xcenter, xcenter), c(0, 0.3),
                                  lwd = .3, col = "grey", type = "-")
                   },bg.border =  NA)
      circos.track(track.index = 2, panel.fun = for(si in get.all.sector.index()[1:total.terms]) {
        xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 2)
        ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 2)
        circos.text(mean(xlim), mean(ylim), si, sector.index = si, track.index = 2,
                    facing = "downward",
                    niceFacing = TRUE,
                    adj = c(0,0.5),
                    cex = cex.terms, 
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
                    cex = cex.genes, #<---------------------
                    font = 2,
                    col = "black")
      }, bg.border = NA)
    }
    
    {
      png(file=paste0("job_KEGG_plots/NeVOmics_Plot_Chord_NULL_downward_",i,".png"),
          width = 30,
          height = 10,
          units = "cm",
          res=1024,bg="white")
      chord.plot.1()
      text(x=.4, y=.9, title.legend, adj = c(0,0), 
           cex = .8,font = 2)
      
      dev.off()}
  }
  
  #............................................................................
  
  # gr\E1fico 2
  
  for (i in -1:0){
    chord.plot.1 = function () {
      chordDiagram(ids.genes.mat, order = c(paste0(nodes$Term, " (", nodes$list_count, ")"), only.entrys),
                   grid.col = col.set.link, scale = TRUE,
                   transparency = 0.25,
                   directional = i,
                   annotationTrack = "grid",
                   #    c(% del datio del circulo, alto del track)
                   annotationTrackHeight = c(0.1, .01),
                   preAllocateTracks = 2)
      circos.track(track.index = 2,ylim = c(0, 1), track.height = uh(5, "cm"),
                   track.margin = c(0, uh(0, "cm")),
                   panel.fun = function(x, y) {
                     xcenter = get.cell.meta.data("xcenter")
                     circos.lines(c(xcenter, xcenter), c(0, 0.3),
                                  lwd = .3, col = "grey", type = "-")
                   },bg.border =  NA)
      circos.track(track.index = 2, panel.fun = for(si in get.all.sector.index()[1:total.terms]) {
        xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 2)
        ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 2)
        circos.text(mean(xlim), mean(ylim), si, sector.index = si, track.index = 2,
                    facing = "clockwise",
                    niceFacing = TRUE,
                    adj = c(0,0.5),
                    cex = cex.terms, 
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
                    cex = cex.genes, #<---------------------
                    font = 2,
                    col = "black")
      }, bg.border = NA)
    }
    
    {
      png(file=paste0("job_KEGG_plots/NeVOmics_Plot_Chord_grid_clockwise",i,".png"),
          width = 30,
          height = 10,
          units = "cm",
          res=1024,bg="white")
      chord.plot.1()
      text(x=.4, y=.9, title.legend, adj = c(0,0), 
           cex = .8,font = 2)
      ## bar
      rasterImage(legend_image, xleft, ybottom, xright, ytop)
      text(x=-1.04, y = seq(ybottom, ytop,l=5),
           labels = round(seq(min(links$values), max(links$values), length.out = 5),digits=1),font = 1,
           pos = 4, cex= 0.37) #<---------------------
      text(x=-1.04,y= .83, barcolortitle, adj = c(0,0) ,font = 2, pos = 3,
           cex = .5) #<---------------------
      
      dev.off()}
  }
  
  # gr\E1fico 2.1
  
  for (i in -1:0){
    chord.plot.1 = function () {
      chordDiagram(ids.genes.mat, order = c(paste0(nodes$Term, " (", nodes$list_count, ")"), only.entrys),
                   grid.col = col.set.link, scale = TRUE,
                   transparency = 0.25,
                   directional = i,
                   annotationTrack = "grid",
                   #    c(% del datio del circulo, alto del track)
                   annotationTrackHeight = c(0.1, .01),
                   preAllocateTracks = 2)
      circos.track(track.index = 2,ylim = c(0, 1), track.height = uh(5, "cm"),
                   track.margin = c(0, uh(0, "cm")),
                   panel.fun = function(x, y) {
                     xcenter = get.cell.meta.data("xcenter")
                     circos.lines(c(xcenter, xcenter), c(0, 0.3),
                                  lwd = .3, col = "grey", type = "-")
                   },bg.border =  NA)
      circos.track(track.index = 2, panel.fun = for(si in get.all.sector.index()[1:total.terms]) {
        xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 2)
        ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 2)
        circos.text(mean(xlim), mean(ylim), si, sector.index = si, track.index = 2,
                    facing = "downward",
                    niceFacing = TRUE,
                    adj = c(0,0.5),
                    cex = cex.terms, 
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
                    cex = cex.genes, #<---------------------
                    font = 2,
                    col = "black")
      }, bg.border = NA)
    }
    
    {
      png(file=paste0("job_KEGG_plots/NeVOmics_Plot_Chord_grid_downward",i,".png"),
          width = 30,
          height = 10,
          units = "cm",
          res=1024,bg="white")
      chord.plot.1()
      text(x=.4, y=.9, title.legend, adj = c(0,0), 
           cex = .8,font = 2)
      ## bar
      rasterImage(legend_image, xleft, ybottom, xright, ytop)
      text(x=-1.04, y = seq(ybottom, ytop,l=5),
           labels = round(seq(min(links$values), max(links$values), length.out = 5),digits=1),font = 1,
           pos = 4, cex= 0.37) #<---------------------
      text(x=-1.04,y= .83, barcolortitle, adj = c(0,0) ,font = 2, pos = 3,
           cex = .5) #<---------------------
      
      dev.off()}
  }
  
  #......................................................................................
  #......................................................................................
  #......................................................................................
  #......................................................................................
  #......................................................................................
  #......................................................................................
  #......................................................................................
  #......................................................................................
  #......................................................................................
  
  # grafico con box legend
  
  rownames(ids.genes.mat)  = indice3$base
  
  #grafico 1 
  
  for (i in -1:0){
    chord.plot.1 = function () {
      chordDiagram(ids.genes.mat, order = c(nodes$base, only.entrys),
                   grid.col = col.set.link, scale = TRUE,
                   transparency = 0.25,
                   directional = i,
                   annotationTrack = "grid",
                   #    c(% del datio del circulo, alto del track)
                   annotationTrackHeight = c(0.1, .01),
                   preAllocateTracks = 2)
      circos.track(track.index = 2,ylim = c(0, 1), track.height = uh(5, "cm"),
                   track.margin = c(0, uh(0, "cm")),
                   panel.fun = function(x, y) {
                     xcenter = get.cell.meta.data("xcenter")
                     circos.lines(c(xcenter, xcenter), c(0, 0.3),
                                  lwd = .3, col = "grey", type = "-")
                   },bg.border =  NA)
      circos.track(track.index = 2, panel.fun = for(si in get.all.sector.index()[1:total.terms]) {
        xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 2)
        ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 2)
        circos.text(mean(xlim), mean(ylim), si, sector.index = si, track.index = 2,
                    facing = "clockwise",
                    niceFacing = TRUE,
                    adj = c(0,0.5),
                    cex = cex.terms, 
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
                    cex = cex.genes, #<---------------------
                    font = 2,
                    col = "black")
      }, bg.border = NA)
    }
    
    {
      png(file=paste0("job_KEGG_plots/NeVOmics_Plot_Chord_Legend_grid_",i,".png"),
          width = 30,
          height = 10,
          units = "cm",
          res=900,bg="white")
      chord.plot.1()
      text(x=1, y=0.8, title.legend, adj = c(0,0), 
           cex = .8,font = 2) #<---------------------
      legend(x=1.1,y=.75, 
             legend=paste0(nodes$Term," (",nodes$list_count,")"),
             cex=cex.legend.box,bg="white", #<---------------------
             border = NA,pt.bg=nodes$term_colors,bty = "n",
             x.intersp=.8,y.intersp=1, #<---------------------
             col=NA,text.font = 1,pch=21,
             pt.cex=punto.legend.box)
      ## bar
      rasterImage(legend_image, xleft, ybottom, xright, ytop)
      text(x=-1.04, y = seq(ybottom, ytop,l=5),
           labels = round(seq(min(links$values), max(links$values), length.out = 5),digits=1),font = 1,
           pos = 4, cex= 0.37) #<---------------------
      text(x=-1.04,y= .83, barcolortitle, adj = c(0,0) ,font = 2, pos = 3,
           cex = .5) #<---------------------
      
      dev.off()
    }
  }
  
  # grafico 2
  
  for (i in -1:0){
    chord.plot.1 = function () {
      chordDiagram(ids.genes.mat, order = c(nodes$base, only.entrys),
                   grid.col = col.set.link, scale = TRUE,
                   transparency = 0.25,
                   directional = i,
                   annotationTrack = NULL,
                   #    c(% del datio del circulo, alto del track)
                   annotationTrackHeight = c(0.1, .01),
                   preAllocateTracks = 2)
      circos.track(track.index = 2,ylim = c(0, 1), track.height = uh(5, "cm"),
                   track.margin = c(0, uh(0, "cm")),
                   panel.fun = function(x, y) {
                     xcenter = get.cell.meta.data("xcenter")
                     circos.lines(c(xcenter, xcenter), c(0, 0.3),
                                  lwd = .3, col = "grey", type = "-")
                   },bg.border =  NA)
      circos.track(track.index = 2, panel.fun = for(si in get.all.sector.index()[1:total.terms]) {
        xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 2)
        ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 2)
        circos.text(mean(xlim), mean(ylim), si, sector.index = si, track.index = 2,
                    facing = "clockwise",
                    niceFacing = TRUE,
                    adj = c(0,0.5),
                    cex = cex.terms, 
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
                    cex = cex.genes, #<---------------------
                    font = 2,
                    col = "black")
      }, bg.border = NA)
    }
    {
      png(file=paste0("job_KEGG_plots/NeVOmics_Plot_Chord_Legend_NULL_",i,".png"),
          width = 30,
          height = 10,
          units = "cm",
          res=900,bg="white")
      chord.plot.1()
      text(x=1, y=0.8, title.legend, adj = c(0,0), 
           cex = .8,font = 2) #<---------------------
      legend(x=1.1,y=.75, 
             legend=paste0(nodes$Term," (",nodes$list_count,")"),
             cex=cex.legend.box,bg="white", #<---------------------
             border = NA,pt.bg=nodes$term_colors,bty = "n",
             x.intersp=.8,y.intersp=1, #<---------------------
             col=NA,text.font = 1,pch=21,
             pt.cex=punto.legend.box)
      
      dev.off()
    }
    
  }


} else if (length(unique(links$entry_colors)) == 1) {
  # datos sin valores asociados a los genes/proteinas
  print('solo hay un color, quiere decir que no ingresaron valores, no crear la rampa')
  
  # aqui agregar el cofigo de R para los graficos
  
  
  circos.clear()
  # gap entre cada cuerda
  circos.par(gap.after = c(rep(2, nrow(ids.genes.mat)-1), 10, rep(0.5, ncol(ids.genes.mat)-1), 10),
             start.degree = angulo, clock.wise = T)
  
  circos.par(track.margin = c(0,0))
  
  
  # grafico 1
  
  for (i in -1:0){
    chord.plot.1 = function () {
      chordDiagram(ids.genes.mat, order = c(paste0(nodes$Term, " (", nodes$list_count, ")"), only.entrys),
                   grid.col = col.set.link, scale = TRUE,
                   transparency = 0.25,
                   directional = i,
                   annotationTrack = NULL,
                   #    c(% del datio del circulo, alto del track)
                   annotationTrackHeight = c(0.1, .01),
                   preAllocateTracks = 2)
      circos.track(track.index = 2,ylim = c(0, 1), track.height = uh(5, "cm"),
                   track.margin = c(0, uh(0, "cm")),
                   panel.fun = function(x, y) {
                     xcenter = get.cell.meta.data("xcenter")
                     circos.lines(c(xcenter, xcenter), c(0, 0.3),
                                  lwd = .3, col = "grey", type = "-")
                   },bg.border =  NA)
      circos.track(track.index = 2, panel.fun = for(si in get.all.sector.index()[1:total.terms]) {
        xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 2)
        ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 2)
        circos.text(mean(xlim), mean(ylim), si, sector.index = si, track.index = 2,
                    facing = "clockwise",
                    niceFacing = TRUE,
                    adj = c(0,0.5),
                    cex = cex.terms, 
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
                    cex = cex.genes, #<---------------------
                    font = 2,
                    col = "black")
      }, bg.border = NA)
    }
    
    {
      png(file=paste0("job_KEGG_plots/NeVOmics_Plot_Chord_NULL_clockwis_",i,".png"),
          width = 30,
          height = 10,
          units = "cm",
          res=1024,bg="white")
      chord.plot.1()
      text(x=.4, y=.9, title.legend, adj = c(0,0), 
           cex = .8,font = 2)
      dev.off()}
  }
  
  # grafico 1.1
  
  for (i in -1:0){
    chord.plot.1 = function () {
      chordDiagram(ids.genes.mat, order = c(paste0(nodes$Term, " (", nodes$list_count, ")"), only.entrys),
                   grid.col = col.set.link, scale = TRUE,
                   transparency = 0.25,
                   directional = i,
                   annotationTrack = NULL,
                   #    c(% del datio del circulo, alto del track)
                   annotationTrackHeight = c(0.1, .01),
                   preAllocateTracks = 2)
      circos.track(track.index = 2,ylim = c(0, 1), track.height = uh(5, "cm"),
                   track.margin = c(0, uh(0, "cm")),
                   panel.fun = function(x, y) {
                     xcenter = get.cell.meta.data("xcenter")
                     circos.lines(c(xcenter, xcenter), c(0, 0.3),
                                  lwd = .3, col = "grey", type = "-")
                   },bg.border =  NA)
      circos.track(track.index = 2, panel.fun = for(si in get.all.sector.index()[1:total.terms]) {
        xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 2)
        ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 2)
        circos.text(mean(xlim), mean(ylim), si, sector.index = si, track.index = 2,
                    facing = "downward",
                    niceFacing = TRUE,
                    adj = c(0,0.5),
                    cex = cex.terms, 
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
                    cex = cex.genes, #<---------------------
                    font = 2,
                    col = "black")
      }, bg.border = NA)
    }
    
    {
      png(file=paste0("job_KEGG_plots/NeVOmics_Plot_Chord_NULL_downward_",i,".png"),
          width = 30,
          height = 10,
          units = "cm",
          res=1024,bg="white")
      chord.plot.1()
      text(x=.4, y=.9, title.legend, adj = c(0,0), 
           cex = .8,font = 2)
      dev.off()}
  }
  
  #............................................................................
  
  # gr\E1fico 2
  
  for (i in -1:0){
    chord.plot.1 = function () {
      chordDiagram(ids.genes.mat, order = c(paste0(nodes$Term, " (", nodes$list_count, ")"), only.entrys),
                   grid.col = col.set.link, scale = TRUE,
                   transparency = 0.25,
                   directional = i,
                   annotationTrack = "grid",
                   #    c(% del datio del circulo, alto del track)
                   annotationTrackHeight = c(0.1, .01),
                   preAllocateTracks = 2)
      circos.track(track.index = 2,ylim = c(0, 1), track.height = uh(5, "cm"),
                   track.margin = c(0, uh(0, "cm")),
                   panel.fun = function(x, y) {
                     xcenter = get.cell.meta.data("xcenter")
                     circos.lines(c(xcenter, xcenter), c(0, 0.3),
                                  lwd = .3, col = "grey", type = "-")
                   },bg.border =  NA)
      circos.track(track.index = 2, panel.fun = for(si in get.all.sector.index()[1:total.terms]) {
        xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 2)
        ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 2)
        circos.text(mean(xlim), mean(ylim), si, sector.index = si, track.index = 2,
                    facing = "clockwise",
                    niceFacing = TRUE,
                    adj = c(0,0.5),
                    cex = cex.terms, 
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
                    cex = cex.genes, #<---------------------
                    font = 2,
                    col = "black")
      }, bg.border = NA)
    }
    
    {
      png(file=paste0("job_KEGG_plots/NeVOmics_Plot_Chord_grid_clockwise",i,".png"),
          width = 30,
          height = 10,
          units = "cm",
          res=1024,bg="white")
      chord.plot.1()
      text(x=.4, y=.9, title.legend, adj = c(0,0), 
           cex = .8,font = 2)
      dev.off()}
  }
  
  # gr\E1fico 2.1
  
  for (i in -1:0){
    chord.plot.1 = function () {
      chordDiagram(ids.genes.mat, order = c(paste0(nodes$Term, " (", nodes$list_count, ")"), only.entrys),
                   grid.col = col.set.link, scale = TRUE,
                   transparency = 0.25,
                   directional = i,
                   annotationTrack = "grid",
                   #    c(% del datio del circulo, alto del track)
                   annotationTrackHeight = c(0.1, .01),
                   preAllocateTracks = 2)
      circos.track(track.index = 2,ylim = c(0, 1), track.height = uh(5, "cm"),
                   track.margin = c(0, uh(0, "cm")),
                   panel.fun = function(x, y) {
                     xcenter = get.cell.meta.data("xcenter")
                     circos.lines(c(xcenter, xcenter), c(0, 0.3),
                                  lwd = .3, col = "grey", type = "-")
                   },bg.border =  NA)
      circos.track(track.index = 2, panel.fun = for(si in get.all.sector.index()[1:total.terms]) {
        xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 2)
        ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 2)
        circos.text(mean(xlim), mean(ylim), si, sector.index = si, track.index = 2,
                    facing = "downward",
                    niceFacing = TRUE,
                    adj = c(0,0.5),
                    cex = cex.terms, 
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
                    cex = cex.genes, #<---------------------
                    font = 2,
                    col = "black")
      }, bg.border = NA)
    }
    
    {
      png(file=paste0("job_KEGG_plots/NeVOmics_Plot_Chord_grid_downward",i,".png"),
          width = 30,
          height = 10,
          units = "cm",
          res=1024,bg="white")
      chord.plot.1()
      text(x=.4, y=.9, title.legend, adj = c(0,0), 
           cex = .8,font = 2)
      dev.off()}
  }
  
  #......................................................................................
  #......................................................................................
  #......................................................................................
  #......................................................................................
  #......................................................................................
  #......................................................................................
  #......................................................................................
  #......................................................................................
  #......................................................................................
  
  # grafico con box legend
  
  rownames(ids.genes.mat)  = indice3$base
  
  #grafico 1 
  
  for (i in -1:0){
    chord.plot.1 = function () {
      chordDiagram(ids.genes.mat, order = c(nodes$base, only.entrys),
                   grid.col = col.set.link, scale = TRUE,
                   transparency = 0.25,
                   directional = i,
                   annotationTrack = "grid",
                   #    c(% del datio del circulo, alto del track)
                   annotationTrackHeight = c(0.1, .01),
                   preAllocateTracks = 2)
      circos.track(track.index = 2,ylim = c(0, 1), track.height = uh(5, "cm"),
                   track.margin = c(0, uh(0, "cm")),
                   panel.fun = function(x, y) {
                     xcenter = get.cell.meta.data("xcenter")
                     circos.lines(c(xcenter, xcenter), c(0, 0.3),
                                  lwd = .3, col = "grey", type = "-")
                   },bg.border =  NA)
      circos.track(track.index = 2, panel.fun = for(si in get.all.sector.index()[1:total.terms]) {
        xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 2)
        ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 2)
        circos.text(mean(xlim), mean(ylim), si, sector.index = si, track.index = 2,
                    facing = "clockwise",
                    niceFacing = TRUE,
                    adj = c(0,0.5),
                    cex = cex.terms, 
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
                    cex = cex.genes, #<---------------------
                    font = 2,
                    col = "black")
      }, bg.border = NA)
    }
    
    {
      png(file=paste0("job_KEGG_plots/NeVOmics_Plot_Chord_Legend_grid_",i,".png"),
          width = 30,
          height = 10,
          units = "cm",
          res=900,bg="white")
      chord.plot.1()
      text(x=1, y=0.8, title.legend, adj = c(0,0), 
           cex = .8,font = 2) #<---------------------
      legend(x=1.1,y=.75, 
             legend=paste0(nodes$Term," (",nodes$list_count,")"),
             cex=cex.legend.box,bg="white", #<---------------------
             border = NA,pt.bg=nodes$term_colors,bty = "n",
             x.intersp=.8,y.intersp=1, #<---------------------
             col=NA,text.font = 1,pch=21,
             pt.cex=punto.legend.box)
      dev.off()
    }
  }
  
  # grafico 2
  
  for (i in -1:0){
    chord.plot.1 = function () {
      chordDiagram(ids.genes.mat, order = c(nodes$base, only.entrys),
                   grid.col = col.set.link, scale = TRUE,
                   transparency = 0.25,
                   directional = i,
                   annotationTrack = NULL,
                   #    c(% del datio del circulo, alto del track)
                   annotationTrackHeight = c(0.1, .01),
                   preAllocateTracks = 2)
      circos.track(track.index = 2,ylim = c(0, 1), track.height = uh(5, "cm"),
                   track.margin = c(0, uh(0, "cm")),
                   panel.fun = function(x, y) {
                     xcenter = get.cell.meta.data("xcenter")
                     circos.lines(c(xcenter, xcenter), c(0, 0.3),
                                  lwd = .3, col = "grey", type = "-")
                   },bg.border =  NA)
      circos.track(track.index = 2, panel.fun = for(si in get.all.sector.index()[1:total.terms]) {
        xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 2)
        ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 2)
        circos.text(mean(xlim), mean(ylim), si, sector.index = si, track.index = 2,
                    facing = "clockwise",
                    niceFacing = TRUE,
                    adj = c(0,0.5),
                    cex = cex.terms, 
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
                    cex = cex.genes, #<---------------------
                    font = 2,
                    col = "black")
      }, bg.border = NA)
    }
    {
      png(file=paste0("job_KEGG_plots/NeVOmics_Plot_Chord_Legend_NULL_",i,".png"),
          width = 30,
          height = 10,
          units = "cm",
          res=900,bg="white")
      chord.plot.1()
      text(x=1, y=0.8, title.legend, adj = c(0,0), 
           cex = .8,font = 2) #<---------------------
      legend(x=1.1,y=.75, 
             legend=paste0(nodes$Term," (",nodes$list_count,")"),
             cex=cex.legend.box,bg="white", #<---------------------
             border = NA,pt.bg=nodes$term_colors,bty = "n",
             x.intersp=.8,y.intersp=1, #<---------------------
             col=NA,text.font = 1,pch=21,
             pt.cex=punto.legend.box)
      dev.off()
    }
    
  }
  
}



