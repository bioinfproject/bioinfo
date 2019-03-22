#
install.packages("xml2", destdir = "./R" , repos='http://cran.us.r-project.org')
install.packages("curl", destdir = "./R" , repos='http://cran.us.r-project.org')
install.packages("httr", destdir = "./R" , repos='http://cran.us.r-project.org')
#
install.packages("GetoptLong", destdir = "./R" , repos='http://cran.us.r-project.org')
install.packages("tidyverse", destdir = "./R" , repos='http://cran.us.r-project.org')
install.packages("tidygraph", destdir = "./R" , repos='http://cran.us.r-project.org')
install.packages("ggraph", destdir = "./R" , repos='http://cran.us.r-project.org')
install.packages("viridis", destdir = "./R" , repos='http://cran.us.r-project.org')
install.packages("circlize", destdir = "./R" , repos='http://cran.us.r-project.org')
install.packages("RColorBrewer", destdir = "./R" , repos='http://cran.us.r-project.org')
install.packages("cowplot", destdir = "./R" , repos='http://cran.us.r-project.org')
install.packages("networkD3", destdir = "./R" , repos='http://cran.us.r-project.org')
install.packages("UpSetR", destdir = "./R" , repos='http://cran.us.r-project.org')
install.packages("gridBase",repos='http://cran.us.r-project.org')
install.packages("reshape2",repos='http://cran.us.r-project.org')
source("http://bioconductor.org/biocLite.R")
biocLite("ComplexHeatmap")
#
cat("\n\n!!! All packages have been installed successfully!!!\n\n")

