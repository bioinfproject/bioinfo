#
install.packages("xml2", lib = "../R" , repos='http://cran.us.r-project.org')
install.packages("curl", lib = "../R" , repos='http://cran.us.r-project.org')
install.packages("httr", lib = "../R" , repos='http://cran.us.r-project.org')
#
install.packages("GetoptLong", lib = "../R" , repos='http://cran.us.r-project.org')
install.packages("tidyverse", lib = "../R" , repos='http://cran.us.r-project.org')
install.packages("tidygraph", lib = "../R" , repos='http://cran.us.r-project.org')
install.packages("ggraph", lib = "../R" , repos='http://cran.us.r-project.org')
install.packages("viridis", lib = "../R" , repos='http://cran.us.r-project.org')
install.packages("circlize", lib = "../R" , repos='http://cran.us.r-project.org')
install.packages("RColorBrewer", lib = "../R" , repos='http://cran.us.r-project.org')
install.packages("cowplot", lib = "../R" , repos='http://cran.us.r-project.org')
install.packages("networkD3", lib = "../R" , repos='http://cran.us.r-project.org')
install.packages("UpSetR", lib = "../R" , repos='http://cran.us.r-project.org')
install.packages("gridBase",lib = "../R" , repos='http://cran.us.r-project.org')
install.packages("reshape2",lib = "../R" , repos='http://cran.us.r-project.org')
source("http://bioconductor.org/biocLite.R")
biocLite("ComplexHeatmap")
#
cat("\n\n!!! All packages have been installed successfully!!!\n\n")

