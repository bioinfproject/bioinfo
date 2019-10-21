
.libPaths('Rlibrary_NeVOmics')
#
install.packages("tidyverse", dependencies = TRUE,repos='http://cran.us.r-project.org')
install.packages("circlize", dependencies = TRUE,repos='http://cran.us.r-project.org')
install.packages("RColorBrewer", dependencies = TRUE,repos='http://cran.us.r-project.org')
install.packages("gridBase",dependencies = TRUE,repos='http://cran.us.r-project.org')
source("http://bioconductor.org/biocLite.R")
biocLite("ComplexHeatmap")
#
