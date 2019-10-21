
.libPaths('Rlibrary_NeVOmics')
#
print("GetoptLong")
install.packages("GetoptLong", dependencies = TRUE,repos='http://cran.us.r-project.org')
print("tidyverse")
install.packages("tidyverse", dependencies = TRUE,repos='http://cran.us.r-project.org')
print("ComplexHeatmap")
source("http://bioconductor.org/biocLite.R")
biocLite("ComplexHeatmap")
print("circlize")
install.packages("circlize", dependencies = TRUE,repos='http://cran.us.r-project.org')
#install.packages("RColorBrewer", dependencies = TRUE,repos='http://cran.us.r-project.org')
#install.packages("gridBase",dependencies = TRUE,repos='http://cran.us.r-project.org')

#
