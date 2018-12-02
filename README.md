# <center> <h1>NeVOmics</h1> </center>

## **Description**
NeVOmics is a functional enrichment analysis tool developed in programming language Python and R that integrates Over-representation analysis (ORA) methodology and network-based visualization. It applies appropriate statistical methods to identify significantly enriched Gene Ontology (GO) terms or pathways in a given list of genes/proteins. It provides several types of graphical visualization to show enrichment results. NeVOmics supports all organisms deposited in UniProt Knowledgebase (UniProtKB) and Kyoto Encyclopedia of Genes and Genomes (KEGG) databases.
# <img src="https://raw.githubusercontent.com/bioinfproject/bioinfo/master/Folder/Figure_1.png" width = 75%>

### **Citation:**<br> **NeVOmics**: An Enrichment Tool for Gene Ontology and Functional Network Analysis and Visualization of Data from OMICs Technologies. [doi.org/10.3390/genes9120569](https://doi.org/10.3390/genes9120569)

### Compatibility with: **Windows** <img src="https://upload.wikimedia.org/wikipedia/sr/thumb/1/14/Windows_logo_-_2006.svg/644px-Windows_logo_-_2006.svg.png" width = 5%> and **Linux** <img src="https://upload.wikimedia.org/wikipedia/commons/thumb/b/b0/NewTux.svg/300px-NewTux.svg.png" width = 5%>
<hr />

### <center> <h1>Run NeVOmics</h1> </center>
# <img src="https://raw.githubusercontent.com/bioinfproject/bioinfo/master/Folder/Before%20run%20NeVOmics.png" width = 85%>
**1**. [**Download NeVOmics**](https://github.com/bioinfproject/bioinfo/blob/master/NeVOmics.zip?raw=true)

**2**. Unzip **NeVOmics.zip**

**3**. Open your terminal window:<br>
**Windows**: press the **Win+R** keys on your keyboard. Then, type **cmd** or **cmd.exe** and press **Enter** or click/tap **OK**.
**Linux**: press **Ctrl+Alt+T** on your keyboard.<br>

**4**. Change directory in **Windows** and **Linux**:
>In Windows
```bash
cd Downloads\NeVOmics
```
>In Linux
```bash
cd Downloads/NeVOmics
```
**5**. Type the next command to run NeVOmics:
```bash 
python NeVOmics.py
```
> - Or double-click the `NeVOmics.exe` file.<br>
**6**. Select one of the three available analyzes:

# <img src="https://raw.githubusercontent.com/bioinfproject/bioinfo/master/Folder/NeVOmics_intro_Windows.png" width = 50%>

**7**. Upload protein list (must be UniProtKB identifiers):

# <img src="https://raw.githubusercontent.com/bioinfproject/bioinfo/master/Folder/Upload_NeVOmics.PNG" width = 50%>

**8**. Selection of a method p-value correction:
```bash
[ Step 2: Choose a correction method (e.g., FDR / Bonferroni) ]
=====> : FDR    # type any method
```
**9**. Choose a significance value for the analysis:
```bash
[ Step 3: Choose a Value (e.g., 0.05) ]
=====> : 0.05   # type a significance value
```
**10**. Optional creation of networks and plots:

# <img src="https://raw.githubusercontent.com/bioinfproject/bioinfo/master/Folder/NeVOmics_save_plots_Windows.png" width = 51%>
The creation of networks and plots may take several minutes (5 - 10 min).

<hr />

### <center> <h1>System requirements</h1> </center>

- **Linux Operating System**<br>
[Ubuntu 14.04](http://releases.ubuntu.com/14.04/?_ga=2.84221786.213236493.1537784541-299536036.1537784541)<br>
[Ubuntu 16.04](http://releases.ubuntu.com/16.04/?_ga=2.84221786.213236493.1537784541-299536036.1537784541)
- **Windows Operating System**<br>
Windows 7<br>
Windows 8<br>
Windows 10
- **Internet Requirements**<br>
Updated Browser Google Chrome<br> 
Minimun Internet speed of 10 MB/sec and ideally 20 MB/sec or more
- **Minimum Hardware Requirements**<br>
Minimun 4GB RAM<br>
Minimum of 500MB of hard-drive<br>
Minimum Screen Resolution of 1024×768

<hr />

### <center> <h1>Operating Systems</h1> </center>


# <img src="https://upload.wikimedia.org/wikipedia/sr/thumb/1/14/Windows_logo_-_2006.svg/644px-Windows_logo_-_2006.svg.png" width = 6%> **Windows** <br>
> `(Tested on Windows 7 and 10)`

#### Install Python 3.6.7

**1**. [Download python 3.6.7 version](https://www.python.org/ftp/python/3.6.7/python-3.6.7-amd64.exe)

**2**. Double-click the `.exe` file.<br>
> - Click on **`Run`**<br> <img src="https://raw.githubusercontent.com/bioinfproject/bioinfo/master/Folder/Python_1.PNG" width = 50%><br>
> - Select the option **`Add Python 3.6 to PATH`** and then **`Install Now`**<br> <img src="https://raw.githubusercontent.com/bioinfproject/bioinfo/master/Folder/Python_2.PNG" width = 50%><br>
> - Click on **`Close`**<br> <img src="https://raw.githubusercontent.com/bioinfproject/bioinfo/master/Folder/Python_3.PNG" width = 50%><br>

**3**. To test your installation, open your terminal window and type:
```bash
python --version
```
> - It should show something like this:<br> <img src="https://raw.githubusercontent.com/bioinfproject/bioinfo/master/Folder/python%20--version.PNG" width = 50%>

More detailed information about the installation of Python on Windows [here](https://www.python.org/downloads/release/python-367/)

**4**. Download [install_python_modules.exe](https://github.com/bioinfproject/bioinfo/blob/master/Folder/install_python_modules.exe?raw=true)

**5**. Go to downloads and to find `install_python_modules.exe`, right click and Run as administrator.

#### Install blast
**1**. Download [BLAST software](http://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.7.1+-win64.exe)

**2**. Double-click the **`.exe`** file to install.

**3**. To test your installation run the command:
```bash
blastp -h
```
#### Install R

**1**. Download [R for Windows](https://cran.r-project.org/bin/windows/base/R-3.5.1-win.exe)

**2**. Double-click the **`.exe`** file.<br>
<span style="color:red">**Important: change destination location to ___`Documents`___ during installing:**</span>

# <img src="https://raw.githubusercontent.com/bioinfproject/bioinfo/master/Folder/path_to_R.PNG" width = 50%><br>
More detailed information about of R installation on Windows [here](https://cran.r-project.org/bin/windows/base/rw-FAQ.html#Does-R-run-under-Windows-Vista_003f)

**3**. Download [add_R_to_Path.exe](https://github.com/bioinfproject/bioinfo/blob/master/Folder/add_R_to_Path.exe?raw=true)
> - Go to downloads and to find `add_R_to_Path.exe`, right click and Run as administrator.

**4**. Download [install_R_packages.exe](https://github.com/bioinfproject/bioinfo/blob/master/Folder/install_R_packages.exe?raw=true
)
> - Go to downloads and to find `install_R_packages.exe`, right click and Run as administrator.

> Installation of R packages required for NeVOmics:
>- [`tidyverse`](https://cran.r-project.org/web/packages/tidyverse/index.html), [`tidygraph`](https://cran.r-project.org/web/packages/tidygraph/index.html), [`ggraph`](https://cran.r-project.org/web/packages/ggraph/index.html), [`viridis`](https://cran.r-project.org/web/packages/viridis/index.html), [`circlize`](https://cran.r-project.org/web/packages/circlize/index.html), [`RColorBrewer`](https://cran.r-project.org/web/packages/RColorBrewer/index.html), [`cowplot`](https://cran.r-project.org/web/packages/cowplot/index.html), [`networkD3`](https://cran.r-project.org/web/packages/networkD3/index.html), [`UpSetR`](https://cran.r-project.org/web/packages/UpSetR/index.html), [`gridBase`](https://cran.r-project.org/web/packages/gridBase/index.html) and [`ComplexHeatmap`](https://bioconductor.org/packages/release/bioc/html/ComplexHeatmap.html)

<hr />

# <img src="https://upload.wikimedia.org/wikipedia/commons/thumb/b/b0/NewTux.svg/300px-NewTux.svg.png" width = 5%> **Linux**<br>
> `(Tested on Ubuntu 14.04 and 16.04)`

#### Install Anaconda (Python)

**1**. [Download python 3.6 version](https://repo.anaconda.com/archive/Anaconda3-5.2.0-Linux-x86_64.sh)

**2**. Open your terminal window (Ctrl+Alt+T) and enter the command:
```bash
bash ~/Downloads/Anaconda3-5.2.0-Linux-x86_64.sh 
```
**3**. Close and open your terminal window for the installation to take effect, or you can enter the command:
```bash
source ~/.bashrc
```
**4**. To test your installation run the command:
```bash
conda list
```
More detailed information about the installation of Anaconda on Linux [here](http://docs.anaconda.com/anaconda/install/linux/)
#### Install blast
**1**. Open your terminal window (Ctrl+Alt+T) and enter the command:
```bash
sudo apt-get install ncbi-blast+
```
**2**. To test your installation run the command:
```bash
blastp -h
```
#### Install R

**1**. [Download R for Linux](https://cran.r-project.org/bin/windows/base/R-3.5.1-win.exe)

**2**. To test your R installation, open your terminal window and run the command:
```bash
R   # R version 3.5.1 (2018-07-02) -- "Feather Spray" 
    # Copyright (C) 2018 The R Foundation for Statistical Computing
    # Platform: x86_64-w64-mingw32/x64 (64-bit)
    # ...
```
> Installation of R packages required for NeVOmics:
>- [`tidyverse`](https://cran.r-project.org/web/packages/tidyverse/index.html)
>- [`tidygraph`](https://cran.r-project.org/web/packages/tidygraph/index.html)
>- [`ggraph`](https://cran.r-project.org/web/packages/ggraph/index.html)
>- [`viridis`](https://cran.r-project.org/web/packages/viridis/index.html)
>- [`circlize`](https://cran.r-project.org/web/packages/circlize/index.html)
>- [`RColorBrewer`](https://cran.r-project.org/web/packages/RColorBrewer/index.html)
>- [`cowplot`](https://cran.r-project.org/web/packages/cowplot/index.html)
>- [`networkD3`](https://cran.r-project.org/web/packages/networkD3/index.html)
>- [`UpSetR`](https://cran.r-project.org/web/packages/UpSetR/index.html)
>- [`gridBase`](https://cran.r-project.org/web/packages/gridBase/index.html)
>- [`ComplexHeatmap`](https://bioconductor.org/packages/release/bioc/html/ComplexHeatmap.html)

**3**. After verifying of R installation, in the same terminal, run each of the following commands to install R packages:
```bash
install.packages("tidyverse", dependencies = TRUE,repos='http://cran.us.r-project.org')
install.packages("tidygraph", dependencies = TRUE,repos='http://cran.us.r-project.org')
install.packages("ggraph", dependencies = TRUE,repos='http://cran.us.r-project.org')
install.packages("viridis", dependencies = TRUE,repos='http://cran.us.r-project.org')
install.packages("circlize", dependencies = TRUE,repos='http://cran.us.r-project.org')
install.packages("RColorBrewer", dependencies = TRUE,repos='http://cran.us.r-project.org')
install.packages("cowplot", dependencies = TRUE,repos='http://cran.us.r-project.org')
install.packages("networkD3", dependencies = TRUE,repos='http://cran.us.r-project.org')
install.packages("UpSetR", dependencies = TRUE,repos='http://cran.us.r-project.org')
install.packages("gridBase",dependencies = TRUE,repos='http://cran.us.r-project.org')
source("http://bioconductor.org/biocLite.R")
biocLite("ComplexHeatmap")
```
