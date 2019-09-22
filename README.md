# Some improvements are being made to the program, if you have any errors please wait for the improved version to be released.
# NeVOmics (`Ne`twork-based `V`isualization for `Omics`) <br>

## **Description**
NeVOmics is a functional enrichment analysis tool developed in programming language Python and R that integrates Over-representation analysis (ORA) methodology and network-based visualization. It applies appropriate statistical methods to identify significantly enriched Gene Ontology (GO) terms or pathways in a given list of genes/proteins. It provides several types of graphical visualization to show enrichment results. NeVOmics supports all organisms deposited in UniProt Knowledgebase (UniProtKB) and Kyoto Encyclopedia of Genes and Genomes (KEGG) databases.
# <img src="https://raw.githubusercontent.com/bioinfproject/bioinfo/master/Folder/Figure_1.png" width = 75%>

### **Citation:**<br> **NeVOmics**: An Enrichment Tool for Gene Ontology and Functional Network Analysis and Visualization of Data from OMICs Technologies. [doi.org/10.3390/genes9120569](https://doi.org/10.3390/genes9120569)

### Compatibility with: **Windows** and **Linux**


> # Download [**NeVOmics**](https://github.com/bioinfproject/bioinfo/blob/master/NeVOmics_Windows.zip?raw=true) for  <img src="https://upload.wikimedia.org/wikipedia/sr/thumb/1/14/Windows_logo_-_2006.svg/644px-Windows_logo_-_2006.svg.png" width = 7%>

> # Download [**NeVOmics**](https://github.com/bioinfproject/bioinfo/blob/master/NeVOmics_Linux.zip?raw=true) for <img src="https://upload.wikimedia.org/wikipedia/commons/thumb/b/b0/NewTux.svg/300px-NewTux.svg.png" width = 7%><br>
<hr />

### <center> <h1>Run NeVOmics</h1> </center>
# <img src="https://raw.githubusercontent.com/bioinfproject/bioinfo/master/Folder/Before%20run%20NeVOmics.png" width = 85%>
**1**. **Download NeVOmics** for **Windows** or **Linux**

**2**. Unzip **NeVOmics_[Windows/Linux].zip**

**3**. Open your terminal window:<br>
**Windows**: press the **Win+R** keys on your keyboard. Then, type **cmd** or **cmd.exe** and press **Enter** or click/tap **OK**.
**Linux**: press **Ctrl+Alt+T** on your keyboard.<br>

**4**. Change directory in **Windows** and **Linux**:
><img src="https://upload.wikimedia.org/wikipedia/sr/thumb/1/14/Windows_logo_-_2006.svg/644px-Windows_logo_-_2006.svg.png" width = 5%>
```bash
cd Downloads\NeVOmics_Windows
```
><img src="https://upload.wikimedia.org/wikipedia/commons/thumb/b/b0/NewTux.svg/300px-NewTux.svg.png" width = 5%>
```bash
cd Downloads/NeVOmics_Linux
```
**5**. Type the next command to run NeVOmics:
><img src="https://upload.wikimedia.org/wikipedia/sr/thumb/1/14/Windows_logo_-_2006.svg/644px-Windows_logo_-_2006.svg.png" width = 5%>
```bash
python NeVOmics
```
> Or double-click on **`NeVOmics.exe`** (**`Only for Windows`**) included in **`NeVOmics.zip`**.<br>

><img src="https://upload.wikimedia.org/wikipedia/commons/thumb/b/b0/NewTux.svg/300px-NewTux.svg.png" width = 5%>
```bash
python3 NeVOmics.py
```

**6**. Select one of the three available analyzes:

# <img src="https://raw.githubusercontent.com/bioinfproject/bioinfo/master/Folder/NeVOmics_intro_Windows.png" width = 50%>

**7**. Upload protein list (must be **UniProtKB** identifiers):

# <img src="https://raw.githubusercontent.com/bioinfproject/bioinfo/master/Folder/Upload_NeVOmics.PNG" width = 50%> <img src="https://raw.githubusercontent.com/bioinfproject/bioinfo/master/Folder/Upload_NeVOmics_Linux.png" width = 48.3%>

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
The creation of networks and plots may take several minutes (5 - 20 min).

<hr />

### <center> <h1>System requirements</h1> </center>

- **Linux (64-bit) Operating System**<br>
[Ubuntu 18.04 LTS (Bionic Beaver)](http://releases.ubuntu.com/18.04/?_ga=2.84221786.213236493.1537784541-299536036.1537784541)
- **Windows Operating System**<br>
Windows 7<br>
Windows 10
- **Internet Requirements**<br> 
Minimun Internet speed of 10 MB/sec and ideally 20 MB/sec or more
- **Minimum Hardware Requirements**<br>
Minimun 4GB RAM<br>
Minimum of 500MB of hard-drive<br>
Minimum Screen Resolution of 1024×768

<hr />

### <center> <h1>Operating Systems</h1> </center>


# <img src="https://upload.wikimedia.org/wikipedia/sr/thumb/1/14/Windows_logo_-_2006.svg/644px-Windows_logo_-_2006.svg.png" width = 6%> **Windows** <br>
> `(Tested on Windows 7 and 10)`

### **Install Python**

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

**4**. Download [`install_python_modules.exe`](https://github.com/bioinfproject/bioinfo/blob/master/Folder/install_python_modules.exe?raw=true)

**5**. Go to downloads and to find `install_python_modules.exe`, right click and Run as administrator.

### **Install blast**
**1**. Download [`BLAST software`](http://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.8.1+-win64.exe)

**2**. Double-click the **`.exe`** file to install.

**3**. To test your installation run the command:
```bash
blastp -h
```
### **Install R**

**1**. Download [R for Windows](https://cran.r-project.org/bin/windows/base/old/3.5.3/R-3.5.3-win.exe)

**2**. Double-click the **`.exe`** file.<br>
### <span style="color:red">**Important: change installation location to **`Documents`**.</span>

# <img src="https://raw.githubusercontent.com/bioinfproject/bioinfo/master/Folder/path_to_R.PNG" width = 50%><br>
More detailed information about of R installation on Windows [here](https://cran.r-project.org/bin/windows/base/rw-FAQ.html#Does-R-run-under-Windows-Vista_003f)

**3**. Install R Packages required by NeVOmics:<br>
Download [install_R_packages.exe](https://github.com/bioinfproject/bioinfo/blob/master/Folder/install_R_packages.exe?raw=true
)
> - Go to downloads and to find `install_R_packages.exe`, right click and Run as administrator.

> R packages needed to build networks and plots with NeVOmics:
>- [`GetoptLong`](https://cran.r-project.org/web/packages/GetoptLong/index.html),[`tidyverse`](https://cran.r-project.org/web/packages/tidyverse/index.html), [`tidygraph`](https://cran.r-project.org/web/packages/tidygraph/index.html), [`ggraph`](https://cran.r-project.org/web/packages/ggraph/index.html), [`viridis`](https://cran.r-project.org/web/packages/viridis/index.html), [`circlize`](https://cran.r-project.org/web/packages/circlize/index.html), [`RColorBrewer`](https://cran.r-project.org/web/packages/RColorBrewer/index.html), [`cowplot`](https://cran.r-project.org/web/packages/cowplot/index.html), [`networkD3`](https://cran.r-project.org/web/packages/networkD3/index.html), [`UpSetR`](https://cran.r-project.org/web/packages/UpSetR/index.html), [`gridBase`](https://cran.r-project.org/web/packages/gridBase/index.html), [`reshape2`](https://cran.r-project.org/web/packages/reshape2/index.html) and [`ComplexHeatmap`](https://bioconductor.org/packages/release/bioc/html/ComplexHeatmap.html)

<hr />

# <img src="https://upload.wikimedia.org/wikipedia/commons/thumb/b/b0/NewTux.svg/300px-NewTux.svg.png" width = 5%> **Linux**<br>
## **NeVOmics was updated to work with Ubuntu 18.04 LTS (Bionic Beaver)**

### **Python**

**1**. Python3 (v3.6.7) is included by default in Ubuntu 18.04<br>

<img src="https://raw.githubusercontent.com/bioinfproject/bioinfo/master/Folder/py3v_ubuntu.png" width = 80%><br>

**2**. NeVOmics requires complementary python3 modules:<br>

Open your terminal and type the next commands to download and run the next python3 script `Python3_mod_Ubuntu.py` to install the modules:
```bash
wget -q https://raw.githubusercontent.com/bioinfproject/bioinfo/master/Folder/Python3_mod_Ubuntu.py 
```
```bash
python3 Python3_mod_Ubuntu.py 
```
<img src="https://raw.githubusercontent.com/bioinfproject/bioinfo/master/Folder/inicio_Python3_mod_Ubuntu.png" width = 50%><br>
<img src="https://raw.githubusercontent.com/bioinfproject/bioinfo/master/Folder/fin_Python3_mod_Ubuntu.png" width = 50%>


**3**. To test modules installation run the next commands:
```bash
python3
```
```bash
help('pandas')
```
<img src="https://raw.githubusercontent.com/bioinfproject/bioinfo/master/Folder/help_pandas.png" width = 50%><br>
<img src="https://raw.githubusercontent.com/bioinfproject/bioinfo/master/Folder/info_pandas.png" width = 50%>

#### **Install blast**
**1**. Open your terminal window (Ctrl+Alt+T) and enter the command:
```bash
sudo apt-get install ncbi-blast+
```
**2**. To test your installation run the command:
```bash
blastp -h
```
<img src="https://raw.githubusercontent.com/bioinfproject/bioinfo/master/Folder/blastp_-h.png" width = 50%>

#### **Install R**

**1**. Open your terminal and type the next command:
```bash
sudo apt-get install r-base
```

**2**. To test your R installation, open your terminal window and run the command:
```bash
R
```
<img src="https://raw.githubusercontent.com/bioinfproject/bioinfo/master/Folder/R.png" width = 50%>

**3**. Preparation of a suitable environment for R

> R packages needed to build networks and plots with NeVOmics:
>- [`GetoptLong`](https://cran.r-project.org/web/packages/GetoptLong/index.html),[`tidyverse`](https://cran.r-project.org/web/packages/tidyverse/index.html), [`tidygraph`](https://cran.r-project.org/web/packages/tidygraph/index.html), [`ggraph`](https://cran.r-project.org/web/packages/ggraph/index.html), [`viridis`](https://cran.r-project.org/web/packages/viridis/index.html), [`circlize`](https://cran.r-project.org/web/packages/circlize/index.html), [`RColorBrewer`](https://cran.r-project.org/web/packages/RColorBrewer/index.html), [`cowplot`](https://cran.r-project.org/web/packages/cowplot/index.html), [`networkD3`](https://cran.r-project.org/web/packages/networkD3/index.html), [`UpSetR`](https://cran.r-project.org/web/packages/UpSetR/index.html), [`gridBase`](https://cran.r-project.org/web/packages/gridBase/index.html), [`reshape2`](https://cran.r-project.org/web/packages/reshape2/index.html) and [`ComplexHeatmap`](https://bioconductor.org/packages/release/bioc/html/ComplexHeatmap.html)

**4**. Installing of Libraries and local R Packages:<br>

Open your terminal and type the next commands to download and run the next python3 script `R-Libs-Packs-Ubuntu.py` to install the libraries and packages:
```bash
wget -q https://raw.githubusercontent.com/bioinfproject/bioinfo/master/Folder/R-Libs-Packs-Ubuntu.py 
```
```bash
python3 R-Libs-Packs-Ubuntu.py 
```
<img src="https://raw.githubusercontent.com/bioinfproject/bioinfo/master/Folder/inicio_python3 R-Libs-Packs-Ubuntu.png" width = 50%><br>
<img src="https://raw.githubusercontent.com/bioinfproject/bioinfo/master/Folder/proceso_python3 R-Libs-Packs-Ubuntu.png" width = 50%><br>
<img src="https://raw.githubusercontent.com/bioinfproject/bioinfo/master/Folder/fin_python3 R-Libs-Packs-Ubuntu.png" width = 50%>

**5**. Setting the local R library path:

Open your terminal and type the next commands:

* To define the R library path create a file .Renviron in your home directory:

* Open the text editor:
```bash
nano ~/.Renviron
```
<img src="https://raw.githubusercontent.com/bioinfproject/bioinfo/master/Folder/1_nano_Renviron.png" width = 50%><br>
* Add the following line to the file:
```bash
R_LIBS=~/R-NeVOmics
```
<img src="https://raw.githubusercontent.com/bioinfproject/bioinfo/master/Folder/2_open_nano.png" width = 50%><br>
* Save file (Crtl+O), then enter and then exit (Ctrl+x) of text editor:

<img src="https://raw.githubusercontent.com/bioinfproject/bioinfo/master/Folder/3_control_O.png" width = 50%><br>

<img src="https://raw.githubusercontent.com/bioinfproject/bioinfo/master/Folder/4_enter.png" width = 50%><br>

<img src="https://raw.githubusercontent.com/bioinfproject/bioinfo/master/Folder/5_final_nano.png" width = 50%><br>

**6**. With this settings the directory ~/R-NeVOmics is added to the list of places to look for R packages. 
Verify the location of the library and packages with the next commands:
```bash
R
```
```bash
.libPaths()
```
<img src="https://raw.githubusercontent.com/bioinfproject/bioinfo/master/Folder/libPaths.png" width = 50%><br>
```bash
library()
```
<img src="https://raw.githubusercontent.com/bioinfproject/bioinfo/master/Folder/library.png" width = 50%><br>

