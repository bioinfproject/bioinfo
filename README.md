## GOmics
#### Description
GOmics is a Functional Enrichment Tool to facilitate the functional characterization of OMICs technologies data such as, proteomic and transcriptomic approaches.

_Please cite the following research paper if you use GOmics in your research_:

##### GOmics: an enrichment tool for gene ontology and visualization of functional networks in data from OMICs technologies

>#### Compatible with: <img src="https://upload.wikimedia.org/wikipedia/commons/thumb/b/b0/NewTux.svg/300px-NewTux.svg.png" width = 5%> <img src="https://upload.wikimedia.org/wikipedia/sr/thumb/1/14/Windows_logo_-_2006.svg/644px-Windows_logo_-_2006.svg.png" width = 6%>

### Run GOmics
1. [**Download GOmics**](https://github.com/bioinfproject/bioinfo/blob/master/GOmics.zip?raw=true)
2. Unzip **GOmics.zip**
3. Open your terminal window:<br>
**Linux**: press **Ctrl+Alt+T** on your keyboard<br>
**Windows**: press the **Win+R** keys on your keyboard. Then, type **cmd** or **cmd.exe** and press **Enter** or click/tap **OK**.
4. Change directory in **Linux** and **Windows**:<br>
```bash
cd Download/GOmics # Linux
```
```bash
cd Download\GOmics # Windows
```
5. Type the next command to run GOmics:
```
python GOmics
```
6. Select one of the three available analyzes: `[1/2/3]`
```bash
[ Functional Enrichment Analysis ]

   1  Gene Ontology Enrichment

   2  KEGG Pathways Enrichment

   3  KEGG Pathways Enrichment using blastp)

***** Select an analysis: [1/2/3]
=====> : 1   #type here any analysis,example: Gene Ontology Enrichment
```
7. Enter a filename with gene/protein list::
```bash
[ Step 1: Enter filename (Uniprot IDs) ]
=====> : filename.tsv   #type here filname, must be in this directory
```
8. Selection of a method p-value correction
```bash
[ Step 2: Choose a correction method (e.g., FDR / Bonferroni) ]
=====> : FDR   # type any method
```
9. Choose a significance value for the analysis
```bash
[ Step 3: Choose a Value (e.g., 0.05) ]
=====> : 0.05   # type a significance value
```
10. Optional creation of network
```bash
[ Step 4: Do you want to create the networks? [y/n]
=====> : y   # the creation of graphics consumes time
```
# _________________________________
### **System requirements**
1.
2.
3.
# _________________________________
<img src="https://upload.wikimedia.org/wikipedia/commons/thumb/b/b0/NewTux.svg/300px-NewTux.svg.png" width = 5%><br>
> ### **For Linux** ``(Tested on Ubuntu 16.04)``

#### Install Anaconda (Python)

1. [Download python 3.6 version](https://repo.anaconda.com/archive/Anaconda3-5.2.0-Linux-x86_64.sh)

2. Open your terminal window (Ctrl+Alt+T) and enter the command:
```bash
bash ~/Downloads/Anaconda3-5.2.0-Linux-x86_64.sh 
```
3. Close and open your terminal window for the installation to take effect, or you can enter the command:
```bash
source ~/.bashrc
```
4. To test your installation run the command:
```bash
conda list
```
More detailed information about the installation of Anaconda on Linux [here](http://docs.anaconda.com/anaconda/install/linux/)
#### Install blast
1. Open your terminal window (Ctrl+Alt+T) and enter the command:
```bash
sudo apt-get install ncbi-blast+
```
2. To test your installation run the command:
```bash
blastp -h
```
# _________________________________
<img src="https://upload.wikimedia.org/wikipedia/sr/thumb/1/14/Windows_logo_-_2006.svg/644px-Windows_logo_-_2006.svg.png" width = 6%><br>
> ### **For Windows** ``(Tested on Windows 7 and 10)``
#### Install R


#### R library usage


#### Install Anaconda (Python)

1. [Download python 3.6 version](https://repo.anaconda.com/archive/Anaconda3-5.2.0-Windows-x86_64.exe)
2. Double-click the `.exe` file.
3. To test your installation, in your Anaconda Prompt, run the command:
```bash
conda list
```
More detailed information about the installation of Anaconda on Windows [here](http://docs.anaconda.com/anaconda/install/windows/)

4. Adding Python to the Windows environment:
>- Open ___`Control Panel`___ > ___`System and Security`___ > ___`System`___
>- ___`Advanced System Settings`___ > ___`Enviroment Variables`___
>- Find Path in ___`System Variables`___ then click on Edit
>- At the end of the row add a semicolon and then copy the python location
>- Click on  

#### Install blast
1. [Download BLAST software](ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.7.1+-win64.exe)
2. Double-click the `.exe` file to install.
3. To test your installation run the command:
```bash
blastp -h
```
#### Install R

#### R library usage
