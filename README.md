## GOmics
#### Description
GOmics is a Functional Enrichment Tool to facilitate the functional characterization of OMICs technologies data such as, proteomic and transcriptomic approaches.

_Please cite the following research paper if you use GOmics in your research_:

##### GOmics: an enrichment tool for gene ontology and visualization of functional networks in data from OMICs technologies

>>>#### Compatible with: <img src="https://upload.wikimedia.org/wikipedia/commons/thumb/b/b0/NewTux.svg/300px-NewTux.svg.png" width = 5%> <img src="https://upload.wikimedia.org/wikipedia/sr/thumb/1/14/Windows_logo_-_2006.svg/644px-Windows_logo_-_2006.svg.png" width = 6%>
<hr />

### Run GOmics

**1**. [**Download GOmics**](https://github.com/bioinfproject/bioinfo/blob/master/GOmics.zip?raw=true)

**2**. Unzip **GOmics.zip**

**3**. Open your terminal window:
>**Linux**: press **Ctrl+Alt+T** on your keyboard<br>
>**Windows**: press the **Win+R** keys on your keyboard. Then, type **cmd** or **cmd.exe** and press **Enter** or click/tap **OK**.<br>

**4**. Change directory in **Linux** and **Windows**:
```bash 
cd Download/GOmics      # On Linux
```
```bash 
cd Download\GOmics      # On Windows
```
**5**. Type the next command to run GOmics:
```bash 
python GOmics
```
**6**. Select one of the three available analyzes: `[1/2/3]`
```python
[ Functional Enrichment Analysis ]

   1  Gene Ontology Enrichment

   2  KEGG Pathways Enrichment

   3  KEGG Pathways Enrichment using blastp)

***** Select an analysis: [1/2/3]
=====> : 1      # type here any analysis,example: Gene Ontology Enrichment
```
**7**. Enter a filename with gene/protein list::
```bash
[ Step 1: Enter filename (Uniprot IDs) ]
=====> : filename.tsv       # type here filname, must be in this directory
```
**8**. Selection of a method p-value correction
```bash
[ Step 2: Choose a correction method (e.g., FDR / Bonferroni) ]
=====> : FDR        # type any method
```
**9**. Choose a significance value for the analysis
```bash
[ Step 3: Choose a Value (e.g., 0.05) ]
=====> : 0.05       # type a significance value
```
**10**. Optional creation of network
```bash
[ Step 4: Do you want to create the networks? [y/n]
=====> : y      # the creation of graphics consumes time
```
<hr />

### **System requirements**
1.
2.
3.
<hr />
