import datetime
start = datetime.datetime.now()
from pandas import Series, DataFrame 
import pandas as pd
from pandas.compat import StringIO
import csv
import pandas
import pathlib
#pd.set_option('max_rows',100000)
#pd.set_option('max_colwidth',100000)
import urllib.request
import webbrowser
import re
import shutil, os
import numpy as np
from urllib.request import urlopen
#from bs4 import BeautifulSoup
import requests
from time import sleep
from subprocess import Popen, PIPE, STDOUT
from subprocess import call
import shlex, subprocess
import subprocess
import sys
import warnings
#warnings.filterwarnings("ignore")
from datetime import datetime 
inicio_total = datetime.now()
import os, fnmatch

urllib.request.urlretrieve('https://raw.githubusercontent.com/bioinfproject/bioinfo/master/Folder/R_packages_for_NeVOmics.R', './R_packages_for_NeVOmics.R')

subprocess.call(["Rscript","--vanilla","R_packages_for_NeVOmics.R"])

if os.path.exists('R_packages_for_NeVOmics.R'): os.remove('R_packages_for_NeVOmics.R')