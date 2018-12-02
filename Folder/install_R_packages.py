import urllib.request
import shutil, os
from urllib.request import urlopen
import requests
from time import sleep
from subprocess import Popen, PIPE, STDOUT
from subprocess import call
import shlex, subprocess
import subprocess
import sys
import warnings
import os, fnmatch

urllib.request.urlretrieve('https://raw.githubusercontent.com/bioinfproject/bioinfo/master/Folder/R_packages_for_NeVOmics.R', './R_packages_for_NeVOmics.R')

subprocess.call(["Rscript","--vanilla","R_packages_for_NeVOmics.R"])

if os.path.exists('R_packages_for_NeVOmics.R'): os.remove('R_packages_for_NeVOmics.R')
