#!/usr/bin/env python3

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


# Preparacion del ambiente de trabajo para R

print("\nPreparing the environment for R\n")

print("***  1. Install Libraries")

print('\nPermissions are needed to install the following libraries\nfor the operation of R in your system:\n')
print("LIBRARY:  libxml2-dev")
print("INFORMATION DETAILS: https://packages.debian.org/es/sid/libxml2-dev [Ctrl+click]\n")
print("LIBRARY:  libcurl4-openssl-dev")
print("INFORMATION DETAILS: https://packages.debian.org/es/sid/libcurl4-openssl-dev [Ctrl+click]\n")
print("LIBRARY:  libssl-dev")
print("INFORMATION DETAILS: https://packages.debian.org/es/sid/libssl-dev [Ctrl+click]")

print("\nType your password")
subprocess.call(['sudo','apt-get','install','libxml2-dev'])
subprocess.call(['sudo','apt-get','install','libcurl4-openssl-dev'])
subprocess.call(['sudo','apt-get','install','libssl-dev'])

# instalacion de paquetes de R

print("\n***  2. Install R Packages\n")

home = os.path.expanduser('~/R-NeVOmics')
os.makedirs(home,exist_ok=True)

urllib.request.urlretrieve('https://raw.githubusercontent.com/bioinfproject/bioinfo/master/Folder/proof_dir.R', './R_packages_for_NeVOmics.R')

print('\nRunning: Rscript\n')
        
subprocess.call(['Rscript','--vanilla','R_packages_for_NeVOmics.R'])

if os.path.exists('R_packages_for_NeVOmics.R'): os.remove('R_packages_for_NeVOmics.R')

#
#
#
#
#
#

