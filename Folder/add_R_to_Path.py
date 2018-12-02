
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

d = os.getcwd() 

a=d.split("\\")[1]

b=d.split("\\")[2]

ruta = "C:/"+a+"/"+b+"/Documents/R-3.5.1/bin;%Path%"

subprocess.call(["setx","Path",ruta])