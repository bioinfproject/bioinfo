#!/usr/bin/env python3

from subprocess import Popen, PIPE, STDOUT
from subprocess import call
import shlex, subprocess
import subprocess

print('\n\n Installation of modules used by NeVOmics\n\n')

#print('MODULES:  pip, tkinter, requests, pandas, scipy, openpyxl, colormap, easydev and networkx')

print("\nType your password")
subprocess.call(['sudo','apt-get','install','python3-pip'])

subprocess.call(['sudo','apt-get','install','python3-tk'])

a=subprocess.Popen(['python3','-m','pip','install','requests'])
a.wait()

a=subprocess.Popen(['python3','-m','pip','install','pandas'])
a.wait()

a=subprocess.Popen(['python3','-m','pip','install','scipy'])
a.wait()

a=subprocess.Popen(['python3','-m','pip','install','openpyxl'])
a.wait()

a=subprocess.Popen(['python3','-m','pip','install','--user','colormap'])
a.wait()

a=subprocess.Popen(['python3','-m','pip','install','--user','easydev'])
a.wait()

a=subprocess.Popen(['python','-m','pip','install','--user','networkx'])
a.wait()

a=subprocess.Popen(['python3','-m','pip','install','--user','googledrivedownloader'])
a.wait()

a=subprocess.Popen(['python3','-m','pip','install','--user','bioservices'])
a.wait()

a=subprocess.Popen(['python3','-m','pip','install','--user','XlsxWriter'])
a.wait()

print('\n\n Modules installed successfully\n\n')
