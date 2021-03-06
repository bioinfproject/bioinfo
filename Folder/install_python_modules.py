from subprocess import Popen, PIPE, STDOUT
from subprocess import call
import shlex, subprocess
import subprocess

print('\n\n Installation of modules used by NeVOmics\n\n')

#print('MODULES:  pip, tkinter, requests, pandas, scipy, openpyxl, colormap, easydev and networkx')

a=subprocess.Popen(['python','-m','pip','install','--upgrade','pip'])
a.wait()

a=subprocess.Popen(['python','-m','pip','install','matplotlib'])
a.wait()

a=subprocess.Popen(['python','-m','pip','install','pandas'])
a.wait()

a=subprocess.Popen(['python','-m','pip','install','scipy'])
a.wait()

a=subprocess.Popen(['python','-m','pip','install','openpyxl'])
a.wait()

a=subprocess.Popen(['python','-m','pip','install','colormap'])
a.wait()

a=subprocess.Popen(['python','-m','pip','install','easydev'])
a.wait()

a=subprocess.Popen(['python','-m','pip','install','--user','networkx'])
a.wait()

a=subprocess.Popen(['python','-m','pip','install','googledrivedownloader'])
a.wait()

a=subprocess.Popen(['python','-m','pip','install','bioservices'])
a.wait()

a=subprocess.Popen(['python','-m','pip','install','XlsxWriter'])
a.wait()

print('\n\n Modules installed successfully\n\n')
