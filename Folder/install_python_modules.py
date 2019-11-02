from subprocess import Popen, PIPE, STDOUT
from subprocess import call
import shlex, subprocess
import subprocess

print('\n\n Installation of modules used by NeVOmics\n\n')

print('MODULE:  pip, tkinter, requests, pandas, scipy, openpyxl and networkx')

a=subprocess.Popen(['python','-m','pip','install','--upgrade','pip'])
a.wait()

a=subprocess.Popen(['python','-m','pip','install','requests'])
a.wait()

a=subprocess.Popen(['python','-m','pip','install','pandas'])
a.wait()

a=subprocess.Popen(['python','-m','pip','install','scipy'])
a.wait()

a=subprocess.Popen(['python','-m','pip','install','openpyxl'])
a.wait()

a=subprocess.Popen(['python','-m','pip','install','--user','networkx'])
a.wait()

print('\n\n Modules installed, Python is ready\n\n')