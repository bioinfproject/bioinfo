#!/usr/bin/env python3

from subprocess import Popen, PIPE, STDOUT
from subprocess import call
import shlex, subprocess
import subprocess

print('\n\n Installation of modules used by NeVOmics\n\n')

print('\nAdministrator permissions are required to install the following modules:\n')

print('MODULE:  pip, tkinter, requests, pandas, scipy and openpyxl')
print('https://docs.python.org/3/installing/index.html [Crtl+click]')

print("\nType your password")
subprocess.call(['sudo','apt-get','install','python3-pip'])

subprocess.call(['sudo','apt-get','install','python3-tk'])

a = subprocess.Popen(['python3','-m','pip','install','requests'])
a.wait()

a = subprocess.Popen(['python3','-m','pip','install','pandas'])
a.wait()

a = subprocess.Popen(['python3','-m','pip','install','scipy'])
a.wait()

a = subprocess.Popen(['python3','-m','pip','install','openpyxl'])
a.wait()

print('\n\n Modules have been installed successfully\n\n')










