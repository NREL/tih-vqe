#!/usr/bin/env python

import os, sys
import argparse 
import numpy as np
import subprocess
from shutil import copy2
import json


def find_key(key_input, tempfile):
    #finds line where key occurs in stored input, last instance
    key_input = str(key_input)
    line = len(tempfile)                  #default to end
    for i in range(0,len(tempfile)):
        if key_input in tempfile[i]:
            line = i
    return line

def find_next_key(key_input, tempfile, startline):
    #finds first line where key occurs in stored input starting from startline
    key_input = str(key_input)
    line = len(tempfile)                  #default to end
    for i in range(startline,len(tempfile)):
        if key_input in tempfile[i]:
            line = i
            break
    return line

def find_first_range_key(key_input, tempfile, startline=0, endline=-1, skip_pound = False):
    #finds all lines that exactly begin with key
    key_input = str(key_input)
    startlen = len(key_input)
    L = []

    if endline == -1:
        endline = len(tempfile)
    for i in range(startline,endline):
        line = tempfile[i]
        if skip_pound == True:
            for j in range(10):  #repeat to make sure no really weird formatting
                line = line.lstrip()
                line = line.lstrip('#')
        line = line[0:startlen]
        if line == key_input:
            L.append(i)
    if not L:
        L = [len(tempfile)]
    return L

def get_dir_files(a_dir = './'):
    #return list of all files within a_dir
    a_dir = a_dir.strip('/') + '/'    #ensure ends in slash
    filenames = [f for f in os.listdir(a_dir) if os.path.isfile(os.path.join(a_dir, f))]
    return filenames

def get_immediate_subdirectories(a_dir = './'):
    #return list of all directories within a_dir
    return [name for name in os.listdir(a_dir) if os.path.isdir(os.path.join(a_dir, name))]

def run_command(command):
    import subprocess
    #run 'command' as a shell command
    #  'command' must be a list of strings
    subprocess.call(command,shell=False)

def run_command_file(command, filename):
    import subprocess
    #run 'command' as a shell command and save output to filename
    #  'command' must be a list of strings
    with open(filename, "w") as outfile:
        subprocess.call(command, stdout = outfile)

def readfile(filename):
    #read a file into a list of strings
    f = open(filename,'r')
    tempfile = f.readlines()
    f.close()
    return tempfile

def writefile(filename,tempfile):
    #write tempfile (list of strings) to filename
    f = open(filename,'w')
    f.writelines(tempfile)
    f.close()

def readjson(filename):
    #read json in dict
    with open(filename) as json_file:
        data = json.load(json_file)
    return data

def writejson(filename,dictdata):
    #write dict data to json
    with open(filename,'w') as outfile:
        json.dump(dictdata, outfile)

def check_file(filename):
    #check if file exists and die if doesn't
    if not os.path.isfile(filename):
        print('\'' + filename + '\' file doesn\'t exist!')
        sys.exit()

def check_file_mild(filename):
    #check if file exists but only return flag, doesn't die
    if not os.path.isfile(filename):
        print('\'' + filename + '\' file doesn\'t exist!')
        return False
    return True


atoms = ['Li','K','Na','Ti']
a = 'Na'
dists = np.arange(0.5,3.51,0.1)
mults = {'Li':[1,3],'Na':[1,3],'K':[1,3],'Ti':[2,4,6]}
basissets = ['STO-3G','3-21G','6-31G','cc-pVDZ','aug-cc-pVQZ']
basisname = ['sto3g','321g','631g','ccd','augccq']

temp = ''
for ib,b in enumerate(basissets):
    print()
    print(b)
    print(mults[a])
    for d in dists:
        E = str(round(d,1))+' '
        for m in mults[a]:
            dirname = a.lower() + 'h-' + str(m) + '-ccsd-'+str(round(d,1))+'-'+ basisname[ib]
            if not os.path.isfile(dirname+'/out'):
                E += '#N/A '
                continue
            temp = readfile(dirname+'/out')
            line = find_key('Wavefunction amplitudes converged.',temp)
            if line == len(temp):
                E += '#N/A '
                continue
            energy = temp[line].split()
            energy = float(energy[4])
            E += str(energy)+ ' '
        print(E)



