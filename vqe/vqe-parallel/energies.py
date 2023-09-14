#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import os, sys
import argparse
import numpy as np
import subprocess
from shutil import copy2
import json


def find_next_key(key_input, tempfile, startline):
    #finds first line where key occurs in stored input starting from startline
    key_input = str(key_input)
    line = len(tempfile)                  #default to end
    for i in range(startline,len(tempfile)):
        if key_input in tempfile[i]:
            line = i
            break
    return line

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

def pygrep(searchstr, tempfile, Blines = 0, Alines = 1, quiet = False):
    line = 0
    count = 0
    output = []
    while line < len(tempfile):
        line = find_next_key(searchstr, tempfile, line)
        if line != len(tempfile):
            lines_to_print = tempfile[line-Blines:line+Alines]
            for x in lines_to_print:
                output.append(x.split('\n')[0])
                if quiet == False:
                    print(x.split('\n')[0])
        line += 1    #so starts search on next line
        if line != len(tempfile):
            count += 1
    return output,count

'''

aggregates all energies from exact, uccs, uccsd runs for inputted op directory

'''

def parse_args():
    #get arguments from command line
    parser = argparse.ArgumentParser(description="Script to ")
    parser.add_argument("op", type=str, help="directory to check")
    args = parser.parse_args()
    return args

def main():

    #store arguments
    args = parse_args()
    op = args.op

    idists = range(31)  #hardcoded idist range

    #initialize
    vals = np.zeros([len(range(31)),3])

    #exact
    c = 0
    for u in ['s','sd']:   #try to get this value from either s or sd files
        for i in idists:
            if os.path.isfile(op + '/parallel.log-'+u+'-'+str(i)):
                out = readfile(op + '/parallel.log-'+u+'-'+str(i))
                output,count = pygrep('reference_energy',out,quiet=True)
                if len(output) == 1:
                    vals[i,c] = output[0].split()[-1]

    #vqe
    for u in ['s','sd']:
        c += 1
        for i in idists:
            if os.path.isfile(op + '/parallel.log-'+u+'-'+str(i)):
                out = readfile(op + '/parallel.log-'+u+'-'+str(i))
                output,count = pygrep('VQE',out,quiet=True)
                if len(output) == 1:
                    vals[i,c] = output[0].split()[-1]
                else:
                    output,count = pygrep('finalE',out,quiet=True)
                    if len(output) > 0:
                        vals[i,c] = output[-1].split()[-1]

    print('Exact UCCS UCCSD')
    vals = vals.tolist()
    vals = [[str(y) for y in x] for x in vals]
    vals = [' '.join(x) for x in vals]
    for x in vals:
        print(x)


if __name__ == '__main__':
    main()