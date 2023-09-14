#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import os, sys
import argparse
import numpy as np


def find_key(key_input, tempfile):
    #finds line where key occurs in stored input, last instance
    key_input = str(key_input)
    line = len(tempfile)                  #default to end
    for i in range(0,len(tempfile)):
        if key_input in tempfile[i]:
            line = i
    return line

def find_prev_key(key_input, tempfile, startline):
    #finds last line where key occurs in stored input starting from startline
    key_input = str(key_input)
    line = -1                  #default to before beginning
    for i in range(startline,-1,-1):
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

'''

aggregates parameters for uccs and uccsd runs for inputted op

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

    if not os.path.isdir(op):
        print('{:} directory not found!'.format(op))
        sys.exit()

    dists = list(range(31))  #hardcoded distance index
    ucc = ['s','sd']

    for u in ucc:
        for d in dists:
            print(u,d)
            if not os.path.isfile(op + '/' + 'parallel.log-'+u+'-'+str(d)):
                print('could not find file for ',u,d)
                continue
            tempfile = readfile(op + '/' + 'parallel.log-'+u+'-'+str(d))
            tempfile = [x for x in tempfile if not 'Starting vqe-wrapper' in x]   #to fix weird bug where this line prints within a param block sometimes
            endline = find_key('finalE',tempfile)
            if endline == len(tempfile):
                print('Could not find data for this job')
                continue
            startline = find_prev_key('params',tempfile,endline)
            if startline == -1:
                print('Could not find data for this job')
                continue
            params = tempfile[startline+1:endline]
            params = [x.split() for x in params]
            finalparams = []
            for p in params:
                finalparams.extend(p)

            finalparams = [x + '\n' for x in finalparams]
            writefile(op + '/' + 'params-'+u+'-'+str(d),finalparams)


if __name__ == '__main__':
    main()