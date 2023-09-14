#!/usr/bin/env python

#have separate Hamilton file for each processor to simplify reading

import os, sys
import numpy as np
import argparse
from shutil import copy2


def check_file_mild(filename):
    #check if file exists but only return flag, doesn't die
    if not os.path.isfile(filename):
        print('\'' + filename + '\' file doesn\'t exist!')
        return False
    return True

def parse_args():
    #get arguments from command line
    parser = argparse.ArgumentParser(description="script to run custom vqe")
    parser.add_argument("op", type=int, help="which operator function to use")
    parser.add_argument('ntasks', type=int, help='total number of tasks in job')
    parser.add_argument('idist', type=int, help='bond distance index')
    args = parser.parse_args()
    return args

def main():

    #store arguments
    args = parse_args()
    op = args.op
    ntasks = args.ntasks
    idist = args.idist

    os.chdir('allH')

    Hnamebase = 'op' + str(op) + 'H-d' + str(idist) + '-'
    Hnameorig = 'op' + str(op) + 'H-d' + str(idist) + '.pkl'
    if not check_file_mild(Hnameorig):
        print('Did not find {:}!!! Quitting...'.format(Hnameorig))
        sys.exit()

    #have separate Hamilton file for each processor to simplify reading
    for i in range(ntasks):
        H = Hnamebase + str(i) + '.pkl'
        copy2(Hnameorig,H)

    os.chdir('../')


if __name__ == '__main__':
    main()