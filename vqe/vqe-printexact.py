#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import os,sys
import numpy as np

import helpers


Nop = 1

dirname = 'opdata/'
filename = 'opdata-bonds.json'

if os.path.isfile(dirname+filename):
    data = helpers.readjson(dirname+filename)
else:
    print('file not found')
    sys.exit()

data = data[str(Nop)]
for i in data:
    print(data[i]['reference_energy'])


