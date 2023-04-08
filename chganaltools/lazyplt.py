#!/usr/bin/env python3
# coding=UTF-8
"""
Created on Fri Sep 6 10:57:09 2019
Name: lazyplt.py
Developer: Ming-Wen Chang
E-mail: ming.wen.c@gmail.com

Read a csc file and plot data
syntax: lazyplt.py filename
output: filename.png
"""

import sys 
import pandas as pd
#import matplotlib.pyplot as plt

if __name__ == "__main__":
    nargs  = len(sys.argv)
    if nargs < 2:
        print ('Syntax: lazyplt.py filename ')
        print ('Output: filename.png')
        exit()
        
    name = sys.argv[1]
    print ('Plotting data in %s file' %(name))
    df = pd.read_csv(name, delim_whitespace=True)
    ax = df.plot(x=df.columns[1], figsize=(4, 9), fontsize=12, legend =True) 
    fig = ax.get_figure()
    fig.savefig('%s.png'%(name))
    print ('Done!')


