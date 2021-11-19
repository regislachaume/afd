#! /usr/bin/env python3.9

import afdtable
import numpy as np
from matplotlib import pylab as plt

def calc(tab):

    cols1 = ['U', 'M', 'S', 'Sp', 'G', 'P']
    cols2 =  ['xi1', 'y1', 'xi2', 'y2',  'xi3', 'y3',  'xi4', 'y4',
              'xi5', 'y5', 'y']

    table = "\\begin{tabular}{l rrrrrr rrrrrrrrrr r rr}\n"
    table += "\\hline\\hline\n"
    table += f"{'University':<28} & "
    table += "\\npup & \\nmaj & \\nprof & \\ngrad  & \\ngrant & P & "
    n = 2021
    for k in range(1, 6):
        table += f"$\\xi_{{i,n,{k}}}$ & $y_{{i,n,{k}}}$ & "
    table += "$y_{i,n}$ &"
    table += "\\multicolumn{2}{c}{5\\% AFD}\\\\\n"
    table += "&" * 18
    table += " [\\%] & [kCLP] \\\\\n"
    table += "\\hline\n"

    for row in tab:

        univ = row['University']
        table += f"{univ:<28} & "
        
        u, m, s, sp, g, p = [row[c] for c in cols1]
        table += f"{u:5.0f} & {m:3.0f} & {s:7.2f} & {sp:7.2f} & {g:3.0f} &{p:7.2f} & "

        for c in cols2:
            table += f"${row[c]:6.3f}$ & "
        table += f"{100*row['%_AFD5%']:6.3f} & {row['AFD5%']:7.0f} \\\\\n"

    table += "\hline\n"
    table += "\\end{tabular}"
    
    with open('../tex/tab-calcdetail.tex', 'w') as out:
        out.write(table)   


tab = afdtable.read(2021)
c = [0.01, 0.15, 0.24, 0.25, 0.35]
y = sum([tab[f'y{k + 1}'] * c[k] for k in range(5)])
tab.add_column(y, name='y')
calc(tab)
