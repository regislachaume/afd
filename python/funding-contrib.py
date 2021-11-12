#! /usr/bin/env python3.9

import afdtable
from matplotlib import pylab as plt
import numpy as np
import re

def scatter(tabs):
    
    fig = plt.figure(2, figsize=(6.5, 7.5))
    fig.clf()
    ax = fig.add_subplot(111)

    tab = tab[:-2]
    
    univ = tab['University']
    nuniv = len(tab)

    # 2021 snapshot

    y = np.array([tab[f'F{i}'] for i in ['h', 1, 2, 3, 4, 5]])
    y_old = y[0] / 1e6 
    
    y_sum = y.sum(axis=0)
    y_science = y[4:6].sum(axis=0) / y[1:6].sum(axis=0)
    y_new = (y_sum - y[0]) / y_sum
   
    y_science_mean = np.median(y_science) #  y[4:6].sum() / y[1:6].sum()
    y_new_mean = 1 - 0.95 ** (2021-2005) # y[0].sum() / y.sum()

    ax.vlines(y_science_mean, 0, 0.32, 'k', '--')
    ax.vlines(y_science_mean, 0.50, 1, 'k', '--')
    ax.hlines(y_new_mean, 0, 1, 'k', '--')
 
    ax.scatter(y_science, y_new, c=[(.5,.5,.5)], 
            s=120*np.sqrt(y_sum/y_sum.max()), edgecolor='k')    
    ax.set_xlabel('teaching $\\rightarrow$ science specialisation --- $\\mathcal{S}_{i} / (\\mathcal{S}_{i} + \\mathcal{T}_{i})$')
    ax.set_xlim(0, 1)

    ax.set_ylabel('established $\\rightarrow$ surging --- $1 - \\mathcal{H}_{i} / F_i$')
    ax.set_ylim(0, 1)

    for ys, yn, un in zip(y_science, y_new, univ):
        if un in ['U. de la Frontera', 'U. de Magallanes', 'U. de Concepción',
                'U. de Antofagasta']:
            ax.text(ys - 0., yn - 0.02, un, fontsize=10, rotation=-70,
                va='top', ha='left')
        elif yn < 0.55:
            ax.text(ys - 0.015, yn - 0.02, un, fontsize=10, rotation=-90,
                va='top', ha='left')
        else:
            ax.text(ys - 0.01, yn + 0.02, un, fontsize=10, rotation=-90,
                va='bottom', ha='left')
    
    fig.tight_layout()
    fig.show()
    
    fig.savefig('../pdf/afd-specialisation.pdf')

def cumulate_afd5(tabs):

    AFD_old = np.zeros((len(tabs[0]),))
    F = np.zeros((5,))
    names = ['F1', 'F2', 'F3', 'F4', 'F5']

    for tab in tabs:
       
        AFD95 = tab['AFD95%'] 
        dampening = tab['AFD95%'] / np.maximum(AFD_old, 0.1)
        AFD = AFD95 + tab['AFD5%']

        F = dampening * F +  np.array([tab[f"f{i}"] for i in range(1, 6)]) 

        tab.add_columns(F, names=names)
        tab.add_column(AFD - F, name='Fh')

        AFD_old = AFD

def plot_afd5(tab):
  
    tab = tab.copy()
    tab.sort(['AFD95%'], reverse=True)
 
    colors = [ (1.0, 1.0, 1.0), (0, 0, 0), (0.4, 0.1, 0.1), 
              (0.2, 0.8, 0.2), (0.4, 0.4, 1.0),
              (0.8, 0.8, 0.8),]

    c = [0.01, 0.15, 0.24, 0.25, 0.35]
    
    fig = plt.figure(1, figsize=(6.5,7.5))
    fig.clf()

    # 5% AFD

    ax = fig.add_subplot(211)

    univ = tab['University']
    label = ['pre-2005 metrics/inital funding', 'students per major', 'students per prof', 
        'graduates per prof', 'grants per prof', 'papers per prof']
    f = [tab[f'f{i + 1}'] / 1e6 for i in range(5)]

    bottom = 0
    for fi, labeli, colori in zip(f, label[1:], colors[1:]):
        ax.bar(univ, fi, bottom=bottom, label=labeli, 
            align='edge', width=-0.8, color=colori, edgecolor='black')
        bottom += fi
    ax.legend(title=f"5\% AFD in {tab.meta['year']}") 
  
    ax.set_xticklabels([]) 
    ax.tick_params(axis='x', bottom=False) 
    
    ax.set_ylabel('funding [$10^9$ CLP]')

    ax2 = ax.twinx()
    fsum = np.sum(f) 
    ysum = np.sum([c[i] * tab[f'y{i+1}'] for i in range(5)])
    factor = ysum / fsum  # lucas -> billions
    ax2.set_ylim(0, ax.get_ylim()[1] * ysum / fsum)

    ax2.set_ylabel('metrics-based score $y$')
    
    # Total 

    ax = fig.add_subplot(212)

    y = [tab[f'F{k}'] / 1e6 for k in ['h', 1, 2, 3, 4, 5]]

    bottom = 0
    for yi, labeli, colori in zip(y, label, colors):
        ax.bar(univ, yi, bottom=bottom, label=labeli, 
            align='edge', width=-0.8, color=colori, edgecolor='black')
        bottom += yi
    ax.legend(title=f"Total AFD in {tab.meta['year']}") 
  
    for tick in ax.get_xticklabels():
        tick.set_rotation(38)
        tick.set_horizontalalignment('right')
    ax.tick_params(axis='x', bottom=False) 
    
    ax.set_ylabel('funding [$10^9$ CLP]')

    fig.show()
    fig.tight_layout()
    fig.subplots_adjust(hspace=0.01)
    
    fig.savefig('../pdf/afd-coefficients.pdf')

tab = afdtable.read_all_years()[-1]
scatter(tab)
plot_afd5(tab)

