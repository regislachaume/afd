#! /usr/bin/env python3

import afdtable
from astropy.table import Table
from matplotlib import pylab as plt
import numpy as np

def plot_afd(tabs):

    afd = [tab['AFD5%'].sum() + tab['AFD95%'].sum() for tab in tabs]
    alum = [tab['U'].sum() for tab in tabs]
    prof = [tab['S'].sum() for tab in tabs]

    tab = Table.read('../src/macro.tsv', format='ascii')
    year = tab['year']
    afd *= (tab['UF'][-1] / tab['UF'])
    alum = alum / alum[-1]
    prof = prof / prof[-1]
    gdp = np.cumprod(1 + tab['growth'] / 100)
    gdp /= gdp[-1]
    wage = np.cumprod(1 + tab['IR'] / 100)
    wage /= wage[-1]
    
    fig = plt.figure(1)
    fig.clf()
    ax = fig.add_subplot(111)
   
    alpha = 1 / (year[-1] - year[0])
    afd_g = ((afd[-1] / afd[0]) ** alpha - 1) * 100
    alum_g = ((alum[-1] / alum[0]) ** alpha - 1) * 100
    prof_g = ((prof[-1] / prof[0]) ** alpha - 1) * 100
    gdp_g = ((gdp[-1] / gdp[0]) ** alpha - 1) * 100
    wage_g = ((wage[-1] / wage[0]) ** alpha - 1) * 100

    ax.plot(year, afd / 1e+6, '-', color=(.4,.4,.4), 
                lw=3, label=f'AFD ({afd_g:+.1f}\\%)')
    m = afd[-1] / 1e+6
    ax.set_ylabel('AFD [$10^9$ CLP]')
    ax.set_ylim(0, 1.1 * m)

    ax2 = ax.twinx()
    ax2.plot(year, 100 * alum, 'k-', label=f'undergraduates ({alum_g:+.1f}\\%)')
    ax2.plot(year, 100 * prof, 'k--', label=f'professors ({prof_g:+.1f}\\%)')
    ax2.plot(year, 100 * gdp,  'k-.', label=f'GDP/capita ({gdp_g:+.1f}\\%)')
    ax2.plot(year, 100 * wage, 'k:', label=f'real wages ({wage_g:+.1f}\\%)')
    ax2.set_ylim(0, 110)
    ax2.set_ylabel('index 100 in 2021')

    fig.legend(bbox_to_anchor=(0.9, 0.4))
 
    fig.show()

    fig.savefig('../pdf/total-afd-timeseries.pdf')

    
years = range(2006, 2022)
tabs = [afdtable.read(year) for year in years]
plot_afd(tabs)
