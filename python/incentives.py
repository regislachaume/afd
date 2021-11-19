#! /usr/bin/env python3.9

import afdtable
import numpy as np
from matplotlib import pylab as plt

def incentives(tab, p=0.02, d=0.05, include_caption=False):
    """Determine the marginal earnings for an additional paper, research
project, or post-grad professor.

Arguments:
    tab:
        AFD table
    p:
        Yearly increase of total AFD (in constant Chilean pesos).
    d:
        Yearly depreciation in constant Chilean pesos 

Returns:
    None

Side effects:
    LaTex table
        
    """

    year = tab.meta['year']
    if not tab.meta['computed']:
        afdtable.compute(tab)
    
    df = np.zeros((4, len(tab)))
    for k in range(len(tab)):
        tab1 = afdtable.change_metrics(tab, k, 'U', increment=1)
        df[0,k] = 1e-3 * (tab1[k]['f'] - tab[k]['f'])
        tab1 = afdtable.change_metrics(tab, k, 'Sp', increment=1)
        df[1,k] = 1e-3 * (tab1[k]['f'] - tab[k]['f'])
        tab1 = afdtable.change_metrics(tab, k, 'G', increment=1)
        df[2,k] = 1e-3 * (tab1[k]['f'] - tab[k]['f'])
        tab1 = afdtable.change_metrics(tab, k, 'P', increment=1)
        df[3,k] = 1e-3 * (tab1[k]['f'] - tab[k]['f'])

    # tex table

    filename = '../tex/tab-incentives.tex'
    nl, tnl = "\n", "\\\\"
    mc = "\\multicolumn{2}{c}"
    mclp = f"{mc}{{[$10^6$ CLP]}}"
    caption = '\\caption'
    label = '\\label'

    with open(filename, 'w') as out:

        if include_caption:
            out.write('\\begin{table}\n')
            out.write(f"{caption}{{Additional earnings in {year} and present value of total earnings if a University has improved the following metrics in {year-1}: one additional undergraduate student, one additional full-time contract for a post-graduate professor during a year, one additional ongoing research grant, and one additional Web of Science (ex-ISI) publication. Assumptions are: the total State funding will continue to grow 2\% per year in real terms; a yearly depreciation of 5\% in real terms ($\\approx 8$% in pesos)) to determine the present value of future earnings.}}{nl}")
            out.write(f"{label}{{tab:incentives:2000}}{nl}")

        out.write('\\begin{tabular}{l crr crr crr crr}\n')
        out.write('\\hline\\hline\n')
        out.write(f"{'university':30} & {nl}")
        out.write(f"& {mc + '{student}':>64} & {nl}")
        out.write(f"& {mc + '{postgrad. staff}':>64} & {nl}")
        out.write(f"& {mc + '{grant}':>98} & {nl}")
        out.write(f"& {mc + '{publication}':>132} {tnl}{nl}")
        out.write(f"{'':30} & ")
        out.write(f"& {year:13} & {'total':15} & ")
        out.write(f"& {year:13} & {'total':15} & ")
        out.write(f"& {year:13} & {'total':15} & ")
        out.write(f"& {year:13} & {'total':15} {tnl}{nl}")
        out.write(f"{'':30} & ")
        out.write(f"& {mclp:31} & ")
        out.write(f"& {mclp:31} & ")
        out.write(f"& {mclp:31} & ")
        out.write(f"& {mclp:31} {tnl}{nl}")
        out.write('\\hline\n')

        f = 1 / (1 - 0.95 * (1+p)/(1+d))

        for k in range(len(tab)):
            df1 = df[0, k]
            df3 = df[1, k]
            df4 = df[2, k]
            df5 = df[3, k]
            u = tab[k]['University']
            out.write(f"{u:30} & ")
            out.write(f"& ${df1:11.3f}$ & ${df1*f:13.3f}$ & ")
            out.write(f"& ${df3:11.3f}$ & ${df3*f:13.3f}$ & ")
            out.write(f"& ${df4:11.3f}$ & ${df4*f:13.3f}$ & ")
            out.write(f"& ${df5:11.3f}$ & ${df5*f:13.3f}$ {tnl}{nl}")

        out.write('\\hline\n')
        out.write('\\end{tabular}\n')

        if include_caption:
            out.write('\\end{table}')
    
    print(f'Incentives written into {filename}')
    
    # plot
    fig = plt.figure(3, figsize=(6.5,7.5))
    fig.clf()
    ax = fig.add_subplot(111)
    univ = tab['University']
    df_tot = df.sum(axis=0)
    indices = np.argsort(df_tot)[::-1]
    univ = univ[indices]
    df = df[:,indices]
    
    label = ['additional student',
             'additional postgrad. prof.',
             'additional grant',
             'additional paper']
    color = [(0, 0, 0), (.6,0,0), (.1, 1, .1), (.4,.4,1)]

    bottom = 0
    for dfi, labeli, colori in zip(df, label, color):
        print(labeli, colori)
        ax.bar(univ, f * dfi, bottom=bottom, label=labeli,
            align='edge', width=-0.8, color=colori, edgecolor='black')
        bottom += f * dfi

    ax.set_ylim(0, ax.get_ylim()[1])
    afdtable.inverse_legend(ax, title='improved 2020 metric')
    ax.set_ylabel('marginal earnings in all years [$10^6$ CLP]')

    ax2 = ax.twinx()
    ax2.set_ylim(0, ax.get_ylim()[1] / f)
    ax2.set_ylabel('2021 marginal earnings [$10^6$ CLP]')
        
    for tick in ax.get_xticklabels():
        tick.set_rotation(38)
        tick.set_horizontalalignment('right')
    ax.tick_params(axis='x', bottom=False)

    fig.tight_layout()

    fig.show()

    fig.savefig('../pdf/marginal-earnings.pdf')

tab = afdtable.read(2021)
incentives(tab)
