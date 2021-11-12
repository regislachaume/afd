#! /usr/bin/env python3.9

import afdtable

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
    
    filename = '../tex/tab-incentives.tex'
    nl, tnl = "\n", "\\\\"
    mc = "\\multicolumn{2}{c}"
    caption = '\\caption'
    label = '\\label'

    with open(filename, 'w') as out:

        if include_caption:
            out.write('\\begin{table}\n')
            out.write(f"{caption}{{Additional earnings in {year} and present value of total earnings if a University has improved the following metrics in {year-1}: one additional full-time contract for a post-graduate professor during a year, one additional ongoing research grant, and one additional Web of Science (ex-ISI) publication. Assumptions are: the total State funding will continue to grow 2\% per year in real terms; a yearly depreciation of 5\% in real terms ($\\approx 8$% in pesos)) to determine the present value of future earnings.}}{nl}")
            out.write(f"{label}{{tab:incentives:2000}}{nl}")

        out.write('\\begin{tabular}{l rr rr rr}\n')
        out.write('\\hline\\hline\n')
        out.write(f"{'university':30} & {nl}")
        out.write(f"{mc + '{postgraduate staff}':>64} & {nl}")
        out.write(f"{mc + '{research grant}':>98} & {nl}")
        out.write(f"{mc + '{WoS publication}':>132} {tnl}{nl}")
        out.write(f"{'':30} & ")
        out.write(f"{year:13} & {'all years':15} & ")
        out.write(f"{year:13} & {'all years':15} & ")
        out.write(f"{year:13} & {'all years':15} {tnl}{nl}")
        out.write(f"{'':30} & ")
        out.write(f"{'':13} & {'[CLP]':15} & ")
        out.write(f"{'':13} & {'[CLP]':15} & ")
        out.write(f"{'':13} & {'[CLP]':15} {tnl}{nl}")
        out.write('\\hline\n')

        f = 1 / (1 - 0.95 * (1+p)/(1+d))

        for k in range(len(tab)):
            tab1 = afdtable.change_metrics(tab, k, 'Sp', increment=1)
            df31 = 1e3 * (tab1[k]['f'] - tab[k]['f'])
            tab1 = afdtable.change_metrics(tab, k, 'G', increment=1)
            df41 = 1e3 * (tab1[k]['f'] - tab[k]['f'])
            tab1 = afdtable.change_metrics(tab, k, 'P', increment=1)
            df51 = 1e3 * (tab1[k]['f'] - tab[k]['f'])
            u = tab[k]['University']
            out.write(f"{u:30} & ")
            out.write(f"{df31:13,.0f} & {df31*f:15,.0f} & ")
            out.write(f"{df41:13,.0f} & {df41*f:15,.0f} & ")
            out.write(f"{df51:13,.0f} & {df51*f:15,.0f} {tnl}{nl}")

        out.write('\\hline\n')
        out.write('\\end{tabular}\n')

        if include_caption:
            out.write('\\end{table}')
    
    print(f'Incentives written into {filename}')

tab = afdtable.read(2021)
incentives(tab)
