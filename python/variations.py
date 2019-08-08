#! /usr/bin/env python3

from matplotlib import pylab as pl
import numpy as np
from numpy import exp

fig = pl.figure(1, figsize=(10,4))
fig.subplots_adjust(left=0.06, right=0.99, top=0.99, wspace=0.22)
fig.clf()

xi1 = np.linspace(-2.776, 3.406)
y1 = (xi1/4 - 7/5)
rel1 = 3/4*y1**2 * (1 - (1 + xi1**2)/27)
abs1 = rel1 * exp(y1 ** 3)

xi2 = np.linspace(-4, 4)
y2 = (xi2/4 - 7/5)
rel2 = 3/4*y2**2 * (1 - (1 + xi2 ** 2)/27)
abs2 = rel2 * exp(y2 ** 3)

ax = fig.add_subplot(121)

ax.plot(xi1, rel1, 'k-', xi2, rel2, 'k:')
ax.set_xlabel('std. dev. from mean')
ax.set_ylabel('relative increase by unit of std. dev.')

ax = fig.add_subplot(122)
ax.semilogy(xi1, abs1, 'k-', xi2, abs2, 'k:')
ax.set_xlabel('std. dev. from mean')
ax.set_ylabel('absolute increase by unit of std. dev.')

fig.show()
fig.savefig('evaluation-coefficient-variation.pdf')
