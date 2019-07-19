#!/usr/bin/env python

import sys
import pandas as pd
from sklearn import linear_model
import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns

F1, F2, TSV, PDF = sys.argv[1:]
df = pd.read_csv(TSV, sep="\t")
print(df.describe())

regr = linear_model.LinearRegression()
X = df[F1].values.reshape(-1, 1)
y = df[F2].values.reshape(-1, 1)
regr.fit(X, y)
a = float(regr.coef_[0])
b = float(regr.intercept_)
eq = "y = {:.3f}x + {:.3f}".format(a, b)
print(eq)

pages = PdfPages(PDF)

g = sns.lmplot(F1, F2, df, fit_reg=True, aspect=1.5, ci=95, scatter_kws={"s": 1})
props = dict(boxstyle='round', alpha=0.5, color=sns.color_palette()[0])
textstr = '$' + eq + '$'
g.ax.text(0.7, 0.05, textstr, transform=g.ax.transAxes, fontsize=14, bbox=props)

pages.savefig()
plt.close()
pages.close()
