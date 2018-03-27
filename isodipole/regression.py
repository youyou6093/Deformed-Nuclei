from scipy import stats
import numpy as np

f = open(raw_input())
xdata = []
ydata = []
for i in f:
	x, y = [float(num) for num in i.split()]
	xdata.append(x)
	ydata.append(y)

print xdata, ydata

print stats.linregress(xdata, ydata)
