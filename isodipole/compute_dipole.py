import numpy as np
import os
from scipy.integrate import simps

fx = np.linspace(0, 20, 801)

def read_data(name):
	denv = []
	f = open(name, 'r')
	for i in f:
		rhov =  float(i.split()[3])
		denv.append(rhov)
	return np.array(denv)


# only take in L = 1 channel
def compute_integral(rhov):
	a = simps( rhov * fx ** 3, fx)
	b = np.sqrt(4 * np.pi / 3)
	return a * b


def main():
	# rhov_list = []
	# result = []
	# for i in names:
	# 	rhov_list.append(read_data(i))
	# 	result.append(compute_integral(read_data(i)))
	# print result
	path = raw_input()
	names= os.listdir(path)
	names = [i for i in names if i.split('.')[-1] == 'txt']
	#print names
	for i in names:
		print i,
		print compute_integral(read_data(path + i))

main()