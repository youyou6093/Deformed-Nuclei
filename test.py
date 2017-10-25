import numpy as np
from scipy.integrate import simps



f1 = open("test_matrix.txt")
f2 = open("matrix_test_python.txt")
M1 = np.zeros(352*352)
M2 = np.zeros(396*396)

pos = 0
for i in f1:
	M1[pos] = float(i)
	pos += 1
pos = 0
for i in f2:
	M2[pos] = float(i)
	pos += 1


M1 = M1.reshape((352,352))

M2 = M2.reshape((396,396))

eig1=np.linalg.eig(M1)  #diagnolize the matrix and get the eigvalue
eigvt1=eig1[1].T                   #get the eigenvector
eig_pair1=zip(eig1[0] * 197.326 - 939,eigvt1)       #pair the eigenvalue with the corresponding eigenvector
energy1=sorted(eig_pair1,key=lambda x:x[0])         #sort the energy from small to big
energy1 = [ i for i in energy1 if (i[0] < 5) & (i[0] > -939)]


eig2=np.linalg.eig(M2)  #diagnolize the matrix and get the eigvalue
eigvt2=eig2[1].T                   #get the eigenvector
eig_pair2=zip(eig2[0] * 197.326 -939 ,eigvt2)       #pair the eigenvalue with the corresponding eigenvector
energy2=sorted(eig_pair2,key=lambda x:x[0])         #sort the energy from small to big
energy2 = [ i for i in energy2 if (i[0] < 5) & (i[0] > -939)]

print energy1[0][0]
print energy2[0][0]