# this program calculate the angular part of the deformation 

from scipy.special import gamma
from scipy.special import hyp1f1
from sympy.physics.quantum.cg import CG
from sympy import S
import numpy as np



#I am not sure whether l1 and l2 has to satsify the selection rule
def Angular_dependence(m,k1,k2,L):
	j1=abs(k1)-0.5
	j2=abs(k2)-0.5
	l1=get_l(k1)[0]
	l2=get_l(k2)[0]
	if check(j1,j2,L)==0:
		#print 1
		return 0
	#if check(l1,l2,L)==0:
	#	#print 2
	#	return 0
	if ((l1+l2+L)%2)!=0:
		#print 3
		return 0
	# print get_m_dependence(m,k1,k2,L)
	# print generate_L_dependence(k1,k2,L)
	return float(get_m_dependence(m,k1,k2,L)*generate_L_dependence(k1,k2,L))   # L_dependence is an dictionary

def generate_L_dependence(k1,k2,L):
	j1=abs(k1)-0.5
	j2=abs(k2)-0.5
	c1=(CG(j2,-0.5,j1,0.5,L,0).doit()) 
	return (1.0/np.sqrt(4*np.pi))*head(j1)*head(j2)*c1/head(L)

def get_m_dependence(m,k1,k2,L):
	j1=abs(k1)-0.5
	j2=abs(k2)-0.5
	# print j2,m,j1,-m,L
	# print 'something',CG(j2,m,j1,-m,L,0).doit()
	return (-1)**(m+0.5)*(CG(j2,m,j1,-m,L,0).doit()) 


def head(a):
	return np.sqrt(2*a+1.0)



def check(l1,l2,L):
	return (l1<=(l2+L)) & (l1>=abs(l2-L))


def get_l(k):          # give kappa and get la,lb as output
	if k>0:
		la=k
		lb=k-1
	else:
		la=-k-1
		lb=-k
	return la,lb

