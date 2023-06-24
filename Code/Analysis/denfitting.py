import pandas as pd
from scipy.optimize import curve_fit
import numpy as np
import matplotlib.pyplot as plt
import glob, os, sys
import regex as re
import math

"""
"""
def func(xvec,x1,x3,s,a,b,c):

	yvec=[]
	for x in xvec:

		if 0<x<=x1:
			y=c
		elif x1<x<=x1+s:
			y=a*math.erfc(-b*(x-x1-s))
		elif x1+s<x<=x3:
			y=a
		elif x3<x<=x3+s:
			y=a*math.erfc(-b*(x3-x))
		else:
			y=c

		yvec.append(y)

	return yvec
"""
"""
prot = sys.argv[1]
seq = sys.argv[2]

sub1=int(sys.argv[3])
sub2=int(sys.argv[4])

x1_1=float(sys.argv[5]) 
x3_1=float(sys.argv[6])
x1_2=float(sys.argv[7])
x3_2=float(sys.argv[8])

outF = open("PDData-%s_%s.txt"%(prot,seq), "w")
fig, axs = plt.subplots(sub1, sub2, sharex=True, sharey=True)
fig.suptitle('Concentration Profile for NVT Slab Simulations for different temperatures')

#axs= axs.flatten()

i=0

for ax in axs.flatten():
    ax.set(xlabel='Z-direction slab cross section', ylabel='Concentration')

for ax in axs.flatten():
    ax.label_outer()

outF.write("Temperature,ConcLow,ConcHigh,block"+"\n")

for file in glob.glob("mean_%s_%s_*-*.density"%(prot,seq)):
#	print("ARCHIVO: ",file)
	match = re.search("mean_%s_%s_(.*)-(.*).density"%(prot,seq), file)
	temp  = match.group(2)
	block = match.group(1)
#	print(temp)
	
	i=i+1

	yarray=pd.read_csv('mean_%s_%s_%s-%03d.density'%(prot,seq,block,int(temp)))
#	print(yarray)
	ydata=yarray.values.flatten()
#	print(ydata, ydata.shape)
	
	xdata=np.linspace(0,1,len(ydata))
#	print(xdata, xdata.shape)
	#Here, write the best guess for the parameters to aid the optimization algorithm

	print(file)
	print(temp)

	popt, pcov = curve_fit(func, xdata, ydata,p0=(0.02,0.6,0.06,1,20,0.),\
		bounds=([x1_1,x3_1,0.05,0.2,1,0.],[x1_2,x3_2,0.5,np.inf,np.inf,0.15]))

	outF.write(str(temp)+" ,  "+str(popt[3])+" ,  "+str(popt[5])+" ,  "+block+"\n")
	#outF.write(str(temp)+" ,  "+str(popt[5])+"\n")

#	print(popt)

	k= int((i-1) % sub1)
	j= int(np.floor((i-1)/ sub2))
	print(i,k,j)
	if j==0 and k==0:
		axs[j,k].plot(xdata, ydata, label='Simulation Data')
		axs[j,k].plot(xdata, func(xdata, *popt), label='Fitting Function')
		axs[j,k].legend()
	else:
		axs[j,k].plot(xdata, ydata)
		axs[j,k].plot(xdata, func(xdata, *popt))

	T=int(temp)
	axs[j,k].set_title('block='+block+' T=%d K, '%T)

plt.savefig("%s_%s.png"%(prot,seq),dpi=200)
#plt.show()
