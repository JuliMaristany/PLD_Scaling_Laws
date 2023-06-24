import pandas as pd
from scipy.optimize import curve_fit
import numpy as np
import matplotlib.pyplot as plt
import os
import matplotlib.colors as mcolors
from scipy import stats
import glob, os, sys
import regex as re
import csv
import seaborn as sns
import argparse
import matplotlib.colors as mcolors
import statistics
import matplotlib.cm as cm

pal = sns.color_palette('Set2')
pal2 = sns.color_palette("husl", 9)
colors = ['black']+pal.as_hex()+pal2.as_hex()

print(colors)

cmap=cm.get_cmap("Spectral")

######

def crit_temp(dx,Tc,a):
    y = Tc - a *( dx **3.06 )
    return y

def binodal(T,Tc,a,s4,cc):
	      R = (-T + Tc)/a
	      S = (-T + Tc)/s4+2.0*cc
	
	      rhoL = (R)**(1.0/3.06)/2.0 + S/2.0
	      rhoV = S - rhoL
	
	      return (rhoL,rhoV)

######

CLI=argparse.ArgumentParser()
CLI.add_argument(
  "--protein",  # name on the CLI - drop the `--` for positional/required parameters
  nargs="*",  # 0 or more values expected => creates a list
  type=str,
  default='FUS',  # default if nothing is provided
)

CLI.add_argument(
  "--variants",
  nargs="*",
  type=str,  
  default='WT',
)

#CLI.add_argument(
#  "--points",
#  nargs="*",
#  type=int,
#  default=[15,15,15,15,15,15,15,15,15,15,15,15,15,15],
#)

CLI.add_argument(
  "--type",
  nargs=1,
  type=str,
  default='test',
)

# parse the command line

args = CLI.parse_args()

######

prot = args.protein[0] #One protein at the time
seqs = args.variants #List of variants to include in the plot

print(prot,seqs)

ctemp_guess = 320

#fpoints= args.points
type=args.type[0]

i=0 

for seq in seqs:

	#####Read them into pandas
	
	df= pd.DataFrame(pd.read_csv('PDData-%s_%s.txt'%(prot,seq)))
	df=df.sort_values('Temperature')
	
	#####Add columns with the data I need
	
	df['DConc']=df['ConcLow']-df['ConcHigh']
	df['SConc']=df['ConcLow']+df['ConcHigh']
	
	#####Find critical temperature

	crit_temps_list=[]

	for fpoints in range(3,len(df)+1):
	
		p0=(ctemp_guess,5.)
		bound=([200,-np.inf],[500,np.inf])
		
		vals=fpoints

		df_trunc = df.tail(int(vals))
		xdata=df_trunc['DConc']
		ydata=df_trunc['Temperature'].astype(float)
			
		popt, pcov = curve_fit(crit_temp, xdata, ydata, p0,bounds=bound)
			
		Tcrit=popt[0]
		
		crit_temps_list.append(Tcrit)

		print(popt[0])

		def crit_den(ds,s4,cc):
		     y=Tcrit-s4*(ds-2*cc)
		     return y
			
		xdata=df_trunc['SConc']
		ydata=df_trunc['Temperature']
			
		p0=(1.,0.5)
		
		#bound=([-np.inf,0],[np.inf,1.])
		poptc, pcovc = curve_fit(crit_den, xdata, ydata, p0)

	
	Tcrit=crit_temps_list[-1]

	if i==0:
		Twt=Tcrit

	Terr=0.5*(max(crit_temps_list)-min(crit_temps_list))

	print(Tcrit, Terr)

	df=pd.melt(df, id_vars=['Temperature'], value_vars=['ConcLow', 'ConcHigh'],  var_name='arm', value_name='Concentration')

	print(df['Temperature'])

	Tscaled=2./175.*Tcrit-209./70.

	if seq[0]=='+' or seq[0]=='W':
		plt.scatter(x=df['Concentration'],y=df['Temperature'],color=cmap(Tscaled),label=seq)
	else:
		plt.scatter(x=df['Concentration'],y=df['Temperature'],color=cmap(Tscaled),label='-'+seq)
	plt.scatter(x=poptc[1],y=Tcrit,label=None,facecolor='none',edgecolor=cmap(Tscaled))
	plt.errorbar(poptc[1], Tcrit,yerr=Terr,ecolor=cmap(Tscaled),capsize=3,fmt='none')
	
	xlin=np.linspace(Tcrit-30,Tcrit,10000001)
		
	plt.plot(binodal(xlin,*popt,*poptc)[0],xlin,color=cmap(Tscaled))
	plt.plot(binodal(xlin,*popt,*poptc)[1],xlin,color=cmap(Tscaled))
	i=i+1
	
#######Write in data file
	op = open("ctemps.dat", "r")
	
	dt = pd.read_csv("ctemps_october.csv")
	df = pd.DataFrame(index=np.arange(len(dt)),columns=["Protein","Variant","Tc","Dc","Tc err","Dc err"])
	
	#print(df)
	###print(dt)
	
	for ind, row in dt.iterrows():
		#split = dt.loc[ind][0].split()
		#print(split)
	    
		df.loc[ind]["Protein"]=dt.loc[ind]["Protein"]
		df.loc[ind]["Variant"]=dt.loc[ind]["Variant"]
		#df.loc[ind]["Tc"]=Tcrit
		#df.loc[ind]["Dc"]=poptc[1]
		#print(ind,row)
		if row["Tc"] is not None:
			df.loc[ind]["Tc"]=dt.loc[ind]["Tc"]
			df.loc[ind]["Dc"]=dt.loc[ind]["Dc"]
			df.loc[ind]["Tc err"]=dt.loc[ind]["Tc err"]
			df.loc[ind]["Dc err"]=dt.loc[ind]["Dc err"]
	
	#print(df)
	indice=df.where((df["Protein"] == prot) & (df["Variant"]==seq)).dropna(how='all').index[0]

	print(np.sqrt(np.diag(pcovc))[1])

	###print(indice)
	df.loc[indice]["Tc"]=Tcrit
	df.loc[indice]["Dc"]=poptc[1]
	df.loc[indice]["Tc err"]=Terr
	df.loc[indice]["Dc err"]=np.sqrt(np.diag(pcovc))[1]

	df.to_csv("ctemps_october.csv")
	
plt.legend(loc='upper right')	
plt.xlabel('Concentration')	
plt.ylabel('Temperature')
plt.title('Binodal for '+prot+' variants.')
	
plt.savefig('binodals_'+prot+'_'+type+'.png', format="png")

plt.show()
	
####up_dt = []
####
####for r in dt:
####    print(r)
####    row = {'Sno': r['Sno'],
####           'Registration Number': r['Registration Number'],
####           'Name': r['Name'],
####           'RollNo': r['RollNo'],
####           'Status': 'P'}
####    up_dt.append(row)
####
####print(up_dt)
####
####op.close()
####op = open("ctemps.dat", "w", newline='')
####
####headers = ['Sno', 'Registration Number', 'Name', 'RollNo', 'Status']
####
####data = csv.DictWriter(op, delimiter=',', fieldnames=headers)
####
####data.writerow(dict((heads, heads) for heads in headers))
####
####data.writerows(up_dt)
####  
####op.close()

