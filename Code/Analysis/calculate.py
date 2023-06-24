import numpy as np      
import sys, os          
import subprocess
import pandas as pd
import glob, os
import regex as re

#os.chdir("new_data")

for file in glob.glob("tail_den-TIA_*-*.*"):
	print("archivo ",file)
	match = re.search(r"tail_den-(.*)_(.*)-(.*)\.(.*)", file)
	print(match)
	prot  = match.group(1)
	seq   = match.group(2)
	temp  = match.group(3)
	block = match.group(4)

	i=0

	#y=250+5*x

	print('Protein: ')
	print(prot)
	print('Sequence: ')
	print(seq)
	print('Temperatura: ')
	print(temp)
	print('Block: ')
	print(block)
	
	b = np.empty([1,100]) # Chunk Size given in density.profile
	c = np.empty([1,100])
	print(b.shape)
	b[:] = 0.0
	c[:] = 0.0
	
	for chunk in pd.read_csv(file,chunksize=101,comment='#', delim_whitespace = True):
	       
	       # iloc should start from 0 since the first line after the 
	       # last comment is not read (probably due to the delim_whitespace varibale)
	       # chunksize=101 first dimenstion of iloc goes from 0 to 100
	       # iloc[0:100] is a slice and goes from 0 to size-1

		a = chunk.iloc[0:100,2]       
		c = c + a.values

		d = chunk.iloc[0:100,0]
		if i==0:
			print(d.values[0])
		if i==1:
			print(d.values[0])
		i = i + 1
	
	j = 0
	
	b = (b + c)/(i+j)
	
	print(b)    
	
	np.savetxt('mean_'+str(prot)+'_'+str(seq)+'_'+str(block)+'-%03d.density' %int(temp),np.transpose(b))
