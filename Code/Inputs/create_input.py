###SCRIPT FOR TIA1 INPUTS

import numpy as np
import pandas as pd

#Defining the sequences numerically

dic={'M':1,'G':2,'K':3,'T':4,'R':5,'A':6,'D':7,'E':8,'Y':9,'V':10,'L':11,'Q':12,'W':13,'F':14,'S':15,'H':16,'N':17,'P':18,'C':19,'I':20}
dicc={'M':0,'G':0,'K':0.75,'T':0,'R':0.75,'A':0,'D':-0.75,'E':0.75,'Y':0,'V':0,'L':0,'Q':0,'W':0,'F':0,'S':0,'H':0.375,'N':0,'P':0,'C':0,'I':0}

TIA_WT     = "MINPVQQQNQIGYPQPYGQWGQWYGNAQQIGQYMPNGWQVPAYGMYGQAWNQQGFNQTQSSAPWMGPNYGVQPPQGQNGSMLPNQPSGYRVAGYETQ"
TIA_YtoF   = "MINPVQQQNQIGFPQPFGQWGQWFGNAQQIGQFMPNGWQVPAFGMFGQAWNQQGFNQTQSSAPWMGPNFGVQPPQGQNGSMLPNQPSGFRVAGFETQ";
TIA_YtoF_2 = "MINPVQQQNQIGFPQPYGQWGQWFGNAQQIGQYMPNGWQVPAFGMYGQAWNQQGFNQTQSSAPWMGPNFGVQPPQGQNGSMLPNQPSGYRVAGYETQ";
TIA_m6Y    = "MINPVQQQNQIGQPQPGGQWGQWSGNAQQIGQYMPNGWQVPAQGMYGQAWNQQGFNQTQSSAPWMGPNGGVQPPQGQNGSMLPNQPSGYRVAGYETQ"
TIA_YtoW   = "MINPVQQQNQIGWPQPWGQWGQWWGNAQQIGQWMPNGWQVPAWGMWGQAWNQQGFNQTQSSAPWMGPNWGVQPPQGQNGSMLPNQPSGWRVAGWETQ";
TIA_YtoW_2 = "MINPVQQQNQIGWPQPYGQWGQWWGNAQQIGQYMPNGWQVPAWGMYGQAWNQQGFNQTQSSAPWMGPNWGVQPPQGQNGSMLPNQPSGYRVAGYETQ";
TIA_WtoF   = "MINPVQQQNQIGYPQPYGQFGQFYGNAQQIGQYMPNGFQVPAYGMYGQAFNQQGFNQTQSSAPFMGPNYGVQPPQGQNGSMLPNQPSGYRVAGYETQ";
TIA_WtoY   = "MINPVQQQNQIGYPQPYGQYGQYYGNAQQIGQYMPNGYQVPAYGMYGQAYNQQGFNQTQSSAPYMGPNYGVQPPQGQNGSMLPNQPSGYRVAGYETQ";
TIA_GtoS   = "MINPVQQQNQISYPQPYSQWSQWYSNAQQISQYMPNSWQVPAYSMYSQAWNQQSFNQTQSSAPWMSPNYSVQPPQSQNSSMLPNQPSSYRVASYETQ";
TIA_GtoS_3 = "MINPVQQQNQISYPQPYSQWGQWYSNAQQISQYMPNGWQVPAYSMYSQAWNQQGFNQTQSSAPWMSPNYSVQPPQGQNSSMLPNQPSSYRVAGYETQ";
TIA_GtoS_9 = "MINPVQQQNQISYPQPYGQWGQWYSNAQQIGQYMPNGWQVPAYSMYGQAWNQQGFNQTQSSAPWMSPNYGVQPPQGQNSSMLPNQPSGYRVAGYETQ";
TIA_StoG   = "MINPVQQQNQIGYPQPYGQWGQWYGNAQQIGQYMPNGWQVPAYGMYGQAWNQQGFNQTQGGAPWMGPNYGVQPPQGQNGGMLPNQPGGYRVAGYETQ";
TIA_QtoN   = "MINPVNNNNNIGYPNPYGNWGNWYGNANNIGNYMPNGWNVPAYGMYGNAWNNNGFNNTNSSAPWMGPNYGVNPPNGNNGSMLPNNPSGYRVAGYETN";
TIA_QtoN_2 = "MINPVNNNNQIGYPQPYGQWGQWYGNANNIGQYMPNGWQVPAYGMYGQAWNQQGFNNTQSSAPWMGPNYGVNPPQGNNGSMLPNNPSGYRVAGYETN";
TIA_NtoQ   = "MIQPVQQQQQIGYPQPYGQWGQWYGQAQQIGQYMPQGWQVPAYGMYGQAWQQQGFQQTQSSAPWMGPQYGVQPPQGQQGSMLPQQPSGYRVAGYETQ";

prot_names=[
'TIA_WT',
'TIA_YtoF',
'TIA_YtoF_2',
'TIA_m6Y',
'TIA_YtoW',
'TIA_YtoW_2',
'TIA_WtoF',
'TIA_WtoY',
'TIA_GtoS',
'TIA_GtoS_3',
'TIA_GtoS_9',
'TIA_StoG',
'TIA_QtoN',
'TIA_QtoN_2',
'TIA_NtoQ']

prot_aa=[
TIA_WT,
TIA_YtoF,
TIA_YtoF_2,
TIA_m6Y,
TIA_YtoW,
TIA_YtoW_2,
TIA_WtoF,
TIA_WtoY,
TIA_GtoS,
TIA_GtoS_3,
TIA_GtoS_9,
TIA_StoG,
TIA_QtoN,
TIA_QtoN_2,
TIA_NtoQ]
seq_list=[[dic[x] for x in list(a)] for a in prot_aa]

seq_c=[[dicc[x] for x in list(a)] for a in prot_aa]

print(seq_list)

d={'Sequence Name':prot_names, 'Sequence AA':prot_aa, 'Sequence Num':seq_list, 'Sequence Charges':seq_c}

df=pd.DataFrame(data=d)

print(df.head())

#Size of the box and positions
system_size=300.0

pos_file=open('positions.dat','r')
positions = pos_file.read().splitlines()

#Write files

for ind in df.index:

	name=df['Sequence Name'][ind]

	seqn=df['Sequence Num'][ind]

	charges=df['Sequence Charges'][ind]

	with open("%s.config" % name,'w') as idata:

		idata.write('Configuration for FUS variant\n\n')

		idata.write('{} atoms\n'.format(len(seqn)))
		idata.write('{} atom types\n\n'.format(40))
		idata.write('{} bonds\n'.format(len(seqn)-1))
		idata.write('{} bond types\n\n'.format(1))


		idata.write('{} {} xlo xhi\n'.format(0.0, system_size))
		idata.write('{} {} ylo yhi\n'.format(0.0, system_size))
		idata.write('{} {} zlo zhi\n'.format(0.0, system_size))
		idata.write('\n')


		with open("masses.dat", "r") as f1:
    			t1 = f1.readlines()

		idata.writelines(t1)

		idata.write('\n')
		idata.write('Atoms\n\n')
		for i,pos in enumerate(positions):
			idata.write('{} 1 {} {} {}\n'.format(i+1,seqn[i],charges[i],pos))
			
			if i==len(seqn)-1:

				break

		idata.write('\n')	
		idata.write('Bonds\n\n')

		for i in range(len(seqn)-1):

			idata.write('{} {} {} {}\n'.format(i+1,1,i+1,i+2))

