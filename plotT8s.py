import BDdb
db=BDdb.get_db('/Users/paigegiorla/Dropbox/BDNYCdb/BDNYC.db')
from matplotlib import pyplot as plt 
import numpy as np
import utilities as u
import astropy.units as q
import scipy.stats as s
import modules as m

def showme():
	filenamelist=['WISEJ225540.75-311842.0_SpeX','WISEJ222623.05+044004.0_SpeX','WISEJ171104.60+350036.8_SpeX','WISEJ165311.05+444422.8_SpeX','WISEJ132233.64-234016.8_SpeX','WISEJ062309.94-045624.6_SpeX','WISEJ050003.04-122343.2_SpeX','WISEJ025409.51+022358.6_SpeX','WISEJ024512.62-345047.8_SpeX']
	spectypelist=[28,28,28,28,28,28,28,28,28]
	for k in range(12):
		with open('/Users/paigegiorla/Code/Python/BDNYC/TDwarfplotting/HD19467B_project/LateTdata/'+'{}'.format(filenamelist[k])+'.txt') as f:
			lines_after_17 = f.readlines()[7:]  
			Wtxt,Ftxt,Utxt=[],[],[]
			for line in lines_after_17:
				columns=line.split()
				wavx=float(columns[0])
				flux=float(columns[1])
				unx=float(columns[2])
				Wtxt.append(wavx) ; Ftxt.append(flux) ; Utxt.append(unx)
			f.close()	
		array=np.array([Wtxt,Ftxt,Utxt])
		plt.plot(array[0],array[1],label=filenamelist[k])
	plt.legend( loc='upper right',fontsize=10,prop={'size':15})
	plt.show()