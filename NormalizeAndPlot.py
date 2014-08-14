import BDdb
db=BDdb.get_db('/Users/paigegiorla/Dropbox/BDNYCdb/BDNYC.db')
from matplotlib import pyplot as plt 
import numpy as np
import utilities as u
import astropy.units as q
import scipy.stats as s

def showme():
	
#!	open HD19467B and read in the arrays for wavelength, flux, and uncertainty
	file=open('/Users/paigegiorla/Code/Python/BDNYC/TDwarfplotting/HD19467B_project/HD19467NewSpect.txt')
	W_obj,F_obj,U_obj=[],[],[]
	for line in file:
		columns=line.split()
		wav=float(columns[0])
		flu=float(columns[1])
		un=float(columns[2])
		W_obj.append(wav) ; F_obj.append(flu) ; U_obj.append(un)
	file.close()	
	W_obj_t=W_obj[4:30]																		#t = trimmed
	F_obj_t=F_obj[4:30]
	U_obj_t=U_obj[4:30]

	myInt=2.3630917239963894e-17
	F_obj_tn = [x/myInt for x in F_obj_t]													#n = normalized
	U_obj_tn = [x/myInt for x in U_obj_t]
	before=np.array([W_obj_t,F_obj_t,U_obj_t])
	after=np.array([W_obj_t,F_obj_tn,U_obj_tn])
	after=u.norm_spec(after,before)																#to make sure all spaces in between points are the same length
	after=[i.value for i in after]
# 	plt.errorbar(before[0],before[1],yerr=before[2])
# 	plt.errorbar(after[0],after[1],yerr=after[2])
 	plt.plot(before[0],before[1])
	plt.show()
	plt.clf()
	plt.plot(after[0],after[1])
	plt.show()