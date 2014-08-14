import logging
import BDdb
db=BDdb.get_db('/Users/paigegiorla/Dropbox/BDNYCdb Backups/BDNYC(Paige_07_03_14).db')
from matplotlib import pyplot as plt 
import numpy as np
import warnings
import utilities as u
import astropy.units as q
import scipy.stats as s
import modules as m
import math
warnings.simplefilter('ignore')
def showme():
	file=open('/Users/paigegiorla/Code/Python/BDNYC/TDwarfplotting/HD19467B_project/HD19467B.txt')
	W_obj,F_obj,U_obj=[],[],[]
	for line in file:
		columns=line.split()
		wav=float(columns[0])
		flu=float(columns[1])
		un=float(columns[2])
		W_obj.append(wav) ; F_obj.append(flu) ; U_obj.append(un)
	file.close()	
	object=np.array([W_obj,F_obj,U_obj])
		
	sources=db.query.execute("select distinct source_id from spectral_types where spectral_type>=20 and spectral_type<=30 and regime='IR'").fetchall()
	sourcelist, ids, idlist=[],[],[]
	index=0
	while index<len(sources):
		sourcelist.append(sources[index][0])
		index=index+1
	for i in sourcelist:
		id=db.query.execute("SELECT distinct id FROM spectra WHERE source_id='{}' AND instrument_id=6".format(i)).fetchone()
		ids.append(id)

	ids=filter(None,ids)	
	index=0
	while index<len(ids):
		idlist.append(ids[index][0])
		index=index+1
# 	idlist=[8470] #7973
	print len(idlist)
	G=m.montecarlo(object,idlist)
	
	sptlist,chilist=[],[]
	for i in G:
		chi,spt=zip(*i)
		
		sptlist+=spt
		chilist+=chi
	print len(sptlist)
# 	sptlist=sptlist[57:]
# 	sptlist=[k-20.0 for k in sptlist]
# 	print sptlist
# 	chilist=chilist[57:]
# 	print chilist
# 	chisquare=np.array([sptlist,chilist])

	plt.scatter(sptlist,chilist,marker='.')
	
 	chilist,sptlist=zip(*sorted(zip(chilist,sptlist),key=lambda x:x[-1]))
	
# 	Ts=['T4','T5','T6','T7','T8','T9']
	
	pfit = np.polyfit(sptlist,chilist, 2)   # Fit a 2nd order polynomial to (x, y) data
	yfit = np.polyval(pfit, sptlist)   # Evaluate the polynomial at x
	
	plt.plot(sptlist, yfit)
# 	ax.set_xlim(4.5,9.5)
# 	ax.set_ylim(12,75)
# 	ax.set_xticklabels(Ts)
	plt.xlabel('Spectral Type')
	plt.ylabel('Chi Squared')
# 	ax.tick_params(labelsize="large")
# 	plt.savefig('/Users/paigegiorla/Code/Python/BDNYC/TDwarfplotting/HD19467B_project/plots/chisquare_paperfigure2.pdf')
# 	plt.clf()	
	plt.show()


# 	h=zip(*G)
# 	y=s.norm.pdf(h[0], loc=0.0, scale=1.0)
# 	new=[y,h[1]]
# 	logging.info(h)
# 	a=plt.subplot(2,1,1)		histogram plot
# 	a=plt.hist2d(h[1],h[0])
# 	b=plt.subplot(2,1,2)
#  	b=plt.scatter(h[1],h[0])
# 	b=plt.xlim(19.5,30)
# 	plt.savefig('/Users/paigegiorla/Code/Python/BDNYC/TDwarfplotting/HD19467B_project/plots/MCplots.pdf')