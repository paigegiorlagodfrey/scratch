import BDdb
db=BDdb.get_db('/Users/paigegiorla/Dropbox/BDNYCdb/BDNYC.db')
from matplotlib import pyplot as plt 
import numpy as np
import utilities as u
import astropy.units as q
import scipy.stats as s
import modules as m
import astrotools as a

def showme(sptlist,chilist):

	sptlist,chilist = (list(x) for x in zip(*sorted(zip(sptlist,chilist))))
	print len(sptlist), len(chilist), type(sptlist), type(chilist)
	sptlist=sptlist[57:]
	sptlist=[k-20.0 for k in sptlist]
	print sptlist
	chilist=chilist[57:]
	print chilist
	chisquare=np.array([sptlist,chilist])
	
 	fig,ax=plt.subplots()
	ax.scatter(chisquare[0],chisquare[1])
	chisquare=zip(*sorted(zip(*chisquare)))
	
	Ts=['T4','T5','T6','T7','T8','T9']
	
	pfit = np.polyfit(chisquare[0],chisquare[1], 2)   # Fit a 2nd order polynomial to (x, y) data
	yfit = np.polyval(pfit, chisquare[0])   # Evaluate the polynomial at x
	
	ax.plot(chisquare[0], yfit)
	ax.set_xlim(4.5,9.5)
	ax.set_ylim(12,75)
	ax.set_xticklabels(Ts)
	ax.set_xlabel('Spectral Type',fontsize="x-large")
	ax.set_ylabel('Chi Squared',fontsize="x-large")
	ax.tick_params(labelsize="large")
	plt.savefig('/Users/paigegiorla/Code/Python/BDNYC/TDwarfplotting/HD19467B_project/plots/chisquare_paperfigure2.pdf')
	plt.clf()	
# 	return p

	
