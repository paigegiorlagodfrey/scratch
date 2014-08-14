from matplotlib import pyplot as plt 
import numpy as np
import BDdb
db=BDdb.get_db('/Users/paigegiorla/Dropbox/BDNYCdb/BDNYC.db')
from matplotlib import pyplot as plt
from matplotlib import gridspec
import numpy as np
import utilities as u
import astropy.units as q
import scipy.stats as s
import modules as m
import astrotools as a
import chisquareplot as ch
import figure_paper as fp

def showme(sptlist,chilist):
	
	plt.figure(figsize=(8,10))
	gs = gridspec.GridSpec(2, 1, height_ratios=[1,2]) 
	ax1 = plt.subplot(gs[0])
# 	ax1=plt.subplot(211)
	ch.showme(sptlist,chilist,ax1)
	ax2 = plt.subplot(gs[1])
# 	ax2=plt.subplot(212)
	fp.showme(ax2)
	plt.savefig('/Users/paigegiorla/Code/Python/BDNYC/TDwarfplotting/HD19467B_project/plots/twoplot.pdf')


	
	
	
	
	
	
	
	