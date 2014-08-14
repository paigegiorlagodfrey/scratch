import BDdb
from matplotlib import pyplot as plt 
import numpy as np
import utilities as u
import astropy.units as q
def showme(source):
	db=BDdb.get_db('/Users/paigegiorla/Dropbox/BDNYCdb/BDNYC.db')
	file=open('/Users/paigegiorla/Code/Python/BDNYC/TDwarfplotting/HD19467B_project/HD19467B.txt')
	wava=[]
	flua=[]
	una=[]
	for line in file:
		columns=line.split()
		wav=float(columns[0])
		flu=float(columns[1])
		un=float(columns[2])
		wava.append(wav)
		flua.append(flu)
		una.append(un)
	file.close()


# 	pick wavelength range between points on object plot
	W,F,unc=db.query.execute("SELECT wavelength,flux, unc from spectra where id='{}'".format(source)).fetchone()
	binned_wl=0
	jump=0.0267097
	start=1.1826776+jump
	times=21
	end=1.73022645
	
# 	tell it for each wavelength range, to add the flux points
	NF=F*1.27/np.interp(1.27,W,F)
	NU=unc*1.27/np.interp(1.27,W,unc)
	flux2=0
	newflux=[]
	unc2=0
	newunc=[]
	w1=start+(jump/2)
	newwl=[]
	wavelength=0
	for i in range(0,20):
		binned_wl=jump + start
		flux2=sum(NF[np.where(np.logical_and(W>start,W<binned_wl))])
		unc2=np.sqrt(sum(NU[np.where(np.logical_and(W>start,W<binned_wl))]**2))	
		start=binned_wl
		
		newflux.append(flux2)
		newunc.append(unc2)
		wavelength=w1+(i*jump)
		newwl.append(wavelength)	

# 	calculate a chi squared
	wava*=q.um
	flua*=q.erg/q.s/q.cm**2/q.AA
	una*=q.erg/q.s/q.cm**2/q.AA
	
	G,C=u.goodness([wava,flua,una],[newwl,newflux,newunc])
	print G

	newflux=C*np.array(newflux)
	newunc=C*np.array(newunc)

	plt.errorbar(newwl,newflux,yerr=newunc,marker="o")
 	plt.errorbar(wava,flua,yerr=una,marker="o")
 	
 	plt.show()
