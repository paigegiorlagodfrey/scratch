import BDdb
from matplotlib import pyplot as plt 
import numpy as np
import utilities as u
import astropy.units as q
import scipy.stats as s

db=BDdb.get_db('/Users/paigegiorla/Dropbox/BDNYCdb/BDNYC.db')

def showme():
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

	wavelength_points=[1.1826776, 1.2093873000000002, 1.236097, 1.2628067, 1.2895164, 1.3162261000000002, 1.3429358, 1.3696455, 1.3963552, 1.4230649000000002, 1.4497746, 1.4764843, 1.503194, 1.5299037000000002, 1.5566134, 1.5833231, 1.6100328, 1.6367425, 1.6634522, 1.6901619, 1.7168716]		
	sources=db.query.execute("select distinct source_id from spectral_types where spectral_type>=20 and spectral_type<=30 and regime='IR'").fetchall()
	
	sourcelist=[]
	
	index=0
	while index<len(sources):
		sourcelist.append(sources[index][0])
		index=index+1
	ids=[]
	for i in sourcelist:
		id=db.query.execute("SELECT distinct id FROM spectra WHERE source_id='{}' AND instrument_id=6".format(i)).fetchone()
		ids.append(id)
	ids=filter(None,ids)
	
	idlist=[]
	
	index=0
	while index<len(ids):
		idlist.append(ids[index][0])
		index=index+1
	
		
	print len(idlist)
	output=[]
	chisquare=[]
	sptlist=[]
	chilist=[]
	for j in idlist:
		name,W,F,unc=db.query.execute("SELECT source_id, wavelength, flux, unc from spectra where id='{}'".format(j)).fetchone()
		binned_wl=0
		jump=0.0267097
		start=1.1826776+(jump/2) #start half way after 1st flux point
		times=21
		end=1.73022645
		spectype=db.query.execute("SELECT spectral_types.spectral_type from spectral_types join spectra on spectral_types.source_id=spectra.source_id where spectral_types.regime='IR' and spectra.id='{}'".format(j)).fetchone()
		
#	normalize the flux incase it's not already
		NF=F*1.27/np.interp(1.27,W,F)
		NU=unc*1.27/np.interp(1.27,W,unc)
		flux2=0
		newflux=[]
		unc2=0
		newunc=[]
		w1=start+(jump/2) #first wl point is at start+the step/2
		newwl=[]
		wavelength=0

		for i in range(0,21):
			binned_wl=jump + start #range is half below and half above w1
			flux2=sum(NF[np.where(np.logical_and(W>start,W<binned_wl))])
			unc2=np.sqrt(sum(NU[np.where(np.logical_and(W>start,W<binned_wl))]**2))	
			start=binned_wl
		
			newflux.append(flux2)
			newunc.append(unc2)
			wavelength=w1+(i*jump)
			newwl.append(wavelength)	
# 
# 		G,C=goodness([wava,flua,una],[newwl,newflux,newunc])
# 		newflux=C*np.array(newflux)
# 		newunc=C*np.array(newunc)
		
# 		plt.errorbar(newwl,newflux,yerr=newunc,marker="o")
# 		plt.errorbar(wava,flua,yerr=una,marker="o")
# 	calculate a chi squared
# 		wava*=q.um
# 		flua*=q.erg/q.s/q.cm**2/q.AA
# 		una*=q.erg/q.s/q.cm**2/q.AA
# # 	

		template=np.array([newwl,newflux,newunc])	
		object=np.array([wava,flua,una])
	
		template=u.norm_spec(template,object)
# 		plt.errorbar(template[0],template[1],yerr=template[2],marker="o")
# 		plt.errorbar(object[0],object[1],yerr=object[2],marker="o")
# 		plt.savefig('/Users/paigegiorla/Code/Python/BDNYC/TDwarfplotting/HD19467B_project/plots/'+'{}'.format(name)+ '.pdf')
# 		plt.clf()

		try:
			chi,p=s.chisquare(template[1],object[1])   #len of template and object dont match
			output.append([chi,name,spectype])
			sptlist.append(spectype)
			chilist.append(chi)

		except ValueError: pass
	
	chisquare=[sptlist,chilist]

	plt.scatter(chisquare[0],chisquare[1])
# 	plt.savefig('/Users/paigegiorla/Code/Python/BDNYC/TDwarfplotting/HD19467B_project/plots/chisquare_distribution.pdf')
# 	plt.clf()	
	plt.show()
# 	answer=([np.argmin(output[0]),name,spectype])
	




		