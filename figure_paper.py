import BDdb
db=BDdb.get_db('/Users/paigegiorla/Dropbox/BDNYCdb/BDNYC.db')
from matplotlib import pyplot as plt 
import numpy as np
import utilities as u
import astropy.units as q
import scipy.stats as s
import modules as m
import astrotools as a 

def showme(ax):
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
	lamb=np.gradient(W_obj)
	xer=lamb/2

# 	T6
	original, template6=m.bin_down(7982,0.0267097,1.1826776,21)	
	source,spectype=db.query.execute("SELECT spectral_types.source_id, spectral_types.spectral_type from spectral_types join spectra on spectral_types.source_id=spectra.source_id where spectral_types.regime='IR' and spectral_types.spectral_type>=20 and spectral_types.spectral_type<=30 and spectra.id=7982").fetchone()
	spectype=a.specType(spectype)
	template_plot26=u.norm_spec(template6,object)												#Joe's code to normalize plots to each other
	template_plot26=[i.value for i in template_plot26]
	ax.plot(template_plot26[0],template_plot26[1],'c:',linewidth=2,label=str(spectype))
	
		
# 	T7
	original, template7=m.bin_down(3121,0.0267097,1.1826776,21)	
	source,spectype=db.query.execute("SELECT spectral_types.source_id, spectral_types.spectral_type from spectral_types join spectra on spectral_types.source_id=spectra.source_id where spectral_types.regime='IR' and spectral_types.spectral_type>=20 and spectral_types.spectral_type<=30 and spectra.id=3121").fetchone()
	spectype=a.specType(spectype)
	template_plot27=u.norm_spec(template7,object)												#Joe's code to normalize plots to each other
	template_plot27=[i.value for i in template_plot27]
	ax.plot(template_plot27[0],template_plot27[1],'g-.',linewidth=2,label=str(spectype))

# 	T7.5
	original, template=m.bin_down(5308,0.0267097,1.1826776,21)	
	source,spectype=db.query.execute("SELECT spectral_types.source_id, spectral_types.spectral_type from spectral_types join spectra on spectral_types.source_id=spectra.source_id where spectral_types.regime='IR' and spectral_types.spectral_type>=20 and spectral_types.spectral_type<=30 and spectra.id=5308").fetchone()
	spectype=a.specType(spectype)
	template_plot2=u.norm_spec(template,object)												#Joe's code to normalize plots to each other
	template_plot2=[i.value for i in template_plot2]
	ax.plot(template_plot2[0],template_plot2[1],'r--',linewidth=2,label=str(spectype))
	
# 	T8
	with open('/Users/paigegiorla/Code/Python/BDNYC/TDwarfplotting/HD19467B_project/LateTdata/WISEJ050003.04-122343.2_SpeX.txt') as f:
		lines_after_17 = f.readlines()[7:]  
		Wtxt,Ftxt,Utxt=[],[],[]
		for line in lines_after_17:
			columns=line.split()
			wavx=float(columns[0])
			flux=float(columns[1])
			unx=float(columns[2])
			Wtxt.append(wavx) ; Ftxt.append(flux) ; Utxt.append(unx)
		f.close()	
	spectype=28
	spectype=a.specType(spectype)
	binned_wl=0																	#defined by the wavelength points of HD19467B data
	start=1.1826776-(0.0267097/2) 															#start half way before 1st flux point
	W,F,U=u.scrub([Wtxt,Ftxt,Utxt])
	W=W.value
	F=F.value
	U=U.value
	original=np.array([Wtxt,Ftxt,Utxt])
	flux2=0
	W_temp, F_temp, U_temp=[],[],[]
	unc2=0																	#first wl point is at start+the step/2
	wavelength=0
	for i in range(21):
		binned_wl=0.0267097 + start 															#range is half below and half above w1
		flux2=sum(F[np.where(np.logical_and(W>start,W<binned_wl))])					#binning flux
		unc2=np.sqrt(sum(U[np.where(np.logical_and(W>start,W<binned_wl))]**2))			#bin uncertainty in quadrature
		start=binned_wl																	#set new range start
		F_temp.append(flux2)
		U_temp.append(unc2)
		wavelength=1.1826776+(i*0.0267097)															#this is where the flux point will go on plot
		W_temp.append(wavelength)	
	template28=np.array([W_temp, F_temp, U_temp])
	template_plot28=u.norm_spec(template28,object)												#Joe's code to normalize plots to each other
	template_plot28=[i.value for i in template_plot28]
	ax.plot(template_plot28[0],template_plot28[1],'b',linewidth=3,label=str(spectype))
	
#	T8.5
	with open('/Users/paigegiorla/Code/Python/BDNYC/TDwarfplotting/HD19467B_project/LateTdata/WISEJ004945.61+215120.0_SpeX.txt') as f:
		lines_after_17 = f.readlines()[7:]  
		Wtxt,Ftxt,Utxt=[],[],[]
		for line in lines_after_17:
			columns=line.split()
			wavx=float(columns[0])
			flux=float(columns[1])
			unx=float(columns[2])
			Wtxt.append(wavx) ; Ftxt.append(flux) ; Utxt.append(unx)
		f.close()	
	spectype=28.5
	spectype=a.specType(spectype)
	binned_wl=0																	#defined by the wavelength points of HD19467B data
	start=1.1826776-(0.0267097/2) 															#start half way before 1st flux point
	W,F,U=u.scrub([Wtxt,Ftxt,Utxt])
	W=W.value
	F=F.value
	U=U.value
	original=np.array([Wtxt,Ftxt,Utxt])
	flux2=0
	W_temp, F_temp, U_temp=[],[],[]
	unc2=0																	#first wl point is at start+the step/2
	wavelength=0
	for i in range(21):
		binned_wl=0.0267097 + start 															#range is half below and half above w1
		flux2=sum(F[np.where(np.logical_and(W>start,W<binned_wl))])					#binning flux
		unc2=np.sqrt(sum(U[np.where(np.logical_and(W>start,W<binned_wl))]**2))			#bin uncertainty in quadrature
		start=binned_wl																	#set new range start
		F_temp.append(flux2)
		U_temp.append(unc2)
		wavelength=1.1826776+(i*0.0267097)															#this is where the flux point will go on plot
		W_temp.append(wavelength)	
	template29=np.array([W_temp, F_temp, U_temp])
	template_plot29=u.norm_spec(template29,object)												#Joe's code to normalize plots to each other
	template_plot29=[i.value for i in template_plot29]
	ax.plot(template_plot29[0],template_plot29[1],'m--',linewidth=2,label=str(spectype))
	
	
# 	T9
	with open('/Users/paigegiorla/Code/Python/BDNYC/TDwarfplotting/HD19467B_project/LateTdata/T9_UGPS0722-0540_STD.txt') as f:
		lines_after_17 = f.readlines()[7:]  
		Wtxt,Ftxt,Utxt=[],[],[]
		for line in lines_after_17:
			columns=line.split()
			wavx=float(columns[0])
			flux=float(columns[1])
			unx=float(columns[2])
			Wtxt.append(wavx) ; Ftxt.append(flux) ; Utxt.append(unx)
		f.close()	
	spectype=29
	spectype=a.specType(spectype)
	binned_wl=0																	#defined by the wavelength points of HD19467B data
	start=1.1826776-(0.0267097/2) 															#start half way before 1st flux point
	W,F,U=u.scrub([Wtxt,Ftxt,Utxt])
	W=W.value
	F=F.value
	U=U.value
	original=np.array([Wtxt,Ftxt,Utxt])
	flux2=0
	W_temp, F_temp, U_temp=[],[],[]
	unc2=0																	#first wl point is at start+the step/2
	wavelength=0
	for i in range(21):
		binned_wl=0.0267097 + start 															#range is half below and half above w1
		flux2=sum(F[np.where(np.logical_and(W>start,W<binned_wl))])					#binning flux
		unc2=np.sqrt(sum(U[np.where(np.logical_and(W>start,W<binned_wl))]**2))			#bin uncertainty in quadrature
		start=binned_wl																	#set new range start
		F_temp.append(flux2)
		U_temp.append(unc2)
		wavelength=1.1826776+(i*0.0267097)															#this is where the flux point will go on plot
		W_temp.append(wavelength)	
	template29=np.array([W_temp, F_temp, U_temp])
	template_plot29=u.norm_spec(template29,object)												#Joe's code to normalize plots to each other
	template_plot29=[i.value for i in template_plot29]
	ax.plot(template_plot29[0],template_plot29[1],'g-.',linewidth=2,label=str(spectype))
	
	
	ax.errorbar(object[0],object[1],xerr=xer,yerr=object[2],fmt=None,ecolor='k',marker='o',label='HD19467B')
	ax.set_ylim(0.001,1.1)
	ax.set_xlim(1.15,1.75)
	ax.legend( loc='upper right',fontsize=10,prop={'size':15}, ncol=2, numpoints=1)
	ax.set_xlabel(r'Wavelength ($\mu$m)', fontsize="x-large")
	ax.set_ylabel('Normalized Flux', fontsize="x-large")
	ax.tick_params(labelsize="large")
	
# 	plt.savefig('/Users/paigegiorla/Code/Python/BDNYC/TDwarfplotting/HD19467B_project/plots/figure.pdf')
# 	plt.clf()
# 	return p2