import BDdb
db=BDdb.get_db('/Users/paigegiorla/Dropbox/BDNYCdb Backups/BDNYC(Paige_07_03_14).db')
from matplotlib import pyplot as plt 
import numpy as np
import utilities as u
import astropy.units as q
import scipy.stats as s
import modules as m

def showme():
	
#!	open HD19467B and read in the arrays for wavelength, flux, and uncertainty
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
	print np.gradient(W_obj)																#to make sure all spaces in between points are the same length


#!	get source_id for all T dwarfs from database, turn it into a list and use list to get spectra_id for each source's spex data	
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
	print len(idlist)

#!	loop through each object
	chisquare, chilist, namelist, specidlist, sptlist=[],[],[],[],[]

	for j in idlist:
		original, template=m.bin_down(j,0.0267097,1.1826776,21)									#create arrays for ease of use		

		source,spectype=db.query.execute("SELECT spectral_types.source_id, spectral_types.spectral_type from spectral_types join spectra on spectral_types.source_id=spectra.source_id where spectral_types.regime='IR' and spectral_types.spectral_type>=20 and spectral_types.spectral_type<=30 and spectra.id='{}'".format(j)).fetchone()
	
		template_plot1=u.norm_spec(template,original)												#Joe's code to normalize plots to each other
		template_plot1=[i.value for i in template_plot1]
# #!		plot template over template's original flux
# 		plt.plot(template_plot1[0],template_plot1[1])
# 		plt.plot(original[0],original[1])
# 		plt.savefig('/Users/paigegiorla/Code/Python/BDNYC/TDwarfplotting/HD19467B_project/plots/binned_flux_with_original/'+'{}'.format(source)+ '.pdf')
# 		plt.clf()  	
		
		template_plot2=u.norm_spec(template,object)												#Joe's code to normalize plots to each other
		template_plot2=[i.value for i in template_plot2]
# #!		plot template over object
# 		plt.errorbar(template_plot2[0],template_plot2[1],yerr=template_plot2[2],marker="o")
# 		plt.errorbar(object[0],object[1],yerr=object[2],marker="o")
# 		plt.savefig('/Users/paigegiorla/Code/Python/BDNYC/TDwarfplotting/HD19467B_project/plots/template_over_object/'+'{}'.format(source)+ '_HD19467.pdf')
# 		plt.clf()

		template=u.norm_spec(template,object)
		template=[i.value for i in template]
		chi=sum(((template[1]-object[1])**2)/((template[2]**2)+(object[2]**2)))     	#both data sets have uncertainties
		namelist.append(source)
		sptlist.append(spectype)
		chilist.append(float(chi))
		specidlist.append(j)

# 	filenamelist=['WISEJ225540.75-311842.0_SpeX','WISEJ222623.05+044004.0_SpeX','WISEJ171104.60+350036.8_SpeX','WISEJ165311.05+444422.8_SpeX','WISEJ132233.64-234016.8_SpeX','WISEJ062309.94-045624.6_SpeX','WISEJ045853.89+643452.5_SpeX','WISEJ050003.04-122343.2_SpeX','WISEJ025409.51+022358.6_SpeX','T9_UGPS0722-0540_STD','WISEJ004945.61+215120.0_SpeX','WISEJ024512.62-345047.8_SpeX']
# 	spectypelist=[28,28,28,28,28,28,28.5,28,28,29,28.5,28]
# 	for k in range(12):
# 		with open('/Users/paigegiorla/Code/Python/BDNYC/TDwarfplotting/HD19467B_project/LateTdata/'+'{}'.format(filenamelist[k])+'.txt') as f:
# 			lines_after_17 = f.readlines()[7:]  
# 			Wtxt,Ftxt,Utxt=[],[],[]
# 			for line in lines_after_17:
# 				columns=line.split()
# 				wavx=float(columns[0])
# 				flux=float(columns[1])
# 				unx=float(columns[2])
# 				Wtxt.append(wavx) ; Ftxt.append(flux) ; Utxt.append(unx)
# 			f.close()	
# 		spectype=spectypelist[k]
# 		
# 		binned_wl=0																	#defined by the wavelength points of HD19467B data
# 		start=1.1826776-(0.0267097/2) 															#start half way before 1st flux point
# 		W,F,U=u.scrub([Wtxt,Ftxt,Utxt])
# 		W=W.value
# 		F=F.value
# 		U=U.value
# 		original=np.array([Wtxt,Ftxt,Utxt])
# 		flux2=0
# 		W_temp, F_temp, U_temp=[],[],[]
# 		unc2=0																	#first wl point is at start+the step/2
# 		wavelength=0
# 
# 		#!		loop through 21 iterations for the amount of wavelength ranges I have
# 		for i in range(21):
# 			binned_wl=0.0267097 + start 															#range is half below and half above w1
# 	
# 			flux2=sum(F[np.where(np.logical_and(W>start,W<binned_wl))])					#binning flux
# 			unc2=np.sqrt(sum(U[np.where(np.logical_and(W>start,W<binned_wl))]**2))			#bin uncertainty in quadrature
# 	
# 			start=binned_wl																	#set new range start
# 
# 			F_temp.append(flux2)
# 			U_temp.append(unc2)
# 			wavelength=1.1826776+(i*0.0267097)															#this is where the flux point will go on plot
# 			W_temp.append(wavelength)	
# 
# 		template=np.array([W_temp, F_temp, U_temp])
# 	
# 		template_plot1=u.norm_spec(template,original)												#Joe's code to normalize plots to each other
# 		template_plot1=[i.value for i in template_plot1]
# # #!		plot template over template's original flux
# 		plt.plot(template_plot1[0],template_plot1[1])
# 		plt.plot(original[0],original[1])
# 		plt.savefig('/Users/paigegiorla/Code/Python/BDNYC/TDwarfplotting/HD19467B_project/plots/binned_flux_with_original/'+'{}'.format(filenamelist[k])+ '.pdf')
# 		plt.clf()  	
# 	
# 		template_plot2=u.norm_spec(template,object)												#Joe's code to normalize plots to each other
# 		template_plot2=[i.value for i in template_plot2]
# # #!		plot template over object
# 		plt.errorbar(template_plot2[0],template_plot2[1],yerr=template_plot2[2],marker="o")
# 		plt.errorbar(object[0],object[1],yerr=object[2],marker="o")
# 		plt.savefig('/Users/paigegiorla/Code/Python/BDNYC/TDwarfplotting/HD19467B_project/plots/template_over_object/'+'{}'.format(filenamelist[k])+ '_HD19467.pdf')
# 		plt.clf()
# 		
# 		template=u.norm_spec(template,object)
# 		template=[i.value for i in template]
# 		chi=sum(((template[1]-object[1])**2)/((template[2]**2)+(object[2]**2)))     	#both data sets have uncertainties
# # 		print filenamelist[k],type(chi), spectype
# 		namelist.append(filenamelist[k])
# 		sptlist.append(spectypelist[k])
# 		chilist.append(float(chi))
# 		specidlist.append('None')

#------------------------------	
	output=np.array([chilist,namelist,specidlist,sptlist]) 	
	chilistfloat=np.array([float(n) for n in output[0]])
	print len(sptlist)
	top5=chilistfloat.argsort()[:10]
	print top5
	for i in top5:
		print output[0][i],output[1][i],output[2][i],output[3][i]
		if output[2][i]!='None':
			k=output[2][i]
			original, template_min=m.bin_down(k,0.0267097,1.1826776,21)									#create arrays for ease of use		
		
			source,spectype=db.query.execute("SELECT spectral_types.source_id, spectral_types.spectral_type from spectral_types join spectra on spectral_types.source_id=spectra.source_id where spectral_types.regime='IR' and spectral_types.spectral_type>=20 and spectral_types.spectral_type<=30 and spectra.id='{}'".format(k)).fetchone()
		
			template_min=u.norm_spec(template_min,object)
			template_min=[i.value for i in template_min]
			plt.errorbar(template_min[0],template_min[1],yerr=template_min[2],marker="o",label=str(spectype))
		else:
			with open('/Users/paigegiorla/Code/Python/BDNYC/TDwarfplotting/HD19467B_project/LateTdata/'+'{}'.format(output[1][i])+'.txt') as f:
				lines_after_17 = f.readlines()[7:]  
				Wtxt,Ftxt,Utxt=[],[],[]
				for line in lines_after_17:
					columns=line.split()
					wavx=float(columns[0])
					flux=float(columns[1])
					unx=float(columns[2])
					Wtxt.append(wavx) ; Ftxt.append(flux) ; Utxt.append(unx)
				f.close()	
			spectype=output[3][i]
		
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

			#!		loop through 21 iterations for the amount of wavelength ranges I have
			for i in range(21):
				binned_wl=0.0267097 + start 															#range is half below and half above w1
	
				flux2=sum(F[np.where(np.logical_and(W>start,W<binned_wl))])					#binning flux
				unc2=np.sqrt(sum(U[np.where(np.logical_and(W>start,W<binned_wl))]**2))			#bin uncertainty in quadrature
	
				start=binned_wl																	#set new range start

				F_temp.append(flux2)
				U_temp.append(unc2)
				wavelength=1.1826776+(i*0.0267097)															#this is where the flux point will go on plot
				W_temp.append(wavelength)	

			template=np.array([W_temp, F_temp, U_temp])
			template_min=u.norm_spec(template,object)
			template_min=[i.value for i in template_min]
			plt.errorbar(template_min[0],template_min[1],yerr=template_min[2],marker="o",label=str(spectype))
	plt.legend( loc='upper right')
	plt.xlim(1.15,1.75)
	plt.ylim(0.001,1.2)
	plt.errorbar(object[0],object[1],yerr=object[2],marker="o",label='HD19467B')
	plt.savefig('/Users/paigegiorla/Code/Python/BDNYC/TDwarfplotting/HD19467B_project/plots/bestfits_w_object.pdf')
 	plt.clf()	
	
 	return sptlist,chilist
	chisquare=[sptlist,chilist]

#!	plot spectral type vs chi square
# 	plt.scatter(chisquare[0],chisquare[1])
# 	chisquare=zip(*sorted(zip(*chisquare)))
# 	pfit = np.polyfit(chisquare[0],chisquare[1], 3)   # Fit a 2nd order polynomial to (x, y) data
# 	yfit = np.polyval(pfit, chisquare[0])   # Evaluate the polynomial at x
# 	plt.plot(chisquare[0], yfit)
# 	plt.savefig('/Users/paigegiorla/Code/Python/BDNYC/TDwarfplotting/HD19467B_project/plots/chisquare_distribution.pdf')
# 	plt.clf()	
