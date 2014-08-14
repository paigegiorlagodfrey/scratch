import BDdb
db=BDdb.get_db('/Users/paigegiorla/Dropbox/BDNYCdb/BDNYC.db')
from matplotlib import pyplot as plt 
import numpy as np
import utilities as u
import astropy.units as q
import scipy.stats as s

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

#! 	pick wavelength range between points on object plot										#not necessary to this code
	wavelength_points=[1.1826776, 1.2093873000000002, 1.236097, 1.2628067, 1.2895164, 1.3162261000000002, 1.3429358, 1.3696455, 1.3963552, 1.4230649000000002, 1.4497746, 1.4764843, 1.503194, 1.5299037000000002, 1.5566134, 1.5833231, 1.6100328, 1.6367425, 1.6634522, 1.6901619, 1.7168716]		
	
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
	chisquare, sptlist, chilist, namelist, specidlist=[],[],[],[],[]

	for j in idlist:
		name,W,F,U=db.query.execute("SELECT source_id, wavelength, flux, unc from spectra where id='{}'".format(j)).fetchone()
		binned_wl=0
		jump=0.0267097																		#defined by the wavelength points of HD19467B data
		start=1.1826776-(jump/2) 															#start half way before 1st flux point
		W,F,U=u.scrub([W,F,U])
		W=W.value
		F=F.value
		U=U.value
		original=np.array([W,F,U])
		source,spectype=db.query.execute("SELECT spectral_types.source_id, spectral_types.spectral_type from spectral_types join spectra on spectral_types.source_id=spectra.source_id where spectral_types.regime='IR' and spectral_types.spectral_type>=20 and spectral_types.spectral_type<=30 and spectra.id='{}'".format(j)).fetchone()
		print name, source
		flux2=0
 		W_temp, F_temp, U_temp=[],[],[]
		unc2=0
		w1=start+(jump/2) 																	#first wl point is at start+the step/2
		wavelength=0
# 		W_temp, F_temp, U_temp=u.rebin_spec([W,F,U],wava)		
#!		loop through 21 iterations for the amount of wavelength ranges I have
		for i in range(21):
			binned_wl=jump + start 															#range is half below and half above w1
			
			flux2=sum(F[np.where(np.logical_and(W>start,W<binned_wl))])					#binning flux
			unc2=np.sqrt(sum(U[np.where(np.logical_and(W>start,W<binned_wl))]**2))			#bin uncertainty in quadrature
			
			start=binned_wl																	#set new range start
		
			F_temp.append(flux2)
			U_temp.append(unc2)
			wavelength=w1+(i*jump)															#this is where the flux point will go on plot
			W_temp.append(wavelength)	

		template=np.array([W_temp, F_temp, U_temp])											#create arrays for ease of use		
		template_plot1=u.norm_spec(template,original)												#Joe's code to normalize plots to each other
		template_plot1=[i.value for i in template_plot1]
#!		plot template over template's original flux
		plt.plot(template_plot1[0],template_plot1[1])
		plt.plot(original[0],original[1])
		plt.savefig('/Users/paigegiorla/Code/Python/BDNYC/TDwarfplotting/HD19467B_project/plots/binned_flux_with_original/'+'{}'.format(source)+ '.pdf')
		plt.clf() 	
		template=np.array([W_temp, F_temp, U_temp])	
		template_plot2=u.norm_spec(template,object)												#Joe's code to normalize plots to each other
		template_plot2=[i.value for i in template_plot2]
#!		plot template over object
		plt.errorbar(template_plot2[0],template_plot2[1],yerr=template_plot2[2],marker="o")
		plt.errorbar(object[0],object[1],yerr=object[2],marker="o")
		plt.savefig('/Users/paigegiorla/Code/Python/BDNYC/TDwarfplotting/HD19467B_project/plots/template_over_object/'+'{}'.format(source)+ '_HD19467.pdf')
		plt.clf()

#-----------------------------------------------------------------------------------------

		template=u.norm_spec(template,object)
		template=[i.value for i in template]
		chi=sum(((template[1]-object[1])**2)/((template[2]**2)+(object[2]**2)))     	#both data sets have uncertainties
		namelist.append(source)
		sptlist.append(spectype)
		chilist.append(chi)
		specidlist.append(j)
# 			output=np.array([name,spectype,chi])		
# 		except ValueError: pass
	
	output=np.array([namelist,specidlist,sptlist,chilist]) 	
	print len(sptlist)
		
#!	plotting the best fits, or minimum chi squared values over the object
	min_vals = np.where(output[3]<25)                                                         
	for i in min_vals:
		print "min_vals = {0}".format(i)
		print output[0][i],output[1][i],output[2][i], output[3][i]

		for k in output[1][i]:
			name,W,F,U=db.query.execute("SELECT source_id, wavelength, flux, unc from spectra where id='{}'".format(k)).fetchone()
			binned_wl=0
			jump=0.0267097																		#defined by the wavelength points of HD19467B data
			start=1.1826776-(jump/2) 															#start half way before 1st flux point

			source,spectype=db.query.execute("SELECT spectral_types.source_id, spectral_types.spectral_type from spectral_types join spectra on spectral_types.source_id=spectra.source_id where spectral_types.regime='IR' and spectral_types.spectral_type>=20 and spectral_types.spectral_type<=30 and spectra.id='{}'".format(k)).fetchone()
			print source, spectype
			flux2=0
			W_temp, F_temp, U_temp,=[],[],[]
			unc2=0
			w1=start-(jump/2) 																	#first wl point is at start+the step/2
			wavelength=0

	#!		loop through 21 iterations for the amount of wavelength ranges I have
			for i in range(21):
				binned_wl=jump + start 															#range is half below and half above w1
				flux2=sum(F[np.where(np.logical_and(W>start,W<binned_wl))])					#binning flux
				unc2=np.sqrt(sum(U[np.where(np.logical_and(W>start,W<binned_wl))]**2))			#bin uncertainty in quadrature
				start=binned_wl																	#set new range start

				F_temp.append(flux2)
				U_temp.append(unc2)
				wavelength=w1+(i*jump)															#this is where the flux point will go on plot
				W_temp.append(wavelength)	

			template=np.array([W_temp, F_temp, U_temp])
			template=u.norm_spec(template,object)				
			template=[i.value for i in template]
			plt.errorbar(template[0],template[1],yerr=template[2],marker="o",label=str(spectype))
			plt.legend( loc='upper right')
	plt.errorbar(object[0],object[1],yerr=object[2],marker="o",label='HD19467B')
	plt.savefig('/Users/paigegiorla/Code/Python/BDNYC/TDwarfplotting/HD19467B_project/plots/bestfits_w_object.pdf')
 	plt.clf()	
	
	chisquare=[sptlist,chilist]

#!	plot spectral type vs chi square
	plt.scatter(chisquare[0],chisquare[1])
	chisquare=zip(*sorted(zip(*chisquare)))
	pfit = np.polyfit(chisquare[0],chisquare[1], 3)   # Fit a 2nd order polynomial to (x, y) data
	yfit = np.polyval(pfit, chisquare[0])   # Evaluate the polynomial at x
	plt.plot(chisquare[0], yfit)
	plt.savefig('/Users/paigegiorla/Code/Python/BDNYC/TDwarfplotting/HD19467B_project/plots/chisquare_distribution.pdf')
	plt.clf()	
