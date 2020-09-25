"""
Created on Fri Dec 20 19:32:18 2019

@author: lokeshmishra
"""

# Standard Imports
import os
import sys
import csv
import time
sys.path.insert(0, '../../')

#get_ipython().run_line_magic('matplotlib', 'inline')
#get_ipython().run_line_magic('matplotlib', 'inline')
import pylab as p
import pandas as pd
import numpy as np
import scipy as sp
import pprint
import matplotlib.pyplot as plt
import matplotlib as mpl
import mpl_toolkits.mplot3d.art3d as art3d
import seaborn as sns
sns.set()
from tqdm import tqdm_notebook as tqdm
from scipy import constants as scipy_constants
from astropy import constants as astropy_constants
from astropy import units as astropy_units
from mpl_toolkits.mplot3d import Axes3D

from mymodules import dataread, analysis, plottingroutine, kobe

from kobe_choices import *
####################################
time_start = time.time()
####################################
#Step 0: Starters
# SI units
mass_sun = 1.9884*10**30
mass_earth = 5.9722*10**24
radius_sun = 6.957*10**8
radius_earth = 6.3781366*10**6

mass_sun = astropy_constants.M_sun
mass_earth = astropy_constants.M_earth
radius_sun = astropy_constants.R_sun
radius_earth = astropy_constants.R_earth

# start seed
np.random.seed(cdpp_seed)

# label the dataset
if n_population_series == 0:
	dataset = 'NGPPS_NG74_020emb'
elif n_population_series == 1:
	dataset = 'NGPPS_NG75_050emb'
elif n_population_series == 2:
	dataset = 'NGPPS_NG76_100emb'
####################################
####################################
####################################
####################################
# Step 1a: Setup for reading Data
####################################
####################################
####################################

# Directory setup for my local computer
if running_machine == 'local':
	projectpath = '/home/lokeshmishra/PaleBlueDot/BernPhD/Projects/'
	rawdatapath = 'data/data_from_horus/NGPPS/'
	outputpath = 'kobe_development/sample_output/'

# Directory setup for horus
if running_machine == 'horus':
	projectpath = '/shares/home0/lokeshmishra/'
	rawdatapath = 'simulated_data/'
	outputpath  = 'kobe/kobe_output/ng74/'

population_series = ['NG74_1Msun_020emb_20Myr_0p3km_alpha2e-3_oldevap/',
					 'NG75_1Msun_050emb_20Myr_0p3km_alpha2e-3_oldevap/',
					 'NG76_1Msun_100emb_20Myr_0p3km_alpha2e-3_oldevap/']

# Create conditional dealing of all vs particular ref_red.
if which_ref_red == "all":
	# Keep a list of all ref_red's and sort them by time
	files = os.listdir(projectpath+rawdatapath+population_series[n_population_series])
	# Keep only 'ref_red' files
	files = [item for item in files if item[:7]=='ref_red']

	# Remove ref_red header file
	del files[files.index('ref_red.hdr')]
	del files[files.index('ref_redP2C.dat')]
	del files[files.index('ref_redtdisk.dat')]
	# let user know how many files being read in current population
	if print_details == True:
		print("Working with %d ref_red files"%(len(files)))

	# teaching order of time
	timestep = [item.lstrip('ref_red').rstrip('.dat') for item in files]
	# sort by time
	time_number = [float(item) for item in timestep]
	files = [item_files for item_time,item_files in sorted(zip(time_number,files))]
	timestep = [item.lstrip('ref_red').rstrip('.dat') for item in files]

	# now we add the removed ref_red files again, but at the end of the queue:
	files.append('ref_redP2C.dat')
	files.append('ref_redtdisk.dat')

elif which_ref_red != "all":
	if print_details == True:
		print('Working with %s file only.'%(which_ref_red))
	files = [which_ref_red]
	timestep = [item.lstrip('ref_red').rstrip('.dat') for item in files]
	

####################################
####################################
####################################
# Step 1b: Read Stellar Data
####################################
####################################
####################################
# part 1 - read noise
df_keplercomplete = pd.read_csv('stellar_data/keplerstellar.csv',sep=',',header=0)
# remove nan from cdpp
na_indices = df_keplercomplete['rrmscdpp06p0'].isna().to_numpy()
na_indices = np.where(na_indices == True)[0]
# na_indices = na_indices[0]
df_keplercomplete = df_keplercomplete.drop(axis=0, index = na_indices, inplace=False)
# we apply cuts to select noise for FGK - Solar Type Stars
condition_radius = (df_keplercomplete['radius']<= 5)
condition_mass = (df_keplercomplete['mass']>= 0.7) &  (df_keplercomplete['mass']<= 1.3) 
condition_temp = (df_keplercomplete['teff'] >= 3880) & (df_keplercomplete['teff']<= 7200)  
#fgk stars - temperatures from Peacaut and Mamajek 2013
#
df_keplercomplete = df_keplercomplete.loc[condition_radius].loc[condition_mass].loc[condition_temp]
#
cdpp3 = df_keplercomplete['rrmscdpp03p0'].to_numpy()
cdpp6 = df_keplercomplete['rrmscdpp06p0'].to_numpy()
cdpp9 = df_keplercomplete['rrmscdpp09p0'].to_numpy()
cdpp12 = df_keplercomplete['rrmscdpp12p0'].to_numpy()

if t_cdpp == 3:
	cdpp = cdpp3
elif t_cdpp == 6:
	cdpp = cdpp6
elif t_cdpp == 9:
	cdpp = cdpp9
elif t_cdpp == 12:
	cdpp = cdpp12
else:
	print('Error with t_cdpp choice.')

####################################
####################################
####################################
# Step 1c: Kepler Completeness
####################################
####################################
####################################
# Since incorporating reliability is still unclear, I will include completeness in two ways
# Path 1 - will include completeness for all planets with a disposition score of 0 and above (i.e. all planets are considered)
# Path 2 - will include completeness for all planets with a disposition score of 0.9 and above. (proxy - for incorporating reliability)

##1 we load the ipac file which contains the kepler robovetter results, and turn into pandas dataframe
df_robovetter = kobe.read_keplercompleteness()

##2 we calculate Kepler completeness
# Explanation of cuts used:
# Period - A min. transit of 2 (eg. Weiss et. al. 2018) for detection implies max 639 day period. Rounded off to 700.
# Radius - Forcing min. injections at 10, shows that there
# Remember - A negative completeness for a bin means either there were no injections in that bin or the number of injects is less than the threshold. (Threshold is kept at 10, see Thompson 2018 or KSCI-19110-001/2) 
kepler_completeness = kobe.calculate_completeness(df_input=df_robovetter,x_col='period [days]',y_col='rp [R_earth]',x_range=completeness_period_range,y_range=completeness_radius_range,injection_threshold=10, spacing = 'linear', dispscore_cutoff=0)
if print_details == True:
	print(kepler_completeness[0])
kepler_completeness = kepler_completeness[3]
# now with disposition score
kepler_completeness_reliable = kobe.calculate_completeness(df_input=df_robovetter,x_col='period [days]',y_col='rp [R_earth]',x_range=completeness_period_range,y_range=completeness_radius_range,injection_threshold=10, spacing = 'linear', dispscore_cutoff=0.9)
if print_details == True:
	print(kepler_completeness_reliable[0])
kepler_completeness_reliable = kepler_completeness_reliable[3]



####################################
####################################
####################################
# Step 2: Loop over files
####################################
####################################
####################################

# Store count for progress bar
iterations_files = len(files)

for index_files in range(iterations_files):
	# Loop 1 
	# This loops over all the ref_red_files. 
	#print progress bar on terminal
	# dataread.printProgressBar(iteration = index_files, total = iterations_files)

	# if print_details == True:
	#   print('Reading file %s'%(files[index_files]))
	print('Reading file %s'%(files[index_files]))

	# dataread processing
	# step 1 read ref red
	df_input = dataread.read_ref_red(projectpath+rawdatapath+population_series[n_population_series]+files[index_files],print_details=print_details)
	
	# step 2 extract selected columns
	df_input = dataread.extract_selected_columns(df_input)

	# step 3 preliminary steps (no masscut)
	df_input, dict_input = dataread.preliminary_steps(df_input,mass_threshold = 0, dataset = dataset)

	# STEP 6 calculate signal for each planet and add to new column (choice - grazing or complete transit)
	# do it here, so it can be taken over after shadow calculation directly
	df_input['transit_signal'] = ((df_input[dataread.find_key(dict_input, 'Radius [Rearth]')]/df_input[dataread.find_key(dict_input, 'Stellar Radius [Rsun]')])*((radius_earth/radius_sun).value))**2
	# step 4 apply cuts : radius, period, mass
	# DECISION - NO CUTS ARE APPLIED.
	# The idea is that the pipeline is desinged to find planets, the cuts should come out naturally.
	# Masscuts
	# The idea is to place a lower limit to get rid of failed embryos (0.1 Mearth). But these should not be detected due
	# to their small size.
	# Radius 
	# The idea is to place upper limit to get rid of unphysical planets (20 Rearth). But these will be detected and
	# should be later dealt with. 
	# Period
	# The idea is to place upper limit to get rid of planets that will not be detected by Kepler due to two conditions:
	# min_transit = 3, time_observation = 3.5 years ==> gives max period of ~420 days
	# But we avoid this because, Weiss et. al. used min_transit = 2. Again, will be dealt in pip analysis and not in 
	# kobe_output.
	
	####################################
	####################################
	####################################
	# Step 1d: Collect Column Numbers
	####################################
	####################################
	####################################
	#get column numbers
	col_sma = dataread.find_key(dict_input,'Semi-Major Axis [AU]')
	col_period = dataread.find_key(dict_input,'Period [days]')
	col_r_star = dataread.find_key(dict_input,'Stellar Radius [Rsun]')
	col_r_planet = dataread.find_key(dict_input, 'Radius [Rearth]')
	col_ecc = dataread.find_key(dict_input,'eccentricity')
	col_inc = dataread.find_key(dict_input,'Inclination [radians]')
	col_long_node = dataread.find_key(dict_input,'long_node')
	col_long_peri = dataread.find_key(dict_input,'long_peri')
	col_system = dataread.find_key(dict_input,'System')
	col_planet = dataread.find_key(dict_input, 'Planet')
	col_planetmultiplicity = dataread.find_key(dict_input, 'Planet Multiplicity')



	####################################
	####################################
	####################################
	# Step 3: KOBE - Loop over systems
	####################################
	####################################
	####################################

	# number of systems in current population
	nsystem = df_input[dataread.find_key(dict_input, 'System')].unique().shape[0]
	# number of views to keep for each system (nsystem*nviews = simulate_nstars)
	nviews = int(simulate_nstars/nsystem)
	# number of views with non-zero transiting planets (Kept for later sanity check on the length of final dataframe)
	# creating two of these for two different purposes
	nsystem_nonzero_views1 = 0
	nsystem_nonzero_views2 = 0
	# total number of transiting planets (for later sanity check)
	nplanet_transiting = 0
	# Start a KOBE output dataframe
	df_kobe_output = pd.DataFrame()
	

	for index_system in range(nsystem):
	# Loop 2
	# This loops over all the systems present in the current ref_red file.

		if print_details == True:
			print('Working on system number %d'%(index_system+1))

		# step 5 call kepler inclination or call kepler shadows
		# df_shadow = kobe.calculate_shadows(df_input=df_input, dict_column=dict_input,system_id=index_system,npoints=1e7,print_probability=print_details,grazing_transits=grazing_transits)
		# keep only nviews
		# np.random.seed(42)
		# df_shadow = df_shadow.loc[np.random.choice(df_shadow.index, size=nviews, replace=False)]

		# Here, we do a small shortcut to cut down computation time. 
		# Instead of calculating shadows at 1e7 point and then randomly selecting few views.
		# We calculate shadows for only required views.
		df_shadow, df_angle_data = kobe.calculate_shadows(df_input=df_input, dict_column=dict_input,system_id=index_system+1,npoints=nviews,print_probability=print_details,grazing_transits=grazing_transits)

		# reject rows where number of transiting planets is zero!
		df_shadow = df_shadow.loc[df_shadow['number_transiting_planet']!=0]
		df_angle_data = df_angle_data.loc[df_shadow.index]

		# check that the two frames have the same shape:
		if df_shadow.shape[0]!= df_angle_data.shape[0]:
			print('BIG PROBLEM BOSS: The number of views in angle dataframe is not the same with the shadow dataframe. WTF?')
		# re-number the indices of this dataframe (used in the loop later on)
		df_shadow.index = np.arange(df_shadow.shape[0])
		df_angle_data.index = np.arange(df_angle_data.shape[0])
		# add number of remaining views to number of views with non-zero transiting planets (for later sanity check)
		nsystem_nonzero_views1 = nsystem_nonzero_views1 + df_shadow.shape[0]
		# add number of transiting planets
		nplanet_transiting = nplanet_transiting + int(df_shadow['number_transiting_planet'].sum())
		
		####################################
		####################################
		####################################
		# Step 4: KOBE -The output file
		####################################
		####################################
		####################################
		# Loop over number of system-views to be added (each row of df_shadow)
		iterations_system_views = df_shadow.shape[0]

		for index_system_views in range(iterations_system_views):
			# Loop 3
			# This loops over all the views which have a transiting planet. 
			# our aim in this loop is to select the rows corresponding to the planets in shadow and 
			# add it to df_kobe_output. 
			nsystem_nonzero_views2 += 1
			# we start by selecting the system number of the current system (add 1 because of python)
			condition1 = (df_input[col_system] == index_system+1)
			# some numpy magic to get planet indices from the current row 
			planet_indices = df_shadow.iloc[[index_system_views],:].where(df_shadow==1).values[0]
			planet_indices = np.delete(planet_indices,-1) # get rid of the last column
			planet_indices = np.delete(planet_indices, [0,1]) # get rid of theta, phi
			planet_indices = list(np.where(planet_indices==1)[0])
			# finally append! 
			df_kobe_output = df_kobe_output.append(df_input.loc[condition1].iloc[planet_indices,:],ignore_index=True,sort=False)
			# reorder the indices of this dataframe
			# df_kobe_output.index = np.arange(df_kobe_output.shape[0])
			#

			# read the inclination seen by kobe for the appended planets
			# numpy magic to extract separate arrays out of (f,alpha) array
			f_alpha_pairs = np.concatenate(df_angle_data.iloc[[index_system_views],2:].to_numpy()[0][planet_indices]).reshape(-1,2)
			f_array = np.take(f_alpha_pairs, indices=0,axis=1)
			alpha_array = np.take(f_alpha_pairs, indices=1,axis=1)

			# note, observer's inclination is same as observer's location given by alpha (along polar)
			df_kobe_output.loc[df_kobe_output.tail(len(planet_indices)).index,'kobe_inclination'] = alpha_array 
			df_kobe_output.loc[df_kobe_output.tail(len(planet_indices)).index,'kobe_observer_azimuth'] = f_array

			# give kobe system numbers and planet numbers
			df_kobe_output.loc[df_kobe_output.tail(len(planet_indices)).index,'kobe_system_id'] = nsystem_nonzero_views2*np.ones(len(planet_indices))
			df_kobe_output.loc[df_kobe_output.tail(len(planet_indices)).index,'kobe_planet_id'] = np.arange(len(planet_indices))+1
			####################################
			####################################
			# KOBE TRANSIT Calculations
			####################################
			####################################
			# step 7 draw CDPP, calculate SNR (for each remaining planet)
			df_kobe_output.loc[df_kobe_output.tail(len(planet_indices)).index,'kobe_cdpp [ppm]'] = np.random.choice(cdpp)*10**-6
			# get out of this loop to continue SNR calculations

		if nsystem_nonzero_views1 != nsystem_nonzero_views2:
			print('Number of systems in output may be wrong. Check system views for system %d'%(index_system+1))

		if index_system == break_after_system:
			break


	# Calculate impact parameter for each planet with respect to the observer
	# we use observer's azimuth as true anomaly --> implies that planet is crossing observer's line of sight
	# we use observer's polar angle as inclination --> identical by definition
	distance_to_star = (df_kobe_output[col_sma]*scipy_constants.astronomical_unit*(1-(df_kobe_output[col_ecc]**2)))/(1 + (df_kobe_output[col_ecc]*np.cos(df_kobe_output['kobe_observer_azimuth'])))
	df_kobe_output['impact_parameter'] = (distance_to_star * np.cos(df_kobe_output['kobe_inclination']))/(df_kobe_output[col_r_star]*radius_sun.value)

	# geometrically, transits occur if impact parameter (sky projected distance of planet from star in stellar radius unit) is less than the sum of planet and stellar radii in stellar radius unit
	# we calculate the this geometric qunatity
	df_kobe_output['rp+rs/rs'] = 1 + ((df_kobe_output[col_r_planet]/df_kobe_output[col_r_star])*(radius_earth/radius_sun))

	## Any observer inside the transit shadow caused by a planet, is bound to see this planet's transit.
	# The geometric condition for transit is that the imapact parameter should be smaller than 'rp+rs/rs'. It should be true for all of the planets in kobe_output
	# We check this by counting how many fail this condition. 
	# first count planets that don't satisfy this strongly condition
	# I'm checking this on absolutue, to ensure that all planets are checked!
	condition_geotransit = abs(df_kobe_output['impact_parameter']) > df_kobe_output['rp+rs/rs']
	failed_planets = df_kobe_output.loc[condition_geotransit].shape[0]
	if failed_planets != 0:
		print('!!!!!!!!!!!!!')
		print('FATAL ERROR : GEOMETRIC TRANSIT CONDITION NOT SATISFIED STRONGLY FOR %d PLANETS !!!!'%(failed_planets))
		print('Strong Condition == b < (rs+rp)/rs')
		print('!!!!!!!!!!!!!')

	# consider a weeker condition as well
	condition_geotransit = abs(df_kobe_output['impact_parameter']) > df_kobe_output['rp+rs/rs']
	failed_planets = df_kobe_output.loc[condition_geotransit].shape[0]
	if failed_planets != 0:
		print('!!!!!!!!!!!!!')
		print('FATAL ERROR : GEOMETRIC TRANSIT CONDITION NOT SATISFIED WEAKLY FOR %d PLANETS !!!!'%(failed_planets))
		print('Weak Condition == b =< (rs+rp)/rs')
		print('!!!!!!!!!!!!!')

	# characteristic transit duration : eq 19, Transits and Occultations by J. Winn in Exoplanets
	df_kobe_output['transit_duration [hours]'] = ((df_kobe_output[col_r_star]*radius_sun.value)*(df_kobe_output[col_period]*24))/((np.pi)*(df_kobe_output[col_sma]*scipy_constants.astronomical_unit))
	# effective cdpp for transit duration. Following, eq 4, Christiansen et. al. 2012.
	df_kobe_output['kobe_cdpp_eff [ppm]'] = df_kobe_output['kobe_cdpp [ppm]']*((t_cdpp/df_kobe_output['transit_duration [hours]'])**0.5)
	# snr for 1 transit, or SES - single event statistic
	df_kobe_output['kobe_ses'] = df_kobe_output['transit_signal']/df_kobe_output['kobe_cdpp_eff [ppm]']

	# calculate number of transits that will be observed by kobe
	df_kobe_output['number_transits'] = t_obs/df_kobe_output[col_period]
	# Finally, snr(multi) or MES - multiple event statistic
	df_kobe_output['kobe_mes'] = df_kobe_output['kobe_ses']*((df_kobe_output['number_transits'])**0.5)


	####################################
	####################################
	# KOBE - tce
	####################################
	####################################
	# those planets which satisfy the following two conditions:
	# number transits >= minimum_transits (currently at 3)
	# kobe_mes >= snr_threshold
	# are flagged as kobe_tce (i.e. kobe_threshold crossing events): 1 means yes tce, and 0 means no tce.
	df_kobe_output['kobe_tce'] = 0
	condition_tce = (df_kobe_output['number_transits']>=minimum_transit) & (df_kobe_output['kobe_mes']>=snr_threshold)
	df_kobe_output.loc[df_kobe_output.loc[condition_tce].index, 'kobe_tce'] = 1

	####################################
	####################################
	# KOBE - Kepler Completeness & Reliability
	####################################
	####################################
	# Two flags are created. 
	# 1. 'flag_completeness' : 'FP' means planet is vetted as False Positive. 'PC' means planet is vetted as Planetary Candidate
	# 2. 'flag_completeness_reliability' : same binary meaning as above
	#
	# This allows us to incorporate and study the effects of incorporating reliability in our studies.
	
	## 0 - create x and y arrays to go through each bin     
	x_array = np.arange(completeness_period_range[0],completeness_period_range[1],completeness_period_range[2])
	y_array = np.arange(completeness_radius_range[0],completeness_radius_range[1],completeness_radius_range[2]) 
	
	## 1 - sanity checks
	drop_c=0
	drop_c_r = 0
	drop1_c=0
	drop1_c_r = 0
	parameter_space_exceeded_c, parameter_space_exceeded_c_r = [],[]

	
	## Create Columns 
	## first create columns and mark each as 'PC'
	df_kobe_output['flag_completeness'] = 'PC'
	df_kobe_output['flag_completeness_reliability'] = 'PC'
	## then all planets that are not tce are marked with NaN
	condition_nan = df_kobe_output['kobe_tce']==0
	df_kobe_output.loc[df_kobe_output.loc[condition_nan].index, 'flag_completeness'] = np.NaN
	df_kobe_output.loc[df_kobe_output.loc[condition_nan].index, 'flag_completeness_reliability'] = np.NaN
	

	for xind in range(x_array.shape[0]):
		for yind in range(y_array.shape[0]):
			## 2 now we are inside a bin
			## we will work two create two simultaneous flags

			completeness_c = kepler_completeness[xind][yind]
			completeness_c_r = kepler_completeness_reliable[xind][yind]
			
			if completeness_c == -0.01:
				completeness_c = 0
				parameter_space_exceeded_c.append(1)

			if completeness_c_r == -0.01:
				completeness_c_r = 0
				parameter_space_exceeded_c_r.append(1)
				
			drop_ratio_c = 1 - completeness_c
			drop_ratio_c_r = 1 - completeness_c_r 
			
			## all planets that fall into current bin and are tce are selected

			condition = (df_kobe_output[col_period]>=x_array[xind]) & (df_kobe_output[col_period]< x_array[xind]+completeness_period_range[2]) & (df_kobe_output[col_r_planet]>=y_array[yind]) & (df_kobe_output[col_r_planet]< y_array[yind]+completeness_radius_range[2]) & (df_kobe_output['kobe_tce']==1)
			planet_index = df_kobe_output.loc[condition].index.to_numpy()

			## estimate number of planets to flag
			drop_planets_size_c = int(np.floor(drop_ratio_c*planet_index.shape[0]))
			drop_planets_size_c_r = int(np.floor(drop_ratio_c_r*planet_index.shape[0]))
			drop_c = drop_c + drop_planets_size_c
			drop_c_r = drop_c_r + drop_planets_size_c_r
		
			## now we randomly select planets within this bin which are flagged as FP
			if planet_index.shape[0] != 0:
				drop_planets_index_c = list(np.random.choice(planet_index,size=drop_planets_size_c,replace=False))
				drop_planets_index_c_r = list(np.random.choice(planet_index,size=drop_planets_size_c_r,replace=False))
				
				df_kobe_output.loc[drop_planets_index_c,'flag_completeness'] = 'FP'
				df_kobe_output.loc[drop_planets_index_c_r, 'flag_completeness_reliability'] = 'FP'

				drop1_c = drop1_c + len(drop_planets_index_c)
				drop1_c_r = drop1_c_r + len(drop_planets_index_c_r)
	## completeness finished
	# step 9 kepler/kobe multiplicity 
	# initalize columns
	df_kobe_output['kobe_cr_multiplicity'] = 0
	df_kobe_output['kobe_c_multiplicity'] = 0
	# now calculate
	pps_system, pps_planet = dataread.planetspersystem(df_input=df_kobe_output.loc[df_kobe_output['flag_completeness_reliability']=='PC'], column_dictionary=dict_input,system_column='kobe_system_id',planet_column='kobe_planet_id')
	df_kobe_output.loc[df_kobe_output.loc[df_kobe_output['flag_completeness_reliability']=='PC'].index,'kobe_cr_multiplicity'] = pps_planet
	# once again but only for completeness
	pps_system, pps_planet = dataread.planetspersystem(df_input=df_kobe_output.loc[df_kobe_output['flag_completeness']=='PC'], column_dictionary=dict_input,system_column='kobe_system_id',planet_column='kobe_planet_id')
	df_kobe_output.loc[df_kobe_output.loc[df_kobe_output['flag_completeness']=='PC'].index,'kobe_c_multiplicity'] = pps_planet

	# SANITY CHECK - number of total systems should be same as nsystem_nonzero_views
	if nsystem_nonzero_views1 != df_kobe_output['kobe_system_id'].unique().shape[0]:
		print('\n \n FATAL ERROR : Mismatch in nsystem_nonzero_views1 and systems in output! \n \n')

	if nsystem_nonzero_views2 != df_kobe_output['kobe_system_id'].unique().shape[0]:
		print('\n \n FATAL ERROR : Mismatch in nsystem_nonzero_views2 and systems in output! \n \n')

	# SANITY Check - number of transiting planets should be same as nplanet_transiting
	if nplanet_transiting != df_kobe_output.shape[0]:
		print('\n \n FATAL ERROR : Mismatch in number of planets from shadow calculations and output \n \n')

	# step 10 print basic info in text output
	f = open(projectpath+outputpath+'kobe_basic_counting_info.txt', "a")
	f.write('\n-------------------')
	f.write('\nReading File : %s'%(files[index_files]))
	f.write('\nDescription of original file: \nDataset = %s \nNumber of Systems = %d \nNumber of Planets =%d'%(dataset,nsystem,df_input.shape[0]))
	f.write('\nKOBE TCE: \nNumber of Systems = %d \nNumber of Planets = %d'%(df_kobe_output.loc[df_kobe_output['kobe_tce']==1]['kobe_system_id'].unique().shape[0],df_kobe_output.loc[df_kobe_output['kobe_tce']==1].shape[0]))
	f.write('\nKOBE Completeness: \nNumber of Systems = %d \nNumber of Planets = %d'%(df_kobe_output.loc[df_kobe_output['flag_completeness']=='PC']['kobe_system_id'].unique().shape[0],df_kobe_output.loc[df_kobe_output['flag_completeness']=='PC'].shape[0]))
	f.write('\nKOBE Completeness+Reliability: \nNumber of Systems = %d \nNumber of Planets = %d'%(df_kobe_output.loc[df_kobe_output['flag_completeness_reliability']=='PC']['kobe_system_id'].unique().shape[0],df_kobe_output.loc[df_kobe_output['flag_completeness_reliability']=='PC'].shape[0]))
	f.write('\n-------------------')
	f.close()

	# step 11 extract additional columns in dictionary
	columns = df_kobe_output.columns
	for key in np.arange(columns.shape[0]):
		if type(columns[key])==str:
			dict_input[key] = columns[key]
	# Rename dataframe columns
	df_kobe_output.columns = np.arange(columns.shape[0])


	# step 12a save dictionary
	w = csv.writer(open(projectpath+outputpath+'dictionary_columns.csv', "w"))
	list_headers = []
	for key, val in dict_input.items():
		if key != 'name':
			if key != 'time':
				w.writerow([key, val])
				list_headers.append(val)

	# step 12b save df_kobe_output
	filename = dataset+'_'+dict_input['time']
	df_kobe_output.to_csv(projectpath+outputpath+filename+'.csv',header=list_headers)
	# break


time_run = time.time() - time_start
print("--- %.5f seconds ---" % (time_run))
print("--- %.5f minutes ---" %(time_run/60))
 
	
