"""
Created on Fri Dec 20 19:32:18 2019

@author: lokeshmishra, pratishtharawat

"""

#-------------------------------------
# Standard Imports
import os
import sys
import csv
import time,datetime
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

from mymodules import dataread, analysis, plottingroutine

import kobe
from kobe_choices import *
from kobe_columns import *

#-------------------------------------
kobe.print_header()
print('KOBE is working on %s.'%(model))
####################################
time_start = time.time()
####################################

#Step 0: Starters
# SI units
mass_sun = astropy_constants.M_sun
mass_earth = astropy_constants.M_earth
radius_sun = astropy_constants.R_sun
radius_earth = astropy_constants.R_earth

# start seed
np.random.seed(cdpp_seed)


####################################
# Step 1a: Setup for reading Data
####################################

# check if output directory exists, otherwise create it
if not os.path.exists(output_directory):
	os.makedirs(output_directory)
	if print_details == True:
		print('Output directory was not found. Output directory created.')

# Create conditional dealing of all vs particular files.
if input_file == "all":
	# Keep a list of all ref_red's and sort them by time
	files = os.listdir(input_directory)
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

elif input_file != "all":
	if print_details == True:
		print('Working with %s file only.'%(input_file))
		print('----------------------------------')
	files = [input_file]
	


####################################
# Step 1b: Read Stellar Data
####################################


if print_details == True:
	print('Started reading Kepler CDPP data.')
	print('----------------------------------')
cdpp = kobe.read_cdpp(cdpp_file=cdpp_file, t_cdpp=t_cdpp)
if print_details == True:
	print('Finished reading Kepler CDPP data.')
	print('----------------------------------')

####################################
# Step 1c: Kepler Completeness
####################################

if print_details ==True:
	print('Started reading Kepler completeness data.')
	print('----------------------------------') 

# To incorporate reliability set dispscore_cutoff      = 0.9 
# To not incorporate reliability set dispscore_cutoff  = 0

# load the ipac file which contains the kepler robovetter results, and turn into pandas dataframe
df_robovetter = kobe.read_keplercompleteness(robovetter_file=robovetter_file)

# A negative completeness for a bin means either there were no injections in that bin or the number of injections is less than the threshold. (Threshold is kept at 10, see Thompson 2018 or KSCI-19110-001/2) 
kepler_completeness = kobe.calculate_completeness(df_input=df_robovetter,x_col='period [days]',y_col='rp [R_earth]',
												x_range=completeness_period_range,y_range=completeness_radius_range,
												injection_threshold=10, spacing = 'linear', dispscore_cutoff=0.9)
if print_fine_details == True:
	print(kepler_completeness[0])
	print('----------------------------------')     
kepler_completeness = kepler_completeness[3]

if print_details ==True:
	print('Finished reading Kepler completeness data.')
	print('----------------------------------') 
	print('----------------------------------')



####################################
# Step 2: Loop over files
####################################


# Store count for progress bar
iterations_files = len(files)

for index_files in range(iterations_files):
	# Loop 1 
	# This loops over all the files. 
	# print progress bar on terminal
	# dataread.printProgressBar(iteration = index_files, total = iterations_files)

	print('KOBE starts with file %s'%(files[index_files]))
	print('----------------------------------')	

	# step 1 read file
	df_input = pd.read_csv(input_directory+files[index_files], delim_whitespace=True,header=None, low_memory=False)


	# step 2 process the input
	print('KOBE is processing input file(s).')
	print('----------------------------------')	
	if model == 'bern_model':
		df_input = kobe.process_input_bern(df_input= df_input,mass_threshold = mass_threshold, 
									calculate_period=calculate_period, mstar=mstar)
		if print_details ==True:
			print('Finished processing input data.')
			print('----------------------------------')
	elif model == 'completo':
		df_input = kobe.process_input_completo(df_input= df_input,mass_threshold = mass_threshold, 
									calculate_period=calculate_period, mstar=mstar)
		if print_details ==True:
			print('Finished processing input data.')
			print('----------------------------------')		

	else:
		print('----------------------------------')
		print('----------------------------------')
		print('----------------------------------')		
		print('Please process input data.')
		print('Ensure: planet numbers are sma-ordered.')
		print('Radius in rearth and inclination in degrees.')
		print('If your data needs no processing, please modify driver script to ensure smooth running.')
		print('----------------------------------')
		print('----------------------------------')
		print('----------------------------------')						
		break


	# STEP 6 calculate signal for each planet and add to new column (choice - grazing or complete transit)
	# do it here, so it can be taken for each observer in kobe_shadows
	df_input['transit_signal'] = ((df_input[col_r_planet]/df_input[col_r_star])*((radius_earth/radius_sun).value))**2

        
	# apply cuts : radius, period, mass
	# if your population has planets which need to be removed by cuts put them here 




	####################################
	####################################
	####################################
	# Step 3: KOBE - Loop over systems
	####################################
	####################################
	####################################
	if print_details ==True:
		print('KOBE Shadows: Started')
		print('----------------------------------') 

	df_kobe_output, sanity_checklist = kobe.kobe_shadows(df_input=df_input,cdpp=cdpp)

	if print_details ==True:
		print('KOBE Shadows: Finished')
		print('----------------------------------')
		print('----------------------------------')     

	if print_details ==True:
		print('KOBE Transits: Started')
		print('----------------------------------') 

	df_kobe_output = kobe.kobe_transits(df_kobe_output=df_kobe_output)

	if print_details ==True:
		print('KOBE Transits: Finished')
		print('----------------------------------')
		print('----------------------------------')

	if print_details ==True:
		print('KOBE Vetter: Started')
		print('----------------------------------') 

	df_kobe_output = kobe.kobe_vetter(df_kobe_output=df_kobe_output, kepler_completeness=kepler_completeness)

	if print_details ==True:
		print('KOBE Vetter: Finished')
		print('----------------------------------')
		print('----------------------------------')


	# SANITY CHECK - number of total systems should be same as nsystem_nonzero_views
	if sanity_checklist[0] != df_kobe_output['kobe_system_id'].unique().shape[0]:
		print('\n \n FATAL ERROR : Mismatch in nsystem_nonzero_views1 and systems in output! \n \n')


	if sanity_checklist[1] != df_kobe_output['kobe_system_id'].unique().shape[0]:
		print('\n \n FATAL ERROR : Mismatch in nsystem_nonzero_views2 and systems in output! \n \n')


	# SANITY Check - number of transiting planets should be same as nplanet_transiting
	if sanity_checklist[2] != df_kobe_output.shape[0]:
		print('\n \n FATAL ERROR : Mismatch in number of planets from shadow calculations and output \n \n')


	# step 10 print basic info in text output
	if break_after_system == None:
		break_after_system = df_input[col_system_pop].unique().shape[0]

	if dataset == None:
		filename = 'kobe_'+str(break_after_system)+'_systems_'+files[index_files]
	else:
		filename = 'kobe_'+dataset+'_'+str(break_after_system)+'_systems_'+files[index_files]

	now = datetime.datetime.now()
	# print ("Current date and time :",now.strftime("%Y-%m-%d %H:%M:%S"))


	f = open(output_directory+'Stats_'+filename+'.txt', "a")
	f.write('\n-------------------')
	f.write('\nCurrent date and time : ')
	f.write(now.strftime("%Y-%m-%d %H:%M:%S"))
	f.write('\nModel = %s; Reading File : %s'%(model,files[index_files]))
	f.write('\nDescription of original file: \nDataset = %s; (minimum mass threshold = %.2f) \nNumber of Systems = %d \nNumber of Planets =%d'%(dataset, mass_threshold,df_input[col_system_pop].unique().shape[0],df_input.shape[0]))
	f.write('\nGrazing Transits Included: %s, Stars simulated = %e, Observation time = %.1f year.'%(str(grazing_transits), simulate_nstars,t_obs/365.25))
	f.write('\nSystems analyzed = %d, Minimum Transits = %d, Minimum SNR (MES) = %.1f'%(break_after_system,minimum_transit,snr_threshold))
	f.write('\nKOBE Shadows: \nNumber of Systems = %d \nNumber of Planets = %d'%(df_kobe_output['kobe_system_id'].unique().shape[0],df_kobe_output.shape[0]))
	f.write('\nKOBE Transits: \nNumber of Systems = %d \nNumber of Planets = %d'%(df_kobe_output.loc[df_kobe_output['kobe_tce']==1]['kobe_system_id'].unique().shape[0],df_kobe_output.loc[df_kobe_output['kobe_tce']==1].shape[0]))
	f.write('\nKOBE Vetter: \nNumber of Systems = %d \nNumber of Planets = %d'%(df_kobe_output.loc[df_kobe_output['flag_completeness']=='PC']['kobe_system_id'].unique().shape[0],df_kobe_output.loc[df_kobe_output['flag_completeness']=='PC'].shape[0]))
	f.write('\n-------------------')
	f.close()


	# step 12b save df_kobe_output     
	df_kobe_output.to_csv(output_directory+filename+'.csv')
	print('KOBE is done with file %s'%(files[index_files]))
	print('----------------------------------')		

print('KOBE is finished. Thanks, bye-bye!')
print('Your output is stored at: %s.'%(output_directory))
time_run = time.time() - time_start
print("---Run time: %.5f seconds ---" % (time_run))
print("---Run time: %.5f minutes ---" %(time_run/60))