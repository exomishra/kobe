# Step 0: Make Choices
# idea: store various choices that are made for the running of the program here 
# so that it can be conveniently be changed by the user afterwards.

print_details = True
# True for debugging 
# False for final running

print_fine_details = False
# Prints even finer details

#-----------------------------------------
# directories, reading, input management/processing

model = 'bern_model'
# specify your model; option: 'bern_model', 'completo'
# used for two purposes:
# reading columns in kobe_columns
# processing input via a suitable function
# please define this function in kobe.py and modify kobe_driver accordingly

dataset = 'ng76'
# string label for the dataset
# used for saving output

input_directory = '/home/lokeshmishra/PaleBlueDot/BernPhD/Projects/datasets/bern_model_data/data_from_horus/NGPPS/playing_with_kobe/'
# absolute path for input directories, ENDING WITH "/"
# examples 
# if on  horus cluster:'/shares/home0/lokeshmishra/simulated_data/NG76_1Msun_100emb_20Myr_0p3km_alpha2e-3_oldevap/'
# local PC : '/home/lokeshmishra/PaleBlueDot/BernPhD/Projects/data/data_from_horus/NGPPS/NG76_1Msun_100emb_20Myr_0p3km_alpha2e-3_oldevap/'

input_file = 'ref_red5e9.dat'
# input your filename
# options: 'all' or name a particular file example 'ref_red4e9.dat'
# 'all' mode is for internal use only
# 'all' mode allows us to read all ref_red files in a directory

time_age = '5e9'
# age of system, or None
# example: 4e9, 5e5, 1000
# if it throws an error, try using an age>2.3e7

auxiliary_file = 'ref_red5e9.dat'
# is none if all information is in input file
# otherwise provide file name (ensure auxiliary file is input_directory)
# this file is read only in your self-writted processing funtion

calculate_period = True
# if period needs to be calculated from sma
# otherwise, False

mstar_column = False
# if mass of star is a column --> True
# else --> False

mstar = 1
# provide mass of star in msun [float]
# if mass of star varies in population, provide column of mstar (string or int)
# note: column should have mstar in msun units. 

output_directory = '/home/lokeshmishra/PaleBlueDot/BernPhD/Projects/datasets/bern_model_data/data_from_horus/NGPPS/playing_with_kobe/kobe_evo_output/'
# absolute path for input directories, ENDING WITH "/"

cdpp_file = '/home/lokeshmishra/PaleBlueDot/BernPhD/Projects/kobe_development/kobe_driver/stellar_data/keplerstellar.csv'
# absolute path to kepler stellar cdpp data
# Download "Robust RMS CDPP" (preferably in csv format) from
# https://exoplanetarchive.ipac.caltech.edu/docs/Kepler_completeness_reliability.html
# Or you can contact me www.lokeshmishra.com

t_cdpp = 6
# CDPP is Combined Differntial Photometric Precision  - empirical estimate for stellar noise
# t_cdpp allows you to chose the transit duration to estimate the noise seen in the flux
# options are: 3,6,9,12
# Weiss et. al. 2018 calculations proceed with t_cdpp = 6 hours. 


robovetter_file = '/home/lokeshmishra/PaleBlueDot/BernPhD/Projects/kobe_development/kobe_driver/stellar_data/kepler_simulated_robovetterresults.txt'
# absolute path to Planet Detection Metrics: Vetting Metrics for Data Release 25
# Obtain from: https://exoplanetarchive.ipac.caltech.edu/docs/KeplerSimulated.html
# Look for 'Robovetter Results', download table corresponding to inj1
# Corresponding document -- J.L.  Coughlin  2017,Planet  Detection  Metrics:Robovetter Completeness and Effectiveness for Data Release 25, KSCI-19114-001

mass_threshold = 0
# planets with mass < mass_threshold will be removed
# reason: bern model starts embryos with 0.01 Me, so mass uncertainaty is high in low mass planets
# we call them failed embryos. they will not be detected by a transit survey like Kepler, anyway

sortsma = True
# sortsma while processing
# False, if your data has planets already sorted from inside out

#---------------------------------------------
# kobe physics

grazing_transits = False
# Reason : Weiss et. al. 2018 had a cut on impact parameter b at 0.9. 
# In the same spiirt, to select transit signals that are of high purity, we will ignore grazing transits.

simulate_nstars = 2e5
# Number of stars to simulate
# Governs the number of times each system is repeated. 
# repetitios = simulate_nstars/number of system in population

t_cdpp = 6
# CDPP is Combined Differntial Photometric Precision  - empirical estimate for stellar noise
# t_cdpp allows you to chose the transit duration to estimate the noise seen in the flux
# options are: 3,6,9,12,15
# Weiss et. al. 2018 calculations proceed with t_cdpp = 6 hours. 

t_obs = 3.5
# Enter number of years for which Kepler made observations
# value same as Weiss. et. al. 2018
# Now we convert this to days
t_obs = t_obs * 365.25

cdpp_seed = 42
# meaning of life
# reproducible draw of random noise

# Kepler Completeness - Choices:
completeness_period_range = [0,700,50] # start, stop, stepsize
completeness_radius_range = [0,20,2] # same

minimum_transit = 3
# planets with >= minimum transits are kobe_tce
# Following Kepler, the minimum number of transits required for detection is 3.
# Following Weiss et. al., the minimum number of transits required for detection is 2.

snr_threshold = 7.1
# planets with >= snr threshold are kobe tce
# Following Kepler, we set the snr_threshold at 7.1
# Following Weiss et. al., we set the snr threshold at 10. 

break_after_system = 10
# break_after = None
# for debugging, we allow the calculations to stop after finishing with certain number of systems
# break_after fixes that number.
# if None, then there is no break!