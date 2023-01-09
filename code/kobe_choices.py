# Step 0: Make Choices
# idea: store various choices that are made for the running of the program here 
# so that it can be conveniently be changed by the user afterwards.

kobe_version = 2 
# Set to 1 for using KOBE version 1, detailed on https://github.com/exomishra/kobe
# Set to 2 for using KOBE version 2, detailed on https://github.com/pratishtha-rawat/kobe

# KOBE 1.0 can be particularly useful when you want to simulate the behavior of a transit survey and:
# Approximate circular orbits 

# KOBE 2.0 can be particularly useful when you want to simulate the behavior of a transit survey and:
# accomodate eccentric orbits
# include grazing transits 
# conduct studies based on transit durations
# include the behavior of the Kepler pipeline detection efficiency with SNR for marking TCEs

# Hereafter, the code can be used to specify choices as per the original first version of KOBE. After the choices for
# KOBE 1.0 are specified, you can make additional choices. These are part of the KOBE upgrade - KOBE 2.0.

################ KOBE 1.0 ################

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

input_directory = 'C:/Users/PR/Desktop/Masters_thesis/'
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
# None if all information is in input file
# otherwise provide absolute file path
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

output_directory = 'C:/Users/PR/Desktop/Masters_thesis/Output/'
# absolute path for input directories, ENDING WITH "/"

cdpp_file = 'C:/Users/PR/Desktop/Masters_thesis/keplerstellar.csv'
# absolute path to kepler stellar cdpp data
# Download "Robust RMS CDPP" (preferably in csv format) from
# https://exoplanetarchive.ipac.caltech.edu/docs/Kepler_completeness_reliability.html
# Or you can contact me www.lokeshmishra.com

t_cdpp = 6
# CDPP is Combined Differntial Photometric Precision  - empirical estimate for stellar noise
# t_cdpp allows you to chose the transit duration to estimate the noise seen in the flux
# options are: 3,6,9,12
# Weiss et. al. 2018 calculations proceed with t_cdpp = 6 hours. 


robovetter_file = 'C:/Users/PR/Desktop/Masters_thesis/kepler_simulated_robovetterresults.txt'
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

################ KOBE 2.0 ################

if kobe_version == 2:
    
    include_ecc = True 
    # Eccentricity of orbit is factored in.
    # the simplification considering circular orbits is removed.
    # False if you wish to run KOBE 1.0 only.
    
    average_transits = False
    # Default is True - calculates the number of transits using an average formula, as in Mishra et al., 2021 (n_transits = t_survey/P)
    # Set to False for calculating the number of transits using the phase of planet for improving the estimate of number of transits
    
    # Make internal choices for using specific expressions of the transit duration 
    t_dur_winn = False
    # set to true for using equations from Transits and Occultations by J. Winn in Exoplanets, 2010
    # Further, make the choice for using one of the specific set of expressions for the transit duration from 
    # Transits and Occultations by J. Winn in Exoplanets. See Figure 2 of the paper for more.
    t_contact_point = 0
    # Specify '0' for taking transit duration from contact point 4 to 1, labelled Total duration in paper
    # Specify '1' for taking transit duration from contact point 3 to 2, labelled Full duration in paper

    # set to true for using eq 15, Investigations of approximate expressions for the transit duration by David M. Kipping, 2010
    t_dur_kipping = True
    
    duty_cycle = True
    # The observing duty cycle is the fraction of t_obs with valid observations. 
    # Incorporate the effect of duty_cycle to calculate the valid number of transits by multiplying with f_duty
    f_duty = 0.88
    # f_duty parameter sets the duty cycle, for your specific telescope
    # For Kepler, set it to 0.88, as per Burke et al., 2015
    # Set to other values based on the observational survey KOBE is being used for
    
    transit_signal_MA2002 = True
    # Define the transit signal that will be calculated precisely based on the impact parameter and size ratio
    # Model transit signal as per eq. 1, Mandel and Agol 2002
    geo_condition_info = True
    # Each planet will satisfy a particular condition from the (1)-(4) in eq. 1, Mandel and Agol 2002
    # Set to true to get information on the condition satisfied, in the Output 'flag_transit_MA2002' column.
    # Flags are: NT (no transit can be observed), GT (grazing transit), FT (full transit can be observed) and 
    # Pdisk>*disk for planet disk larger than star's disk
    
    detection_efficiency = True
    # Marks potential transiting planets as TCEs based on detection efficiency, as defined 
    # and calculated in Christiansen, 2017
    # Flag is flag_det_eff with outputs '1' and '-1' for TCE and not TCE, respectively
    
elif kobe_version == 1:
    
    # Resets inputs given above to ensure only KOBE version 1.0 is run
    
    include_ecc = False
    average_transits = True
    t_dur_winn = False
    t_dur_kipping = False
    duty_cycle = False
    transit_signal_MA2002 = False
    geo_condition_info = False
    detection_efficiency = False
    
