# Step 0: Make Choices
# idea: store various choices that are made for the running of the program here 
# so that it can be conveniently be changed by the user afterwards.

grazing_transits = False
# Reason : Weiss et. al. 2018 had a cut on impact parameter b at 0.9. 
# In the same spiirt, to select transit signals that are of high purity, we will ignore grazing transits.

print_details = False
# True for debugging 
# False for final running

sortsma = True
# sortsma while still at preliminary steps
# Should we sortsma at that step?

running_machine = 'local'
# Which machine are we running on?
# 'local' - for my own pc
# 'horus' - for horus cluster

which_ref_red = 'ref_red4e9.dat'
# options: 'all' or name a particular file example 'ref_red4e9.dat'

simulate_nstars = 2e5
# Number of stars to simulate
# Governs the number of times each system is repeated. 
# repetitios = simulate_nstars/number of system in population

n_population_series = 0
# which population to run? 
# 0 -- NG74 20  embryos
# 1 -- NG75 50  embryos
# 2 -- NG76 100 embryos

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

minimum_transit = 2
# planets with >= minimum transits are kobe_tce
# Following Kepler, the minimum number of transits required for detection is 3.
# Following Weiss et. al., the minimum number of transits required for detection is 2.

snr_threshold = 10
# planets with >= snr threshold are kobe tce
# Following Kepler, we set the snr_threshold at 7.1
# Following Weiss et. al., we set the snr threshold at 10. 

break_after_system = 100
# break_after = None
# for debugging, we allow the calculations to stop after finishing with certain number of systems
# break_after fixes that number.
# if None, then there is no break!
