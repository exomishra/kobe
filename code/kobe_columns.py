# Some columns are essential for running KOBE:
# note: -1 is coming because python indexing starts at 0

# use the following as is for running with Bern Model

# essential columns:
## star columns
col_r_star = 51 - 1 
col_m_star = 'mstar'
# do not change mstar column here
# instead provide the column name in kobe_choice:mstar

## planet columns
col_r_planet = 'radius' # [rearth], processed in kobe.process_input

## system info to write
col_system = 'system'
col_planet = 'planet'
col_planetmultiplicity = 'multiplicity'

## system info to read
col_system_pop = 127
col_planet_pop = 128

##  orbital elements
col_sma = 19-1
col_period = 'period [days]'  #calculated in kobe.process_input
col_ecc = 63-1
col_inc = 64-1
col_long_node = 110-1
col_long_peri = 111-1

# following are needed for running with bern model
col_mass = 5-1 
col_status = 68-1 
col_core_radius = 10-1
col_total_radius = 15-1
col_transit_radius = 96-1 # radius where chord optical depth is 2/3 (light ray grazing terminator)
col_p1bar_radius = 92-1   # radius where pressure is 1 bar


