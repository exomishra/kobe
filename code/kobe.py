"""
Created on Fri Nov 29 14:38:18 2019

@author: lokeshmishra
"""

# Standard Imports
import os
import sys

#get_ipython().run_line_magic('matplotlib', 'inline')
#get_ipython().run_line_magic('matplotlib', 'inline')
import pylab as p
import pandas as pd
import numpy as np
import scipy as sp
import pprint
import matplotlib.pyplot as plt
import matplotlib as mpl

import seaborn as sns
sns.set()

from scipy import constants as sp_constants
from astropy import constants as astropy_constants
from astropy import units as astropy_units


# from mymodules import dataread, analysis, plottingroutine
from kobe_choices import *
from kobe_columns import *

# SI units
mass_sun = 1.9884*10**30
mass_earth = 5.9722*10**24
radius_sun = 6.957*10**8
radius_earth = 6.3781366*10**6

mass_sun = astropy_constants.M_sun
mass_earth = astropy_constants.M_earth
radius_sun = astropy_constants.R_sun
radius_earth = astropy_constants.R_earth

def rotation_matrix(axis, angle, invert=False):
    """
    This function returns a 3D rotation matrix.
    
    Parameters
    ----------
    axis : string
        rotation axis. example: 'x', 'y', 'z'
    angle : float
        angle of rotation in radians
    invert : bool (default - False)
        if True - inverse of rotation matrix (by angle --> - angle)

    Examples:
    a) Rotate a vector along x, by 90 degrees, about z:
        i/o: rotation_matrix('z',np.pi/2, invert=False) @ np.array((1,0,0))
        o/p: array([0,1,0])

    b) Rotating arrays of vectors:
        input:
        v = np.array((np.array((0,0,1)),np.array((0,1,0)),np.array((1,0,0)), np.array((1,1,1))))
        v.reshape(3,-1).reshape(-1,3)       # reshape
        np.einsum('bc,ac',rotation_matrix('x',np.pi/2), v)  # use numpy's einsum notation format to carry out rotation
        output:
        array([ [ 0., -1.,  0.],
                [ 0.,  0.,  1.],
                [ 1.,  0.,  0.],
                [ 1., -1.,  1.]])

    """
    if invert == True:
        angle = - angle
    
    if axis == 'x':
        rotation_matrix = np.array([[1,0,0],[0, np.cos(angle), - np.sin(angle)],[0, np.sin(angle), np.cos(angle)]])
    
    if axis == 'y':
        rotation_matrix = np.array([[np.cos(angle), 0,  np.sin(angle)], [0,1,0], [-np.sin(angle),0, np.cos(angle)]])
        
    if axis == 'z':
        rotation_matrix = np.array([[np.cos(angle), - np.sin(angle) , 0], [ np.sin(angle), np.cos(angle), 0], [0,0,1]])
    
    # replace the small number with 0 
    small_number = 6.5e-17
    rotation_matrix[abs(rotation_matrix) < small_number] = 0.0
    
    return rotation_matrix

def calculate_shadows(df_input, system_id, npoints=1e7, 
 stellar_radius=None,seed=None, print_probability=True, grazing_transits=False):
    """
    This function calculates the transit shadow bands for all planets in a chosen system.
    
    Parameters
    ----------
    
    df_input : pandas dataframe
        dataframe containing planets with their orbital elements
    npoints : int/float
        number of grid points on sphere (default is 1e7)
    system_id : int
        system number for which shadows have to be calculated
    stellar_radius : float, optional
        Stellar radius in solar units.  
    seed : int
        seed the random number generator for grid points on sphere
    print_probability : bool, default True
        if true, prints individual transit probabilities for all planets
    grazing_transits : bool, default False
        if true, include grazing transits. if False only complete transits are calculated
    
    Returns
    -------
    df_shadow : pandas dataframe
        dataframe - column 1 is theta, the azimuthal angle range [0,2 pi]
        column 2 is phi, the polar angle range [0,pi]
        next column is for planet 1, and so on.
        Entry is 1 if (theta, phi) in shadow, else 0.

    df_angle_data : pandas dataframe
        (f,alpha) poisition of observer (wrt to each planet i.e. after rotations)
        useful in calculating impact parameter
                    
    Notes
    -----
    Timing for npoints per planet:
    1e7 points  - 1.2s
    1e8 points  - 7s
    1e9 points  - memory allocation error
    Results unreliable for less than 1e7 points.
    
    Individual transit probability of planet i is given by:
    df_shadow_output['planet i'].nonzero()[0].shape[0]/df_shadow_output.shape[0]
    
    """
    from astropy import units as astropy_units
    from astropy import constants as astropy_constants
    
    npoints = int(npoints)
        
    # Get system info: number of planets, stellar radius, etc..
    nplanets = df_input.loc[df_input[col_system]==system_id][col_planetmultiplicity].iat[0]
    if stellar_radius == None:
        r_star = df_input.loc[df_input[col_system]==system_id][col_r_star].mean()   * astropy_constants.R_sun.value
    else:
        r_star = stellar_radius                                                     * astropy_constants.R_sun.value
    series_sma = df_input.loc[df_input[col_system]==system_id][col_sma]             * astropy_constants.au.value
    series_r_planet = df_input.loc[df_input[col_system]==system_id][col_r_planet]   * astropy_constants.R_earth.value
    series_ecc = df_input.loc[df_input[col_system]==system_id][col_ecc]             #* astropy_units.one   
    series_inc = df_input.loc[df_input[col_system]==system_id][col_inc]             #* astropy_units.rad
    series_long_node = df_input.loc[df_input[col_system]==system_id][col_long_node] #* astropy_units.rad
    series_long_peri = df_input.loc[df_input[col_system]==system_id][col_long_peri] #* astropy_units.rad
    series_arg_peri = series_long_peri - series_long_node
    # series_transitprob = df_input.loc[df_input[col_system]==system_id][col_transitprob]
    
    # Step 1 - create grid on sphere
    r = 1                       # unit sphere
    if seed == None:
        seed = 42               # meaning of life
    np.random.seed(seed)
    u1 = np.random.rand(npoints) # a random number generator
    np.random.seed(seed+1)
    u2 = np.random.rand(npoints) # a random number generator
    theta = 2*np.pi*u1        #* astropy_units.rad
    phi = np.arccos(1-2*u2)   #* astropy_units.rad
    X = r*np.sin(phi)*np.cos(theta)
    Y = r*np.sin(phi)*np.sin(theta)
    Z = r*np.cos(phi)
    # Create an array of vectors, V
    V = np.column_stack((X,Y,Z))
    V_norm = np.linalg.norm(V, axis=1)
    ### SANITY TEST 1 ####
    #sanity check: norm of all vectors should be 1 (1 is True in bool)
    # check that the difference between norm of each vector and 1 is less than 1e-15. 
    # If not replace it with 0 (0 is flase).
    check_norm = np.where(abs(V_norm-1)<1e-15, V_norm, 0)
    if check_norm.all() == False:
        print('Sanity test 1 - FAILED. Not all vectors are unit normed. ')
    
    
    # Initiate output dataframe: first column is (theta) and second column is phi
    df_shadow_output = pd.DataFrame()
    df_shadow_output['theta'] = theta
    df_shadow_output['phi'] = phi

    df_angle_data = pd.DataFrame()
    df_angle_data['theta'] = theta
    df_angle_data['phi'] = phi

    # Step 2 - start loop over planets
    for index_planet in range(nplanets):
        # Step 2a - align orbit of current planet
        # inverse rotate along z by long_node
        # inverse rotate along x by inclination
        # finally inverse rotate along z by arg_peri
            
        v_prime = np.copy(V)
        v_prime = np.einsum('bc,ac',rotation_matrix('z',series_long_node.iat[index_planet],invert=True),v_prime)
        v_prime = np.einsum('bc,ac',rotation_matrix('x',series_inc.iat[index_planet],invert=True),v_prime)
        v_prime = np.einsum('bc,ac',rotation_matrix('z',series_arg_peri.iat[index_planet],invert=True),v_prime)
        v_prime_norm = np.linalg.norm(v_prime, axis=1)
        
        ### SANITY TEST 2 ####
        #sanity check: norm of all vectors should be 1 (1 is True in bool)
        # check that the difference between norm of each vector and 1 is less than 1e-15. 
        # If not replace it with 0 (0 is flase).
        check_norm = np.where(abs(v_prime_norm-1)<1e-15, v_prime_norm, 0)
        if check_norm.all() == False:
            print('Sanity Test 2 - FAILED. Not all vectors are unit normed. ')

        check_norm = np.where((abs(v_prime_norm)-1)<1e-2, v_prime_norm,0)
        if check_norm.all() == False:
            print('Sanity Test 3 - FAILED. Not all vectors are unit normed. ')


        
        ### Get cartesian of primed (rotated) vectors
        x_prime, y_prime, z_prime = np.hsplit(v_prime,3)
        x_prime = x_prime.reshape(-1,)
        y_prime = y_prime.reshape(-1,)
        z_prime = z_prime.reshape(-1,)
        r_xy_prime = np.sqrt(1 - z_prime*z_prime)
        #r_xy is norm of x & y
        
        ### Get inclination angle for each rotated observer
        # i.e. angle between Z and v_prime 
        # since both vectors have unit norm, and since one vector points only along z, the dot product simplifies
        # inclination_rad = np.arccos(z_prime)

        # no output because this is same as alpha of observer, defined below


        ### Get polar angles from rotated vectors (f- true anomaly is azimuthal, alpha is polar)
        alpha = np.arccos(z_prime)         #* astropy_units.rad        # in radians
        f = np.arccos(x_prime/r_xy_prime)  #* astropy_units.rad        # in radians
        
        # Step 2b - Calculate pairs of (f, alpha) inside shadow
        if grazing_transits == True:
            height_shadow = ((r_star - series_r_planet.iat[index_planet])*(1 + series_ecc.iat[index_planet] * np.cos(f)))/((series_sma.iat[index_planet])*(1 - series_ecc.iat[index_planet]**2))
        elif grazing_transits == False:
            height_shadow = ((r_star + series_r_planet.iat[index_planet])*(1 + series_ecc.iat[index_planet] * np.cos(f)))/((series_sma.iat[index_planet])*(1 - series_ecc.iat[index_planet]**2))
        
        alpha_rad = np.arcsin(height_shadow)
        # for each f (true anomaly), we calculate the shadow band by calculating the min and max of shadow
        # shadow_min = np.pi/2 * astropy_units.rad - alpha_rad
        # shadow_max = np.pi/2 * astropy_units.rad + alpha_rad
        shadow_min = np.pi/2 - alpha_rad
        shadow_max = np.pi/2 + alpha_rad        
        
        # For a simple probability calculation
        shadow = np.where((alpha >= shadow_min) & (alpha <= shadow_max), alpha,0)
        shadow_points = shadow.nonzero()[0].shape[0]
        #         print('Probability planet %d = %.3f %%'%(index_planet+1,100*shadow_points/npoints))

        # Step 2c - mark (theta,phi) pair of current planet inside shadow 
        # zeroes everywhere, except at shadow where it is 1
        # For storing location of shadow band
        shadow_index = np.where((alpha >= shadow_min) & (alpha <= shadow_max))[0]
        boolean_result = np.zeros((npoints))
        boolean_result[shadow_index] = 1 
        ## OUTPUT!
        df_shadow_output['planet %d'%(index_planet+1)] = boolean_result
        df_angle_data['(f,alpha) planet %d'%(index_planet+1)] = list(np.column_stack((f,alpha)))
        # sanity check - number of f,alpha pairs = theta phi pairs
        if shadow_points != boolean_result.nonzero()[0].shape[0]:
            print('DANGER: Number of shadow points is not same in output! HELP!')
            print(shadow_points)
            print(boolean_result.nonzero()[0].shape[0])
        
        # Step 2d - repeat for next planet
    # Step 3 - calculate number of transiting planets for each point on the grid
    df_shadow_output['number_transiting_planet'] = df_shadow_output[list(df_shadow_output.columns)[2:]].sum(axis=1)
    
    if print_probability == True:
        for index_planet in range(nplanets):
            probability = df_shadow_output.iloc[:,2+index_planet].to_numpy().nonzero()[0].shape[0]/df_shadow_output.shape[0]
            print('Transit probability of planet %d = %.4f %%'%(index_planet+1,100*probability))
        print('Multi Transit Probability for %d planets = %.4f %%'%(nplanets,100*df_shadow_output.loc[df_shadow_output['number_transiting_planet']==nplanets].shape[0]/df_shadow_output.shape[0]))
        if nplanets >= 6:
            print('Multi Transit Probability for %d planets = %.4f %%'%(nplanets-5,100*df_shadow_output.loc[df_shadow_output['number_transiting_planet']==nplanets-5].shape[0]/df_shadow_output.shape[0]))
        if nplanets >= 11:
            print('Multi Transit Probability for %d planets = %.4f %%'%(nplanets-10,100*df_shadow_output.loc[df_shadow_output['number_transiting_planet']==nplanets-10].shape[0]/df_shadow_output.shape[0]))
        if nplanets >= 16:
            print('Multi Transit Probability for %d planets = %.4f %%'%(nplanets-15,100*df_shadow_output.loc[df_shadow_output['number_transiting_planet']==nplanets-15].shape[0]/df_shadow_output.shape[0]))
        if nplanets >= 21:
            print('Multi Transit Probability for %d planets = %.4f %%'%(nplanets-20,100*df_shadow_output.loc[df_shadow_output['number_transiting_planet']==nplanets-20].shape[0]/df_shadow_output.shape[0]))

        # print % of observers that see atleast one transit.
        print('Shadow calculation: %d planets in system %d in %s population with %.1e grid points.\n%.2f %% observers see atleast one transiting planet.'%(nplanets, system_id,dataset, npoints,100*df_shadow_output.loc[df_shadow_output['number_transiting_planet']!=0].shape[0]/npoints))
        print('----------------------------------')    
    
    return df_shadow_output, df_angle_data

def read_keplercompleteness(robovetter_file):
    """
    This function reads the IPAC file containing Keplar Simulated Robovetter Results.
    
    Parameter
    ---------
    robovetter_file : string
        absolute path to robovetter file
        set/see detains in kobe_choices

    """
    #Step 1 - read kepler stellar file
    # Reading the file

    f = open(robovetter_file,'r')

    # Read and ignore header lines
    for index in range(75):
        line = f.readline()

    # Create container lists
    tceid, kicid, disp, dispscore = [],[],[],[]
    ntl, StellarEclipse,CentroidOffset,EphemerisMatch = [],[],[],[]
    period, epoch, expect_mes, detect_mes = [],[],[],[]
    ntran, depth_ppm, t_duration = [],[],[]
    rp, rs, teff, logg = [],[],[],[]
    sma, snr_dv, fit_prov, rp_over_rs, a_over_rs = [],[],[],[],[]
    impact = []

    # Loop over lines and extract variables of interest
    for line in f:
        data = line.strip().split() 
        tceid.append(data[0])
        kicid.append(data[1])
        disp .append(data[2])
        dispscore.append(data[3])
        ntl .append(data[4])
        StellarEclipse.append(data[5])
        CentroidOffset.append(data[6])
        EphemerisMatch.append(data[7])
        period.append(data[8])
        epoch.append(data[9])
        expect_mes.append(data[10])
        detect_mes.append(data[11])
        ntran.append(data[12])
        depth_ppm.append(data[13])
        t_duration.append(data[14])
        rp.append(data[15]) # earth radii
        rs.append(data[16]) # solar radii
        teff.append(data[17])
        logg.append(data[18])
        sma.append(data[19])
        rp_over_rs.append(data[20])
        a_over_rs.append(data[21])
        impact.append(data[22])
        snr_dv.append(data[23])
        # data[24] is insolation on planet
        fit_prov.append(data[25])

    # Step 2 - create dataframe out of it
    df_robovetter = pd.DataFrame()

    df_robovetter['tceid'] = tceid
    df_robovetter['kicid'] = kicid
    df_robovetter['disp'] = disp
    df_robovetter['dispscore'] = dispscore
    df_robovetter['Flag - Not Transit Like'] = ntl
    df_robovetter['Flag - StellarEclipse'] =  StellarEclipse
    df_robovetter['Flag - CentroidOffset'] = CentroidOffset
    df_robovetter['Flag - EphemerisMatch'] = EphemerisMatch
    df_robovetter['period [days]'] = period
    df_robovetter['epoch [barycentric kepler julian date]'] = epoch
    df_robovetter['expect_mes'] = expect_mes
    df_robovetter['detect_mes'] = detect_mes
    df_robovetter['ntran'] = ntran
    df_robovetter['depth_ppm'] = depth_ppm
    df_robovetter['t_duration [hours]'] = t_duration
    df_robovetter['rp [R_earth]'] = rp
    df_robovetter['rs [R_sun]'] = rs
    df_robovetter['t_star [Kelvin]'] = teff
    df_robovetter['logg [cm/s2]'] = logg
    df_robovetter['sma [AU]'] = sma
    df_robovetter['rp_over_rs'] = rp_over_rs
    df_robovetter['a_over_rs'] = a_over_rs
    df_robovetter['impact'] = impact
    df_robovetter['snr_dv'] = snr_dv
    df_robovetter['fit_prov'] = fit_prov


    # Step 3 - Turn robovetter output to boolean
    # Fix robovetting output to boolean:
    # 'PC' (planetary candidate) => 1
    # 'FP' (false positive) => 0
    bool_disp = np.zeros(df_robovetter.shape[0])
    for index in range(df_robovetter.shape[0]):
        if disp[index] == 'PC':
            bool_disp[index] = 1
    df_robovetter['vetting outcome'] = bool_disp

    # Change datatype to float

    cols = ['kicid', 'dispscore', 'Flag - Not Transit Like',
           'Flag - StellarEclipse', 'Flag - CentroidOffset',
           'Flag - EphemerisMatch', 'period [days]',
           'epoch [barycentric kepler julian date]', 'expect_mes', 'detect_mes',
           'ntran', 'depth_ppm', 't_duration [hours]', 'rp [R_earth]',
           'rs [R_sun]', 't_star [Kelvin]', 'logg [cm/s2]', 'sma [AU]',
           'rp_over_rs', 'a_over_rs', 'impact', 'snr_dv', 'fit_prov',
           'vetting outcome']

    for index in cols:
        df_robovetter[index] = df_robovetter[index].astype(float)

    return df_robovetter

def calculate_completeness(df_input, x_col, y_col, x_range, y_range, injection_threshold,spacing,dispscore_cutoff):
    """
    This function calculates completeness via counting number of injected TCE's vetted as PC.

    Parameters
    ----------

    df_input : pandas dataframe
        dataframe with Kepler Simulated results
        output of function 'read_keplercompleteness()'

    x_col : string
        provide one axis on which completeness is projected
        usually : 'period [days]'

    y_col : string
        usually : 'rp [R_earth]'

    x_range : list
        format : [min, max, step size]

    y_range : list
        format : [min, max, step size]

    injection_threshold : int
        the minimum number of injected TCE's per bin (inclusive) for which completeness is calculated.
        if number of injected TCS in a bin falls below this threshold, 
        then completeness for this bin is fixed to -0.01.

    spacing : string
        options are: 'linear' or 'log'

    dispscore_cutoff : float (between 0 and 1)
        used to incorporate Kepler 'reliability' (see sec: 7.3.4 in Thompson et. al. 2018, or Mulders et. al. 2018 or Bryson et. al. 2019)
        use : 0 to not invoke do anything
        use : 0.9 to follow Mulder et. al.


    Note
    ----
    At each bin lower point in included, and upper point is excluded.
    """
    # spacing
    if spacing == 'linear':
        x_array = np.arange(x_range[0],x_range[1],x_range[2])
        y_array = np.arange(y_range[0],y_range[1],y_range[2])
        
    elif spacing =='log':
        xrange = [x_range[0],x_range[1]]
        yrange = [y_range[0],y_range[1]]
        xbins = x_range[2]
        ybins = y_range[2]
        xdeltalog = np.log10(xrange[1]/xrange[0])/(xbins-1)
        ydeltalog = np.log10(yrange[1]/yrange[0])/(ybins-1)

        x_array, y_array = np.zeros(xbins),np.zeros(ybins)
        for index in range(xbins):
            x_array[index] = 10**(np.log10(xrange[0]) + index*xdeltalog)
        for index in range(ybins):
            y_array[index] = 10**(np.log10(yrange[0]) + index*ydeltalog)
    
    total_injections = np.zeros((x_array.shape[0],y_array.shape[0]))
    pc_injections = np.zeros((x_array.shape[0],y_array.shape[0]))
    completeness = np.zeros((x_array.shape[0],y_array.shape[0]))
    thresholdbins = 0
    emptybins = 0
    totalbins = 0

    
    for xind in range(x_array.shape[0]):
        for yind in range(y_array.shape[0]):
            if spacing == 'linear':
                condition = (df_input[x_col]>=x_array[xind]) & (df_input[x_col]< x_array[xind]+x_range[2]) & (df_input[y_col]>=y_array[yind]) & (df_input[y_col]< y_array[yind]+y_range[2])
            elif spacing == 'log':    
                condition = (df_input[x_col]>=10**np.log10(x_array[xind])) & (df_input[x_col]< 10**(np.log10(x_array[xind])+xdeltalog)) & (df_input[y_col]>=10**(np.log10(y_array[yind]))) & (df_input[y_col]< 10**(np.log10(y_array[yind])+ydeltalog))
            
            total_injections[xind][yind] = int(df_input.loc[condition].shape[0])           
            pc_injections[xind][yind] = int(df_input.loc[condition & (df_input['vetting outcome']==1) & (df_input['dispscore'] >= dispscore_cutoff)].shape[0])
            ## threshold check
            if total_injections[xind][yind] == 0:
                completeness[xind][yind] = -0.01
                emptybins += 1
            elif total_injections[xind][yind] < injection_threshold:
                completeness[xind][yind] = -0.01
                thresholdbins += 1
            else:
                completeness[xind][yind] = pc_injections[xind][yind]/total_injections[xind][yind]
            totalbins += 1
    if print_details == True:
        print('Bins along x and y: ', total_injections.shape) 
        print('Number of bins with less than %d injections = %d'%(injection_threshold,thresholdbins))
        print('Number of empty bins = %d'%(emptybins))
        print('Number of total bins = %d'%(totalbins))
    metadata = 'Completeness Calculation Information: \n Column along x is %s with min %.2f, max %.2f and step size %d. \n Column along y is %s with min %.2f, max %.2f and step size %.2f. \n Order of output: this text, total injections, pc injections, completeness(over threshold), \n x_range(mix, max, step), y_range(min, max, step). \n Note: Completeness ranges from 0-1. \n Completeness of -1 flags when total injected TCEs are either 0 or less than threshold injections of %d. \n Spacing is %s. \n Disposition Score cutoff is at %.1f'%(x_col,x_range[0],x_range[1],x_range[2],y_col,y_range[0],y_range[1],y_range[2],injection_threshold, spacing, dispscore_cutoff)
    
    output =  np.array([metadata,total_injections, pc_injections, completeness, x_range, y_range,dispscore_cutoff], dtype=object)
    
    return output

def read_cdpp(cdpp_file, t_cdpp):
    """
    This function read the KEPLER noise data and returns cdpp.
    
    Parameters
    ----------

    cdpp_file : string
        absolute path to stellar noise data from Kepler
        set path in kobe_choices

    t_cdpp : int
        approximate duration of transit
        options are: 3,6,9,12

    Output
    ------

    cdpp : numpy array
        a distribution of cdpp 

    """
    # part 1 - read noise
    df_keplercomplete = pd.read_csv(cdpp_file,sep=',',header=0,low_memory=False)
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

    return cdpp

def process_input_bern(df_input, mass_threshold, calculate_period, mstar):
    """
    Process the input a bit. This function has been developed for Bern Model.
    You may need to do your input processing. 

    Parameters
    ----------
    
    df_input : pandas dataframe
        an input of data

    mass_threshold : float/int
        minimum mass of a planet which is kept (planets with mass < mass_threshold are dropped)

    calculate_period : boolean
        whether to calculate period from sma, mass of star

    mstar : float or string
        constant mass of star in msun
    
    Output
    ------

    df_output : pandas dataframe
        an output of data

    Notes (for processing bern model output)
    ----
    - place mass cut
    - remove ejected/accreeted/collided planets
    - deal with radius definition
    - fix negative sma
    - calculate periods
    - convert Rj into Rearth
    - convert inclination to degrees
    - sort planets by sma
    - correct system numbers, planet numbers and multiplicity    
    - correct index of final output 

    """

    # masscut
    df_input = df_input.loc[df_input[col_mass]>= mass_threshold]

    # remove planets that did not end well
    # in bern model this is checked via column 'emps_status'
    # status codes: 0 fine, negative accreted by planet, 2 ejected, 3 accreted by star
    # 5 couldn't start, 6 not started, 8 core evaporated, 9 tidally spiralled
    df_input = df_input.loc[df_input[col_status]==0]


    # convert Rj to Rearth
    # astropy_constants.R_jup/astropy_constants.R_earth = 11.208981
    df_input[col_core_radius] = df_input[col_core_radius]*11.208981
    df_input[col_total_radius] = df_input[col_total_radius]*11.208981
    df_input[col_transit_radius] = df_input[col_transit_radius]*11.208981
    df_input[col_p1bar_radius] = df_input[col_p1bar_radius]*11.208981 

    # fixing radius
    # if time_age is < 2.3e7 then use max of core or p1bar radius (planete)
    # if time_age is > 2.3e7 then use max of core or transit radius (completo)
    # all final radius values are stored in col_total_radius:
    if time_age != None:
        if float(time_age) < 2.3e7:
            df_input['radius'] = df_input[[col_core_radius,col_p1bar_radius]].max(axis=1)
            
        elif float(time_age) > 2.3e7:
            df_input['radius'] = df_input[[col_core_radius,col_transit_radius]].max(axis=1) 
  

    # fix planets with negative sma (planets with sma less than 1e3 are put at 1e2)
    # note: 4e3 au is solar radii
    if df_input.loc[df_input[col_sma]<0].shape[0]>0:
        print('Error: %d Planets found with SMA < 0 au. SMA replaced to 0.01 au.'%(df_input.loc[df_in[col_sma]<0][col_planet].unique().shape[0]))
    df_input = df_input.where(df_input[col_sma]>=0,other=10**-2)

    
    # if mstar is not a column then set constant value:
    # else copy column
    if mstar_column == False:
        print('mstar is not a column; reading value mstar = %.2f msun.'%(mstar))        
        mstar = mstar * mass_sun.value
        df_input[col_m_star] = mstar
    elif mstar_column == True:
        df_input[col_m_star] = df_input[mstar]
        print('mstar column found.')

    # calculate period (assuming circular orbit)
    if calculate_period == True:
        kepler_orbit_constant = 4*(np.pi**2)*(sp_constants.gravitational_constant**-1)*(df_input[col_m_star]**-1)        
        df_input[col_period] = ((((df_input[col_sma]*sp_constants.astronomical_unit)**3)*kepler_orbit_constant)**0.5)*(sp_constants.day**-1)    
        


    # convert inclination from radians to degrees
    # df_input[col_inc] = df_input[col_inc]*(180/np.pi)    

    # sort by sma
    df_input = sort_by_column(df_input, which_column=col_sma, col_system=col_system_pop,col_planet=col_planet_pop)

    # correct system, planet numbers and add multiplicity
    df_input = correct_system_planet_numbers(df_input, col_system=col_system_pop, col_planet=col_planet_pop,
                                            new_system_column=col_system, new_planet_column=col_planet,
                                            col_planetmultiplicity=col_planetmultiplicity, correct_index=True)

    return df_input

def process_input_completo(df_input, mass_threshold, calculate_period, mstar):
    """
    Process the input a bit. This function has been developed for Bern Model.
    You may need to do your input processing. 

    Parameters
    ----------
    
    df_input : pandas dataframe
        an input of data

    mass_threshold : float/int
        minimum mass of a planet which is kept (planets with mass < mass_threshold are dropped)

    calculate_period : boolean
        whether to calculate period from sma, mass of star

    mstar : float or string
        constant mass of star in msun
    
    Output
    ------

    df_output : pandas dataframe
        an output of data

    Notes (for processing bern model output)
    ----
    - floor emps status to make integer
    - remove ejected/accreeted/collided planets
    - deal with radius definition
    - fix negative sma
    - calculate periods
    - read orbital elements from ref_red
    - convert inclination to degrees
    - sort planets by sma
    - correct system numbers, planet numbers and multiplicity    
    - correct index of final output 

    """

    # masscut
    df_input = df_input.loc[df_input[col_mass]>= mass_threshold]

    # remove planets that did not end well
    # in bern model this is checked via column 'emps_status'
    # status codes: 0 fine, negative accreted by planet, 2 ejected, 3 accreted by star
    # 5 couldn't start, 6 not started, 8 core evaporated, 9 tidally spiralled
    df_input[col_status] = np.floor(df_input[col_status])
    df_input = df_input.loc[df_input[col_status]==0]

    df_input.index = np.arange(df_input.shape[0])    

    # fixing radius
    # if time_age is < 2.3e7 then use max of core or p1bar radius (planete)
    # if time_age is > 2.3e7 then use max of core or transit radius (completo)
    # all final radius values are stored in col_total_radius:
    if time_age != None:
        if float(time_age) < 2.3e7:
            df_input['radius'] = df_input[[col_core_radius,col_p1bar_radius]].max(axis=1)
            
        elif float(time_age) > 2.3e7:
            df_input['radius'] = df_input[[col_core_radius,col_transit_radius]].max(axis=1) 
  

    # fix planets with negative sma (planets with sma less than 1e3 are put at 1e2)
    # note: 4e3 au is solar radii
    if df_input.loc[df_input[col_sma]<0].shape[0]>0:
        print('Error: %d Planets found with SMA < 0 au. SMA replaced to 0.01 au.'%(df_input.loc[df_in[col_sma]<0][col_planet].unique().shape[0]))
    df_input = df_input.where(df_input[col_sma]>=0,other=10**-2)

    
    # if mstar is not a column then set constant value:
    # else copy column
    if mstar_column == False:
        print('mstar is not a column; reading value mstar = %.2f msun.'%(mstar))        
        mstar = mstar * mass_sun.value
        df_input[col_m_star] = mstar
    elif mstar_column == True:
        df_input[col_m_star] = df_input[mstar] * mass_sun.value
        print('mstar column found.')

    # calculate period (assuming circular orbit)
    if calculate_period == True:
        kepler_orbit_constant = 4*(np.pi**2)*(sp_constants.gravitational_constant**-1)*(df_input[col_m_star]**-1)        
        df_input[col_period] = ((((df_input[col_sma]*sp_constants.astronomical_unit)**3)*kepler_orbit_constant)**0.5)*(sp_constants.day**-1)    
        
    # convert Rj to Rearth
    # astropy_constants.R_jup/astropy_constants.R_earth = 11.208981
    # df_input[col_core_radius] = df_input[col_core_radius]*11.208981
    # df_input[col_total_radius] = df_input[col_total_radius]*11.208981
    # df_input[col_transit_radius] = df_input[col_transit_radius]*11.208981
    # df_input[col_p1bar_radius] = df_input[col_p1bar_radius]*11.208981  

    # read orbital elements from auxiliary file:
    df_refred = pd.read_csv(input_directory+auxiliary_file, delim_whitespace=True, header=None, low_memory=False)
    col_refred_system = 128-1
    col_refred_planet = 129-1

    # read from corresponding ref_red 
    col_refred_ecc = 63-1
    col_refred_inc = 64-1
    col_refred_long_node = 110-1
    col_refred_long_peri = 111-1

    # loop over input dataframe
    if print_details == True:
        print('Copying data from auxiliary file.')
    for index_row in range(df_input.shape[0]):
        planet_id = df_input.iloc[index_row,col_planet_pop]
        system_id = df_input.iloc[index_row,col_system_pop]
        condition = (df_refred[col_refred_system]==system_id) & (df_refred[col_refred_planet]==planet_id)
        df_input.loc[index_row, col_ecc] = df_refred.loc[condition][col_refred_ecc].to_numpy()[0]
        df_input.loc[index_row, col_inc] = df_refred.loc[condition][col_refred_inc].to_numpy()[0]
        df_input.loc[index_row, col_long_node] = df_refred.loc[condition][col_refred_long_node].to_numpy()[0]
        df_input.loc[index_row, col_long_peri] = df_refred.loc[condition][col_refred_long_peri].to_numpy()[0] 
        # print(index_row, planet_id, system_id,df_refred.loc[condition][col_refred_ecc].to_numpy()[0])


    if print_details == True:
        print('Finished copying data.')

    # convert inclination from radians to degrees
    # df_input[col_inc] = df_input[col_inc]*(180/np.pi)    

    # sort by sma
    df_input = sort_by_column(df_input, which_column=col_sma, col_system=col_system_pop,col_planet=col_planet_pop)

    # correct system, planet numbers and add multiplicity
    df_input = correct_system_planet_numbers(df_input, col_system=col_system_pop, col_planet=col_planet_pop,
                                            new_system_column=col_system, new_planet_column=col_planet,
                                            col_planetmultiplicity=col_planetmultiplicity, correct_index=True)

    return df_input

# Count number of planets per system
def planetspersystem(df_input, column_dictionary=None,system_column=None,planet_column=None):
    """
    Given a data frame containing planetary system with system number and planet number, 
    this function gives two lists. 
    List 1 - system_wise -- entries correspond to planetary systems
    List 2 - planet_wise -- entries correspond to each planet in the the planetary system
    """
    if system_column==None:
        system_column = dataread.find_key(column_dictionary,'System')
    
    if planet_column==None:
        planet_column = dataread.find_key(column_dictionary,'Planet')
    
    pps_system=[]
    ind=0
    while ind < df_input.shape[0]:
        nsys = df_input.iloc[ind][system_column]
        temp=[]
        while df_input.iloc[ind][system_column]==nsys:
            temp.append(int((df_input.iloc[ind][planet_column])))
            ind+=1
            pmax=len(temp)
            if ind==df_input.shape[0]:
                break
        pps_system.append(len(temp))
    
    # Now we write the number of planets per system for each planet
    pps_planet = []
    for ind1 in np.arange(len(pps_system)):
        nsys = pps_system[ind1]
        ind2 = 0
        while ind2 < nsys:
            pps_planet.append(nsys)
            ind2 += 1
    # nsystems = len(pps_system)
    # nplanets1 = sum(pps_system)
    # nplanets2 = len(pps_planet)
    # check = nplanets1 - nplanets2
    #     print('Total number of systems = %d'%(nsystems))
    #     print('Total number of planets = %d and check is %d (0 means ok)'%(nplanets1,check))
    return pps_system,pps_planet

# Function to sort dataframe by a column
def sort_by_column(df_input,which_column, col_system, col_planet, ascending_order=True):
    """  
    This function is used to sort planets within a system according to the sorting column.

    Parameters
    ----------

    df_inputv: dataframe

    which_column : column info
        which column to sort on

    col_system : column_info
        system column to read

    col_planet : column info
        planet column to read
  
    acending_order : (boolean)
        True means low to high
        False means high to low

    Output
    ------
    df_sorted : dataframe sorted by column

    """
    

    # Correction for index list
    df_input.index = np.arange(df_input.shape[0])
    
    # Total number of systems
    iterations = int(df_input.iloc[-1][col_system_pop])
    system_array = df_input[col_system_pop]
    
    # Create new empty dataframe
    df_sorted = pd.DataFrame()
    # failed_systems = []    
    # Loop over each system
    for current_system in range(iterations):
        current_system+=1 #correction because for loop indices start from 0!
        
        current_indices = np.where(system_array==current_system)
        entries = df_input.iloc[current_indices[0]][which_column]
        
        # Sanity Check 1
        # multiplicity_check = df_input.iloc[current_indices[0][0]][col_planetmultiplicity]
        # if multiplicity_check != len(entries):
        #     failed_systems.append(current_system)
        
        #Sort entries:
        entries = entries.sort_values(ascending=ascending_order)
        sorted_indices = entries.index
        
        #Place sorted entries in the new dataframe
        df_sorted = df_sorted.append(df_input.iloc[sorted_indices][:], ignore_index=True)
    # exit loop
    #Correct planet numbers in sorted dataframe
    # df_sorted = correct_system_planet_numbers(df_sorted, column_dictionary)
        
    # if len(failed_systems)>0:
    #     print('Multiplicity mismatch found in %d systems!'%(len(failed_systems)))
        
    return df_sorted

# Function to correct the system and planet number columns.
def correct_system_planet_numbers(df_input, column_dictionary=None,col_system=None, col_planet=None,new_system_column=None, new_planet_column = None, col_planetmultiplicity=None,correct_index=True):
    """
    This function corrects system and planet numbers to ensure continuity.
    Can also correct multiplicity if user provides col_planetmultiplicity. 
    """ 
    if col_system==None:
        col_system = dataread.find_key(column_dictionary,'System')
    
    if col_planet==None:
        col_planet = dataread.find_key(column_dictionary,'Planet')
    
    pps_system, pps_planet = planetspersystem(df_input,column_dictionary,col_system, col_planet)
    # Correct system numbers to ensure continuity
    # creating system list
    syslist = []
    nsys = 1
    for ind in range(len(pps_system)):
        iterations = pps_system[ind]
        ind2 = 0
        while ind2 < iterations:
            syslist.append(nsys)
            ind2 += 1
        nsys += 1
    if new_system_column == None:
        df_input[col_system] = syslist
    else:
        df_input[new_system_column] = syslist
    
    # Correct planet numbers to ensure continuity
    # We use pps_system(planets per system count - system wise) to create a list of planet index
    planetlist = []
    for ind in range(len(pps_system)):
        iterations = pps_system[ind]
        ind2 = 1
        plist = []
        while ind2 <= iterations:
            plist.append(ind2)
            ind2 += 1
        planetlist.append(plist)
    # Flatten this list
    planetlist = [item for sublist in planetlist for item in sublist]
    if new_planet_column == None:
        df_input[col_planet] = planetlist
    else:
        df_input[new_planet_column] = planetlist


    # Correct multiplicity column is provided:
    if col_planetmultiplicity != None:
        df_input[col_planetmultiplicity] = pps_planet

    #     print('Total Number of Systems and Planets is %d and %d'%(len(pps_system),len(pps_planet)))

    if correct_index == True:
        df_input.index = np.arange(df_input.shape[0])

    return df_input

def kobe_shadows(df_input,cdpp):
    """
    The KOBE Shadows module: calculate transit shadow bands for all planets

    NOTE
    ----
    1. the transit signal (rplanet/rstar)**2 is calculated in driver for all planets
        because this remains same for all planets irrespective of observers

    2. CDPP for each KOBE shadow system is sampled at end of this function
    
    Parameters
    ----------
    df_input : pandas dataframe
        input dataframe processed suitably

    Output 
    ------
    df_kobe_output : pandas dataframe
        dataframe of planets that transit (potentially detectable)
        columns: all from input, some are added

    sanity_checklist : list
        list containing values for sanity checks in driver

    """

    # number of systems in current population
    nsystem = df_input[col_system].unique().shape[0]
    # number of views to keep for each system (nsystem*nviews = simulate_nstars)
    nviews = int(simulate_nstars/nsystem)
    # number of views with non-zero transiting planets (Kept for later sanity check on the length of final dataframe)
    # creating two of these for two different purposes 
    nsystem_nonzero_views1 = 0  # via number of observers who see atleast a planet
    nsystem_nonzero_views2 = 0  # via loop appending all such potentially observable planets
    # total number of transiting planets (for later sanity check)
    nplanet_transiting = 0
    # Start a KOBE output dataframe
    df_kobe_output = pd.DataFrame()
    

    for index_system in range(nsystem):
    # Loop 2
    # This loops over all the systems present in the current ref_red file.

        if print_details == True:
            print('Calculating Transit Shadow Band for system - %d'%(index_system+1))

        # step 5 call kepler inclination or call kepler shadows
        # df_shadow = kobe.calculate_shadows(df_input=df_input, dict_column=dict_input,system_id=index_system,npoints=1e7,print_probability=print_details,grazing_transits=grazing_transits)
        # keep only nviews
        # np.random.seed(42)
        # df_shadow = df_shadow.loc[np.random.choice(df_shadow.index, size=nviews, replace=False)]

        # Here, we do a small shortcut to cut down computation time. 
        # Instead of calculating shadows at 1e7 point and then randomly selecting few views.
        # We calculate shadows for only required views.
        df_shadow, df_angle_data = calculate_shadows(df_input=df_input, system_id=index_system+1,npoints=nviews,print_probability=print_fine_details,grazing_transits=grazing_transits)

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
    
    return df_kobe_output, [nsystem_nonzero_views1,nsystem_nonzero_views2,nplanet_transiting]

def kobe_transits(df_kobe_output):
    """
    The KOBE Transits module: calculate transit related parameters for all planets

    Parameters
    ----------
    df_input : pandas dataframe
        input dataframe processed suitably

    """

    # Calculate impact parameter for each planet with respect to the observer
    # we use observer's azimuth as true anomaly --> implies that planet is crossing observer's azimuth
    # we use observer's polar angle as inclination --> identical by definition
    distance_to_star = (df_kobe_output[col_sma]*sp_constants.astronomical_unit*(1-(df_kobe_output[col_ecc]**2)))/(1 + (df_kobe_output[col_ecc]*np.cos(df_kobe_output['kobe_observer_azimuth'])))
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
    df_kobe_output['transit_duration [hours]'] = ((df_kobe_output[col_r_star]*radius_sun.value)*(df_kobe_output[col_period]*24))/((np.pi)*(df_kobe_output[col_sma]*sp_constants.astronomical_unit))
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

    return df_kobe_output

def kobe_vetter(df_kobe_output,kepler_completeness):
    """
    The KOBE Vetter module: applies Kepler completenes.

    Parameters
    ----------
    df_input : pandas dataframe
        input dataframe processed suitably    
    
    kepler_completeness : binned information
        an output of calculate_completeness

    Note
    ----
    1. flag_completeness is
        'PC' for planet candidates
        'FP' for false positives
        '-1' for planets which are not tces

    """
    ####################################
    ####################################
    # KOBE - Kepler Completeness & Reliability
    ####################################
    ####################################
    # Two flags are created. 
    # 'flag_completeness' : 'FP' means planet is vetted as False Positive. 'PC' means planet is vetted as Planetary Candidate
    
    ## 0 - create x and y arrays to go through each bin     
    x_array = np.arange(completeness_period_range[0],completeness_period_range[1],completeness_period_range[2])
    y_array = np.arange(completeness_radius_range[0],completeness_radius_range[1],completeness_radius_range[2]) 
    
    ## 1 - sanity checks
    drop_c=0
    drop1_c=0
    parameter_space_exceeded_c= []

    
    ## Create Columns 
    ## first create columns and mark each as 'PC'
    df_kobe_output['flag_completeness'] = 'PC'
    ## then all planets that are not tce are marked with NaN
    condition_nan = df_kobe_output['kobe_tce']==0
    df_kobe_output.loc[df_kobe_output.loc[condition_nan].index, 'flag_completeness'] = -1
   

    for xind in range(x_array.shape[0]):
        for yind in range(y_array.shape[0]):
            ## 2 now we are inside a bin
            ## we will create a flag

            completeness_c = kepler_completeness[xind][yind]
            
            if completeness_c == -0.01:
                completeness_c = 0
                parameter_space_exceeded_c.append(1)
                
            drop_ratio_c = 1 - completeness_c
            
            ## all planets that fall into current bin and are tce are selected

            condition = (df_kobe_output[col_period]>=x_array[xind]) & (df_kobe_output[col_period]< x_array[xind]+completeness_period_range[2]) & (df_kobe_output[col_r_planet]>=y_array[yind]) & (df_kobe_output[col_r_planet]< y_array[yind]+completeness_radius_range[2]) & (df_kobe_output['kobe_tce']==1)
            planet_index = df_kobe_output.loc[condition].index.to_numpy()

            ## estimate number of planets to flag
            drop_planets_size_c = int(np.floor(drop_ratio_c*planet_index.shape[0]))
            drop_c = drop_c + drop_planets_size_c

        
            ## now we randomly select planets within this bin which are flagged as FP
            if planet_index.shape[0] != 0:
                drop_planets_index_c = list(np.random.choice(planet_index,size=drop_planets_size_c,replace=False))
                df_kobe_output.loc[drop_planets_index_c,'flag_completeness'] = 'FP'
                drop1_c = drop1_c + len(drop_planets_index_c)

    ## completeness finished
    # step 9 kepler/kobe multiplicity 
    # initalize columns
    df_kobe_output['kobe_c_multiplicity'] = 0
    # now calculate
    pps_system, pps_planet = planetspersystem(df_input=df_kobe_output.loc[df_kobe_output['flag_completeness']=='PC'], system_column='kobe_system_id',planet_column='kobe_planet_id')
    df_kobe_output.loc[df_kobe_output.loc[df_kobe_output['flag_completeness']=='PC'].index,'kobe_c_multiplicity'] = pps_planet

    return df_kobe_output

def print_header():
    """
    """
    header1 ="""
     -----------------------------------------------------------------------------------------           
     -----------------------------------------------------------------------------------------                                                                                 
        KKKKKKKKK    KKKKKKK     OOOOOOOOO     BBBBBBBBBBBBBBBBB   EEEEEEEEEEEEEEEEEEEEEE
        K:::::::K    K:::::K   OO:::::::::OO   B::::::::::::::::B  E::::::::::::::::::::E
        K:::::::K    K:::::K OO:::::::::::::OO B::::::BBBBBB:::::B E::::::::::::::::::::E
        K:::::::K   K::::::KO:::::::OOO:::::::OBB:::::B     B:::::BEE::::::EEEEEEEEE::::E
        KK::::::K  K:::::KKKO::::::O   O::::::O  B::::B     B:::::B  E:::::E       EEEEEE
          K:::::K K:::::K   O:::::O     O:::::O  B::::B     B:::::B  E:::::E             
          K::::::K:::::K    O:::::O     O:::::O  B::::BBBBBB:::::B   E::::::EEEEEEEEEE   
          K:::::::::::K     O:::::O     O:::::O  B:::::::::::::BB    E:::::::::::::::E   
          K:::::::::::K     O:::::O     O:::::O  B::::BBBBBB:::::B   E:::::::::::::::E   
          K::::::K:::::K    O:::::O     O:::::O  B::::B     B:::::B  E::::::EEEEEEEEEE   
          K:::::K K:::::K   O:::::O     O:::::O  B::::B     B:::::B  E:::::E             
        KK::::::K  K:::::KKKO::::::O   O::::::O  B::::B     B:::::B  E:::::E       EEEEEE
        K:::::::K   K::::::KO:::::::OOO:::::::OBB:::::BBBBBB::::::BEE::::::EEEEEEEE:::::E
        K:::::::K    K:::::K OO:::::::::::::OO B:::::::::::::::::B E::::::::::::::::::::E
        K:::::::K    K:::::K   OO:::::::::OO   B::::::::::::::::B  E::::::::::::::::::::E
        KKKKKKKKK    KKKKKKK     OOOOOOOOO     BBBBBBBBBBBBBBBBB   EEEEEEEEEEEEEEEEEEEEEE
            K E P L E R       O B S E R V E S       B E R N          E X O P L A N E T S
     -----------------------------------------------------------------------------------------
     Code: Kepler Observes Bern Exoplanets
     Starting date: 25.10.2019
     Last edit date: 09.12.2020
     Author: Lokesh Mishra (University of Bern & Geneva Observatory)
     Contact: www.lokeshmishra.com
     -----------------------------------------------------------------------------------------
     -----------------------------------------------------------------------------------------"""

    print(header1)


