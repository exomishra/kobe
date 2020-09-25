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
import mpl_toolkits.mplot3d.art3d as art3d
import seaborn as sns
sns.set()
from tqdm import tqdm_notebook as tqdm
from scipy import constants as scipy_constants
from astropy import constants as astropy_constants
from astropy import units as astropy_units
from mpl_toolkits.mplot3d import Axes3D

from mymodules import dataread, analysis, plottingroutine, kobe

# SI units
mass_sun = 1.9884*10**30
mass_earth = 5.9722*10**24
radius_sun = 6.957*10**8
radius_earth = 6.3781366*10**6

mass_sun = astropy_constants.M_sun
mass_earth = astropy_constants.M_earth
radius_sun = astropy_constants.R_sun
radius_earth = astropy_constants.R_earth



def calculate_shadows(df_input, dict_column, system_id, npoints=1e7, 
 stellar_radius=None,seed=None, print_probability=True, grazing_transits=False):
    """
    This function calculates the transit shadow bands for all planets in a chosen system.
    
    Parameters
    ----------
    
    df_input : pandas dataframe
        dataframe containing planets with their orbital elements
    dict_column : dictionary 
        dictionary of column names (key: number, value: column name)
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
        
    #get column numbers
    col_sma = dataread.find_key(dict_column,'Semi-Major Axis [AU]')
    col_r_star = dataread.find_key(dict_column,'Stellar Radius [Rsun]')
    col_r_planet = dataread.find_key(dict_column, 'Radius [Rearth]')
    col_ecc = dataread.find_key(dict_column,'eccentricity')
    col_inc = dataread.find_key(dict_column,'Inclination [radians]')
    col_long_node = dataread.find_key(dict_column,'long_node')
    col_long_peri = dataread.find_key(dict_column,'long_peri')
    col_system = dataread.find_key(dict_column,'System')
    col_planet = dataread.find_key(dict_column, 'Planet')
    col_planetmultiplicity = dataread.find_key(dict_column, 'Planet Multiplicity')
    col_transitprob = dataread.find_key(dict_column, 'Transit Probability')
    
    # Get system info: number of planets, stellar radius, etc..
    nplanets = df_input.loc[df_input[col_system]==system_id][col_planetmultiplicity].iat[0]
    if stellar_radius == None:
        r_star = df_input.loc[df_input[col_system]==system_id][col_r_star].mean()   * astropy_constants.R_sun
    else:
        r_star = stellar_radius * astropy_constants.R_sun
    series_sma = df_input.loc[df_input[col_system]==system_id][col_sma]             * astropy_constants.au
    series_r_planet = df_input.loc[df_input[col_system]==system_id][col_r_planet]   * astropy_constants.R_earth
    series_ecc = df_input.loc[df_input[col_system]==system_id][col_ecc]             * astropy_units.one   
    series_inc = df_input.loc[df_input[col_system]==system_id][col_inc]             * astropy_units.rad
    series_long_node = df_input.loc[df_input[col_system]==system_id][col_long_node] * astropy_units.rad
    series_long_peri = df_input.loc[df_input[col_system]==system_id][col_long_peri] * astropy_units.rad
    series_arg_peri = series_long_peri - series_long_node
    series_transitprob = df_input.loc[df_input[col_system]==system_id][col_transitprob]
    
    # Step 1 - create grid on sphere
    r = 1                       # unit sphere
    if seed == None:
        seed = 42               # meaning of life
    np.random.seed(seed)
    u1 = np.random.rand(npoints) # a random number generator
    np.random.seed(seed+1)
    u2 = np.random.rand(npoints) # a random number generator
    theta = 2*np.pi*u1                                                               * astropy_units.rad
    phi = np.arccos(1-2*u2)                                                          * astropy_units.rad
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
        v_prime = np.einsum('bc,ac',analysis.rotation_matrix('z',series_long_node.iat[index_planet],invert=True),v_prime)
        v_prime = np.einsum('bc,ac',analysis.rotation_matrix('x',series_inc.iat[index_planet],invert=True),v_prime)
        v_prime = np.einsum('bc,ac',analysis.rotation_matrix('z',series_arg_peri.iat[index_planet],invert=True),v_prime)
        v_prime_norm = np.linalg.norm(v_prime, axis=1)
        
        ### SANITY TEST 2 ####
        #sanity check: norm of all vectors should be 1 (1 is True in bool)
        # check that the difference between norm of each vector and 1 is less than 1e-15. 
        # If not replace it with 0 (0 is flase).
        check_norm = np.where(abs(v_prime_norm-1)<1e-15, v_prime_norm, 0)
        if check_norm.all() == False:
            print('Sanity Test 2 - FAILED. Not all vectors are unit normed. ')
        
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
        alpha = np.arccos(z_prime)         * astropy_units.rad        # in radians
        f = np.arccos(x_prime/r_xy_prime)  * astropy_units.rad        # in radians
        
        # Step 2b - Calculate pairs of (f, alpha) inside shadow
        if grazing_transits == True:
            height_shadow = ((r_star - series_r_planet.iat[index_planet])*(1 + series_ecc.iat[index_planet] * np.cos(f)))/((series_sma.iat[index_planet])*(1 - series_ecc.iat[index_planet]**2))
        elif grazing_transits == False:
            height_shadow = ((r_star + series_r_planet.iat[index_planet])*(1 + series_ecc.iat[index_planet] * np.cos(f)))/((series_sma.iat[index_planet])*(1 - series_ecc.iat[index_planet]**2))
        
        alpha_rad = np.arcsin(height_shadow)
        # for each f (true anomaly), we calculate the shadow band by calculating the min and max of shadow
        shadow_min = np.pi/2 * astropy_units.rad - alpha_rad
        shadow_max = np.pi/2 * astropy_units.rad + alpha_rad
        
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
        df_angle_data['(f,alpha) planet %d'%(index_planet+1)] = list(np.column_stack((f.value,alpha.value)))
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
    print('Shadow calculation: %d planets in system %d in %s population with %.1e grid points.\n%.2f %% observers see atleast one transiting planet.'%(nplanets, system_id,dict_column['name'], npoints,100*df_shadow_output.loc[df_shadow_output['number_transiting_planet']!=0].shape[0]/npoints))
    
    return df_shadow_output, df_angle_data

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
    print('Bins along x and y: ', total_injections.shape) 
    print('Number of bins with less than %d injections = %d'%(injection_threshold,thresholdbins))
    print('Number of empty bins = %d'%(emptybins))
    print('Number of total bins = %d'%(totalbins))
    metadata = 'Completeness Calculation Information: \n Column along x is %s with min %.2f, max %.2f and step size %d. \n Column along y is %s with min %.2f, max %.2f and step size %.2f. \n Order of output: this text, total injections, pc injections, completeness(over threshold), \n x_range(mix, max, step), y_range(min, max, step). \n Note: Completeness ranges from 0-1. \n Completeness of -1 flags when total injected TCEs are either 0 or less than threshold injections of %d. \n Spacing is %s. \n Disposition Score cutoff is at %.1f'%(x_col,x_range[0],x_range[1],x_range[2],y_col,y_range[0],y_range[1],y_range[2],injection_threshold, spacing, dispscore_cutoff)
    
    output =  np.array([metadata,total_injections, pc_injections, completeness, x_range, y_range,dispscore_cutoff], dtype=object)
    
    return output
