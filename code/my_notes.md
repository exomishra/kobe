Steps in the running of KOBE
----------------------------

In this document, I will list the steps and choices implemented in KOBE (17.02.2020).

----------------------------

Preliminary Notes:

1. The choices that should be taken as input from a user are now put into a separate python file called kobe_choices.py.

2. Two files are necessary for the running of KOBE:
    a) Kepler stellar data from DR25. We read CDPP from this file.
    b) Kepler robovetter simulated data results. We read the results and calculate completeness and reliability from this data.

3. These files are in the folder 'stellar_data' which is in the same location as the kobe_driver.

4. The number of 'views' or the number of times the shadows of a system is calculated is decided by:
    - number of stars to simulate (a choice in kobe_choices)
    - number of systems in a given ref_red file (after removing planets that did not end well)
    - Thus, number of views = number of stars to simulate/number of systems


Main Steps in KOBE:
-------------------

1. We read constants from the astropy package. (radius and mass of sun and earth)

2. A seed is given to ensure reproducibility for each of the random choices made in the program.

3. Datasets are labelled, depending on the choice of population in kobe_choices.

4. step 1a - setup for reading data
    a) directory setup for my local machine and horus is defined
    b) the program learns to read name of ref_red files and also get their time from the filename
    c) a list of filenames, sorted by time is prepared. (the program will follow this order of working)
    
5. step 1b - read stellar data (for cdpp calculation)
    a) read file and remove entries where cdpp is NaN.
    b) Place cuts on type of star: 
        Mass    - 0.7 to 1.3 Solar masses
        Radius  - 0 to 5 Solar Radius
        Temp    - 3880 to 7200 K (selection for FGK stars  - temperatures from Peacaut and Mamajek 2013)
    c) create numpy arrays of cdpp for t= 3, 6, 9, 12 hours
    
6. step 1c - Kepler completeness
    Note - since we are unsure of the best way to  incorporate Kepler's reliability, we create two pathways.
    Path 1 - will include completeness for all planets with a disposition score of 0 and above 
    (i.e. all planets are considered)
    Path 2 - will include completeness for all planets with a disposition score of 0.9 and above. 
    (proxy - for incorporating reliability, see Fig 31 in Bryson 2019, or Mulders et. al. 2018)
    
    a) User decides the following:
        -- period range : range for which completeness is calculated
        recommended settings [0,700] days with 50 days in 1 bin
        
        -- radius range : range for which completeness is calculated
        recommended settings [0, 20] earth radii with 2 rearth in 1 bin
        
        injection threshold : if number of injections in a bin are below this number, then this bin is flagged with value = -0.01 (so, in percent it gives -1%)
        recommended settings = 10
        
        spacing : 'linear' or 'log' distribution of bins
        
        disposition_cutoff : number below which planets are ignored
        
    b) completeness is calculated
    
7. Step 2: Loop 1 -- loop over files
    a) We use this loop to loop over all the ref_red files in a given population.
    b) A file is read. (only planets with emps_status = fine are kept)
    c) selected columns are extracted
    d) Preliminary steps are performed: 
        - correct system and planet number
        - add column for planet multiplicity
        - add column for mass in Mjupiters
        - add column for radius in Rearth
        - add column for inclination in degrees
        - add column for period in days
        - add column for density in g/cm3
        - add column for transit probability
        - sort planets by sma
        
    e) Calculate the Transit Signal for each planet (Rplanet/Rstar)^2

8. Step 1d: Get column numbers of some commonly used columns

9. Step 3: Loop 2 -- loop over systems
    a) In this loop, we go through all the systems in the current ref_red file. 
    b) Inside the loop, we call kobe.calculate_shadows on the current system:
        - user can specify, if grazing transits are included or not. (Excluded in current run)
    c) The output of shadows is a dataframe:
        - each row is one view of the system from a random point in space
        - each columns is for a planet: 
            -- 1 if they cast a shadow at the current viewing position (meaning they transit)
            -- 0 if they don't cast a shadow
    d) Empty rows (meaning no planet transit from this viewing angle) are removed
    e) The total number of transiting planets from all views is collected (in nplanet_transiting) for sanity check
    
10. Step 4: KOBE - Output file
    Loop 3 -- loop over all the views for each systems
    
    a) planets which transit in current view (current row in shadow output) are added to kobe output
    b) these planets form "a planetary system" -- same kobe_system_id
    c) For each view, a CDPP is randomly picked from the cdpp array created earlier.
    
11. Loop 3 is now over. This means our ouput has all transiting planetary systems from all views for the first system in the input file. We make a quick sanity check on the number of transiting planetary systems found, so far.

12. Loop 2 is now over. This means our output has all transiting planetary systems from all views for all the systems in the input file. 

13. For all planets we calculate:
    a) characteristic transit duration (hours) (eq 19, Transits and Occultations by J. Winn in Exoplanets)
    b) effective CDPP seen by each planet (eq 4, Christiansen et. al. 2012.)
    c) single event statistic for one transit
    d) number of transits
    e) multiple event statistic or mes for all transits. (Final, SNR)
    
14. KOBE - tce (threshold crossing event)
    a) Those planets which transit at least 3 times (specified in kobe_choices), and
    b) Have an SNR signal of at least 7.1 (threshold from kobe_choices)
    c) are marked as kobe_tce: 1 - means they are a tce.

15. KOBE - Completeness and reliability
    a) Two flags are created: flag_completeness and flag_completeness_reliability 
    b) the flags are: 'FP' is for False Positive and 'PC' is for Planetary Candidate (for TCE's) and NaN for not TCEs.
    c) Another loop to access and apply completeness on kobe's output. For each bin, the value of completeness calculated above is read. 
    - If it is 0.5, say. Then, 50% of planets in this bin (and are also a tce) are randomly vetted as 'PC' and the rest as 'FP'.
    - If it is 0.1, say. Then, 10% of planets in this bin (and are also a tce) are randomly vetted as 'PC' and the rest as 'FP'.
    - If it is -0.01. Then completeness is put to 0. Thus, all planets in this bin are 'FP'.
    d) kobe's multiplicity is calculated/

16. Sanity Checks:
    a) the number of systems with transiting planets in the output should match the number of systems counted from views
    b) same but for the number of transiting planets.
    
17. The following counts are output to file "kobe_basic_counting_info.txt"
    a) file name
    b) dataset name, number of system in original, number of planets
    c) kobe_tce - number of systems and number of planets
    d) kobe_completeness - number of systems and number of planets
    e) kobe_completeness_reliability - number of systems and number of planets
    
18. Additional columns are extracted from the dataframe and put back into the dictionary. The dictionary is saved.

19. The dataframe is saved as .csv file. 

20. Loop 1 now over. That means we have run over all the files in the population. 

21. Run times are printed.
