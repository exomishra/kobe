# How does KOBE work?

## The physical picture
The physical aspect of KOBE is explained in detail in my paper [Mishra et. al. 2020](https://ui.adsabs.harvard.edu/abs/2021arXiv210512745M). Please see this paper for obtaining a physical understanding of how KOBE works.

## The computational aspect

This part may not be very user friendly. Consider contacting me or see my_notes!

The computational part is sub-divided into three parts:
1. kobe.py  
Contains the definitions of functions which do the BIG calculations.

2. kobe_driver.py  
The main running script. It assembles everything that's necessary for running KOBE.

3. kobe_choices.py  
There are several choices that a user has to specify. These are put together in this file for convenience. 

In addition, two other files are necessary for running KOBE. 
    * Kepler stellar data from DR25. We read stellar noice CDPP from this file.
    * Kepler robovetter simulated data results. We read the results and calculate completeness and reliability from this data.
    
***
