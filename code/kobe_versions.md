The first version of KOBE (version 1.0) and the second version of KOBE (version 2.0) have the following key differences:

## KOBE version 1.0: 
1) Assumes a circular orbit for calculating the impact parameter. 
2) The transit duration is calculated assuming circular, coplanar orbits. 
3) The number of transits are calculated using an average estimation.
4) The threshold crossing events are identified by placing a cut on the SNR (>7.1). 
5) Grazing transits can be included in version 1.0, with the grazing transit signal estimated using the full-transit approximation of [Seager & Mallén-Ornelas, 2003](https://iopscience.iop.org/article/10.1086/346105/fulltext/). 
6) The duty cycle, which accounts for planned gaps for data download and spacecraft operations, is assumed to be 100%. 

## KOBE version 2.0: 
1) Accounts for eccentric orbits for calculating the impact parameter. 
2) The transit durations can be calculated using expressions from [Winn (2010)](https://arxiv.org/abs/1001.2010) or [Kipping (2010)](https://arxiv.org/abs/1004.3819). These transit durations capture the effect of eccentricity, impact parameter and the planet's varying speed in eccentric orbits. [Winn (2010)](https://arxiv.org/abs/1001.2010) equations also capture the effect of the radius of the planet (significant for planet radius > 15 earth radii).
3) The number of transits are calculated in a monte-carlo manner, which further improves the SNR estimate.
4) The detection efficiency is modelled on the lines of [Christiansen (2017)](https://exoplanetarchive.ipac.caltech.edu/docs/KSCI-19110-001.pdf), improving the identification of threshold crossing events with the SNR.
5) The transit signal is modelled using exact, analytic formulae as in [Mandel and Agol (2002)](https://iopscience.iop.org/article/10.1086/345520/pdf), improving the treatment of grazing transits.
6) The duty cycle is incorporated using a factor that accounts for the fraction of valid observations, as in [Burke et al. (2015)](https://iopscience.iop.org/article/10.1088/0004-637X/809/1/8). In case you are unsure about the duty cycle of the transit survey you are using KOBE for, turn duty_cycle to False. The calculations will then be performed taking the duty cycle by default as 100%. 

Both of the versions ultimately help you place the biases and limitations of a transit survey on your synthetic population of planets, thus facilitating a comparison with observations. You can choose between the versions based on the requirements of your study, and make internal choices for using/not using each of the specific features in each version.

