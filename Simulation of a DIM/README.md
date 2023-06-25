## DIM code

The disc instability code which models the time-dependent evolution of the accretion disk (**Fig.4**)

## References

J.-M. Hameury, K. Menou, G. Dubus, J.-P. Lasota, J.-M. Hure, 1998, MNRAS 298, 1048.

Dubus, G., Hameury, J.M., Lasota, J.P., 2001, A&A 373, 251.

## How to use it

1. compile:

```
gfortran code.f -L /usr/X11R6/lib -lX11 -L /data/pgplot/ -lpgplot 
```

2. run: 

```
./a.out
```



`Data.txt` and `Data_nowind.txt` are the outputs by the code running with and without wind, which is used to explain the delayed optical emission and generate **Figure 4** (using the column "t" and "mv"). 

The only difference between the wind and nowind code is the value of the parameter f and g .

The time t=15222(days) in the column "t" corresponds to MJD 58380 in this work.

mv - 14.85 gives the absolute magnitude Mv and the simulated luminosity by the code running, i.e. 

$L_V=3.85\times10^{33}\times10^{-0.4(-4.81+m_V-14.85)} \rm{erg \ s^{-1}}$.
