### Determine the the lags between flares of different wavelengths:
- `piccf_mc.c`: the script to determine the time-delay
- `iccf_xr.py`: is used for generate the correlation figure between x-ray and radio data.
- `iccf_xv.py`: is used for generate the correlation figure between x-ray and V-band data.
- `corona_58380_1sigma.txt`: Compton luminosity after MJD 58380 in **Fig. 2A**.
- `radio_58380_296.txt`: Radio luminosity after MJD 58380 in **Fig. 2B**.
- `V_58380_0620.txt`: V-band luminosity after MJD 58380 in **Fig. 2C** (orange points).
### How-to-use:

1. Compile: 

   ```
    gcc -O2 piccf_mc.c -o libpiccf_mc.out -lm -lgsl -lplplot -lgslcblas 
   ```

2. Run: 

   ```
    ./libpiccf_mc.out
   ```

  Note that, in the script `piccf_mc.c` , after the "continuum_name" and "line_name" are specified with the respective filenames for two light currves (see above), the running will return the resulting time-delay.

### How-to-plot (Figure S1):
 `python iccf.py` Note that, after the "data1" and "data2" are specified with the respective filenames for two light currves (see above), the script will plot the resulting ICCF analysis.

  `piccf_mc.py`: the script invoked by  `iccf.py`, to use this script, you should first change the cdll = ctypes.CDLL parameter to the path which stored the  `libpiccf_mc.so` file.
