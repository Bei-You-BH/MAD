### Determine the the lags between flares of different wavelengths:
- `piccf_mc.c`: the script to determine the time-delay
- `iccf_xr.py`: is used for generate the correlation figure between x-ray and radio data (**Fig. S1A**).
- `iccf_xv.py`: is used for generate the correlation figure between x-ray and V-band data (**Fig. S1B**).
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

  Note that, in the script `piccf_mc.c` , after the "continuum_name" and "line_name" are specified with the respective filenames for two light curves (see above), the running will return the resulting time-delay.

### How-to-plot (Figure S1):

​    First change the `cdll = ctypes.CDLL` parameter in `piccf_mc.py` to the absolute path which stored the  `libpiccf_mc.so` file.

​    `python iccf_xr.py`: to generate the correlation figure between x-ray and radio data (**Fig. S1A**).

​    `python iccf_xv.py`: to generate the correlation figure between x-ray and V-band data (**Fig. S1B**).

  - Notice: the "data1" and "data2" parameters in `iccf_xr.py` and `iccf_xv.py` are specified with the respective filenames for two light curves.