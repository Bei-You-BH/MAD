## Figure S3-S6

- The script `code.ipynb`  is used to plot **Figure S3-S6**
- In the block `Figure S3`, the file `r_tr2.dat`  is used to plot **Figure S3**. This file is calculated by `adaf_tr.for` in which `lum_x.dat` is the input file contains the Compton, disk and total X-ray luminosity displayed in **Figure 2A**.
- In the block `Figure S4`,  `mad_bz_br32p8rm1p17w0p92c.dat` and  `mad_bz_br6p79rm0w0p92c.dat`  are used to plot **Figure S4**, which are calculated by the code `mad_b_config_h.for` with parameters r_tr=32.8r_isco, r_m=1.17, f_m=0.92 and r_tr=6.79r_isco; r_m=0, f_m=0.92,  respectively. For both two filed, the columns are arranged as Radius, B_z and B_r, 
- In the block  `Figure S5`,  `b_t.dat`  is used to plot **Figure S5**. The columns are arranged as: 
  (1) MJD; 
  (2) the simulated magnetic strength advected by the ADAF; 
  (3) the truncation radius; 
  (4)  the MAD criterion.
- In the block  `Figure S6`, the script is used to plot **Figure S6**, which calculates BZ power, $P_{\rm BZ}$, as a function of the BH spin parameter $a$ for a BH with 8.5 solar mass and the field strength at the horizon {$B_{\rm h}=7\times 10^7$}~Gauss is adopted.