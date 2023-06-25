/* This program is written for PICCF calculation and Monte-Carlo simulation.
 * Author: Pu Du 
 * Date: 2015.01.20                                                         */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_statistics.h>
#include <plplot/plplot.h>

/* if lc arrays are longer than this value, please increase this number */
#define LC_NMAX 50000

int piccf_mc(int nc, double *jdc, double *fc, double *efc, 
        int nl, double *jdl, double *fl, double *efl, 
        int nt, double tbeg, double tend, 
        int *num_mc, double *tau_cent_tot, double *tau_peak_tot);
double cal_mean(int n, double *x);
double cal_std(int n, double *x);
double cal_r(int n1, double *f1, int n2, double *f2);
void interpolate_lc(int n1, double *jd1, double *f1, int n2, double *jd2, double *f2, 
        double dt, int *ni1, double *jdi1, double *fi1, int *ni2, double *jdi2, double *fi2);
void piccf(int nc, double *jdc, double *fc, int nl, double *jdl, double *fl, 
    int nt, double tbeg, double tend, double *ccf_t, double *ccf_r, double *rmax, 
    double *tau_cent, double *tau_peak);
double cal_lag_peak_r(int nt, double *t, double *r);
double cal_lag_peak_t(int nt, double *t, double *r);
double cal_lag_centroid_t(int nt, double *t, double *r);
void plotrange(int num, double *x, double *x1, double *x2);
void plotrange0(int num, double *x, double *x1, double *x2);
void plotrange_jd(int nc, double *xc, int nl, double *xl, double *x1, double *x2);
void norm_bin_data(int n, double *x, int nbin, double bin1, double bin2, 
        double *binx, double *biny);
void plspec(int num, double *wave, double *flux);

int main()
{
    /* input lc names */
    char continuum_name[] = "../data_use/corona_58380_1sigma.txt";
    char line_name[] = "../data_use/radio_58380_296.txt";

    int nc = 0, nl = 0;
    double jdc[LC_NMAX], fc[LC_NMAX], efc[LC_NMAX];
    double jdl[LC_NMAX], fl[LC_NMAX], efl[LC_NMAX];

    /* open and read lc files */
    FILE *fptr = NULL;
    fptr = fopen(continuum_name, "r"); // continuum lc
    while (!feof(fptr))
    {
        fscanf(fptr, "%lf %lf %lf\n", &jdc[nc], &fc[nc], &efc[nc]);
        //printf("%i %f %e %e\n", nc, jdc[nc], fc[nc], efc[nc]);
        nc++;
    }
    fclose(fptr);

    fptr = fopen(line_name, "r"); // line lc
    while (!feof(fptr))
    {
        fscanf(fptr, "%lf %lf %lf\n", &jdl[nl], &fl[nl], &efl[nl]);
        //printf("%i %f %e %e\n", nl, jdl[nl], fl[nl], efl[nl]);
        nl++;
    }
    fclose(fptr);

    fptr = NULL;

    /* define array to store output ccf and lag distribution */

    int num_mc = 10000, nt = 451;
    double tbeg = -5.0, tend = 20.0;
    int nbin = 101; // bin number of lag distribution
    double binx[nbin], biny[nbin];
    double binxp[nbin], binyp[nbin];

    double rmax, tau_cent, tau_peak;
    double t[nt], r[nt], tau_cent_tot[num_mc], tau_peak_tot[num_mc], a[nt];

    /* print informations */
    printf("continuum LC - length:%i\n", nc);
    printf("line LC - length:%i\n", nl);
    printf("output time array - begin:%.2f end:%.2f length:%i step:%.2f\n", 
            tbeg, tend, nt, (tend - tbeg) / ((float) nt - 1.0));

    /* acf */
    piccf(nc, jdc, fc, nc, jdc, fc, nt, tbeg, tend, t, a, &rmax, &tau_cent, &tau_peak);
    //printf("rmax:%f tau_cent:%f tau_peak:%f\n", rmax, tau_cent, tau_peak);

    /* ccf */
    piccf(nc, jdc, fc, nl, jdl, fl, nt, tbeg, tend, t, r, &rmax, &tau_cent, &tau_peak);
    printf("rmax:%f tau_cent:%f tau_peak:%f\n", rmax, tau_cent, tau_peak);

    /* Monte-Carlo simulation */
    piccf_mc(nc, jdc, fc, efc, nl, jdl, fl, efl, nt, tbeg, tend, &num_mc, tau_cent_tot, tau_peak_tot);
    //norm_bin_data(num_mc, tau_cent_tot, nbin, tbeg, tend, binx, biny);

    /* output to file */
    int i;
    char output[] = "piccf_out.txt";
    fptr = fopen(output, "w");
    for (i = 0; i < nt; i++)
    {
        fprintf(fptr, "%.3f  %.3f\n", t[i], r[i]);
    }
    fclose(fptr);

    norm_bin_data(num_mc, tau_peak_tot, nbin, tbeg, tend, binxp, binyp);
    fptr = fopen("peak_dist.txt", "w");
    for (i = 0; i < nbin; i++)
    {
        fprintf(fptr, "%.3f  %.3f\n", binxp[i], binyp[i]);
    }
    fclose(fptr);

    norm_bin_data(num_mc, tau_cent_tot, nbin, tbeg, tend, binx, biny);
    fptr = fopen("cent_dist.txt", "w");
    for (i = 0; i < nbin; i++)
    {
        fprintf(fptr, "%.3f  %.3f\n", binx[i], biny[i]);
    }
    fclose(fptr);

    fptr = fopen("cent_mc.txt", "w");
    for (i = 0; i < num_mc; i++)
    {
        fprintf(fptr, "%.3f\n", tau_cent_tot[i]);
    }
    fclose(fptr);

    fptr = fopen("peak_mc.txt", "w");
    for (i = 0; i < num_mc; i++)
    {
        fprintf(fptr, "%.3f\n", tau_peak_tot[i]);
    }
    fclose(fptr);

    /* percentile */
    gsl_sort(tau_cent_tot, 1, num_mc);
    double median, upper_lim, lower_lim;
    median = gsl_stats_median_from_sorted_data(tau_cent_tot, 1, num_mc);
    lower_lim = gsl_stats_quantile_from_sorted_data(tau_cent_tot, 1, num_mc, 0.1585);
    upper_lim = gsl_stats_quantile_from_sorted_data(tau_cent_tot, 1, num_mc, 0.8415);
    printf("valid monte carlo simulation:%i\n", num_mc);
    printf("percentile(15.85, 50, 84.15) of tau centroid: %.2f %.2f %.2f\n", lower_lim, median, upper_lim);

    gsl_sort(tau_peak_tot, 1, num_mc);
    median = gsl_stats_median_from_sorted_data(tau_peak_tot, 1, num_mc);
    lower_lim = gsl_stats_quantile_from_sorted_data(tau_peak_tot, 1, num_mc, 0.1585);
    upper_lim = gsl_stats_quantile_from_sorted_data(tau_peak_tot, 1, num_mc, 0.8415);
    printf("valid monte carlo simulation:%i\n", num_mc);
    printf("percentile(15.85, 50, 84.15) of tau peak: %.2f %.2f %.2f\n", lower_lim, median, upper_lim);

    /* plot */
    plsdev("xcairo");
    plspage(50, 50, 900, 600, 400, 200);
    plinit();
    pladv(0);
    plschr(0, 0.8);
    plvpor(0.10, 0.65, 0.51, 0.91);
    double xlim1, xlim2, ylim1, ylim2;
    plotrange_jd(nc, jdc, nl, jdl, &xlim1, &xlim2);
    plotrange(nc, fc, &ylim1, &ylim2);
    plwind(xlim1, xlim2, ylim1, ylim2);
    plpoin(nc, jdc, fc, 17);
    double fcmin[nc], fcmax[nc];
    for (i = 0; i < nc; i++)
    {
        fcmin[i] = fc[i] - efc[i];
        fcmax[i] = fc[i] + efc[i];
    }
    plerry(nc, jdc, fcmin, fcmax);
    plbox("bcst", 0.0, 0, "bcnst", 0.0, 0);
    pllab("", "F#dcontinuum#u", "");

    plvpor(0.10, 0.65, 0.10, 0.50);
    plotrange(nl, fl, &ylim1, &ylim2);
    plwind(xlim1, xlim2, ylim1, ylim2);
    plpoin(nl, jdl, fl, 17);
    double flmin[nl], flmax[nl];
    for (i = 0; i < nl; i++)
    {
        flmin[i] = fl[i] - efl[i];
        flmax[i] = fl[i] + efl[i];
    }
    plerry(nl, jdl, flmin, flmax);
    plbox("bncst", 0.0, 0, "bcnst", 0.0, 0);
    pllab("MJD-58000", "F#dline#u", "");

    plvpor(0.66, 0.93, 0.51, 0.91);
    plwind(tbeg, tend, -0.99, 1.0);
    plline(nt, t, a);
    plbox("bcst", 0.0, 0, "bcmst", 0.0, 0);
    plmtex("r", 3.5, 0.5, 0.5, "ACF");

    plvpor(0.66, 0.93, 0.10, 0.50);
    plwind(tbeg, tend, -0.99, 1.0);
    plcol0(2);
    pllsty(2);
    //pljoin(tau_peak, -0.99, tau_peak, 1.0);
    pllsty(3);
    plcol0(3);
    pljoin(tau_cent, -0.99, tau_cent, 1.0);
    pllsty(1);
    plcol0(1);
    plline(nt, t, r);
    plbox("bcnst", 0.0, 0, "bcmst", 0.0, 0);
    pllab("Lag", "", "");
    plmtex("r", 3.5, 0.5, 0.5, "CCF");

    plvpor(0.66, 0.93, 0.10, 0.50);
    plotrange0(nbin, biny, &ylim1, &ylim2);
    plwind(tbeg, tend, 0.0, ylim2 / rmax);
    plcol0(2);
    //norm_bin_data(num_mc, tau_peak_tot, nbin, tbeg, tend, binx, biny, rmax);
    plspec(nbin, binxp, binyp);
    plcol0(3);
    plspec(nbin, binx, biny);
    plcol0(1);

    plend();

    return 0;
}

void norm_bin_data(int n, double *x, int nbin, double bin1, double bin2, 
        double *binx, double *biny)
{
    int i, tempn;
    //double maxbin = 0.0;
    double delta = 0.0;
    delta = (bin2 - bin1) / ((float) nbin - 1.0);
    for (i = 0; i < nbin; i++)
    {
        binx[i] = bin1 + delta * i;
        biny[i] = 0.0;
    }
    for (i = 0; i < n; i++)
    {
        tempn = (x[i] - bin1) / delta;
        if ((tempn >= 0) && (tempn <= nbin - 1))
        {
            biny[tempn] += 1;
        }
    }
    //maxbin = biny[0];
    //for (i = 0; i < nbin; i++)
    //{
    //    if (biny[i] >= maxbin)
    //    {
    //        maxbin = biny[i];
    //    }
    //}
    //printf("%f %f\n", maxbin, norm);
    for (i = 0; i < nbin; i++) biny[i] = biny[i] / ((double) n);
    //for (i = 0; i < nbin; i++) printf("%f\n", biny[i]);
}

int piccf_mc(int nc, double *jdc, double *fc, double *efc, 
        int nl, double *jdl, double *fl, double *efl, 
        int nt, double tbeg, double tend, 
        int *num_mc, double *tau_cent_tot, double *tau_peak_tot)
{
    /* Calculate PICCF of two light curves and apply Monte-Carlo simulation to
     * obtain the error bar of time lag.  nc is the length of continuum light
     * curves, and nl is the length of line light curves. fc and efc are the
     * array of continuum lc and the error bars, fl and efl are the array of
     * line lc and the error bars. nt is the length of output time lag array,
     * tbeg is the beginning of time lag array and tend is the end of that
     * array. */

    /* check if the order of light curve array is ascending */
    int i = 0;
    for (i = 0; i < nc - 1; i++)
    {
        if (jdc[i + 1] <= jdc[i])
        {
            printf("check if the order of continuum light curve array is ascending - no\n");
            printf("%f(%i) %f(%i)\n", jdc[i], i, jdc[i + 1], i);
            exit(-1);
        }
    }
    for (i = 0; i < nl - 1; i++)
    {
        if (jdl[i + 1] <= jdl[i])
        {
            printf("check if the order of line light curve array is ascending - no\n");
            printf("%f(%i) %f(%i)\n", jdl[i], i, jdl[i + 1], i);
            exit(-1);
        }
    }
    //printf("check if the order of light curve arrays is ascending - ok\n");

    /* calculate step of output time lag array */
    double dt;
    double t[nt], r[nt];
    dt = (tend - tbeg) / ((float) nt - 1.0);
    for (i = 0; i < nt; i++) t[i] = tbeg + dt * i;  /* initialization */ 

    /* Monte-Carlo simulation */
    const gsl_rng_type *T;
    gsl_rng *R_unic, *R_unil, *R_norm;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    R_unic = gsl_rng_alloc(T);
    gsl_rng_set(R_unic, (unsigned long) time(NULL));

    R_unil = gsl_rng_alloc(T);
    gsl_rng_set(R_unil, (unsigned long) (time(NULL) - 2));

    R_norm = gsl_rng_alloc(T);
    gsl_rng_set(R_norm, (unsigned long) (time(NULL) - 4));

    /* create input light curves */
    double rmax_tot[*num_mc];
    double rmax = 0.0, tau_cent = 0.0, tau_peak = 0.0;
    double jdmc[nc], fmc[nc], efmc[nc];
    double jdml[nl], fml[nl], efml[nl];
    double ccf_t[nt], ccf_r[nt];

    int nmc, nml;
    int cindex[nc], lindex[nl];

    int j = 0, k = 0;
    for (j = 0; j < *num_mc; j++)
    {
        /* initialize */
        nmc = 0;
        nml = 0;
        for (i = 0; i < nc; i++) cindex[i] = 0;
        for (i = 0; i < nl; i++) lindex[i] = 0;

        /* sampling continuum light curves */
        for (i = 0; i < nc; i++) cindex[gsl_rng_uniform_int(R_unic, nc)] += 1;
        for (i = 0; i < nc; i++)
        {
            //printf("%i %i ", i, cindex[i]);
            if (cindex[i] != 0)
            {
                jdmc[nmc] = jdc[i];
                fmc[nmc] = fc[i];
                efmc[nmc] = efc[i] / pow((float) cindex[i], 0.5);
                //printf("%i %f %e %e", nmc, jdmc[nmc], fmc[nmc], efmc[nmc]);
                nmc++;
            }
            //printf("\n");
        }
        for (i = 0; i < nmc; i++)
        {
            //printf("%e %e ", fmc[i], efmc[i]);
            fmc[i] = fmc[i] + efmc[i] * gsl_ran_gaussian_ziggurat(R_norm, 1.0);
            //printf("%e\n", fmc[i]);
        }

        /* sampling line light curves */
        for (i = 0; i < nl; i++) lindex[gsl_rng_uniform_int(R_unil, nl)] += 1;
        for (i = 0; i < nl; i++)
        {
            //printf("%i %i ", i, lindex[i]);
            if (lindex[i] != 0)
            {
                jdml[nml] = jdl[i];
                fml[nml] = fl[i];
                efml[nml] = efl[i] / pow((float) lindex[i], 0.5);
                //printf("%i %f %e %e", nml, jdml[nml], fml[nml], efml[nml]);
                nml++;
            }
            //printf("\n");
        }
        for (i = 0; i < nml; i++)
        {
            fml[i] = fml[i] + efml[i] * gsl_ran_gaussian_ziggurat(R_norm, 1.0);
        }

        piccf(nmc, jdmc, fmc, nml, jdml, fml, nt, tbeg, tend, ccf_t, ccf_r, &rmax, &tau_cent, &tau_peak);

        //printf("%i rmax:%f tau_cent:%f tau_peak:%f\n", j + 1, rmax, tau_cent, tau_peak);

        if (rmax >= 0.0)
        {
            rmax_tot[k] = rmax;
            tau_cent_tot[k] = tau_cent;
            tau_peak_tot[k] = tau_peak;
            k++;
        }

        //plsdev("xcairo");
        //plinit();
        //pladv(0);
        //plvpor(0.15, 0.95, 0.55, 0.95);
        //plwind(jdc[0] - 10.0, jdc[nc - 1] + 10.0, 5.8e-16, 8.3e-16);
        //plpoin(nc, jdc, fc, 3);
        //double clim1[nc], clim2[nc];
        //for (i = 0; i < nc; i++)
        //{
        //    clim1[i] = fc[i] - efc[i];
        //    clim2[i] = fc[i] + efc[i];
        //}
        //plerry(nc, jdc, clim1, clim2);
        //plcol0(3);
        //plpoin(nmc, jdmc, fmc, 3);
        //plcol0(1);
        //plbox("bncst", 0.0, 0, "bcnst", 0.0, 0);
        //plend();
    }
    *num_mc = k;

    //printf ("generator type: %s\n", gsl_rng_name (R_uni));
    //printf ("seed = %lu\n", gsl_rng_default_seed);
    //printf ("first value = %lu\n", gsl_rng_get (R_uni));

    /* test */
    //double x[1000000];
    //double y[1000000];
    //////double x[100000];
    //for (i = 0; i < 1000000; i++)
    //{
    //    y[i] = 1.0 + 0.1 * gsl_ran_gaussian_ziggurat(R_norm, 1.0);
    //    x[i] = gsl_rng_uniform_int(R_uni, 11);
    //    //printf("%f\n", x[i]);
    //}

    //plsdev("xcairo");
    //plinit();
    //plhist(1000000, x, -1.0, 11.0, 100, PL_HIST_DEFAULT);
    //plend();

    gsl_rng_free(R_unic);
    gsl_rng_free(R_unil);
    gsl_rng_free(R_norm);
    return 1;
}

void piccf(int nc, double *jdc, double *fc, int nl, double *jdl, double *fl, 
    int nt, double tbeg, double tend, double *t, double *r, double *rmax, 
    double *tau_cent, double *tau_peak)
{
    /* Calculate PICCF of two light curves.  nc is the length of continuum
     * light curves, and nl is the length of line light curves. fc is the array
     * of continuum lc, and fl is the array of line lc. nt is the length of 
     * output time lag array, tbeg is the beginning of time lag array and tend 
     * is the end of that array. */

    /* calculate step of output time lag array */
    int i = 0;
    double dt;
    double rcl[nt], rlc[nt];
    dt = (tend - tbeg) / ((float) nt - 1.0);
    for (i = 0; i < nt; i++) t[i] = tbeg + dt * i;  /* initialization */ 

    /* calculate CCF of interpolated continuum light curve vs. line light curve */
    int nic = 0, nil = 0;
    double jdic[LC_NMAX], fic[LC_NMAX], jdil[LC_NMAX], fil[LC_NMAX];
    for (i = 0; i < nt; i++)
    {
        interpolate_lc(nc, jdc, fc, nl, jdl, fl, t[i], &nic, jdic, fic, &nil, jdil, fil);
        rcl[i] = cal_r(nic, fic, nil, fil);
    }

    /* calculate CCF of continuum light curve vs. interpolated line light curve */
    for (i = 0; i < nt; i++)
    {
        interpolate_lc(nl, jdl, fl, nc, jdc, fc, -t[i], &nil, jdil, fil, &nic, jdic, fic);
        rlc[i] = cal_r(nil, fil, nic, fic);
    }

    /* calculate combined partially integrated CCF */
    for (i = 0; i < nt; i++) r[i] = 0.5 * (rcl[i] + rlc[i]);

    *tau_peak = cal_lag_peak_t(nt, t, r);
    *rmax = cal_lag_peak_r(nt, t, r);
    *tau_cent = cal_lag_centroid_t(nt, t, r);

    /* plot to check CCF */
    //plsdev("xcairo");
    //plinit();
    //pladv(0);
    //plvpor(0.15, 0.95, 0.15, 0.95);
    //plwind(t[0], t[nt - 1], -1.0, 1.0);
    //plline(nt, t, r);
    ////double tau_peak, tau_cent, rmax;
    ////tau_peak = cal_lag_peak_t(nt, t, r);
    ////rmax = cal_lag_peak_r(nt, t, r);
    ////tau_cent = cal_lag_centroid_t(nt, t, r);
    ////printf("%f %f\n", tau_peak, tau_cent);
    //pllsty(2);
    //plcol0(5);
    ////pljoin(tau_peak, -1.0, tau_peak, 1.0);
    //pllsty(3);
    //plcol0(6);
    ////pljoin(tau_cent, -1.0, tau_cent, 1.0);
    //pllsty(4);
    //plcol0(7);
    ////pljoin(t[0], rmax, t[nt - 1], rmax);
    ////pljoin(t[0], 0.8 * rmax, t[nt - 1], 0.8 * rmax);
    //pllsty(1);
    //plcol0(1);
    //plbox("bcnst", 0.0, 0, "bcnst", 0.0, 0);
    //pllab("time", "CCF", "");
    //plend();
}

double cal_lag_centroid_t(int nt, double *t, double *r)
{
    double rmax, rlim;
    rmax = cal_lag_peak_r(nt, t, r);
    rlim = 0.8 * rmax;
    //printf("%f %f\n", rmax, rlim);

    int i;
    double dt, area = 0.0, total = 0.0;
    dt = t[1] - t[0];
    for (i = 0; i < nt; i++)
    {
        if (r[i] >= rlim)
        {
            area += dt * r[i];
            total += dt * t[i] * r[i];
        }
    }

    return total / area;
}

double cal_lag_peak_t(int nt, double *t, double *r)
{
    int i;
    double maxr, maxt;
    maxt = t[0];
    maxr = r[0];
    for (i = 0; i < nt; i++)
    {
        if (r[i] >= maxr)
        {
            maxt = t[i];
            maxr = r[i];
        }
    }
    return maxt;
}

double cal_lag_peak_r(int nt, double *t, double *r)
{
    int i;
    double maxr, maxt;
    maxt = t[0];
    maxr = r[0];
    for (i = 0; i < nt; i++)
    {
        if (r[i] >= maxr)
        {
            maxt = t[i];
            maxr = r[i];
        }
    }
    return maxr;
}

double cal_r(int n1, double *f1, int n2, double *f2)
{
    /* calculate cross-correlation function */
    double mean1 = 0.0, mean2 = 0.0;
    double std1 = 0.0, std2 = 0.0;
    mean1 = cal_mean(n1, f1);
    mean2 = cal_mean(n2, f2);
    std1 = cal_std(n1, f1);
    std2 = cal_std(n2, f2);
    //printf("%e %e %e %e\n", mean1, std1, mean2, std2);

    int i = 0;
    double t[n2];
    for (i = 0; i < n2; i++)
    {
        t[i] = (f1[i] - mean1) * (f2[i] - mean2);
    }

    return cal_mean(n2, t) / std1 / std2;
}

double cal_mean(int n, double *x)
{
    /* calculate mean first */
    int i = 0;
    double mean = 0.0;
    for (i = 0; i < n; i++)
    {
        mean += x[i];
    }
    mean = mean / ((float) n);
    return mean;
}

double cal_std(int n, double *x)
{
    /* calculate mean first */
    int i = 0;
    double mean = 0.0;
    for (i = 0; i < n; i++)
    {
        mean += x[i];
    }
    mean = mean / ((float) n);

    /* calculate std */
    double std = 0.0;
    for (i = 0; i < n; i++)
    {
        std += (x[i] - mean) * (x[i] - mean);
    }
    std = pow(std / ((float) n), 0.5);
    return std;
}

void interpolate_lc(int n1, double *jd1, double *f1, int n2, double *jd2, double *f2, 
        double dt, int *ni1, double *jdi1, double *fi1, int *ni2, double *jdi2, double *fi2)
{
    /* interpolate light curve 1 and obtain corresponding points to (jd2 - dt)
     * in light curve 2. Output interpolated light curve are jdi1, fi1, jdi2,
     * fi2 and its length of arrays are ni1 and ni2.  ni1 must be equal to ni2.*/

    int i1 = 0, i2 = 0, ii1 = 0, ii2 = 0; /* index for searching interpolated points */
    double temp = 0.0;

    for (i2 = 0; i2 < n2; i2++)
    {
        temp = jd2[i2] - dt;
        if (temp == jd1[0])
        {
            jdi1[ii1] = temp;
            fi1[ii1] = f1[0];
            jdi2[ii2] = jd2[i2];
            fi2[ii2] = f2[i2];
            ii1++;
            ii2++;
        }
        else if (temp == jd1[n1 - 1])
        {
            jdi1[ii1] = temp;
            fi1[ii1] = f1[n1 - 1];
            jdi2[ii2] = jd2[i2];
            fi2[ii2] = f2[i2];
            ii1++;
            ii2++;
        }
        else if ((temp > jd1[0]) && (temp < jd1[n1 - 1]))
        {
            while (temp > jd1[i1]) i1++;
            jdi1[ii1] = temp;
            fi1[ii1] = f1[i1 - 1] + (f1[i1] - f1[i1 - 1]) / (jd1[i1] - jd1[i1 - 1]) 
                * (temp - jd1[i1 - 1]);
            jdi2[ii2] = jd2[i2];
            fi2[ii2] = f2[i2];
            ii1++;
            ii2++;
        }
    }
    *ni1 = ii1;
    *ni2 = ii2;
    //printf("ori con:%i lin:%i inte con:%i lin:%i\n", n1, n2, ni1, ni2);

    /* plot to check */
//    int i = 0;
//    double xlim1 = 0.0, xlim2 = 0.0;
//    double xlim1a, xlim2a;
//    double xlim1b, xlim2b;
//    double ylim1, ylim2;
//    plotrange(n1, jd1, &xlim1a, &xlim2a);
//    plotrange(n2, jd2, &xlim1b, &xlim2b);
//    if (xlim1a >= xlim1b) xlim1 = xlim1b;
//    else xlim1 = xlim1a;
//    if (xlim2a >= xlim2b) xlim2 = xlim2a;
//    else xlim2 = xlim2b;
//    plsdev("xcairo");
//    plinit();
//    pladv(0);
//    plvpor(0.15, 0.9, 0.55, 0.9);
//    plotrange(n1, f1, &ylim1, &ylim2);
//    plwind(xlim1, xlim2, ylim1, ylim2);
//    pllsty(2);
//    plcol0(9);
//    for (i = 0; i < n2; i++) pljoin(jd2[i] - dt, ylim1, jd2[i] - dt, ylim2);
//    plcol0(1);
//    pllsty(1);
//    pljoin(xlim1, f1[0], jd1[0], f1[0]);
//    pljoin(xlim2, f1[n1 - 1], jd1[n1 - 1], f1[n1 - 1]);
//    plline(n1, jd1, f1);
//    plcol0(3);
//    plpoin(n1, jd1, f1, 2);
//    plcol0(2);
//    plpoin(ni1, jdi1, fi1, 21);
//    plcol0(1);
//    plbox("bcst", 0.0, 0, "bcnst", 0.0, 0);
//    pllab("", "F1", "");
//
//    plvpor(0.15, 0.9, 0.15, 0.5);
//    plotrange(n2, f2, &ylim1, &ylim2);
//    plwind(xlim1, xlim2, ylim1, ylim2);
//    pllsty(2);
//    plcol0(9);
//    for (i = 0; i < i2; i++) pljoin(jd2[i], ylim1, jd2[i], ylim2);
//    plcol0(1);
//    pllsty(1);
//    plline(n2, jd2, f2);
//    plcol0(3);
//    plpoin(n2, jd2, f2, 2);
//    plcol0(2);
//    plpoin(ni2, jdi2, fi2, 21);
//    plcol0(1);
//    plbox("bcnst", 0.0, 0, "bcnst", 0.0, 0);
//    pllab("time", "F2", "");
//    plend();
}

void plotrange_jd(int nc, double *xc, int nl, double *xl, double *x1, double *x2)
{
    // calculate plot jd range of light c and l.
    double xlimc1, xlimc2, xliml1, xliml2;
    plotrange(nc, xc, &xlimc1, &xlimc2);
    plotrange(nl, xl, &xliml1, &xliml2);
    
    if (xlimc1 >= xliml1) *x1 = xliml1;
    else *x1 = xlimc1;

    if (xlimc2 >= xliml2) *x2 = xlimc2;
    else *x2 = xliml2;
}

void plotrange(int num, double *x, double *x1, double *x2)
{
    // calculate the x or y plot range.
    int i;
    double dx;
    *x1 = x[0];
    *x2 = x[0];
    for (i = 0; i < num; i++)
    {
        if (x[i] <= *x1)
        {
            *x1 = x[i];
        }
        if (x[i] >= *x2)
        {
            *x2 = x[i];
        }
    }
    dx = *x2 - *x1;
    *x1 = *x1 - 0.1 * dx;
    *x2 = *x2 + 0.1 * dx;
}

void plotrange0(int num, double *x, double *x1, double *x2)
{
    // calculate the x or y plot range.
    int i;
    double dx;
    *x1 = x[0];
    *x2 = x[0];
    for (i = 0; i < num; i++)
    {
        if (x[i] <= *x1)
        {
            *x1 = x[i];
        }
        if (x[i] >= *x2)
        {
            *x2 = x[i];
        }
    }
    dx = *x2 - *x1;
}

void plspec(int num, double *wave, double *flux)
{
    /* plot spectrum using plplot */
    int i;
    double penx1, peny1, penx2, peny2;
    penx1 = wave[0];
    peny1 = flux[0];
    penx2 = 0.5 * (wave[1] - wave[0]) + wave[0];
    peny2 = flux[0];
    pljoin(penx1, peny1, penx2, peny2);
    for (i = 1; i < num - 1; i++)
    {
        penx1 = penx2;
        peny1 = flux[i];
        pljoin(penx2, peny2, penx1, peny1);
        peny2 = peny1;
        penx2 = penx1 + 0.5 * (wave[i] - wave[i - 1]) + 0.5 * (wave[i + 1] - wave[i]);
        pljoin(penx1, peny1, penx2, peny2);
    }
    penx1 = penx2;
    peny1 = flux[num - 1];
    pljoin(penx2, peny2, penx1, peny1);
    penx2 = wave[num - 1];
    peny2 = flux[num - 1];
    pljoin(penx1, peny1, penx2, peny2);
}
