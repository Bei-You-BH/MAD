#!/usr/bin/env python

"""
A interface to piccf_mc.c. It's used for PICCF calculation and
Monte-Carlo simulations.
"""

import ctypes
import numpy as np

def piccf(jdc, fc, jdl, fl, nt, tbeg, tend):
    """
    Calculate PICCF of two light curves.  jdc, fc are the arrays of continuum
    lc.  jdl, fl are the arrays of line lc. nt is the length of output time lag
    array, tbeg is the beginning of time lag array and tend is the end of that
    array.
    """
    cnc = ctypes.c_int(len(jdc))
    cjdc = (ctypes.c_double * len(jdc))()
    cfc = (ctypes.c_double * len(jdc))()

    for i in range(len(jdc)):
        cjdc[i] = jdc[i]
        cfc[i] = fc[i]

    cnl = ctypes.c_int(len(jdl))
    cjdl = (ctypes.c_double * len(jdl))()
    cfl = (ctypes.c_double * len(jdl))()

    for i in range(len(jdl)):
        cjdl[i] = jdl[i]
        cfl[i] = fl[i]

    cnt = ctypes.c_int(nt)
    ctbeg = ctypes.c_double(tbeg)
    ctend = ctypes.c_double(tend)
    ct = (ctypes.c_double * nt)()
    cr = (ctypes.c_double * nt)()
    crmax = ctypes.c_double(0.0)
    ctau_cent = ctypes.c_double(0.0)
    ctau_peak = ctypes.c_double(0.0)

    cdll = ctypes.CDLL("/data/users/yuc/0624_ccf_test/libpiccf_mc.so")
    cdll.piccf(cnc, cjdc, cfc, cnl, cjdl, cfl, cnt, ctbeg, 
            ctend, ct, cr, ctypes.byref(crmax), ctypes.byref(ctau_cent), ctypes.byref(ctau_peak))

    t = []
    r = []
    for i in range(nt):
        t.append(ct[i])
        r.append(cr[i])
    rmax = crmax.value
    tau_cent = ctau_cent.value
    tau_peak = ctau_peak.value

    (t, r) = np.array((t, r))

    return t, r, rmax, tau_cent, tau_peak

def piccf_mc(jdc, fc, efc, jdl, fl, efl, nt, tbeg, tend, num_mc):
    """
    Calculate PICCF of two light curves and apply Monte-Carlo simulation to
    obtain the error bar of time lag. fc and efc are the array of continuum lc
    and the error bars, fl and efl are the array of line lc and the error bars.
    nt is the length of output time lag array, tbeg is the beginning of time
    lag array and tend is the end of that array. num_mc is the number of
    Monte-Carlo simulation.
    """
    cnc = ctypes.c_int(len(jdc))
    cjdc = (ctypes.c_double * len(jdc))()
    cfc = (ctypes.c_double * len(jdc))()
    cefc = (ctypes.c_double * len(jdc))()

    for i in range(len(jdc)):
        cjdc[i] = jdc[i]
        cfc[i] = fc[i]
        cefc[i] = efc[i]

    cnl = ctypes.c_int(len(jdl))
    cjdl = (ctypes.c_double * len(jdl))()
    cfl = (ctypes.c_double * len(jdl))()
    cefl = (ctypes.c_double * len(jdl))()

    for i in range(len(jdl)):
        cjdl[i] = jdl[i]
        cfl[i] = fl[i]
        cefl[i] = efl[i]

    cnt = ctypes.c_int(nt)
    ctbeg = ctypes.c_double(tbeg)
    ctend = ctypes.c_double(tend)
    cnum_mc = ctypes.c_int(num_mc)
    ctau_cent_tot = (ctypes.c_double * num_mc)()
    ctau_peak_tot = (ctypes.c_double * num_mc)()

    cdll = ctypes.CDLL("/data/users/yuc/0624_ccf_test/libpiccf_mc.so")
    cdll.piccf_mc(cnc, cjdc, cfc, cefc, cnl, cjdl, cfl, cefl, cnt, 
            ctbeg, ctend, ctypes.byref(cnum_mc), ctau_cent_tot, ctau_peak_tot)

    tau_cent_tot = []
    tau_peak_tot = []
    for i in range(num_mc):
        tau_cent_tot.append(ctau_cent_tot[i])
        tau_peak_tot.append(ctau_peak_tot[i])

    (tau_cent_tot, tau_peak_tot) = np.array((tau_cent_tot, tau_peak_tot))

    return tau_cent_tot, tau_peak_tot

def main():

    import sys
    import matplotlib.pyplot as plt
    data = np.loadtxt(sys.argv[1])
    data = data.tolist()
    data.sort(key=lambda a:a[0])
    data = np.array(data)
    jdc = data[:,0]
    jdl = data[:,0]
    fc = data[:,1]
    efc = (data[:,3]-data[:,2])/2
    fl = data[:,4]
    efl = (data[:,6]-data[:,5])/2
#     f = open(sys.argv[1])
#     l = f.readlines()
#     f.close()
#     jdc = [float(i.split()[0]) for i in l]
#     fc = [float(i.split()[1]) for i in l]
#     efc = [float(i.split()[2]) for i in l]
#     jdl = [float(i.split()[0]) for i in l]
#     fl = [float(i.split()[4]) for i in l]
#     efl = [float(i.split()[5]) for i in l]

    t, r, rmax, tau_cent, tau_peak = piccf(jdc, fc, jdl, fl, 1001, -50.0, 50.0)
    print('ccf rmax:%.3f tau_cent:%.3f tau_peak:%.3f' % (rmax, tau_cent, tau_peak))

    #output = open('pypiccf_ccf_out.txt', 'w')
    #output.write('t  r\n')
    #for i in xrange(len(t)):
    #    text = '%.3f %.3f\n' % (t[i], r[i])
    #    output.write(text)
    #output.close()

    tau_cent_mc, tau_peak_mc = piccf_mc(jdc, fc, efc, jdl, fl, efl, 1001, -50.0, 50.0, 5000)
    print('Monte-Carlo simulation:%.3f %.3f %.3f' % (np.percentile(tau_cent_mc, 15.85), 
            np.percentile(tau_cent_mc, 50.0), np.percentile(tau_cent_mc, 84.15)))

    output = open('pypiccf_mc_out.txt', 'w')
    output.write('tau_cent  tau_peak\n')
    for i in range(len(tau_cent_mc)):
        text = '%.3f %.3f\n' % (tau_cent_mc[i], tau_peak_mc[i])
        output.write(text)
    output.close()

if __name__ == "__main__":
    main()

