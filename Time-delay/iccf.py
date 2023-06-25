#!/usr/bin/env python

"""
A interface to piccf_mc.c. It's used for PICCF calculation and
Monte-Carlo simulations.
"""
import numpy as np
import matplotlib.pyplot as plt
import piccf_mc


data2 = np.loadtxt("04_data/radio_58380.txt")
data1 = np.loadtxt("05_data/corona_58380_cor.txt")


data1 = data1.tolist()
data1.sort(key=lambda a:a[0])
data1 = np.array(data1)
data2 = data2.tolist()
data2.sort(key=lambda a:a[0])
data2 = np.array(data2)


jdc = data1[:,0]
fc = data1[:,1]
efc = data1[:,2]

jdl = data2[:,0]
fl = data2[:,1]
efl = data2[:,2]

fig = plt.figure(figsize=(8, 6))

#ax = fig.add_subplot(311)
#ax.errorbar(jdc, fc, yerr=efc, ls='none')
#ax.errorbar(jdl, fl, yerr=efl, ls='none')

#ax.set_xlabel("mjd")
plt.rcParams['font.family'] = 'sans-serif'
ax1 = fig.add_subplot(111)
t, r, rmax, tau_cent, tau_peak = piccf_mc.piccf(jdc, fc, jdl, fl, 1001, -20.0, 20.0)


ax1.plot(t, r,color='red')
ax1.set_xlabel("lag (days)",fontsize=15)
ax1.set_ylabel("CCF",fontsize=15,color='red')

#ax = fig.add_subplot(212)

ax2 = ax1.twinx()
tau_cent_mc, tau_peak_mc = piccf_mc.piccf_mc(jdc, fc, efc, jdl, fl, efl, 1001, -20.0, 20.0, 5000)

tau_cent_mc = np.array(tau_cent_mc)
tau_cent_mc = np.delete(tau_cent_mc,np.where(tau_cent_mc==0.))
tau_peak_mc = np.array(tau_peak_mc)
tau_peak_mc = np.delete(tau_peak_mc,np.where(tau_peak_mc==0.))


bins = np.arange(-15.0,15.0,1.25)
#print(bins)
ax2.hist(tau_cent_mc, density=True, bins=bins,color='blue',linestyle='--',label='centroid',histtype='step')
ax2.hist(tau_peak_mc, density=True, bins=bins,color='blue',linestyle='dotted',stacked = True,alpha=0.5, label='peak',histtype = 'step')

#plt.hist(tau_peak_mc, density=True, bins=300,alpha=0.5, label='peak')

#ax.set_xlabel("lag")
ax2.legend()


#ax2.set_xlim(-10,10)
ax2.set_ylabel("density", fontsize=15,color='blue')
ax2.set_xlabel("lag (days)", fontsize=15)
ax2.tick_params(axis='both',direction='in',which='both',top=True, bottom=True, left=False, right=True, width=1.0)
ax2.tick_params(axis='both',direction='in',which='major',top=True, bottom=True, left=False, right=True, width=1.0, length=5)
ax1.tick_params(axis='both',direction='in',which='both',top=True, bottom=True, left=True, right=False, width=1.0)
ax1.tick_params(axis='both',direction='in',which='major',top=True, bottom=True, left=True, right=False, width=1.0, length=5)




ax2.spines['left'].set_color('red')
ax2.spines['left'].set_linewidth(1)
ax1.tick_params(axis='y',colors='red')
#ax1.yaxis.set_ticks_position('right')

ax2.spines['right'].set_color('blue')
ax2.spines['right'].set_linewidth(1)
ax2.tick_params(axis='y',colors='blue')
#ax2.yaxis.set_ticks_position('left')
ax1.set_xlim(-5,15)
ax2.set_xlim(-5,15)
ax1.set_ylim(-1.01,1.01)
ax1.text(-2.5,0.75,"A",fontsize=15,fontweight='bold')
#plt.show()
plt.title("Radio  vs.  X-ray",fontsize=20)
plt.savefig("xr.png",dpi=800,bbox_inches='tight')
plt.show()
#print(tau_cent_mc)


H = open("para.txt",'w')
for i in range(len(tau_cent_mc)):
    H.write(str(tau_cent_mc[i])+" ")
#H.write(tau_peak_mc)
H.close()
