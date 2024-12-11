import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import eigh_tridiagonal, eigvalsh
import scipy.linalg.lapack
from matplotlib import rc
from scipy import special
from scipy.special import gamma
import scipy.special as sc
import mpmath
from mpmath import *
import matplotlib
from matplotlib.ticker import AutoMinorLocator, MultipleLocator
import matplotlib.ticker
from scipy import interpolate

rc('text', usetex=True)
font = {"family" : "Latin Modern Roman", "size" : 12}
rc("font", **font)

'Loading data from telescopes'

data_l = open("data/Obs.csv")
plot_data_l=np.genfromtxt(data_l, delimiter=",", dtype=float)
Plot_data_lX=np.transpose(plot_data_l)[0]
Plot_data_lY=np.transpose(plot_data_l)[1]*4*np.pi*1e3
Plot_data_lX=np.delete(Plot_data_lX, 63)
Plot_data_lY=np.delete(Plot_data_lY, 63)
Interpol_l= interpolate.interp1d(Plot_data_lX,Plot_data_lY,kind='cubic')
data_l2 = open("data/Obs-2.csv")
plot_data_l2=np.genfromtxt(data_l2, delimiter=",", dtype=float)
Plot_data_l2X=np.transpose(plot_data_l2)[0]
Plot_data_l2Y=np.transpose(plot_data_l2)[1]*4*np.pi*1e3
Plot_data_l2X=np.delete(Plot_data_l2X, 63)
Plot_data_l2Y=np.delete(Plot_data_l2Y, 63)
Interpol_l= interpolate.interp1d(Plot_data_lX,Plot_data_lY,kind='cubic')
Interpol_l2= interpolate.interp1d(Plot_data_l2X,Plot_data_l2Y,kind='cubic')
#Interpol2= interpolate.interp1d(Plot_dataX,Plot_dataY,kind='cubic')
#

Plot_FluxM1 = np.loadtxt("data/Plot_FluxM1.txt")
Plot_FluxM2 = np.loadtxt("data/Plot_FluxM2.txt")

#####
folder = "data/"
#folder = ""
data_old = open(folder+"ITO23.csv")
plot_data_old=np.genfromtxt(data_old, delimiter=",", dtype=float)
Plot_data_oldX=np.transpose(plot_data_old)[0]
Plot_data_oldY=np.transpose(plot_data_old)[1]
#####
data_ALP = open(folder+"ALP.csv")
plot_data_ALP=np.genfromtxt(data_ALP, delimiter=",", dtype=float)
Plot_data_ALPX=np.transpose(plot_data_ALP)[0]
Plot_data_ALPY=np.transpose(plot_data_ALP)[1]
#Interpol_ALP= interpolate.interp1d(Plot_data_ALPX,Plot_data_ALPY,kind='cubic')
#####
data_OSC = open(folder+"OSC.csv")
plot_data_OSC=np.genfromtxt(data_OSC, delimiter=",", dtype=float)
Plot_data_OSCX=np.transpose(plot_data_OSC)[0]
Plot_data_OSCY=np.transpose(plot_data_OSC)[1]
######
data_CAST = open(folder+"CAST.csv")
plot_data_CAST=np.genfromtxt(data_CAST, delimiter=",", dtype=float)
Plot_data_CASTX=np.transpose(plot_data_CAST)[0]
Plot_data_CASTY=np.transpose(plot_data_CAST)[1]
######
data_ARC = open(folder+"ARC.csv")
plot_data_ARC=np.genfromtxt(data_ARC, delimiter=",", dtype=float)
Plot_data_ARCX=np.transpose(plot_data_ARC)[0]
Plot_data_ARCY=np.transpose(plot_data_ARC)[1]
######
data_Geo = open(folder+"Geomagnetic.csv")
plot_data_Geo=np.genfromtxt(data_Geo, delimiter=",", dtype=float)
Plot_data_GeoX=np.transpose(plot_data_Geo)[0]
Plot_data_GeoY=np.transpose(plot_data_Geo)[1]
######
data_Int = open(folder+"intergalectic.csv")
plot_data_Int=np.genfromtxt(data_Int, delimiter=",", dtype=float)
Plot_data_IntX=np.transpose(plot_data_Int)[0]
Plot_data_IntY=np.transpose(plot_data_Int)[1]
#####
######
data_BBN = open(folder+"BBN.csv")
plot_data_BBN=np.genfromtxt(data_BBN, delimiter=",", dtype=float)
Plot_data_BBNX=np.transpose(plot_data_BBN)[0]
Plot_data_BBNY=np.transpose(plot_data_BBN)[1]

#####
#####


f_mat=np.logspace(7.8,26,200)
f_mat_long=np.logspace(9,25,200)

colors = plt.cm.hot(np.linspace(0,1,15))
colors2 = plt.cm.Greens(np.linspace(0,1,7))
colors3 = plt.cm.PuBu(np.linspace(0,1,7))
colors4 = plt.cm.Reds(np.linspace(0,1,5))




fig=plt.figure(figsize=(10,6))
ax = fig.add_subplot()

ax.plot(f_mat,np.sqrt(2*np.pi)*np.sqrt(Interpol_l(f_mat)/(Plot_FluxM1)),'-',color=colors3[5],linewidth=3)
ax.plot(f_mat,np.sqrt(2*np.pi)*np.sqrt(Interpol_l(f_mat)/(Plot_FluxM2)),'-',color=colors2[2],linewidth=3)

const_B = np.vstack((f_mat,np.sqrt(2*np.pi)*np.sqrt(Interpol_l(f_mat)/(Plot_FluxM1)))).transpose()
decay_B = np.vstack((f_mat,np.sqrt(2*np.pi)*np.sqrt(Interpol_l(f_mat)/(Plot_FluxM2)))).transpose()

np.savetxt("bound_const_B.dat",const_B)
np.savetxt("bound_decay_B.dat",decay_B)


ax.scatter(Plot_data_oldX,Plot_data_oldY,marker='x',color=colors4[4],label='Single Pulsar, Ito 23', zorder = 11)

ax.plot(Plot_data_GeoX,Plot_data_GeoY,'--',color=colors4[2],linewidth=1,label='Geomagnetic')
ax.plot(Plot_data_IntX,Plot_data_IntY,'--',color=colors4[1],linewidth=1,label='Intergalactic')
ax.plot(Plot_data_OSCX,Plot_data_OSCY,'-',color='#01295f',linewidth=1.5,label='OSCAR')
ax.plot(Plot_data_ALPX,Plot_data_ALPY,'-',color='#437f97',linewidth=1.5,label='ALP')
ax.plot(Plot_data_ARCX,Plot_data_ARCY*1e12,'-',color='#849324',linewidth=1.5,label='ARCADE')
ax.plot(Plot_data_ARCX,Plot_data_ARCY*1e3,'-',color='#849324',linewidth=1.5,label='ARCADE')
ax.plot(Plot_data_CASTX,Plot_data_CASTY,'-',color='#ffb30f',linewidth=1.5,label='CAST')
ax.plot(pow(10,Plot_data_BBNX),pow(10,Plot_data_BBNY),'-.',color='black',linewidth=1.5,label='BBN')

ax.fill_between(Plot_data_OSCX,Plot_data_OSCY, 1e-5, color='#01295f',alpha=0.1)
ax.fill_between(Plot_data_ALPX,Plot_data_ALPY, 1e-5, color='#437f97',alpha=0.1)
ax.fill_between(Plot_data_ARCX,Plot_data_ARCY*1e12, 1e-5, color='#849324',alpha=0.1)
ax.fill_between(Plot_data_ARCX,Plot_data_ARCY*1e3, 1e-5, color='#849324',alpha=0.05)
ax.fill_between(Plot_data_CASTX,Plot_data_CASTY, 1e-5, color='#ffb30f',alpha=0.1)

ax.fill_between(f_mat,np.sqrt(2*np.pi*Interpol_l(f_mat)/(Plot_FluxM1)),1e-5, color=colors3[5],alpha=0.1)
ax.fill_between(f_mat,np.sqrt(2*np.pi*Interpol_l(f_mat)/(Plot_FluxM2)),1e-5, color=colors2[2],alpha=0.15)

ax.text(1e17/3, 1e-26/4, 'CAST ', fontsize = 14,color='#ffb30f',rotation=0)
ax.text(1e13, 4e-27, 'OSCQAR', fontsize = 14,color='#01295f',rotation=0)
ax.text(1e14/2.5, 1e-23/4, 'ALP', fontsize = 14,color='#437f97',rotation=0)
ax.text(1e10/1.5, 1e-25, 'ARCADE', fontsize = 14,color='#849324',rotation=00)
ax.text(1e13, 1e-36, 'BBN', fontsize = 16,color='black',rotation=-20)
# ax.text(1e20, 1e-12, 'Galactic NSs', fontsize = 14,color=colors2[3],rotation=0)
ax.text(1e20, 3e-25, r'\textbf{Decaying} $\vec{\mathbf{B}}$', fontsize = 18,color=colors2[5],rotation=-18, weight = 'bold', zorder = 10)
ax.text(1e20, 1.5e-28, r'\textbf{Constant} $\vec{\mathbf{B}}$', fontsize = 18,color=colors3[5],rotation=-20, weight = 'bold', zorder = 10)
ax.text(1e23, 2e-32, r'Geomagnetic', fontsize = 16,color=colors4[3],rotation=-23, zorder = 10)
ax.text(1e22, 1e-37, r'Galactic', fontsize = 16,color=colors4[2],rotation=-22, zorder = 10)
ax.text(2e20, 1e-18, r'Single-pulsar', fontsize = 16,color=colors4[4],rotation=0, zorder = 10)



ax.set_ylim(1e-38,1e-10)
ax.set_xlim(1e10/4,1.2e25)
ax.set_xlim(0.5e8,1e26)
#ax.set_xlim(1e10,1e17)
ax.set_xscale('log')
ax.set_yscale('log')
ax.tick_params(labelsize=22)
ax.tick_params( which='minor',size=3,labelsize=22)
ax.set_xlabel(r'$f\, \rm [Hz]$ ',fontsize=22)
ax.set_ylabel(r'$h_c$  ',fontsize=22)
#ax.legend(fontsize=15)
x_major = matplotlib.ticker.LogLocator(base = 10.0, numticks = 5)
ax.xaxis.set_major_locator(x_major)
x_minor = matplotlib.ticker.LogLocator(base = 10.0, subs = np.arange(1.0, 10.0) * 0.1, numticks = 10)
ax.xaxis.set_minor_locator(x_minor)
ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
ax.set_xticks([1e11,1e13,1e15,1e17, 1e19,1e21,1e23,1e25])
#ax.grid(linestyle='--', linewidth=0.1)
x_major = matplotlib.ticker.LogLocator(base = 10.0, numticks = 10)
ax.xaxis.set_major_locator(x_major)
x_minor = matplotlib.ticker.LogLocator(base = 10.0, numticks = 50)
ax.xaxis.set_minor_locator(x_minor)
ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
#
y_minor = matplotlib.ticker.LogLocator(base = 10.0, numticks = 50)
ax.yaxis.set_minor_locator(y_minor)
ax.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())


for axis in ['top', 'bottom', 'left', 'right']:
    ax.spines[axis].set_linewidth(2.5)

plt.tight_layout()

#plt.show()

plt.savefig('Bounds.pdf')

