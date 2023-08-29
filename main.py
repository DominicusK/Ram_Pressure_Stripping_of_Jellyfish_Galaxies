# from Evolution_Class import Particle_Class_3D_3_1_1 as p3d_1
# from Evolution_Class import Particle_Class_3D_3_1 as p3d
from Evolution_Class import Particle_Class_3D_3_2 as p3d
from Graphing import RAM_Plots as rmplt
import pickle
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import math as m
import mpmath as mp
from matplotlib import gridspec
from mpl_toolkits.axes_grid1 import make_axes_locatable
from galpy.potential import plotRotcurve, HernquistPotential, BurkertPotential, MiyamotoNagaiPotential, MN3ExponentialDiskPotential,NFWPotential
from astropy import units
import scipy as sp
import scipy.integrate as integrate
import sys
from Analysis_Code import Analysis as al

#np.set_printoptions(threshold=sys.maxsize)

### EXAMPLe OF HOW DATA IS LOADED INTO MAIN.
## loading files.
# file1=open("RP_Runs/JO201_Bellhouse/Mass_disk_RAM_0_72Gyr_0_5_Beta.p","rb")
# file2=open("RP_Runs/JO201_Bellhouse/Mass_disk_data_0_72Gyr_0_5_Beta.p","rb")
# file3=open("RP_Runs/JO201_Bellhouse/Mass_disk_RAM_0_91Gyr_0_5.p","rb")
# file4=open("RP_Runs/JO201_Bellhouse/Mass_disk_data_0_91Gyr_0_5.p","rb")
#
# RAMs=pickle.load(file1)
# data=pickle.load(file2)
#
# RAMs2=pickle.load(file3)
# data2=pickle.load(file4)

## Loadind Disk Setup data.
# fname3 = 'Disk_Setup/JO201_24kpc_V1_1.txt'
# part_data3 = np.loadtxt(fname3)
# xpart_JO= part_data3[:,0]
# ypart_JO =part_data3[:,1]
# zpart_JO =part_data3[:,2]
# thetapart_JO=part_data3[:, 3]
# rpart2_JO=part_data3[:, 4]
# surf_rho_JO=part_data3[:, 5] #*1.05026504E-34
# mass_JO=part_data3[:, 6]
# p_area_elm_JO=part_data3[:,7]
# #
# MO2kg=1.989e+30
# fname4 = 'Disk_Setup/JO201_24kpc_V1_2.txt'
# part_data4= np.loadtxt(fname4)
# rpart1_JO= part_data4[:, 0]
# mass_enc_JO =part_data4[:, 1] #/(1e10*MO2kg)

## Loading CLuster data
# cluster_data=open("Clusters/ST2019_Mock_Cluster/MOCK_cluster_TVST2019.p", "rb")
# data_CL=pickle.load(cluster_data)
# wind_velocity=data_CL[:,4]
# radial_infall=data_CL[:,3]

## Loading Analytic Stripping radius data
# strip_radius_1=open("Analysis_Code/cw_0_5_10kpc_JO201_72_1890_va.p", "rb")
# strip_radius_2=open("Analysis_Code/cw_0_94_10kpc_JO201_72_1594_va.p", "rb")
# strip_radius_3=open("Analysis_Code/cw_1_5_10kpc_JO201_72_1299_va.p", "rb")
#
# R_st_1=pickle.load(strip_radius_1)
# R_st_2=pickle.load(strip_radius_2)
# R_st_3=pickle.load(strip_radius_3)
###

#Constants.
GRAV = 6.674E-20               # Gravitational constant (km^3/kg*s^2)
pc2cm = 3.086E18
cm2kpc = 3.2408E-22
pc2km = 3.086E+13
kg2M0= 5.0279E-31
npart_in_seg=20
MO2kg=1.989e+30
