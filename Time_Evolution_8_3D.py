import os
import matplotlib.pyplot as plt
import math as m
import numpy as np
import mpmath as mp
from matplotlib import gridspec
from mpl_toolkits.axes_grid1 import make_axes_locatable
from tqdm import tqdm
import pickle
from Evolution_Class import Particle_Class_3D_3_2 as p3d

# Constants-not all are used but useful for plots.
GRAV = 6.674E-20               # Gravitational constant (km^3/kg*s^2)
pc2cm = 3.086E18
cm2kpc = 3.2408E-22
pc2km = 3.086E+13
kg2M0= 5.0279E-31
MO2kg=1.989e+30

# loading text file from disk setup and defining data .
fname1 = 'Disk_Setup/JO201_24kpc_V1_1.txt'
part_data1 = np.loadtxt(fname1)
xpart= part_data1[:,0]
ypart =part_data1[:,1]
zpart =part_data1[:,2]
thetapart=part_data1[:, 3]
rpart2=part_data1[:, 4]
surf_rho_at_r=part_data1[:, 5]
mass=part_data1[:, 6]
p_area_elm=part_data1[:,7]

fname2 = 'Disk_Setup/JO201_24kpc_V1_2.txt'
part_data2 = np.loadtxt(fname2)
rpart1= part_data2[:, 0]
m_enc =part_data2[:, 1]

# Loading pickle file Cluster Infall class and defining parameters.
cluster_data=open("Clusters/Abell_85_corrected_a/Bellhouse_GASPXV_0_72Gyr_2279.p", "rb")
data_CL=pickle.load(cluster_data)
wind_velocity=data_CL[:,4]
radial_infall=data_CL[:,3]

# Properties of the DM halo
p3d.DM_type='Mori'                    #Type of potential Hernquist 'Hernquist' or Mori-Burkert 'Mori'.

# Properties of the DM halo (Mori & Burkert 2000).
p3d.r0DM_pc = 23E3  # 18.5E3 #        # scale radius (in pc).
p3d.r0DM = p3d.r0DM_pc * pc2km        # scale radius (in km).
p3d.rho0DM = ((3.084E-24 * (p3d.r0DM_pc / 1E3) ** (-2 / 3)) * (1E-3 / 1E-15))   # central density (kg/km^3).

# Properties of the DM halo (Hernquist).
p3d.M_DM=0.5E12*MO2kg                 # Mass of DM halo (kg).
p3d.b_DM=18.7E3*pc2km                 # Scaling length (km).

# Properties of Bulge (Hernquist potential).
p3d.M_B = 0.89E10*MO2kg         # Mass of Bulge (kg).
p3d.b_B = 700*pc2km             # Scaling length (km).

# Properties of Stellar Disk (Kuzmin-Plummer Disk).
p3d.M_SD = 4E10*MO2kg #1.15E11*MO2kg  #          # Mass of Stellar Disk (kg).
p3d.a_SD = 4E3*pc2km #3620*pc2km#     #          # scaling lengths (km).
p3d.b_SD=  300*pc2km# 700*pc2km # 700*pc2km

# Particle Parameters.
p3d.npart_in_seg=20              # Total number of segments.
p3d.npart =200 #200 #199         # Total number of (gas cloud) particles in galaxy.

# Ram Pressure Constants.

p3d.v_cwind=1500                 # Constant wind.
p3d.rho_wind=5.9422E-16          # Constant density for ICM wind kg/km^3 (used for constant wind)-

p3d.v_vwind=wind_velocity*-1     # Time varying wind list.
p3d.r_cluster=radial_infall      # The radial infall of galaxy into cluster centre list.
p3d.rho0_cl=2.6E-14#2.6E-14#1.29E-14               # Central density for cluster used for beta profile.
p3d.rc= 86E3*pc2km #  82E3*pc2km #                 # Coefficient for beta profile (km).
p3d.beta=0.62 #0.532#                              # Constant for beta profile.

p3d.dt_rp=1E13                                     # Time increment for evolution of Ram_Pressure() function in Particle class.


p3d.wind_type="varying"                            # ICM wind type either 'constant' or 'varying'.
p3d.cw=1.2                                         # Drag coefficient.

# Calling class into list.
RAMs = []
t=10#2279  #2650 #9465# 6310                          # Duration of RP in terms of dt.
dt=1E13                                            # Time Step seconds
rp_i= 0 #2019 #3155                                # Start of RP (in terms of varaible t).
rp_f= 5 #2279#2650 #5174 #5805                        # End of RP (in terms of varaible t)

data_shape = (p3d.npart, t, 14)                    # Defining the shape of data
data = np.ndarray(data_shape)                      # ndarray is arbitrary array of <shape> dimension: data is a 3D array.

for i in tqdm(range(0, p3d.npart)):                # Iterates over the number of particles; tqdm gives a loading bar.
    RAMs.append(p3d.Particle(x=xpart[i]*pc2km, y=ypart[i]*pc2km, z=0, theta=thetapart[i],surf_rho=surf_rho_at_r[i]*1E-27,mass=mass[i])) # Calling class into list with values from data. Creating a list of class objects one for each particle.
    RAMs[-1].Cvel(t,dt,[rp_i*dt,rp_f*dt])                     # For dt=1E13 seconds: 3155.76 (approx 3155) is 1Gyr, approx 3.155 for 1MYR.
                                                              #  Conversion for RP units: 1 dyne/cm^2 =0.1 pascals =100kg/kms^2.
    data[i] = RAMs[i].history                                 # data array is the list history from every class object in RAMs.


# Creating directory to a text file that records the paramters for each simulation run.
def Saving_data_to_file(parent_directory,galaxy,file_name_1,file_name_2,parameter_list):
    directory= f"{galaxy}"                 # The name of the directory (folder) the simulation is saved to.
    parent_dir = f"{parent_directory}"     # The parent directory (folder) for directory, in this case, called RP_Runs.

    # Path
    path = os.path.join(parent_dir, directory)   #  Joins the path components: os.path.join concatenates various path components with a one directory separation.
    if os.path.exists(path):                                                         # If path exists to directory pickle RAMs and data in the folder.
        pickle.dump(RAMs,open(f"{parent_directory}/{galaxy}/{file_name_1}.p", "wb"))
        pickle.dump(data,open(f"{parent_directory}/{galaxy}/{file_name_2}.p", "wb"))
    else:
        os.makedirs(path)                                                            # If path dosen't exist create the directorys and path. Then pickle RAMs and data.
        pickle.dump(RAMs,open(f"{parent_directory}/{galaxy}/{file_name_1}.p", "wb"))
        pickle.dump(data,open(f"{parent_directory}/{galaxy}/{file_name_2}.p", "wb"))
                                                                                     # These next lines are for writing a text file in the directory location
    if os.path.exists(f'{parent_directory}/{galaxy}/{parameter_list}.txt')==False:   # If  text file doesn't exist then -
        fh=open(f'{parent_directory}/{galaxy}/{parameter_list}.txt', 'w')            # - create file and wite
        fh.write(f'{file_name_1}\nType: {p3d.wind_type}\nWind Vel: {p3d.v_cwind} (km/s)\nRun t: {t/3155} (Gyr)\nWind t: {rp_i/3155} to {rp_f/3155} (Gyr)\ncw: {p3d.cw}\nICM rho: {p3d.rho_wind*1e-12} (g/cm^3)\nCluster rho_0: {p3d.rho0_cl*1e-12} (g/cm^3)\nrc: {p3d.rc/(pc2km*1E3)} (Kpc)\nBeta: {p3d.beta}\nDM Halo: {p3d.DM_type}\nMori: {p3d.r0DM_pc/1E3} (Kpc), {p3d.rho0DM*1e-12} (g/cm^3)\nHernquist: {p3d.M_DM/(MO2kg*1e12)}^12 (M0), {p3d.b_DM/1000} (Kpc)\nBulge: {p3d.M_B/(MO2kg*1e10)}^10 (M0), {p3d.b_B/(1E3*pc2km)}(Kpc)\nStellar disk: {p3d.M_SD/(MO2kg*1e10)}^10 (M0), a;{p3d.a_SD/(1E3*pc2km)} (Kpc), b;{p3d.b_SD/(1E3*pc2km)} (Kpc)')  # paramters to in text file
        fh.close()

    else:
        fh=open(f'{parent_directory}/{galaxy}/{parameter_list}.txt', 'a')  # Text file does exist open file and alter.
        fh.write(f'\n\n{file_name_1}\nType: {p3d.wind_type}\nWind Vel: {p3d.v_cwind} (km/s)\nRun t: {t/3155} (Gyr)\nWind t: {rp_i/3155} to {rp_f/3155} (Gyr)\ncw: {p3d.cw}\nICM rho: {p3d.rho_wind*1e-12} (g/cm^3)\nCluster rho_0: {p3d.rho0_cl*1e-12} (g/cm^3)\nrc: {p3d.rc/(pc2km*1E3)} (Kpc)\nBeta: {p3d.beta}\nDM Halo: {p3d.DM_type}\nMori: {p3d.r0DM_pc/1E3} (Kpc), {p3d.rho0DM*1e-12} (g/cm^3)\nHernquist: {p3d.M_DM/(MO2kg*1e12)}^12 (M0), {p3d.b_DM/1000} (Kpc)\nBulge: {p3d.M_B/(MO2kg*1e10)}^10 (M0), {p3d.b_B/(1E3*pc2km)}(Kpc)\nStellar disk: {p3d.M_SD/(MO2kg*1e10)}^10 (M0), a; {p3d.a_SD/(1E3*pc2km)} (Kpc), b; {p3d.b_SD/(1E3*pc2km)} (Kpc)')
        fh.close()
    return

Save_data=Saving_data_to_file("RP_Runs","D100","RAM_cw1_2", "Data_cw1_2","Parameters") # The arguement are as follows: 1st the parent directory, 2nd-3rd the RAM and data list/array to save,-
                                                                                        # -4th the text file that lists the parameters used to create RAM and Data.
