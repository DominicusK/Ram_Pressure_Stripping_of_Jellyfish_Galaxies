import pickle
import numpy as np
import math as m
import mpmath as mp
import os

# Find previous directory.
curdire = os.path.dirname(__file__)
predire= os.path.dirname(curdire)

# Loading Infall data: velocity and radial postion.
cluster_data=open(f"{predire}/Clusters/Abell_85_corrected_a/Bellhouse_GASPXV_0_72Gyr_2279.p", "rb")
data_CL=pickle.load(cluster_data)
wind_velocity=data_CL[:,4]
radial_infall=data_CL[:,3]

# Constants
GRAV = 6.674E-20               # Gravitational constant (km^3/kg*s^2).
pc2cm = 3.086E18
cm2kpc = 3.2408E-22
kg2M0= 5.0279E-31
pc2km = 3.086E+13
MO2kg=1.989e+30

# Disk Sufrace density constant (Softened exponential).
Mgas_tot=0.6e10*MO2kg          # Mass value used before intergration in Disk_Setup.
a=7500*pc2km                   # Scaling constant.

# DM halo Properties
# Mori & Burkert DM potential.
r0DM_pc = 23E3                 # scale radius (in pc).
r0DM = r0DM_pc * pc2km         # scale radius (in km).
rho0DM = ((3.084E-24 * (r0DM_pc / 1E3) ** (-2 / 3)) * (1E-3 / 1E-15))   # Central density (kg/km^3).

# Hernquist DM Potential.
M_DM=0.5E12*MO2kg              # Mass of DM halo (kg).
b_DM=18.7E3*pc2km              # Scaling length (km).
#print(rho0DM*1e-12)

#Bulge (Hernquist potential).
M_b = 0.89E10*MO2kg            # Mass of Bulge (kg).
b_B = 700*pc2km                # Scaling length (km).

# Properties of Stellar Disk (Kuzmin-Plummer Disk).
M_sd = 4E10*MO2kg              # Mass of Stellar Disk (kg).
a_SD = 4E3*pc2km               # scaling lengths (km).
b_SD=  300*pc2km

# ICM wind paramters.
# Constant wind.
v_wind_c=1500                # Constant wind velocity.
icm_density=5.9422E-16       # Constant ICM density.

# Varying wind.
v_wind_tv=wind_velocity*-1   # Time varying wind list.
r_cl=radial_infall           # The radial infall of galaxy into cluster centre list.
rho_0_cl=2.6E-14             # Central density for cluster used for beta profile.
rc=82E3*pc2km                # Coefficient for beta profile (km).
beta=0.532                   # Constant for beta profile.

cw=0.5                       # Drag coefficient.
wind_type='varying'          # ICM wind type either 'constant' or 'varying'.

class Stripping_Radius:
    def __init__(self,maxr,r_step,z_height):    # Intial Conditions.
        self.maxr=maxr
        self.r_step=r_step

        self.stripping=[]
        self.stripping_tv=[]
        self.pot_pressure=[]

        self.z_height=z_height                  # Height above the galaxtic plane.
        self.z=self.z_height*pc2km

    def Potential_z(self,r):
            self.r=r
            #Disk Surface density (Softended Expontenital).
            self.surf_rho=((Mgas_tot*0.25/(4*(a**2)))*(mp.sech((self.r*pc2km)/a)))

            #3D radius.
            self.radius=m.sqrt((self.r*pc2km)**2+(self.z)**2)

            # The z components of the Gravitational potential accerlations.
            # DM Halo
            self.Rcontrib_DM = -GRAV*M_DM/((self.radius+b_DM)**2)
            self.zcontrib_DM = self.Rcontrib_DM * self.z / self.radius
            # Bulge
            self.Rcontrib_B = -GRAV*M_b/((self.radius+b_B)**2)
            self.zcontrib_B = self.Rcontrib_B * self.z / self.radius
            # Stellar disk (dphi/dz).
            self.zcontrib_SD=-(GRAV*M_sd*self.z*(a_SD+m.sqrt(self.z**2+b_SD**2)))/(m.sqrt(self.z**2+b_SD**2)*((self.r*pc2km)**2+(a_SD+m.sqrt(self.z**2+b_SD**2))**2)**(3/2))

            # Total contribution of z components:
            self.total_zcontrib=self.zcontrib_DM+self.zcontrib_B+self.zcontrib_SD
            self.pot_pressure = abs(self.total_zcontrib*self.surf_rho)  # Absolute value to satisfy the Gunn & Gott criterion.

    def RP(self,icm_density,v_wind_c,v_wind_tv,r_cl,rho_0_cl,rc,beta,wind):  # Determing the RPs and the accerlation due to ICM wind.

        if wind=="constant":
            self.pram_original=icm_density*v_wind_c**2
            self.pram_w=0.5*cw*icm_density*v_wind_c**2
        if wind=="varying":
            self.rho_vw = rho_0_cl * (1 + (r_cl / rc) ** 2) ** (-3 * beta / 2)
            self.pram_original = self.rho_vw * v_wind_tv ** 2
            self.pram_w=0.5*cw*self.rho_vw*v_wind_tv**2

    def GG_criterion(self,RP_start,RP_finish):                               # A function that determines the Gun and Gott criterion.
        self.stripping.clear()
        self.stripping_tv.clear()
        self.RP_duration=RP_finish-RP_start                                  # Duration of the RP.

        if wind_type=="constant":                                            # Criterion for constant case.
            self.RP(icm_density,v_wind_c,v_wind_tv,0,rho_0_cl,rc,beta,wind=wind_type)
            for i in reversed(range(1,self.maxr,self.r_step)):               # Reversing lists to read outer most radii first.
                self.Potential_z(i)
                if self.pram_w  < self.pot_pressure:                         # The criterion itself. If the potentials' 'pressure' > RP-
                    self.stripping.append(self.r)                            # - then append the radii to list.

        if wind_type=="varying":                                             # Criterion for varying case.
            for t in range(0,self.RP_duration):
                self.RP(icm_density,v_wind_c,v_wind_tv[t],r_cl[t],rho_0_cl,rc,beta,wind=wind_type)
                for i in reversed(range(1,self.maxr,self.r_step)):
                    self.Potential_z(i)
                    if self.pram_w  < self.pot_pressure:
                        self.stripping.append(i)                             # All the unstripped particles for index (time step)-
                                                                             # - in v_wind_tv[t] & r_cl[t].
                self.stripping_tv.append(max(self.stripping))                # The maximum value i.e., the stripping radius  for each index of v_wind_tv[t] & r_cl[t].
                self.stripping.clear()                                       # Make sure list continaully appended.

RP_st=1         # Start of RP.
RP_fin=1891     # End of RP.

StrRad=[]
StrRad.append(Stripping_Radius(maxr=24600,r_step=100,z_height=10000))        # In aguemnents in PARSECS! (Only exception in code).
StrRad[-1].GG_criterion(RP_st,RP_fin)

#print(StrRad[0].stripping)                                                  # Check Strppinf radii and correct wind type is used.
#print(StrRad[0].stripping_tv)

R_stripped=StrRad[0].stripping_tv                                            # This list is uncommented if wind type is varying.
#R_stripped=StrRad[0].stripping                                              # This list is uncommented if wind type is constant.
pickle.dump(R_stripped,open(f"cw_0_5_10kpc_JO201_72_1890_va.p", "wb"))       # Saving data with pickle dump.
