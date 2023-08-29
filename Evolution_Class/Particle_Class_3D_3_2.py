import numpy as np
import matplotlib.pyplot as plt
import math as m
import mpmath as mp
import pickle

# Constants-not all are used but useful for quick tests in class file.
GRAV = 6.674E-20               # Gravitational constant (km^3/kg*s^2)
pc2cm = 3.086E18
cm2kpc = 3.2408E-22
pc2km = 3.086E+13
kg2M0= 5.0279E-31
MO2kg=1.989e+30

# Properties of the DM halo
DM_type='Hernquist' #Type of potential Hernquist 'Hernquist' or Mori-Burkert 'Mori'

# DM halo (Mori & Burkert 2000).
r0DM_pc = 23E3                # scale radius (in pc).
r0DM = r0DM_pc * pc2km        # scale radius (in km).
rho0DM = ((3.084E-24 * (r0DM_pc / 1E3) ** (-2 / 3)) * (1E-3 / 1E-15))   # central density (kg/km^3).

# DM halo (Hernquist).
M_DM=1E12*MO2kg
b_DM=18.7E3*pc2km

# Properties of Bulge (Herquist potential  potential).
M_B = 1E40                    # Mass of Bulge (kg)
b_B = 600*pc2km               # scaling length (km).

# Properties of Stellar Disk (Kuzmin-Plummer Disk).
M_SD = 1.15E11*MO2kg          # Mass of Stellar Disk (kg)
a_SD = 3500*pc2km             # scaling lengths (km).
b_SD= 700*pc2km

# Particle Parameters.
npart_in_seg=20              # Total number of segments.
npart = 200                  # Total number of particles in galaxy (gas).

# Ram Pressure Constants.
v_cwind=1500                 # Constant wind.
rho_wind=5.9422E-16          # Constant density for ICM wind kg/km^3 (used for constant wind)

v_vwind= 1000                # Time varying wind variable is a list in Time_Evolution_8_3D.py.
r_cluster=1                  # The radial infall of galaxy into cluster centre. In Time_Evolution_8_3D.py the input is a list.
rho0_cl=1.069597918E-14      # Central density for cluster used for beta profile.
rc= 22.10546894964873E3*pc2km    # Coefficient for beta profile (km).
beta=0.666                       # Constant for beta profile.

cw=0.5                       # Drag coefficient.

wind_type="constant"         # ICM wind type either 'constant' or 'varying'

dt_rp=1E13                   # Time increment for evolution of Ram_Pressure() function. Usally the same as time dt in Velocity_Step() & Cvel()

# Class System (base units for ENTIRE CODE: Km,s,Kg i.e a surface density will be Kg/Km^2)
class Particle:
    def __init__(self, x=0.0, y=0.0, z=0.0,theta=0.0,surf_rho=0.0, mass=0.0):  # Initial parameters.
        self.x = x
        self.y = y
        self.z = z
        self.theta = theta
        self.surf_rho=surf_rho
        self.mass=mass

        self.r0=m.sqrt(self.x**2+self.y**2)                                    # Initial radius (x-y plane) of the particle used for stability tests.

        self.pram_orig=0
        self.pram_max=0
        self.pram_tv=0
        self.zram_vel=0
        self.rho_vw=0

        self.ramhistory = []

        self.scatter = None

    def Calc_Radius(self):
        self.r= m.sqrt(self.x**2 + self.y**2)                   # x-y plane radius.
        self.radius =m.sqrt(self.r**2 + self.z**2)              # Spherical coordinates radius.

        return self.radius,self.r, self.x, self.y

    def Darkmatter_Contribtuion(self):                          # Acceleration due to Dark Matter Halo (km/s^2).

        if DM_type=='Mori':
            self.contrib_DM = -2. * (-(r0DM / (self.radius ** 2)) * m.atan(self.radius / r0DM) + (1. + (r0DM / self.radius)) * (r0DM / (r0DM ** 2 + self.radius ** 2))) + 2. * (((-r0DM / (self.radius ** 2)) * m.log(1. + (self.radius / r0DM))) + (1 + (r0DM / self.radius)) / (r0DM + self.radius)) - (((r0DM / self.radius ** 2) * m.log(1. + ((self.radius ** 2) / r0DM ** 2))) + (1. - (r0DM / self.radius)) * 2. * self.radius / (r0DM ** 2 + self.radius ** 2))
            self.Rcontrib_DM = self.contrib_DM * m.pi * GRAV * rho0DM * r0DM ** 2

        if DM_type=='Hernquist':
            self.Rcontrib_DM = -GRAV*M_DM/((self.radius+b_DM)**2)

        # Components
        self.rcontrib_DM = self.Rcontrib_DM * self.r / self.radius
        self.zcontrib_DM = self.Rcontrib_DM * self.z / self.radius
        self.xcontrib_DM = self.Rcontrib_DM * self.x / self.radius
        self.ycontrib_DM = self.Rcontrib_DM * self.y / self.radius

        return self.r, self.radius, self.rcontrib_DM,  self.Rcontrib_DM, self.xcontrib_DM, self.ycontrib_DM,self.zcontrib_DM

    def Bulge_Contribution(self):                              # Acceleration due to Bulge (km/s^2).

       self.Rcontrib_B = -GRAV*M_B/((self.radius+b_B)**2)

       # Components.
       self.rcontrib_B = self.Rcontrib_B * self.r / self.radius
       self.zcontrib_B = self.Rcontrib_B * self.z / self.radius
       self.xcontrib_B = self.Rcontrib_B * self.x / self.radius
       self.ycontrib_B = self.Rcontrib_B * self.y / self.radius

       return self.r, self.radius, self.Rcontrib_B, self.rcontrib_B, self.xcontrib_B, self.ycontrib_B,self.zcontrib_B

    def Stellar_Disk(self):                                    # Acceleration due to Stellar Disk (km/s^2).

       self.rcontrib_SD= -GRAV *M_SD * self.r/ ((self.r**2+(a_SD+m.sqrt((self.z)**2+b_SD**2))**2)**(3/2))
       self.zcontrib_SD=-(GRAV*M_SD*self.z*(a_SD+m.sqrt(self.z**2+b_SD**2)))/(m.sqrt(self.z**2+b_SD**2)*(self.r**2+(a_SD+m.sqrt(self.z**2+b_SD**2))**2)**(3/2))

       # Components.
       self.xcontrib_SD = self.rcontrib_SD * self.x / self.r
       self.ycontrib_SD = self.rcontrib_SD * self.y / self.r

       return self.rcontrib_SD, self.xcontrib_SD, self.ycontrib_SD,self.zcontrib_SD

    def Total_Acc(self):                                  # The total acceleration of all the components x,y,z,r from the potentials.
        self.Calc_Radius()
        self.Bulge_Contribution()
        self.Darkmatter_Contribtuion()
        self.Stellar_Disk()

        self.xtot_acc=self.xcontrib_DM + self.xcontrib_B + self.xcontrib_SD
        self.ytot_acc=self.ycontrib_DM + self.ycontrib_B + self.ycontrib_SD
        self.ztot_acc=self.zcontrib_DM + self.zcontrib_B + self.zcontrib_SD
        self.rtot_acc=self.rcontrib_DM + self.rcontrib_B + self.rcontrib_SD


        return self.rtot_acc, self.xtot_acc, self.ytot_acc, self.zcontrib_DM,self.zcontrib_B

    def Initial_Velocity(self):                           # Initial Velocity Setup (km/s).
        self.Total_Acc()
        # Individual circular velocities calculated from mass profile (Individual needed to find rotational curve).

        # Circular velocity of DM Halo: Mori-Burkert.
        if DM_type=='Mori':
            self.Menc_DM =m.pi * rho0DM * (r0DM ** 3) * (( -2 * m.atan(self.r / r0DM) + 2 * m.log(1 + (self.r / r0DM)) + m.log(1 + (self.r / r0DM) ** 2)))
            self.Vcirc_DM =m.sqrt(GRAV * self.Menc_DM / (self.r))
        # Circular velocity of DM Halo :Hernquist.
        if DM_type=='Hernquist':
            self.Menc_DM =M_DM * (self.r**2 / ((self.r + b_DM)**2))
            self.Vcirc_DM =m.sqrt(GRAV * self.Menc_DM / (self.r))

        # Circular velocity of Bulge: Hernquist.
        self.Menc_BP =M_B * (self.r**2 / ((self.r + b_B)**2))
        self.Vcirc_B =m.sqrt(GRAV * self.Menc_BP / (self.r))

        # Circular velocity of Stellar disk: Plummer-Kuzmin.
        self.Vcirc_SD =m.sqrt((GRAV *M_SD * self.r**2/ ((self.r**2+(a_SD+b_SD)**2)**(3/2))))

        # Circular velocities given to galaxy at t=0, with x-y components calculated.
        self.Vcirc_sph =GRAV * (self.Menc_DM + self.Menc_BP) / self.r
        self.Vcirc_tot =m.sqrt(self.Vcirc_sph + self.Vcirc_SD**2)


        self.xvel =-1.0 * (self.Vcirc_tot) * m.sin(self.theta)
        self.yvel =(self.Vcirc_tot) * m.cos(self.theta)
        self.zvel= 0

        self.rplot= 0                                      # Parameter for stability test of orbit.

    def Velocity_Step(self,dt):                            # Velocity of particles over time step dt (in seconds).
        self.Total_Acc()
        self.xvel =self.xvel+self.xtot_acc*dt
        self.yvel =self.yvel+self.ytot_acc*dt
        self.zvel =self.zvel+ self.zram_vel +self.ztot_acc*dt

        self.x=self.x+self.xvel*dt
        self.y=self.y+self.yvel*dt
        self.z = self.z + self.zvel * dt

        self.r=m.sqrt(self.x**2+self.y**2)
        self.radius=m.sqrt(self.x**2+self.y**2+self.z**2)
        self.theta =m.atan(self.y/self.x)

        self.rplot =((self.r-self.r0)/self.r0)*100         # Parameter for stability test of orbit: The radial varation in orbit of each particle as a %.


    def Cvel(self, t, dt, ram):                            # All functions combined to evolve galaxy overtime and apply RP.
        self.ramhistory.clear()
        self.Initial_Velocity()
        history_shape=(t,14)
        self.history=np.ndarray(history_shape) # list contains parameters wanted from class (before time evolution).
        self.history[0]=[self.x, self.y, self.r, self.radius, self.Vcirc_DM, self.Vcirc_B, self.Vcirc_SD, self.z, self.rplot,self.zvel, self.Rcontrib_DM,self.xcontrib_DM,self.ztot_acc,self.Vcirc_tot]

        for I in range(0,int(t)*int(dt),int(dt)):  # Evolving the galaxy in time increments dt.
            self.Velocity_Step(dt)                 # Calling Velocity_Step().

            if ram[0] < I < ram[1]:                # This parameter (list) defines time the RP is 'turned on' and 'turned off'.

                I_develop=int((I-ram[0])/dt)       # Index for lists tracking galaxys cluster postions and velocity, when using time varing RP.
                I_infall=int(I/dt)

                if wind_type=="constant":
                    self.Ram_Pressure(I-ram[0], v_cwind, rho_wind, r_cluster, dt_rp, wind=wind_type)   # ICM density (Kg/Km^3) (conversion from g/cm^3 is 1E12).
                if wind_type=="varying":
                    self.Ram_Pressure(I-ram[0], v_vwind[I_develop], rho_wind, r_cluster[I_develop], dt_rp,wind=wind_type)
                if wind_type=="test":
                    self.Ram_Pressure(I-ram[0], v_cwind, rho_wind, r_cluster, dt_rp, wind=wind_type)

            J=int(I/dt)   # index to update list self.history
            self.history[J]=[self.x, self.y, self.r, self.radius, self.Vcirc_DM, self.Vcirc_B, self.Vcirc_SD, self.z, self.rplot,self.zvel,self.Rcontrib_DM,self.xcontrib_DM,self.ztot_acc,self.Vcirc_tot] # list contains parameters wanted from class (after time evolution).

    def Ram_Pressure(self,t_rp,v_w,rho_w,r_cl,dt,wind):            # Ram Pressure function (from galaxies reference frame). RP units kg/kms^2.

            if wind== "constant":                                  # Constant RP.
            # Constant RP Equation.
                self.pram_orig=rho_w*v_w**2                        # The RP as defined by Gott and Gunn in the constant case.
                self.pram_max=0.5*cw*rho_w*v_w**2                  # The maximum RP the galaxy can experinece (derived from drag force equation). Constant RP.

                self.zram_acc=(0.5*cw*rho_w*(v_w-self.zvel)**2)/self.surf_rho # The relative acceleration due to the RP.

                if self.zvel>=v_w:                                            # This 'if' & 'else' statement is to ensure the that the limit on the relative velocity is 0.
                    self.zram_acc =0
                else:
                    self.zram_acc=(0.5*cw*rho_w*(v_w-self.zvel)**2)/self.surf_rho
#                self.f_ram=0.5*cw*rho_w*(v_w**2-self.zvel**2)                # Force due to Constant RP.
                self.zram_vel = self.zram_acc * dt                            # The Velocity due to RP.

            if wind=="varying":                                               # Varying RP.
                self.rho_vw = rho0_cl * (1 + (r_cl / rc) ** 2) ** (-3 * beta / 2)   # ICM density function: beta profile.

                self.pram_orig = self.rho_vw * v_w ** 2                       # The RP as defined by Gott and Gunn in the varying case.
                self.pram_tv=0.5*cw*self.rho_vw*v_w**2                        # The maximum RP the galaxy can experinece (derived from drag force equation). Varying RP.

                self.zram_acc = (0.5 * cw * self.rho_vw * (v_w-self.zvel) ** 2) / self.surf_rho  # The relative acceleration due to the RP.

                if self.zvel >=v_w:                                           # This 'if' statement is to ensure the that the limit on the relative velocity is 0.
                    self.zram_acc =0

                self.zram_vel = self.zram_acc * dt

            if wind== 'test':                                                 # This is 'wind type' is to ensure the class is behaving as intended and a test to see if the constant RP obeys analytic equations.
                self.z=200
                #self.zram_vel=200
                #self.zvel=self.zram_vel

                # Analytic equation to compare constant RP profile.
                # c=0.5*cw*rho_w/self.surf_rho
                # self.zveltest=v_w-1/(c*t_rp+(1/v_w))
            self.ramhistory.append(self.pram_orig)
