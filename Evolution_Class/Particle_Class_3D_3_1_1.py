import numpy as np
import matplotlib.pyplot as plt
import math as m
import mpmath as mp
import pickle

# Constants-not all are used but useful for plots.
GRAV = 6.674E-20               # Gravitational constant (km^3/kg*s^2)
pc2cm = 3.086E18
cm2kpc = 3.2408E-22
pc2km = 3.086E+13
kg2M0= 5.0279E-31
MO2kg=1.989e+30

# Galaxy Properties.
# Radius.
maxr=26000
maxr_Kpc=int(maxr/1E3)
Mgas_tot =1.59E41

# Properties of the DM halo (Mori & Burkert 2000).
r0DM_pc = 23E3                # scale radius (in pc).
r0DM = r0DM_pc * pc2km        # scale radius (in km).
rho0DM = ((3.084E-24 * (r0DM_pc / 1E3) ** (-2 / 3)) * (1E-3 / 1E-15))   # central density (kg/km^3).

# Properties of the DM halo (Hernquist).
M_DM=1E12*MO2kg
b_DM=18.7E3*pc2km

# Properties of Bulge (Herquist potential  potential).
M_b = 1E40                    # Mass of Bulge (kg) (equivalent to 2.0E10 Solar Masses).
b_B = 600*pc2km               # scaling length b   (in km).

# Properties of Stellar Disk (Kuzmin-Plummer Disk).
M_sd = 1.15E11*MO2kg              # Mass of Stellar Disk in Km (equivalent to 2.0E11 Solar Masses).
a_SD = 3500*pc2km             # scaling length a km.
b_SD= 700*pc2km

# Gas disk/cloud properties.
rs=8500                       # scaling factors for surface density
c=1.47

# Particle Parameters.
npart_in_seg=20              # Total number of segments.
npart = 200                 # Total number of particles in galaxy (gas).

# Ram Pressure Constants.
c_rpv=900
v_cw=1500
v_wind= 1000                # For vw code this can also be a list
rho_wind=5.9422E-16
dt_rp=1E13
cw=0.5                     # Drag coefficient.
rc= 22.10546894964873E3*pc2km              # Coefficient for beta profile , for vw code
r_cluster=1                 # radius of cluster infall, for vw code
rho0_cl=1.069597918E-14       # rho0 (kg/km^3) if dont have then use n0 and u
                                # Central number density for cluster (beta profile)
                       # molecular weight used to convert to density instead of number density for beat profile
beta=0.666                   # beta constant for beat profile
wind_type="constant"


# Class System (base units for ENTIRE CODE: Km,s,Kg i.e a density will be Kg/Km^2)
class Particle:
    def __init__(self, x=0.0, y=0.0, z=0.0,theta=0.0,surf_rho=0.0, mass=0.0):  # Initial parameters.
        self.x = x
        self.y = y
        self.z = z
        self.theta = theta
        self.surf_rho=surf_rho
        self.mass=mass

        self.r0=m.sqrt(self.x**2+self.y**2)

        self.ffglimit = 10E3 * pc2km
        self.height = m.sqrt(self.r0 ** 2 + self.z**2)

        self.zram_vel=0
        self.pram_tv=0
        self.pram_tvrel=0
        self.pram_orig=0
        self.pram_c=0
        self.rho_vw=0

        self.ramhistory = []

        self.scatter = None

    def Calc_Radius(self):
        self.r= m.sqrt(self.x**2 + self.y**2)                        # x-y plane radius.
        self.radius =m.sqrt(self.r**2 + self.z**2)              # Spherical coordinates radius.

        return self.radius,self.r, self.x, self.y

    def Darkmatter_Contribtuion(self):                               # Acceleration due to Dark Matter Halo (km/s^2) (Mori & Burkert 2000).
 #       self.Calc_Radius()

        self.Rcontrib_DM = -GRAV*M_DM/((self.radius+b_DM)**2)

#        self.contrib_DMtest=-2*((-r0DM/self.radius**2)* m.atan(self.radius / r0DM)+(1+(r0DM/self.radius))*(r0DM/(r0DM**2+self.radius**2)))+2*((-r0DM/self.radius**2)* np.log(1+(self.radius / r0DM))+(1+(r0DM/self.radius))*(1/(r0DM+self.radius)))-((r0DM/self.radius**2)* np.log(1+(self.radius**2 / r0DM**2))+(1-(r0DM/self.radius))*(2*self.radius/(r0DM**2+self.radius**2)))
#        self.Rcontrib_DMtest=self.contrib_DMtest*m.pi * GRAV * rho0DM * r0DM ** 2
#        self.zcontrib_DMtest=self.Rcontrib_DMtest*(self.z / self.radius)

        # Components
        self.rcontrib_DM = self.Rcontrib_DM * self.r / self.radius
        self.zcontrib_DM = self.Rcontrib_DM * self.z / self.radius
        self.xcontrib_DM = self.Rcontrib_DM * self.x / self.radius
        self.ycontrib_DM = self.Rcontrib_DM * self.y / self.radius


        return self.r, self.radius, self.rcontrib_DM,  self.Rcontrib_DM, self.xcontrib_DM, self.ycontrib_DM,self.zcontrib_DM

    def Bulge_Contribution(self):                                    # Acceleration due to Bulge (km/s^2) (Plummer Model).
 #     self.Calc_Radius()

      self.Rcontrib_B = -GRAV*M_b/((self.radius+b_B)**2)
          # -(GRAV * M_b * self.radius / ((self.radius ** 2 + b_B ** 2) ** (3 / 2)))

      # component
      self.rcontrib_B = self.Rcontrib_B * self.r / self.radius
      self.zcontrib_B = self.Rcontrib_B * self.z / self.radius
      self.xcontrib_B = self.Rcontrib_B * self.x / self.radius
      self.ycontrib_B = self.Rcontrib_B * self.y / self.radius


      self.Rcontrib_Bffg = -GRAV * M_b / ((self.height + b_B) ** 2)
      self.zcontrib_Bffg = self.Rcontrib_B * self.ffglimit / self.height
      return self.r, self.radius, self.Rcontrib_B, self.rcontrib_B, self.xcontrib_B, self.ycontrib_B,self.zcontrib_B

    def Stellar_Disk(self):                                          # Acceleration due to Stellar Disk (km/s^2) (Kuzmin Disk).
#       self. Calc_Radius()

       # self.rcontrib_SD = -GRAV *M_sd * self.r/ ((self.r ** 2 +(abs(self.z)+a_SD) ** 2) ** (3 / 2))
       # self.zcontrib_SD = (-GRAV *M_sd * self.r/ ((self.r ** 2 +(abs(self.z)+a_SD) ** 2) ** (3 / 2)))* self.z/self.radius
       # self.xcontrib_SD = self.rcontrib_SD*self.x / self.r
       # self.ycontrib_SD = self.rcontrib_SD*self.y / self.r

       self.rcontrib_SD= -GRAV *M_sd * self.r/ ((self.r**2+(a_SD+m.sqrt((self.z)**2+b_SD**2))**2)**(3/2)) # may need to change self.r to self.radius
      # self.zcontrib_SD= (self.z*(a_SD+m.sqrt(self.z**2+b_SD**2)))/(m.sqrt(self.z**2+b_SD**2)*(self.radius**2+(a_SD+m.sqrt(self.z**2+b_SD**2))**2))
       self.zcontrib_SD=-(GRAV*M_sd*self.z*(a_SD+m.sqrt(self.z**2+b_SD**2)))/(m.sqrt(self.z**2+b_SD**2)*(self.r**2+(a_SD+m.sqrt(self.z**2+b_SD**2))**2)**(3/2))
       self.xcontrib_SD = self.rcontrib_SD * self.x / self.r
       self.ycontrib_SD = self.rcontrib_SD * self.y / self.r


       self.zcontrib_SDffg = (self.ffglimit * (a_SD + m.sqrt(self.ffglimit ** 2 + b_SD ** 2))) / (m.sqrt(self.ffglimit ** 2 + b_SD ** 2) * (self.height ** 2 + (a_SD + m.sqrt(self.ffglimit ** 2 + b_SD ** 2)) ** 2))
       return self.rcontrib_SD, self.xcontrib_SD, self.ycontrib_SD,self.zcontrib_SD

    def Total_Acc(self):                                             # The acceleration components in x,y,z,r from all potentials.
        self.Calc_Radius()
        self.Bulge_Contribution()
        self.Darkmatter_Contribtuion()
        self.Stellar_Disk()

        self.xtot_acc=self.xcontrib_DM + self.xcontrib_B + self.xcontrib_SD
        self.ytot_acc=self.ycontrib_DM + self.ycontrib_B + self.ycontrib_SD
        self.ztot_acc=self.zcontrib_DM + self.zcontrib_B + self.zcontrib_SD
        self.rtot_acc=self.rcontrib_DM + self.rcontrib_B + self.rcontrib_SD


        self.f_gg = 0
                # The force from the gravitational potentials as stated by Gott & Gun (2D).

        return self.rtot_acc, self.xtot_acc, self.ytot_acc, self.zcontrib_DM,self.zcontrib_B,

    def Initial_Velocity(self):                                      # Initial Velocity Setup (km/s).
        self.Total_Acc()

        # Individual circular velocities calculated from mass profile (Individual needed to find rotational curve).
        # self.Menc_DM =m.pi * rho0DM * (r0DM ** 3) * (( -2 * m.atan(self.r / r0DM) + 2 * m.log(1 + (self.r / r0DM)) + m.log(1 + (self.r / r0DM) ** 2)))
        # self.Vcirc_DM =m.sqrt(GRAV * self.Menc_DM / (self.r))

        self.Menc_DM =M_DM * (self.r**2 / ((self.r + b_DM)**2))
        self.Vcirc_DM =m.sqrt(GRAV * self.Menc_DM / (self.r))

        self.Menc_BP =M_b * (self.r**2 / ((self.r + b_B)**2))
        self.Vcirc_B =m.sqrt(GRAV * self.Menc_BP / (self.r))


        self.Menc_SD=M_sd-(M_sd*a_SD/m.sqrt(self.r**2+a_SD**2))                                         # potentially wrong
        self.Vcirc_SDM =m.sqrt(GRAV *M_sd * self.r**2/ ((self.r**2+(a_SD+m.sqrt(b_SD**2))**2)**(3/2))) # potentially wrong

        # Circular velocities given to galaxy at t=0, with x-y components calculated.
        self.Vcirc_SD =m.sqrt((GRAV *M_sd * self.r**2/ ((self.r**2+(a_SD+b_SD)**2)**(3/2))))
        self.Vcirc_SP =GRAV * (self.Menc_DM + self.Menc_BP) / self.r
        self.Vcirc_tot =m.sqrt(self.Vcirc_SP + self.Vcirc_SD**2)


        self.xvel =-1.0 * (self.Vcirc_tot) * m.sin(self.theta)
        self.yvel =(self.Vcirc_tot) * m.cos(self.theta)
        self.zvel= 0

        self.rplot= 0                                                # Parameter for stability tests.

    def Velocity_Step(self,dt):                                      # Velocity over time step dt (in seconds).
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

        self.rplot =((self.r-self.r0)/self.r0)*100                  # Parameter for stability tests.


    def Cvel(self, t, dt, ram):                                      # All functions combined to evolve galaxy overtime and apply ram pressure.
        self.ramhistory.clear()
        self.Initial_Velocity()
        history_shape=(t,14)
        self.history=np.ndarray(history_shape)
        self.history[0]=[self.x, self.y, self.r, self.radius, self.Vcirc_DM, self.Vcirc_B, self.Vcirc_SD, self.z, self.rplot,self.zvel, self.f_gg,self.xcontrib_DM,self.ztot_acc,self.Vcirc_tot] # list contains parameters wanted from class (before time evolution).

        for I in range(0,int(t)*int(dt),int(dt)):
            self.Velocity_Step(dt)
 #           self.Calc_Radius()      # added this term later
            if ram[0] < I < ram[1]:
                h_rp=(ram[1]-ram[0])/2
#               F = ((ram[1] - ram[0]) / (dt))
                I_infall=int(I/dt)
                I_develop=int((I-ram[0])/dt)
                if wind_type=="constant":
                    self.Ram_Pressure(I-ram[0],h_rp, c_rpv, v_cw, rho_wind, r_cluster, dt_rp, wind=wind_type)   # ICM density (Kg/Km^3) (conversion from g/cm^3 is 1E12).
                if wind_type=="varying":
                    self.Ram_Pressure(I - ram[0], h_rp, c_rpv, v_wind[I_develop], rho_wind, r_cluster[I_develop], dt_rp,wind=wind_type)
                if wind_type=="test":
                    self.Ram_Pressure(I-ram[0],h_rp, c_rpv, v_cw, rho_wind, r_cluster, dt_rp, wind=wind_type)
            # if int(F) > 0:
 #                   ramhistory_shape = (int(F), 1)
 #                   self.ramhistory = np.ndarray(ramhistory_shape)
 #                   self.ramhistory[int(F-1)] = [self.pram_max]
            J=int(I/dt)
            self.history[J]=[self.x, self.y, self.r, self.radius, self.Vcirc_DM, self.Vcirc_B, self.Vcirc_SD, self.z, self.rplot,self.zvel,self.f_gg,self.xcontrib_DM,self.ztot_acc,self.Vcirc_tot] # list contains parameters wanted from class (after time evolution).

    def Ram_Pressure(self,t_rp,t_c,c_p,v_w,rho_w,r_cl,dt,wind):            # Ram Pressure function (from galaxies reference frame).

            self.pram_max=0.5*cw*rho_w*v_w**2                          # RP units kg/kms^2
            if wind=="varying":
                self.rho_vw = rho0_cl * (1 + (r_cl / rc) ** 2) ** (-3 * beta / 2)

                self.pram_orig = self.rho_vw * v_w ** 2
                self.pram_tv=0.5*cw*self.rho_vw*v_w**2
                self.pramzvelrel=0.5 * cw *self.rho_vw * (self.zvel) ** 2
#                self.pram_tvrel=(self.pram_tv - 0.5 * cw *self.rho_vw * (self.zvel) ** 2)
                self.pram_tvrel = (0.5 * cw * self.rho_vw * (v_w-self.zvel) ** 2)
#                self.zram_acc = (self.pram_tv - 0.5 * cw * self.rho_vw * (self.zvel) ** 2) / self.surf_rho
                self.zram_acc = (0.5 * cw * self.rho_vw * (v_w-self.zvel) ** 2) / self.surf_rho
#                if self.pramzvelrel > self.pram_tv:
                if self.zvel >=v_w:
                    self.zram_acc =0
                    self.pram_tvrel=0
                self.zram_vel = self.zram_acc * dt

            if wind== "Volmer_varying":
             # Time varying RP - Volmer 2001.
                self.pram_tv=self.pram_max*c_p/(((t_rp-t_c)/1E13)**2+c_p)
                self.zram_acc =(self.pram_tv - 0.5*cw*rho_w*self.zvel ** 2) / self.surf_rho
                self.zram_vel = self.zram_acc * dt

            if wind== "constant":
            # Constant RP Equation.
                self.pram_orig=rho_w*v_w**2
                self.pram_c = 0.5 * cw * rho_w * (v_w-self.zvel) ** 2
                self.zram_max = (0.5 * cw * rho_w * v_w ** 2) / self.surf_rho
                self.zram_acc=(0.5*cw*rho_w*(v_w-self.zvel)**2)/self.surf_rho
                #self.zram_vel = self.zram_acc * dt
                if self.zvel>=v_w:
                    self.zram_acc =0
                    self.pram_c = 0
                else:
                    self.zram_acc=(0.5*cw*rho_w*(v_w-self.zvel)**2)/self.surf_rho
#                self.f_ram=0.5*cw*rho_w*(v_w**2-self.zvel**2)         # Force due to Constant RP.
                self.zram_vel = self.zram_acc * dt
            if wind== 'test':
                self.z=2E3*pc2km
                #self.zram_vel=200
                #self.zvel=self.zram_vel
            # Analytic equation to compare Constant RP profile to.
            # c=0.5*cw*rho_w/self.surf_rho
            # self.zveltest=v_w-1/(c*t_rp+(1/v_w))
            self.ramhistory.append(self.pram_orig)
