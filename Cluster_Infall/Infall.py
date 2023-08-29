import numpy as np
import matplotlib.pyplot as plt
import math as m
import mpmath as mp
import pickle
import os

# Constants-not all are used but useful for plots.
GRAV = 6.674E-20               # Gravitational constant (km^3/kg*s^2)
pc2cm = 3.086E18
cm2kpc = 3.2408E-22
pc2km = 3.086E+13
M02kg= 1.989e+30

# Galaxy Properties
maxr=26000
maxr_Kpc=int(maxr/1E3)
M_G=1.5E12*M02kg
Mgas_tot =1.59E41


# Galaxy Position
x_pos=0
y_pos=0
z_pos=2.4E6 # 1.94E6 #2.4E6   #3.0E6

# Galaxy intail velocity:
Infall_v0=-2000

# Time elements
dt=1E13
t=6310

# Cluster properties (Herquist potential).
m_CL =1.58E15*M02kg #1.23E15*M02kg #1.83E15*M02kg# 4.10E14*M02kg                # Mass of cluster (kg)
b_CL =902E3*pc2km #902E3*pc2km  #1140E3*pc2km #663E3*pc2km  #                 # Scaling length b (km).
radius_CL=1940E3#2400E3                      # Cluster radius

rho0_cl= 2.6E-14
beta=0.532
rc= 82E3*pc2km#82E3*pc2km

# Class System (base units for ENTIRE CODE: Km,s,Kg)
class CLuster_Infall:
    def __init__(self, x=0.0, y=0.0, z=0.0,R_vel=0.0,theta_vel=0.0,phi_vel=0.0):  # Initial parameters.
        self.x = x
        self.y = y
        self.z = z

        self.R_vel = R_vel
        self.theta_vel = theta_vel
        self.phi_vel = phi_vel

        self.mass=mass

        self.icm_zvel= 0
        self.scatter = None


    def Calc_Radius(self):
        self.r= m.sqrt(self.x**2 + self.y**2)                        # x-y plane radius.
        self.radius = m.sqrt(self.r ** 2 + self.z ** 2)              # Spherical coordinates radius.

        return self.radius,self.r,self.x, self.y, self.z


    def Cluster_Contribution(self):                                  # Acceleration due to cluster (km/s^2) (Hernquist Model).
      self.Calc_Radius()

      self.Rcontrib_CL = -GRAV*m_CL/((self.radius+b_CL)**2)

      # component
      self.rcontrib_CL = self.Rcontrib_CL * self.r / self.radius
      self.zcontrib_CL = self.Rcontrib_CL * self.z / self.radius
      self.xcontrib_CL = self.Rcontrib_CL * self.x / self.radius
      self.ycontrib_CL = self.Rcontrib_CL * self.y / self.radius

      return self.r, self.radius, self.Rcontrib_CL, self.rcontrib_CL, self.xcontrib_CL, self.ycontrib_CL,self.zcontrib_CL


    def Total_Acc(self):                                             # The Acceleration components in x,y,z,r from all potentials.
        self.Cluster_Contribution()

        self.xtot_acc=self.xcontrib_CL
        self.ytot_acc=self.ycontrib_CL
        self.ztot_acc=self.zcontrib_CL
        self.rtot_acc=self.rcontrib_CL

        return self.rtot_acc, self.xtot_acc, self.ytot_acc

    def Initial_Velocity_Conversion(self):                                      # Initial Velocity Setup (km/s).
        self.Total_Acc()

        self.xvel = self.R_vel*np.cos(self.theta_vel)*np.sin(self.phi_vel)
        self.yvel = self.R_vel * np.sin(self.theta_vel) * np.sin(self.phi_vel)
        self.zvel = self.R_vel #* np.cos(self.phi_vel)
        return  self.zvel

    def Velocity_Step(self,dt):                                      # Velocity over time step dt (in seconds).
        self.Total_Acc()

        self.xvel =self.xvel+self.xtot_acc*dt
        self.yvel =self.yvel+self.ytot_acc*dt
        self.zvel =self.zvel+self.ztot_acc*dt

        #self.zvelrel=self.zvel+self.icm_zvel

        self.x=self.x+self.xvel*dt
        self.y=self.y+self.yvel*dt
        self.z=self.z+self.zvel*dt

        #self.radius=m.sqrt(self.x**2+self.y**2+self.z)


    def Vel(self, t, dt):                                      # All functions combined to evolve galaxy overtime and apply ram pressure.
        self.Initial_Velocity_Conversion()

        history_shape=(t,6)
        self.history_CL=np.ndarray(history_shape)
        self.history_CL[0]=[self.x, self.y, self.z, self.radius,self.zvel,self.R_vel]    # list contains parameters wanted from class (before time evolution).

        for I in range(0,int(t)*int(dt),int(dt)):
            self.Velocity_Step(dt)
            self.Calc_Radius()

            J=int(I/dt)
            self.history_CL[J]=[self.x, self.y, self.z, self.radius,self.zvel,self.R_vel]      # list contains parameters wanted from class (after time evolution).

INFALL = []
data_CL = []
INFALL.append(CLuster_Infall(x=x_pos*pc2km, y=y_pos*pc2km, z=z_pos*pc2km,R_vel=Infall_v0,theta_vel=0.0,phi_vel=0.0)) # Calling class into list with values from data and converting into code base units.
INFALL[-1].Vel(t,dt)
for i in INFALL:
     data_CL= i.history_CL


def Saving_cl_data_to_file(parent_directory,cluster,file_name_1,parameter_list):
    curdire = os.path.dirname(__file__)
    predire= os.path.dirname(curdire)

    back_dir=f"{predire}"
    directory= f"{cluster}"

    parent_dir = f"{parent_directory}"

    # Path
    path = os.path.join(back_dir, parent_dir, directory)
    if os.path.exists(path):
        pickle.dump(data_CL,open(f"{predire}/{parent_directory}/{cluster}/{file_name_1}.p", "wb"))
    else:
        os.makedirs(path)
        pickle.dump(data_CL,open(f"{predire}/{parent_directory}/{cluster}/{file_name_1}.p", "wb"))

    if os.path.exists(f'{predire}/{parent_directory}/{cluster}/{parameter_list}.txt')==False:
        fh=open(f'{predire}/{parent_directory}/{cluster}/{parameter_list}.txt', 'w')

        fh.write(f'{file_name_1}\nM_200: {m_CL/M02kg} (M0)\nR_200:{radius_CL/1E3} (kpc)\nRho_0: {rho0_cl*1E-12} (g/cm^3)\nR_core (r_c): {rc/(1E3*pc2km)} (kpc)\nBeta: {beta}\nHernquist scaling paramater b: {b_CL/(1E3*pc2km)} (kpc)\nIntial velocity: {Infall_v0} (Km/s)\nIntial Galaxy Postion: {z_pos/1E3} (Kpc)')
        fh.close()

    else:
        fh=open(f'{predire}/{parent_directory}/{cluster}/{parameter_list}.txt', 'a')
        fh.write(f'\n\n{file_name_1}\nM_200: {m_CL/M02kg} (M0)\nR_200:{radius_CL/1E3} (kpc)\nRho_0: {rho0_cl*1E-12} (g/cm^3)\nR_core (r_c): {rc/(1E3*pc2km)} (kpc)\nBeta: {beta}\nHernquist scaling paramater b: {b_CL/(1E3*pc2km)} (kpc)\nIntial velocity: {Infall_v0} (Km/s)\nIntial Galaxy Postion: {z_pos/1E3} (Kpc)')
        fh.close()
    return

# Save_cl_data=Saving_cl_data_to_file("Clusters","Abell_85_corrected_a","Bellhouse_GASPXV_0_91Gyr_2853_HM","Cl_Parameters")

def closest_value(input_list, input_value):      # A function that finds the closest value in a list to the input value.

  arr = np.asarray(input_list)

  i = (np.abs(arr - input_value)).argmin()

  return arr[i]


gal_z=data_CL[:,2]   # List of z positions.
gal_r=data_CL[:,3]   # List of the radii.
gal_v=data_CL[:,4]   # List of the velocity.

gal_z_kpc=[i/(pc2km*1E3) for i in gal_z]   # Z position in kpc.
gal_v_list=[i*-1 for i in gal_v]           # Flips the sign of the negative infall velocities to positive.

z_value=int(360)                           # Desired final postions of galaxy (Kpc).
v_value=int(3364)                          # Desired final velocity of galaxy (km/s).

val_z=closest_value(gal_z_kpc,z_value)     # Finding closest infall position to desired cluster centric distance.
index_z = gal_z_kpc.index(val_z)           # The index of the closest z postion.

val_v=closest_value(gal_v_list,v_value)    # Finding closest velocity to desired final velocity.
index_v = gal_v_list.index(val_v)          # The index for val_v.

vmax=max(gal_v_list)                       # The maxmimum velocity experienced during infall.
index_vmax=gal_v_list.index(vmax)

print('Closest radius to the cluster centre', val_z,index_z)      # Closest infall position to desired cluster centric distance.
print('Closest velocity to desired velocity', val_v,index_v)      # Closest velocity to desired final velocity.
print('Velocity at the radius index', gal_v[index_z]*-1)          # The velocity at cluster centric distance, val_z.
print('Time taken to reach radius index (Gyr)', index_z/3155)     # The Time (Gyr) to reach cluster centric distance, val_z.
print('Maximum velocity galaxy infalling', vmax,index_vmax)       # The maximum velocity experienced by the galaxy.








# print(gal_z[2650])
# print(gal_v[2650])

#print(gal_z[2019])
#print(gal_v[2019])


rho0_cl1= 2.6E-14
beta1=0.62
rc1= 86E3*pc2km #22.1E3*pc2km

rho0_cl2= 2.6E-14
beta2=0.532
rc2=82E3*pc2km#82E3*pc2km

RP=[]
for i in range(0,index_z):
    r_s=(gal_z[i]/rc2)**2
    RP.append((rho0_cl2*(1+r_s)**(-3*beta2/2))*(gal_v[i])**2)
RP1=[]
for i in range(0,index_z):
    r_s1=(gal_z[i]/rc1)**2
    RP1.append((rho0_cl1*(1+r_s1)**(-3*beta1/2))*(gal_v[i])**2)

print(RP[2581])

raise SystemExit(0)
time=[]
for i in range(0,int(index_z)):
     time.append((i)/3155)

fig2, ax = plt.subplots()
ticksize=14
xyfontsize=14
ax.grid(linestyle='dotted')
ax.tick_params(direction='in')
ax.tick_params(top=True, right=True, labelsize=ticksize)
ax.set_xlabel("Distance from Cluster Centre (Mpc)",fontsize = ticksize)
ax.set_ylabel("Density (10$^{-26}$ g/cm$^{3}$)",fontsize = xyfontsize)
ax.set_title('Ram Pressure Profile: Time Varying' )
ax.plot(time,RP,linewidth=1,c='lightcoral',label='Abell 85')
ax.plot(time,RP1,linewidth=1,c='black',label='Abell 85')
#ax.legend(prop={'size': 11})

print(max(RP))




 # 1.0694E-14 #


# u=0.6
# n0=(13.0E-3)*(100000**3)


#print(gal_z[6310],gal_v[6310])

distance_step=100e3*pc2km
distance_step2=100e3*pc2km

distance_cl=[]
distance_cl2=[]
for i in range(0,int(2.4E6*pc2km),int(distance_step)):
    distance_cl.append(i)
for i in range(0,int(2.4E6*pc2km),int(distance_step2)):
    distance_cl2.append(i)
rho_vw=[]
rho_vw2=[]

for i in distance_cl:
   rho_vw.append((rho0_cl1*(1+(i/rc1)**2)**(-3*beta1/2))*1E-12*1E26)
for i in distance_cl2:
   rho_vw2.append((rho0_cl2 * (1 + (i / rc2) ** 2) ** (-3 * beta2 / 2))*1E-12*1E26)
#   rho_vw3.append((rho0_cl * (1 + (i / rc3) ** 2) ** (-3 * beta / 2))*1E-12)
#   rho_vw4.append((rho0_cl * (1 + (i / rc4) ** 2) ** (-3 * beta / 2)) * 1E-12)

distance_clMyr=[]
for i in distance_cl:
    distance_clMyr.append(i/(1e6*pc2km))

distance_clMyr2=[]
for i in distance_cl2:
    distance_clMyr2.append(i/(1e6*pc2km))
# pram=[]
# for i in range(0,len(gal_z)):
#     pram.append((rho0_cl*(1+((gal_z[i]*pc2km*1e6)/rc)**2)**(-3*beta/2))) #*gal_v[i]**2

time=range(0,len(gal_z))
fig2, ax = plt.subplots()
ticksize=14
xyfontsize=14
ax.grid(linestyle='dotted')
ax.tick_params(direction='in')
ax.tick_params(top=True, right=True, labelsize=ticksize)
ax.set_xlabel("Distance from Cluster Centre (Mpc)",fontsize = ticksize)
ax.set_ylabel("Density (10$^{-26}$ g/cm$^{3}$)",fontsize = xyfontsize)
ax.set_title('Ram Pressure Profile: Time Varying' )
ax.plot(distance_clMyr,rho_vw,linewidth=2,c='lightcoral',label='Abell 85')
ax.plot(distance_clMyr2,rho_vw2,linewidth=2,linestyle='--',c='black', label='Mock Cluster')
ax.legend(prop={'size': 11})
plt.show()

raise SystemExit(0)
#pram=rho_vw[2650]*gal_v[2650]**2
#print(pram)

raise SystemExit(0)





# time=list(range(t))
#
#
timegyr=[]
for i in range(0,t):
    timegyr.append(i/3155)


fig2, ax = plt.subplots()

ax.grid(linestyle='dotted')
ax.tick_params(direction='in')
ax.tick_params(top=True, right=True)
ax.set_xlabel("Distance (Mpc)")
ax.set_ylabel(r"$\rho$$_{Cluster}$ (10$^{-13}$ Kgkm$^{-3}$)")
ax.set_title('Cluster ICM Density against Distance from Cluster Centre' )
ax.plot(distance_clMyr,rho_vw,c='r',label="r$_{Core}$=22kpc")
#ax.plot(distance_clMyr,rho_vw2,c='b',label="r$_{Core}$=442kpc")
#ax.plot(distance_clMyr,rho_vw3,c='g',label="r$_{Core}$=79kpc")
#ax.plot(distance_clMyr,rho_vw4,c='purple',label="r$_{Core}$=135kpc")
ax.set_yscale("log")
ax.legend()
plt.show()

fig3, ax1 = plt.subplots()

ax1.grid(linestyle='dotted')
ax1.tick_params(direction='in')
ax1.tick_params(top=True, right=True)
ax1.set_xlabel("Distance in Cluster (Mpc)")
ax1.set_ylabel("Galaxy Velocity Km s$^{-1}$")
ax1.set_title("Position of Galaxy in Cluster against Galaxy Velocity")
ax1.plot(gal_z,gal_v,c='r')
ax1.scatter(gal_z[2650],gal_v[2650],s=18, label="Velocity$_{Gal}$ at 2 Gyr")
ax1.legend()

fig4, ax2= plt.subplots()

ax2.grid(linestyle='dotted')
ax2.tick_params(direction='in')
ax2.tick_params(top=True, right=True)
ax2.set_xlabel("Time (Gyr)")
ax2.set_ylabel("Distance (Mpc)")
ax2.set_title("Position of Galaxy in Cluster against Infall Time")
ax2.plot(timegyr,gal_z,c='r')
ax2.scatter(2650/3155,gal_z[2650],s=18,label="Postion$_{Gal}$ at 2 Gyr")
ax2.legend()


plt.show()
