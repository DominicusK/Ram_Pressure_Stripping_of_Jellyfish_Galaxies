import numpy as np
import matplotlib.pyplot as plt
import math as m
import mpmath as mp
from matplotlib import gridspec
from mpl_toolkits.axes_grid1 import make_axes_locatable

GRAV = 6.674E-20               # Gravitational constant (km^3/kg*s^2)
pc2cm = 3.086E18
cm2kpc = 3.2408E-22
pc2km = 3.086E+13
kg2M0= 5.0279E-31
MO2kg=1.989e+30

# Creates 2 plots showing the evolution of the galaxy in the x-y plane with no ram pressue.
def Galaxy_Evolution_NoRP(data,t_0,t_1,x,y):                                     # data - is the 3D array created when calling Particle in Time_Evolution_8_3D.
                                                                                 # t_0, t_1 - the times that the snapshots of the galaxy ,in the x-y plane, are taken.
    plt.figure(figsize=(8,8))                                                    # x - the index in data that correspondeds to the x coordinate.
                                                                                 # y - the index in data that correspondeds to the y coordinate.
    gs1 = gridspec.GridSpec(2, 1)
    ax1 = plt.subplot(gs1[0],)
    ax1.grid(linestyle='dotted')
    ax1.tick_params(direction='in')
    ax1.tick_params(top=True, right=True)
    ax1.set_ylabel("Y-coordinate (Kpc)")


    ax2 = plt.subplot(gs1[1],sharex=ax1)
    ax2.grid(linestyle='dotted')
    ax2.tick_params(direction='in')
    ax2.tick_params(top=True, right=True)
    ax2.set_xlabel('X-coordinate (Kpc)')
    ax2.set_ylabel('Y-coordinate (Kpc)')


    ax1.set_box_aspect(1)
    ax2.set_box_aspect(1)

    ax1.scatter(data[:,t_0,x] / (pc2km * 1E3), data[:,t_0,y] / (pc2km * 1E3), s=10)
    ax2.scatter(data[:,t_1,x] / (pc2km * 1E3), data[:,t_1,y] / (pc2km * 1E3), s=10)

    plt.subplots_adjust(hspace=.0)
    return ax1,ax2

# Plotting An Individuals Particles Orbit in the x-y plabe, no RP.
def Particle_Orbit(data,particle,title_time,title_radius):                       # data - the 3D array created when calling Particle in Time_Evolution_8_3D.
    fig2, ax = plt.subplots()                                                    # particle - the index of the desired particle from data.
                                                                                 # title_time - the time in Gyr of the particel orbit for the title,
    ax.grid(linestyle='dotted')                                                  # title_radius - the radius of the particle for the title.
    ax.tick_params(direction='in')
    ax.tick_params(top=True, right=True)
    ax.set_xlabel('X-coordinate (Kpc)')
    ax.set_ylabel('Y-coordinate (Kpc)')
    ax.set_title(f'Cloud Orbit Over {title_time} Gyr: Radius {title_radius} Kpc')


    ax.scatter(data[particle,:,0] /(pc2km*1E3), data[particle,:,1] /(pc2km*1E3), s=10)
    plt.gca().set_aspect('equal', adjustable='box')
    return ax

# PLotting radial and z stabilities with respect to time.
def RZ_Stability_Tests(data,t, t_scale,particle,r_st, z):                        # data - the 3D array created when calling Particle in Time_Evolution_8_3D.
    timescale = list(range(t))                                                   # t - the time range plotting the z & r stability, units dt.
    num = t_scale                                                                # t_scale - determines the units of the time scale in plot i.e, /3155 gives Gyr.
    newtime = [i / num for i in timescale]                                       # particle  - the particle for which the orbital stabilities will be plotted plots.
                                                                                 # r_st - the index of radial stability from data.
                                                                                 # z - the index of z position from data.
    plt.figure()
    gs = gridspec.GridSpec(2, 1)

    ax = plt.subplot(gs[0])
    ax.plot(newtime, data[particle,:t,r_st], c='green')
    ax.tick_params(top=True, right=True)
    ax.grid(linestyle='dotted')
    ax.tick_params(direction='in')
    ax.set_ylabel("Difference in radial orbit (%)", labelpad=5)

    ax1 = plt.subplot(gs[1], sharex=ax)
    ax1.plot(newtime, data[particle,:t,z]/(pc2km * 1E3), c='deepskyblue')
    ax1.grid(linestyle='dotted')
    ax1.tick_params(direction='in')
    ax1.tick_params(top=True, right=True)
    ax1.set_xlabel('Time (Gyr)')
    ax1.set_ylabel('Z coordinate (Kpc)')
    plt.subplots_adjust(hspace=.0)
    return ax, ax1


# Quick Overview plot of the galaxy undergoing RP.
def RP_Quick_Overview(data,x,y,z):
    fig, ((ax, ax1), (ax2, ax3)) = plt.subplots(2, 2)                            # data - the 3D array created when calling Particle in Time_Evolution_8_3D.
                                                                                 # x - the index of x position from data.
    ax.set_box_aspect(1)                                                         # y - the index of y position from data.
    ax1.set_box_aspect(1)                                                        # z - the index of z position from data.
    ax2.set_box_aspect(1)
    ax3.set_box_aspect(1)

    ax.scatter(data[:,0,x], data[:,0,y], s=10)
    ax1.scatter(data[:,-1,x], data[:,-1,y], s=10)
    ax2.scatter(data[:,0,x], data[:,0,z])
    ax3.scatter(data[:,-1,x], data[:,-1,z])
    return ax, ax1, ax2, ax3


# Neat RP plot at a specific time during the simulation.
def RP_Plot(data,t,x,y,z):                                                       # data - the 3D array created when calling Particle in Time_Evolution_8_3D.
    plt.figure(figsize=(8,8))                                                    # t - the desired tine for a snapshot fo the galaxy, units of dt.
                                                                                 # x - the index of x position from data.
    gs1 = gridspec.GridSpec(2, 1)                                                # y - the index of y position from data.
    ax = plt.subplot(gs1[0])                                                     # z - the index of z position from data.
    ax.grid(linestyle='dotted')
    ax.tick_params(direction='in')
    ax.tick_params(top=True, right=True)
    ax.set_ylabel("Z-coordinate (Kpc)")
    #ax.set_xticks(np.arange(-25, 25, 5.0))

    #ax.set_ylim(((-6000,-5000)))

    ax1 = plt.subplot(gs1[1])
    ax1.grid(linestyle='dotted')
    ax1.tick_params(direction='in')
    ax1.tick_params(top=True, right=True)
    ax1.set_xlabel('X-coordinate (Kpc)')
    ax1.set_ylabel('Y-coordinate (Kpc)')


    ax.set_box_aspect(1)
    ax1.set_box_aspect(1)

    for P in range(0, len(data)):
        ax.scatter(data[P][t][x]/(pc2km*1E3), data[P][t][z]/(pc2km*1E3), s=10)
        ax1.scatter(data[P][t][x]/(pc2km*1E3), data[P][t][y]/(pc2km*1E3), s=10)

    #plt.subplots_adjust(hspace=.0)
    return ax, ax1


# Plotting the gas mass within in the disk
def Mass_Disk_Plot(t,scale,mass_disk,ticksize,xyfontsize):                       # t - the time range for plotting the gas mass with the disk, units dt.
                                                                                 # scale - determines the units of the time scale in plot i.e, /3155 gives Gyr.
    timescale=list(range(t))                                                     # mass_disk - mass within disk over the time range. Extracted using the function mass_tail.
    num=scale                                                                    # ticksize - tick size for plot.
    newtime= [i/num for i in timescale]                                          # xyfontsize - font size for plot.


    fig2, ax3 = plt.subplots()
    ax3.set(ylim=(0, 8.5))
    ax3.grid(linestyle='dotted')
    ax3.tick_params(direction='in')
    ax3.tick_params(top=True, right=True)
    ax3.set_ylabel('Gas Mass within ±10kpc (10\N{SUPERSCRIPT ONE}\N{SUPERSCRIPT ZERO} M$_{☉}$) ', fontsize = xyfontsize)
    ax3.set_xlabel('Time (Gyr)',fontsize = xyfontsize)
    ax3.set_title('Cumulative Gas Mass in Disk (Constant Ram Pressure)')
    ax3.tick_params(axis='both', which='major', labelsize=ticksize)

    ax3.plot(newtime, mass_disk , c='orange', label='Particle Model: c$_{w}$=0.5')
#    ax3.plot(newtime, mass_disk2, c='blue',label='Particle Model: c$_{w}$=0.94')
#    ax3.plot(newtime, mass_disk3, c='crimson', label='Particle Model: c$_{w}$=1.5')
#    ax3.plot(t_com, mass_com,c='black',label='S.Tonnesen 2019: NCFOVW')
    ax3.legend()
#    plt.gca().set_aspect('equal', adjustable='box')

    return ax3

# A basic 2D histogram plotting the z displacement above the disk wrt to the z velocity (at a given time).
def Velocity_Dist_Hist(z,z_v,ticksize,xyfontsize):                               # z-  list of z postions of each particle.
    fig4, ax5 = plt.subplots()                                                   # z_v - list of z velocities of each particle.
                                                                                 # ticksize - tick size for plot.
    ax5.grid(linestyle='dotted')                                                 # xyfontsize - font size for plot.
    ax5.tick_params(direction='in')
    ax5.tick_params(top=True, right=True)
    ax5.set_ylabel('Z Velocity (Km/s)',fontsize = xyfontsize)
    ax5.set_xlabel('Z height above disk (kpc)',fontsize = xyfontsize)

#    ax5.set_title('Gas Mass Distribution as a Function of Z Coordinate & Wind direction')
    ax5.tick_params(axis='both', which='major', labelsize=ticksize)
    h1=ax5.hist2d(z,z_v, bins=(30, 30), cmap=plt.cm.Blues, label='1.5Gyr')
    #
    divider1= make_axes_locatable(ax5)
    cax1= divider1.append_axes('right',size='6%', pad=0.5)
    fig4.colorbar(h1[3],cax=cax1,ax=ax5,orientation='vertical',)
    ax5.set_xlim([0, 210])
    ax5.set_ylim([-5, 1250])
    #plt.gca().set_aspect('equal', adjustable='box')

# # Animation Plot of galaxy in the x-z or y-z plane
# time_dis=3155                                                                  # time starting animation, in units of dt.
# x_y=0                                                                          # index in 3D array data for the x or y coordinate.
# z=7                                                                            # index in 3D array data for the z coordinate.
# r=2                                                                            # index in 3D array data for the radial coordinate (sets a unique colour particles at a given radius).
#
# fig, ax = plt.subplots()
# ax.set(ylim=(-10, 70))
# ax.grid(linestyle='dotted')
# ax.tick_params(direction='in')
# ax.tick_params(top=True, right=True)
# ax.set_xlabel('X-coordinate (Kpc)')
# ax.set_ylabel('Z-coordinate (Kpc)')
# ax.set_title('Particle Model: Ram Pressure Stripping')
# scatter=ax.scatter(data[:,0,x_y]/(pc2km*1E3), data[:,0,z]/(pc2km*1E3), c=data[:,0,r],cmap='rainbow', s=10)
#
# def frame(i):
#         data_x_y=np.array(data[:,time_dis+i,x_y]/(pc2km*1E3))
#         data_z =np.array(data[:,time_dis+i,z]/(pc2km*1E3))
#         scatter.set_offsets(np.c_[data_x_y, data_z])
#         return scatter
#
#
# anim = FuncAnimation(fig, func=frame,frames =4000, interval=10)
# plt.show(
