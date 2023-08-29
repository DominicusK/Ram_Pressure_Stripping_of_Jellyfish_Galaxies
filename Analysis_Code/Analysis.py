import numpy as np
import matplotlib.pyplot as plt
import math as m
import mpmath as mp
import pickle

# Functions to manipulate class data.

# Mass retained in the disk within a given z height.
def mass_disk(RAMs, time, limit):                     # RAMs - list of Particle class objects which contains parameter list history.
    massleft = 0                                      # time - the gas mass at a certain time step, dt.
    for particle in RAMs:                             # limit - z axis height limit to measure the mass contained within the disk, km.
        z = particle.history[time][7]
        if abs(z) < limit:
            massleft += particle.mass
    return massleft

# The culmaltive mass along the galaxy's tail.
def mass_tails(RAMs,time,boundary, limit):            # RAMs - list of Parricle class objects which contains parameter list history.
    masstails=0                                       # time - the mass in the tails at a certian time step, dt.
    for particle in RAMs:                             # boundary - the z height at which the tail starts, km.
        z=particle.history[time][7]                   # limit - z axis upper limit to measure tail length, km.
        if boundary < z < limit:
            masstails +=particle.mass
    return  masstails

# Returns the z postion and z velocity at a given time: useful for 2D histogram plots.
def velocity_distribution(RAMs,time,boundary):        # RAMs - list of Particle class objects which contains parameter list history.
    z_position=[]                                     # time - the z velocity and z position at a certian time step, dt.
    z_velocity=[]                                     # limit - z axis upper limit to particle (tail) measurements, km.
    for particle in RAMs:
        z=particle.history[time][7]
        if abs(z)<boundary:
            z_position.append(particle.history[time][7])
            z_velocity.append(particle.history[time][9])
    return z_position, z_velocity

# Returns the maxr radius (i.e the stripping radius) at a given height above the disk and time.
def Model_Stripping_Radius(RAMs, time, limit):        # RAMs - list of Particle class objects which contains parameter list history.
    r=[]                                              # time - records the radius, at a given height, at time step, dt.
    rstrip=0                                          # limit - The radii contained with a z axis limit.
    for particle in RAMs:
        z=particle.history[time][7]
        if z < limit:
            r.append(particle.history[time][2])
            rstrip=max(r)/pc2km
            r.clear()
    return rstrip

# Mean average.
def Average(lst):
    return sum(lst) / len(lst)
