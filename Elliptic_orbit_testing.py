# -*- coding: utf-8 -
"""
Created on Sat Jun 24 22:17:30 2023

@author: Nate Whitsett
"""

import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
import math as math
from scipy import special as sp
from scipy import integrate as it


class Planet:
    def __init__(self, radius, a, e, B, stellar_mass=1, arg_periastron=0, orbit_resolution=0.1,inclination=90):
        self.period = np.sqrt(a**3*3.1536e7)
        self.radius = radius*7.149e9
        self.e = e
        self.a = a
        self.B = B
        self.phase = np.arange(0,2*np.pi, orbit_resolution)
        self.orbit = a*(1-e**2)/(1+e*np.cos(self.phase))
        self.orbit_resolution = orbit_resolution
        self.arg_periastron = arg_periastron*0.0174533
        if arg_periastron > 90:
            self.true_anomaly = np.pi*2 + np.pi/2 - arg_periastron*0.0174533
        elif arg_periastron <= 90:
            self.true_anomaly = np.pi/2 - arg_periastron*0.0174533
        self.inclination = inclination*0.0174533
        self.periastron = None
        self.orbit2 = a*(1-e**2)/(1+e*np.cos(self.phase+self.true_anomaly+np.pi))
        circum = self.a*4*sp.ellipe(self.e)
        time_steps = []
        for distances in self.orbit2:
            time_steps.append(circum/np.sqrt(1.989e30*stellar_mass*6.67e-11*(2/(distances*1.496e+11)-(1/(self.a*1.496e+11)))))
        self.time = np.cumsum(time_steps)
        self.time = 2*self.time/(max(self.time)) - 1
        
            
        
        

planet = Planet(1, 1, .7, 100, arg_periastron=90, orbit_resolution = 0.0002)
plt.plot(planet.time, planet.orbit2)
print(planet.period)