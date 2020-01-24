''' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * File name  :  RocketFlightCharacteristics.java
 * Purpose    :  Holds utilities useful for calculating rocket performance
 * @author    :  Harrison Leece, John, Brian, Alejandro
 * Date       :  2020-01-23
 * Notes      :  None
 * Warnings   :  None
 *  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Revision History
 * ================
 *   Ver      Date     Modified by:  Reason for change or modification
 *  -----  ----------  ------------  ---------------------------------------------------------------------
 *
 *  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ '''

import math

import numpy as np
import matplotlib.pyplot as plt

#compThrust and compThrust2 from thrustCurve


def pressurant_thrust(self, init_thrust, init_m_dot, rocket_mass):
    thrust = init_thrust
    mDot = init_m_dot
    if rocket_mass < .1:
        thrust = 0
        mDot = 0
    return thrust, mDot

def comp_mach(self, rocket_velocity, temp):
    #convert v from ft/s to m/s
    v = rocket_velocity * .3048
    #R is the universal gas constant in kJ/(kg*K)
    R = .287
    #Specific heat ratio (gamma) is approx 1.4 for whole domain of temperatures
    gamma = 1.4
    #speed of sound = a, meters/s
    a = (gamma * R * 1000 * temp)**.5
    mach = v / a
    return mach

'''
Curve fits below this line.
@Author Harrison Leece harryleecemail@gmail.com
'''

#compDensity from atmosphereDensity
def compute_atm_density(alt):
    #altitude is in feet
    s = -.2457039080 * (alt / 10**4)**1 - .0351850842 * (alt / 10**4)**2 + 0.0021044026 * (alt / 10**4)**3 - 0.0000390562 * (alt / 10**4)**4
    return 23.77 * 10**(-4) * np.exp(s)
def compute_atm_density_test():
    altitude = np.linspace(0,400000,10000)
    density_list = []
    for alt in altitude:
        density = compute_atm_density(alt)
        density_list.append(density)
    plt.plot(altitude,density_list)
    plt.show()


def comp_temp(self, altitude):
    #convert altiude from feet to meters
    x = altitude * .3048
    if(x > 0 and x < 15354.016):
        C = 2.5574 * 10**(-7) * x**2 - 8.8754 * 10**(-3) * x + 18.1061
    elif(x > 15354.016 and x < 49697.695):
        C = 3.6838 * 10**(-8) * x**2 -7.6707 * 10**(-4) * x - 54.7841
    elif(x > 49697.695 and x < 120000):
        C = -2.4347 * 10**(-3) * x + 119.078
    else:
        C = -172
    temp = C + 273.15
    return temp

def wave_drag(self, mach, cd):
    adj = cd - .3
    if(mach < .5085):
        ncd = .6827 * mach**3 - .4297 * mach**2 - .0358 * mach + .3 + adj
    elif(mach > .5085 and mach < 1.3618):
        ncd = -0.7809 * mach**4 + 2.324 * mach**3 - 2.0189 * mach**2 + 0.6793 * mach + 0.1837 + adj
    elif(mach > 1.3618 and mach < 4):
        ncd = -.003495 * mach**6 + .07004 * mach**5 -.583 * mach**4 + 2.564 * mach**3 -6.186 * self.mach**2 + 7.466 * self.mach -2.923 + adj
    else:
        ncd = .2184 + adj
    return ncd

#F_wr is wind in r direction
#import drag, gravity, time, mass, thrust, gamma, alpha, and wind in all directions
class SIXDOF:
    def __init__(self, x_1, x_2, x_3, x_4, x_5, x_6, fw, alt, x, step_size, nthrust):
        self.state = [x_1, x_2, x_3, x_4, x_5, x_6]
        self.drag = 0
        self.gravity = 32.174 #ft / s^2
        self.time = 0
        self.thrust = 0
        self.gamma = 88
        self.alphap = 2
        self.stateZ = [0, 0, 0, 0, 0, 0]
        self.fw = fw
        self.alt = alt
        self.x = x
        self.step_size = step_size
        self.nthrust = nthrust

    def sixeqdiff(self):
        x_1 = self.state[0]
        x_2 = self.state[1]
        x_3 = self.state[2]
        x_4 = self.state[3]
        x_5 = self.state[4]
        x_6 = self.state[5]

        F_wr = self.fw
        F_wtheta = 0
        F_wphi = 0

        x_dot_1 = x_2
        x_dot_2 = -self.gravity + (-self.drag * sin(self.gamma) + (self.thrust * sin(self.gamma + self.alphap)) + F_wr) / mass + x_1 * x_4**2 * sin(x_5)**2 + x_1 * x_5**2
        x_dot_3 = x_4
        x_dot_4 = ((-self.drag * cos(self.gamma) + self.thrust * cos(self.gamma + self.alphap) + F_wtheta) / m - 2 * x_2 * x_4 * sin(x_5) - 2 * x_1 * x_2 * x_6 * cos(x_5)) * (1 / (x_1 * sin(x_5)))
        x_dot_5 = x_6
        x_dot_6 = ((F_wphi / mass) - 2 * x_2 * x_6 + x_1 * x_4**2 * cos(x_5) * sin(x_5)) * (1 / x_1)
        return  x_dot_1, x_dot_2, x_dot_3, x_dot_4, x_dot_5, x_dot_6

    def windForce(self):
        #determines wind velocity value for certain altitude (alt)
        altitude = self.alt
        if altitude <= 1 and altitude >= 0:
            windVel = 5 * (altitude - 1) + 8
        elif altitude <= 2 and altitude > 1:
            windVel = 8
        elif altitude <= 3 and altitude > 2:
            windVel = -6 * altitude + 20
        elif  altitude <= 11 and altitude > 3:
            windVel = ((43 * altitude) - 113.004) / 8
        elif altitude <= 15 and altitude > 11:
            windVel = 45
        elif altitude < 22 and altitude > 15:
            windVel = -6 * altitude + 135
        elif altitude <= 25 and altitude >= 22:
            windVel = 4
        elif altitude <= 30 and altitude > 25:
            windVel = 2 * altitude - 46
        else:
            windVel = 0

        return self.alt, windVel


    def gusts(self):
        num_gusts = random.randrange(0, 20)
        alt_range = random.randrange(0, 30000, 1)
        mag_gust = random.randrange(0, 12, 1)
        direction = random.choice((-1, 1))

        #creates an array of values and an empty list for all gust altitudes
        gust_matrix = np.array([0, alt_range, mag_gust, direction])
        gust_altitudes = []

        for i in range(num_gusts-1):
            alt_range = random.randrange(0, 30000, 1)
            mag_gust = random.randrange(0, 12, 1)
            direction = random.choice((-1, 1))
            temp_array = np.array([i+1, alt_range, mag_gust, direction])
            gust_matrix = np.vstack((temp_array, gust_matrix))
            gust_altitudes.append(alt_range / 1000)
        gust_matrix = np.flipud(gust_matrix)

        return gust_matrix

    def sixdofsolver(self):
        init_state = [x_1, x_2, x_3, x_4, x_5, x_6]
        sixeq = sixeqdiff(self, fw)
        state = odeint(sixeq, init_state)
        return state

if __name__ == '__main__':
    compute_atm_density_test()
