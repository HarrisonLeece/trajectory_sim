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

def eom(thrust, theta, omega, mass, density, area, velocity, c_drag):
	'''check this EoM for correct multiplcation to gravity'''
	a_prime = thrust/mass - (1/2)*c_drag*area*density*(v_mag**2)/mass
	a_t = a_prime - GRAV*np.sin(theta)
	a_x = a_prime*np.cos(theta)*np.cos(omega)
	a_y = a_prime*np.cos(theta)*np.sin(omega)
	a_z = a_prime*np.sin(theta) - GRAV

	a_array = np.array([a_x, a_y, a_z, a_t])

	return a_array

def eom_vector(thrust, mass, density, area, velocity, c_drag):
	v_mag = np.linalg.norm(velocity[0:3])
	v_unit = velocity[0:3]/v_mag
	grav_vec = np.array((0,0,32.174))
	a = thrust/mass * v_unit - v_unit*(1/2)*c_drag*area*density*(v_mag**2)/mass - grav_vec
	a_t = np.linalg.norm(a)
	a = np.hstack((a, a_t))
	return a


def eom_tester():
	thrust, mass, density, area, velocity, c_drag = 3000, 22 , .0025, .065, np.array([30,40,1000,0]), .25
	a = eom_vector(thrust, mass, density, area, velocity, c_drag)
	print(a)


def integrate_eom(a_array, step_size, d_matrix, v_matrix):
	v_x = v_matrix[-1,0] + a_array[0]*step_size
	v_y = v_matrix[-1,1] + a_array[1]*step_size
	v_z = v_matrix[-1,2] + a_array[2]*step_size
	v_t = v_matrix[-1,3] + a_array[3]*step_size

	v_array = np.array([v_x, v_y, v_z, v_t])

	d_x = d_matrix[-1,0] + v_matrix[-1,0]*step_size + .5*a_array[0]*step_size**2
	d_y = d_matrix[-1,1] + v_matrix[-1,1]*step_size + .5*a_array[1]*step_size**2
	d_z = d_matrix[-1,2] + v_matrix[-1,2]*step_size + .5*a_array[2]*step_size**2
	d_t = d_matrix[-1,3] + v_matrix[-1,3]*step_size + .5*a_array[3]*step_size**2

	d_array = np.array([d_x, d_y, d_z, d_t])

	return v_array, d_array

def integrate_eom_vector(a_array, step_size, d_matrix, v_matrix):
	v_x = v_matrix[-1,0] + a_array[0]*step_size
	v_y = v_matrix[-1,1] + a_array[1]*step_size
	v_z = v_matrix[-1,2] + a_array[2]*step_size
	v_array = np.array([v_x, v_y, v_z,v_t])
	v_t = np.linalg.norm(v_array)
	v_array = np.hstack((v_array, v_t))

	d_x = d_matrix[-1,0] + v_matrix[-1,0]*step_size + .5*a_array[0]*step_size**2
	d_y = d_matrix[-1,1] + v_matrix[-1,1]*step_size + .5*a_array[1]*step_size**2
	d_z = d_matrix[-1,2] + v_matrix[-1,2]*step_size + .5*a_array[2]*step_size**2
	d_array = np.array([d_x, d_y, d_z, d_t])
	d_t = np.linalg.norm(d_array)
	d_array = np.hstack((d_array, d_t))

	return v_array, d_array


def comp_mach(rocket_velocity, temp):
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

def calc_thrust(mdot, optimal_thrust, nozzle_exit_pressure, area, alt, propellantmass):
	'''
	Note:  test must be written for this function.
	Hold m_dot, exit pressure, area, nozzle_exit_velocity constant
	provide an initial sea level thrust, and plot thrust vs altitude
	'''
	atm_press = calt_atm_press(alt)
	F = mdot * nozzle_exit_velocity + (nozzle_exit_pressure - press1) * area
	if propellantmass < .1:
		thrust = 0
		mDot = 0
	else:
		thrust = F
		#Redundant,but explicit
		mdot = mdot
	return thrust, mDot

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


def comp_temp(altitude):
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

def compute_atm_temp_test():
	altitude = np.linspace(1,400000,10000)
	temp_list = []
	for alt in altitude:
		temp = comp_temp(alt)
		temp_list.append(temp)
	plt.plot(altitude,temp_list)
	plt.show()

def wave_drag(mach, cd):
	adj = cd - .3
	if(mach < .5085):
		ncd = .6827 * mach**3 - .4297 * mach**2 - .0358 * mach + .3 + adj
	elif(mach > .5085 and mach < 1.3618):
		ncd = -0.7809 * mach**4 + 2.324 * mach**3 - 2.0189 * mach**2 + 0.6793 * mach + 0.1837 + adj
	elif(mach > 1.3618 and mach < 4):
		ncd = -.003495 * mach**6 + .07004 * mach**5 -.583 * mach**4 + 2.564 * mach**3 -6.186 * mach**2 + 7.466 * mach -2.923 + adj
	else:
		ncd = .2184 + adj
	return ncd

def compute_atm_wave_drag_test():
	# set coeff of drag to .3 and vary mach between 0 and 6
	cd = 0.3
	mach = np.linspace(0,6,100)
	ncd_list = []
	for i in mach:
		ncd = wave_drag(i, cd)
		ncd_list.append(ncd)
	plt.plot(mach,ncd_list)
	plt.show()


def calc_atm_press(alt):
	'''@author Hannah Lane
	@param alt, The altitude the rocket is at
	@return pressure, the absolute pressure of the atmosphere at the
			rocket's altitude in psia
	'''
	pressure = 14.7*np.exp((-0.366713637542122*10**-4)*alt+(-0.000001623765497*10**(-4)*alt**2+(0.000000000007584*10**(-4)*alt**3)))
	return pressure

def atm_press_test():
	'''@author Hannah Lane
	Testing block for the pressure function
	'''
	alts = np.linspace(0,150000, 10000)
	pressure_vector = np.array([])
	for alt in alts:
		pressure = calc_atm_press(alt)
		pressure_vector = np.append(pressure_vector, pressure)


	plt.plot(alts, pressure_vector)
	plt.xlabel("Pressure")
	plt.ylabel("Altitude")
	plt.show()

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
	#compute_atm_density_test()
	#compute_atm_temp_test()
	#compute_atm_wave_drag_test()
	#tm_press_test()
	eom_tester()
