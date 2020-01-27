# eom_solver

# tangential acceleration*mass = thrust(~) + drag(~) - gravity(~)*mass
#
# vectors

import numpy as np
from .utilites import compute_atm_density as comp_dense


class Rocket():

	def __init__(self, sl_thrust, i_sp, total_mass, dry_mass, front_area, cd):
		#Initialize rocket parameters at sea level, in lists
		self.thrust = [sl_thrust]
		self.i_sp = [i_sp]
		self.total_mass = [total_mass]
		self.cd = [cd]
		#initialize rocket parameters which do not change
		self.dry_mass = dry_mass
		self.front_area = front_area
		self.m_dot = 1
		self.
		#initialize some non-rocket data, importantfor analysis
		self.time = []

	def three_dof_sim(self):
		#initiate matricies
		a_matrix = np.zeros(4)
		v_matrix = np.zeros(4)
		d_matrix = np.zeros(4)
		d_matrix = np.vstack((d_matrix, d_matrix))

		print('Change!')

		while (d_matrix[-1,2] > -.01 ):
			'''
			1) Get current atomspheric state
			'''
			#atmospheric density is a function of altitude only (d_z)
			atm_rho = compute_atm_density(d_matrix[-1,2])
			#atmospheric temperature is a function of altitude only (d_z)
			atm_t = compute_temp(d_matrix[-1,2])
			#Here would be the code to implement the gust model, magnitude
			#depending on current altitude
			'''
			2) Get current rocket state at time, altitude and etc...
			'''
			#Calculate coefficeient of drag, accounting for wave drag

			#Calculate thrust, accounting for nozzle expansion ratio and
			#aerostatic back pressure loss

			'''
			3) Run Equations of Motion, using gathered state information
			'''

			'''
			4) Organize/Parse data.  Update instance varaibles
			'''

			#Retart the loop

		print('Simulation done!')


def tangent_eom(thrust, theta, omega, mass, density, area, velocity, c_drag):
	'''check this EoM for correct multiplcation to gravity'''
	a_t = (thrust/mass + (1/2)*c_drag*density*(velocity**2)) / mass - gravity*np.sin(theta)

	a_x = a_t*np.cos(theta)*np.cos(omega)
	a_y = a_t*np.cos(theta)*np.sin(omega)
	a_z = a_t*np.sin(theta)

	a_array = np.array([a_x, a_y, a_z, a_t])

	return a_array

def integrate_eom(a_list, step_size, d_matrix, v_matrix):
	v_x = v_matrix[-1,0] + a_list[0]*step_size
	v_y = v_matrix[-1,1] + a_list[1]*step_size
	v_z = v_matrix[-1,2] + a_list[2]*step_size
	v_t = v_matrix[-1,3] + a_list[3]*step_size

	v_list = np.array([v_x, v_y, v_z, v_t])

	d_x = d_matrix[-1,0] + v_list[0]*step_size + .5*a_list[0]*step_size**2
	d_y = d_matrix[-1,1] + v_list[1]*step_size + .5*a_list[1]*step_size**2
	d_z = d_matrix[-1,2] + v_list[2]*step_size + .5*a_list[2]*step_size**2
	d_t = d_matrix[-1,3] + v_list[3]*step_size + .5*a_list[3]*step_size**2

	d_list = np.array([d_x, d_y, d_z, d_t])

	return v_list, d_list

if __name__ == '__main__':
	thedude = Rocket(3000, 230, 730, 300, np.pi*(11.2/12)^2, .25)
	thedude.three_dof_sim()
	print('Running the main!')
