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

	def three_dof_sim(self):
		#initiate matricies
		a_matrix = np.zeros(4)
		v_matrix = np.zeros(4)
		d_matrix = np.zeros(4)
		d_matrix = np.vstack((d_matrix, d_matrix))

		print('Change!')

		while (d_matrix[-1,-2] > -.01 ):
			'''
			1) Get current atomspheric state
			'''

			'''
			2) Get current rocket state at time, altitude and etc...
			'''

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

def integrate_eom(a_list, step_size):
	v_x = a_list[0]*step_size
	v_y = a_list[1]*step_size
	v_z = a_list[2]*step_size
	v_t = a_list[3]*step_size

	v_list = np.array([v_x, v_y, v_z, v_t])

	d_x = v_list[0]*step_size
	d_y = v_list[1]*step_size
	d_z = v_list[2]*step_size
	d_t = v_list[3]*step_size

	d_list = np.array([d_x, d_y, d_z, d_t])

	return v_list, d_list

if __name__ == '__main__':
	thedude = Rocket(3000, 230, 730, 300, np.pi*(11.2/12)^2, .25)
	thedude.three_dof_sim()
	print('Running the main!')
