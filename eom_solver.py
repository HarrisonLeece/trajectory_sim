# eom_solver

# tangential acceleration*mass = thrust(~) + drag(~) - gravity(~)*mass
#
# vectors

import numpy as np

class Rocket():
	def __init__(self, thrust, i_sp, total_mass, dry_mass, front_area, cd):
		#Initialize rocket parameters at sea level, in lists
		self.thrust = [thrust]
		self.i_sp = [i_sp]
		self.total_mass = [total_mass]
		self.cd = [cd]
		#initialize rocket parameters which do not change
		self.dry_mass = dry_mass
		self.front_area = front_area





def tangent_eom(thrust, theta, omega, mass, density, area, velocity, c_drag):
	a_t = (thrust/mass + (1/2)*c_drag*density*(velocity**2)) / mass - gravity

	a_x = a_t*np.cos(theta)*np.cos(omega)
	a_y = a_t*np.cos(theta)*np.sin(omega)
	a_z = a_t*np.sin(theta)

	a_list = np.array([a_x, a_y, a_z, a_t])

	return a_list

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
	#initiate matricies
	a_matrix = np.zeros(4)
	v_matrix = np.zeros(4)
	d_matrix = np.zeros(4)
	d_matrix = np.vstack((d_matrix, d_matrix))
	thrust = [3000]
	density = [14]
	print('Change!')



	while (d_matrix[-1,-2] > -.01 ):
		d = [-1,-1,-1,-1]
		d_matrix = np.vstack((d_matrix, d))
		print(d_matrix)
