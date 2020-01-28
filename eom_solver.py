# eom_solver

# tangential acceleration*mass = thrust(~) + drag(~) - gravity(~)*mass
#
# vectors

import numpy as np
from utilities import compute_atm_density as comp_dens
from utilities import comp_temp
from utilities import calc_thrust
from utilities import comp_mach
from utilities import wave_drag
import matplotlib.pyplot as plt
#The below allows matplotlibt to plot 3d axis
from mpl_toolkits import mplot3d


GRAV = 32.174

class Rocket():
	#acceleration due to gravity in ft/s^2


	def __init__(self, sl_thrust, optimal_thrust, nozzle_exit_p, nozzle_area, i_sp_sl, total_weight, dry_weight, front_area, cd):
		#Initialize rocket parameters at sea level, in lists
		self.thrust = [sl_thrust]
		self.i_sp = [i_sp_sl]
		self.total_mass = [total_weight/GRAV]
		self.cd = [cd]
		#initialize rocket parameters which do not change
		self.dry_mass = dry_weight/GRAV
		self.front_area = front_area
		self.optimal_thrust = optimal_thrust
		self.nozzle_exit_p = nozzle_exit_p
		self.nozzle_area = nozzle_area
		self.mdot = sl_thrust/(i_sp_sl * GRAV)

		#initialize some non-rocket data which is important for analysis
		self.time = [0]
		self.mach = [0]
		#initiate matricies
		self.a_matrix = np.zeros(4)
		self.v_matrix = np.zeros(4)
		self.d_matrix = np.zeros(4)
		self.v_matrix = np.vstack((self.v_matrix, self.v_matrix))
		self.d_matrix = np.vstack((self.d_matrix, self.d_matrix))

	def three_dof_sim(self, step):
		print('Starting up the simulation:')
		while (self.d_matrix[-1,2] > -.01 ):
			#Get altitude and tangential velocity (we use these a lot)
			current_alt = self.d_matrix[-1,2]
			print('Current Altitude: {}'.format(current_alt))
			v_t = self.v_matrix[-1,3]
			theta = 85 * np.pi/180
			omega = 45 * np.pi/180
			'''
			1) Get current atomspheric state
			'''
			#atmospheric density is a function of altitude only (d_z)
			atm_rho = comp_dens(current_alt)
			#atmospheric temperature is a function of altitude only (d_z)
			atm_temp = comp_temp(current_alt)
			#mach number rocket is at
			mach = comp_mach(v_t, atm_temp)
			#Here would be the code to implement the gust model, magnitude
			#depending on current altitude
			'''
			2) Get current rocket state at time, altitude and etc...
			'''
			#Calculate coefficeient of drag, accounting for wave drag
			cd = wave_drag(mach, self.cd[0])
			#Calculate thrust, accounting for nozzle expansion ratio and
			#aerostatic back pressure loss
			thrust, self.mdot = calc_thrust(self.mdot, self.optimal_thrust, self.nozzle_exit_p, self.nozzle_area, current_alt, (self.total_mass[-1] - self.dry_mass))
			'''
			3) Run Equations of Motion, using gathered state information
			'''
			a_array = tangent_eom(thrust, theta, omega, self.total_mass[-1], atm_rho, self.front_area, v_t, cd)
			#very crude numerical integration method
			v_array, d_array = integrate_eom(a_array,step, self.d_matrix, self.v_matrix)
			'''
			4) Update values on rocket which change only as a function of time
			Eg.  Mass, staging, parachute if these values are modeled
			'''
			self.total_mass.append(self.total_mass[-1] - self.mdot *step)
			'''
			5) Organize/Parse data.  Update instance varaibles
			'''
			#update time
			self.time.append(self.time[-1] + step)
			#Update positional data
			self.a_matrix = np.vstack((self.a_matrix, a_array))
			self.v_matrix = np.vstack((self.v_matrix, v_array))
			self.d_matrix = np.vstack((self.d_matrix, d_array))
			#other data
			self.thrust.append(thrust)
			self.mach.append(mach)
			self.cd.append(cd)
			#Retart the loop
		'''
		6) Do a little, post simulation data clean-up on the hacks used to circumvent matrix out of bounds
		'''
		self.v_matrix = np.delete(self.v_matrix,0,0)
		self.d_matrix = np.delete(self.d_matrix,0,0)
		print('Simulation done!')

	def plot_xyz(self):
		ax = plt.axes(projection='3d')
		ax.set_xlabel('x')
		ax.set_ylabel('y')
		ax.set_zlabel('z');
		ax.plot3D(self.d_matrix[:,0], self.d_matrix[:,1], self.d_matrix[:,2], c=self.d_matrix[:,2], cmap='Greens')
		plt.show()
	def plot_mains(self):
		plt.figure()
		plt.subplot(2,1,1)
		plt.title('Displacement vs time')
		plt.xlabel('time (s)')
		plt.ylabel('displacement (feet)')
		plt.plot(self.time, self.d_matrix[:,3])
		plt.subplot(2,1,2)
		plt.title('Displacement vs time')
		plt.xlabel('time (s)')
		plt.ylabel('tangential velocity (feet/s)')
		plt.plot(self.time, self.v_matrix[:,3])
		plt.show()
def tangent_eom(thrust, theta, omega, mass, density, area, velocity, c_drag):
	'''check this EoM for correct multiplcation to gravity'''
	a_t = thrust/mass - (1/2)*c_drag*density*(velocity**2)/mass - GRAV*np.sin(theta)
	a_x = a_t*np.cos(theta)*np.cos(omega)
	a_y = a_t*np.cos(theta)*np.sin(omega)
	a_z = a_t*np.sin(theta)

	a_array = np.array([a_x, a_y, a_z, a_t])

	return a_array

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

if __name__ == '__main__':
	#enter parameters in this order:
	#sl_thrust(lbs), optimal_thrust, nozzle_exit_p (psi), nozzle_area (in^2), i_sp_sl, total_weight(lbf), dry_weight, front_area(ft), cd
	thedude = Rocket(3038.868, 3250, 8.103, 43.78772307748711, 230, 730, 300, np.pi*(11.2/12)**2, .25)
	thedude.three_dof_sim(.001)
	thedude.plot_mains()
