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
from utilities import integrate_eom
from utilities import eom
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
		#might use these to convert to vector form later
		#v_hat_0 = np.array(.996,.7071,.7071)
		#v_hat = v_hat_0 / linalg.norm(v_hat_0)
		rail = True
		theta = 85 * np.pi/180
		omega = 45 * np.pi/180
		while (self.d_matrix[-1,2] > -.01 ):
			if (rail == True):
				self.v_matrix[-1,2] = .000001*np.sin(theta)
				self.v_matrix[-1,0] = .000001*np.cos(theta)*np.cos(omega)
				self.v_matrix[-1,1] = .000001*np.cos(theta)*np.sin(omega)
				rail = False


			#Get altitude and tangential velocity (we use these a lot)
			current_alt = self.d_matrix[-1,2]
			#print(current_alt)
			#print('Current Altitude: {}'.format(current_alt))
			v_t = self.v_matrix[-1,3]
			'''
			1) Get current atomspheric state
			'''
			#atmospheric density is a function of altitude only (d_z)
			atm_rho = comp_dens(current_alt)
			#atmospheric temperature is a function of altitude only (d_z)
			atm_temp = comp_temp(current_alt)
			#mach number rocket is at
			mach = comp_mach(abs(v_t), atm_temp)
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
			a_array = eom_vector(thrust, self.total_mass[-1], atm_rho, self.front_area, self.v_matrix[-1,:], cd)
			#very crude numerical integration method
			v_array, d_array = integrate_eom_vector(a_array,step, self.d_matrix, self.v_matrix)
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
		ax.set_zlabel('z')
		ax.plot3D(self.d_matrix[:,0], self.d_matrix[:,1], self.d_matrix[:,2], 'gray')
		plt.show()
	def plot_all(self):
		for i in range (0,2):
			plt.figure()
			plt.subplot(3,1,1)
			plt.plot(self.time, self.d_matrix[:,i])
			plt.subplot(3,1,2)
			plt.plot(self.time, self.v_matrix[:,i])
			plt.subplot(3,1,3)
			plt.plot(self.time, self.a_matrix[:,i])
		plt.figure()
		plt.plot(self.time,self.thrust)
		plt.show()
	def plot_accel(self):
		plt.figure()
		plt.subplot(2,1,1)
		plt.title('Displacement vs time')
		plt.xlabel('time (s)')
		plt.ylabel('displacement (feet)')
		plt.plot(self.time, self.a_matrix[:,3])
		plt.subplot(2,1,2)
		plt.title('Displacement vs time')
		plt.xlabel('time (s)')
		plt.ylabel('tangential velocity (feet/s)')
		plt.plot(self.time, self.a_matrix[:,0])
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


if __name__ == '__main__':
	#enter parameters in this order:
	#sl_thrust(lbs), optimal_thrust, nozzle_exit_p (psi), nozzle_area (in^2), i_sp_sl, total_weight(lbf), dry_weight, front_area(ft), cd
	thedude = Rocket(3038.868, 3250, 8.103, 43.78772307748711, 230, 730, 330, np.pi/4*(11.2/12)**2, .25)
	thedude.three_dof_sim(.05)
	#thedude.plot_mains()
	thedude.plot_xyz()
	thedude.plot_all()
