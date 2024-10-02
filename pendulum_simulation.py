# pendulum_simulation.py

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import os

# =========================
# Pendulum Simulation Class
# =========================

class PendulumSimulation:
    def __init__(self, m=1.0, L=1.0, g=9.81, c=0.05,
                 theta0=0.2, omega0=0.0,
                 dt=0.01, t_max=20.0):
        """
        Initialize the pendulum simulation parameters.

        Parameters:
            m (float): Mass of the pendulum bob (kg).
            L (float): Length of the pendulum string (m).
            g (float): Acceleration due to gravity (m/s^2).
            c (float): Damping coefficient (kg*m^2/s).
            theta0 (float): Initial angular displacement (radians).
            omega0 (float): Initial angular velocity (rad/s).
            dt (float): Time step for simulation (s).
            t_max (float): Total simulation time (s).
        """
        self.m = m
        self.L = L
        self.g = g
        self.c = c
        self.theta0 = theta0
        self.omega0 = omega0
        self.dt = dt
        self.t_max = t_max
        self.num_steps = int(t_max / dt)
        self.time = np.linspace(0, t_max, self.num_steps)
        
        # Initialize arrays for Euler's Method
        self.theta_euler = np.zeros(self.num_steps)
        self.omega_euler = np.zeros(self.num_steps)
        self.theta_euler[0] = theta0
        self.omega_euler[0] = omega0
        
        # Initialize arrays for RK4 Method
        self.theta_rk4 = np.zeros(self.num_steps)
        self.omega_rk4 = np.zeros(self.num_steps)
        self.theta_rk4[0] = theta0
        self.omega_rk4[0] = omega0
        
        # Ensure images directory exists
        if not os.path.exists('images'):
            os.makedirs('images')

    def derivatives(self, theta, omega):
        """
        Compute the derivatives for the pendulum system.

        Parameters:
            theta (float): Current angular displacement (radians).
            omega (float): Current angular velocity (rad/s).

        Returns:
            dtheta_dt (float): Derivative of theta.
            domega_dt (float): Derivative of omega.
        """
        dtheta_dt = omega
        domega_dt = - (self.g / self.L) * np.sin(theta) - (self.c / (self.m * self.L)) * omega
        return dtheta_dt, domega_dt

    def run_euler(self):
        """
        Execute Euler's Method for the pendulum simulation.
        """
        for i in range(1, self.num_steps):
            theta = self.theta_euler[i-1]
            omega = self.omega_euler[i-1]
            
            dtheta_dt, domega_dt = self.derivatives(theta, omega)
            
            self.theta_euler[i] = theta + self.dt * dtheta_dt
            self.omega_euler[i] = omega + self.dt * domega_dt

    def run_rk4(self):
        """
        Execute the 4th-Order Runge-Kutta (RK4) Method for the pendulum simulation.
        """
        for i in range(1, self.num_steps):
            theta = self.theta_rk4[i-1]
            omega = self.omega_rk4[i-1]
            
            # First estimate
            k1_theta, k1_omega = self.derivatives(theta, omega)
            
            # Second estimate
            k2_theta, k2_omega = self.derivatives(theta + 0.5 * self.dt * k1_theta,
                                                 omega + 0.5 * self.dt * k1_omega)
            
            # Third estimate
            k3_theta, k3_omega = self.derivatives(theta + 0.5 * self.dt * k2_theta,
                                                 omega + 0.5 * self.dt * k2_omega)
            
            # Fourth estimate
            k4_theta, k4_omega = self.derivatives(theta + self.dt * k3_theta,
                                                 omega + self.dt * k3_omega)
            
            # Update states
            self.theta_rk4[i] = theta + (self.dt / 6.0) * (k1_theta + 2*k2_theta + 2*k3_theta + k4_theta)
            self.omega_rk4[i] = omega + (self.dt / 6.0) * (k1_omega + 2*k2_omega + 2*k3_omega + k4_omega)

    def calculate_energy(self, theta, omega):
        """
        Calculate the total mechanical energy of the pendulum.

        Parameters:
            theta (ndarray): Angular displacement array (radians).
            omega (ndarray): Angular velocity array (rad/s).

        Returns:
            energy (ndarray): Total mechanical energy array (Joules).
        """
        potential = self.m * self.g * self.L * (1 - np.cos(theta))
        kinetic = 0.5 * self.m * (self.L * omega)**2
        energy = potential + kinetic
        return energy

    def plot_displacement(self):
        """
        Plot Angular Displacement vs Time for both Euler's Method and RK4.
        """
        plt.figure(figsize=(12, 6))
        plt.plot(self.time, self.theta_euler, label="Euler's Method", alpha=0.7)
        plt.plot(self.time, self.theta_rk4, label='RK4 Method', alpha=0.7)
        plt.title('Angular Displacement vs Time')
        plt.xlabel('Time (s)')
        plt.ylabel('Angle (radians)')
        plt.legend()
        plt.grid(True)
        plt.savefig('images/angular_displacement.png')
        plt.show()

    def plot_phase_space(self):
        """
        Plot Phase Space (Angular Velocity vs Angular Displacement) for both Euler's Method and RK4.
        """
        plt.figure(figsize=(12, 6))
        plt.plot(self.theta_euler, self.omega_euler, label="Euler's Method", alpha=0.7)
        plt.plot(self.theta_rk4, self.omega_rk4, label='RK4 Method', alpha=0.7)
        plt.title('Phase Space Plot')
        plt.xlabel('Angle (radians)')
        plt.ylabel('Angular Velocity (rad/s)')
        plt.legend()
        plt.grid(True)
        plt.savefig('images/phase_space.png')
        plt.show()

    def plot_energy(self):
        """
        Plot Total Mechanical Energy vs Time for both Euler's Method and RK4.
        """
        energy_euler = self.calculate_energy(self.theta_euler, self.omega_euler)
        energy_rk4 = self.calculate_energy(self.theta_rk4, self.omega_rk4)
        
        plt.figure(figsize=(12, 6))
        plt.plot(self.time, energy_euler, label="Euler's Method", alpha=0.7)
        plt.plot(self.time, energy_rk4, label='RK4 Method', alpha=0.7)
        plt.title('Total Mechanical Energy vs Time')
        plt.xlabel('Time (s)')
        plt.ylabel('Energy (Joules)')
        plt.legend()
        plt.grid(True)
        plt.savefig('images/energy_conservation.png')
        plt.show()

    def create_animation(self, method='RK4', save=False):
        """
        Create an animation of the pendulum's motion.

        Parameters:
            method (str): Numerical method used ('Euler' or 'RK4').
            save (bool): Whether to save the animation as a GIF.
        """
        if method == 'Euler':
            theta = self.theta_euler
        elif method == 'RK4':
            theta = self.theta_rk4
        else:
            raise ValueError("Method must be 'Euler' or 'RK4'")
        
        x = self.L * np.sin(theta)
        y = -self.L * np.cos(theta)
        
        fig, ax = plt.subplots(figsize=(6,6))
        ax.set_xlim(-self.L-0.5, self.L+0.5)
        ax.set_ylim(-self.L-0.5, self.L+0.5)
        ax.set_aspect('equal')
        ax.grid()
        ax.set_title(f'Pendulum Simulation using {method} Method')
        
        line, = ax.plot([], [], 'o-', lw=2)
        time_text = ax.text(0.05, 0.95, '', transform=ax.transAxes)
        
        def init():
            line.set_data([], [])
            time_text.set_text('')
            return line, time_text
        
        def animate(i):
            line.set_data([0, x[i]], [0, y[i]])
            time_text.set_text(f'Time = {self.time[i]:.2f}s')
            return line, time_text
        
        ani = animation.FuncAnimation(fig, animate, frames=range(0, self.num_steps, 10),
                                      init_func=init, blit=True, interval=20)
        
        if save:
            ani.save(f'images/pendulum_{method.lower()}_animation.gif', writer='imagemagick')
            print(f'Animation saved as images/pendulum_{method.lower()}_animation.gif')
        
        plt.show()

    def run_all(self):
        """
        Execute the entire simulation workflow.
        """
        print("Running Euler's Method...")
        self.run_euler()
        print("Euler's Method completed.\n")
        
        print("Running RK4 Method...")
        self.run_rk4()
        print("RK4 Method completed.\n")
        
        print("Generating plots...")
        self.plot_displacement()
        self.plot_phase_space()
        self.plot_energy()
        print("Plots saved in the 'images/' directory.\n")
        
        print("Creating animation for RK4 Method...")
        self.create_animation(method='RK4', save=True)
        print("Animation completed.\n")
        
        print("Creating animation for Euler's Method...")
        self.create_animation(method='Euler', save=True)
        print("Animation completed.\n")
        
        print("Simulation completed successfully.")

# =========================
# Main Execution
# =========================

if __name__ == "__main__":
    # Create a PendulumSimulation instance with default parameters
    pendulum = PendulumSimulation()
    
    # Run the simulation
    pendulum.run_all()
