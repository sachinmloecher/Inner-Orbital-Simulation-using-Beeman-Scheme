import numpy as np
from numpy.linalg import norm
import math
import copy
import matplotlib.pyplot as plt
from matplotlib import rc
from matplotlib.animation import FuncAnimation
import imp
import os
import matplotlib.image as mpimg
import pandas as pd


# File reader helper function
def getVarFromFile(filename):
    import imp
    f = open(filename)
    global parameters
    param_path = os.path.join(os.path.dirname(os.path.abspath(filename)), 'SimParameters.txt')
    parameters = imp.load_source('parameters', param_path,f)
    f.close()

class Body(object):
    def __init__(self, name, color, mass, position, velocity, diameter):
        self.name = name
        self.color = color
        self.mass = mass
        self.position = position
        self.velocity = velocity
        self.curr_accel = None
        self.prev_accel = None
        self.diameter = diameter
        self.KE = (1/2)*(self.mass)*(norm(self.velocity)**2)
    def update_position(self, timestep):
        # Beeman Method
        next_position = self.position + self.velocity*timestep + (1/6)*(4*self.curr_accel - self.prev_accel)*(timestep**2)
        self.position = next_position
    def update_velocity(self, timestep, next_acceleration):
        # Beeman Method
        next_velocity = self.velocity + (1/6)*(2*next_acceleration + 5*self.curr_accel - self.prev_accel)*timestep
        self.velocity = next_velocity
        self.KE = (1/2)*(self.mass)*(norm(self.velocity)**2)
        # Update accelerations
        self.prev_accel = self.curr_accel
        self.curr_accel = next_acceleration
        

class Simulation(object):
    def __init__(self, filename, n_bodies):
        # Read simulation parameters from file
        self.filename = filename
        self.n_bodies = n_bodies
        getVarFromFile(filename)
        self.timestep = parameters.timestep
        self.iterations = parameters.iterations
        self.body_list = []
        self.patch_list = []
        self.orbital_periods = {}
        self.total_energy = []
        self.iteration = []
        self.KE = []
        self.PE = []
        # Read Bodies from file
        for i in range(1, n_bodies + 1):
            var = "parameters.body" + str(i)
            body_params = eval(var)
            name = body_params[0]
            color = body_params[1]
            mass = body_params[2]
            position = np.array(body_params[3])
            velocity = np.array(body_params[4])
            diameter = body_params[5]
            body = Body(name, color, mass, position, velocity, diameter)
            self.add_body(body)
        self.G = (6.67430e-11)
        # Set Initial Accelerations
        for body in self.body_list:
            body.curr_accel = self.compute_next_acceleration(body)
            body.prev_accel = body.curr_accel
    
    def add_body(self, body):
        # Add body object to body list and create/add corresponding circle for display
        self.body_list.append(body)
        patch = plt.Circle(tuple(body.position) , body.diameter, color = body.color , animated = True, label=body.name)
        self.patch_list.append(patch)
    
    def compute_next_acceleration(self, body):
        accel_components = []
        # Loop over all other bodies
        for p in (self.body_list):
            if p == body:
                continue
            else:
                # Calculate acceleration contributed by p and add to accel_components
                r_ji = body.position - p.position
                accel_components.append((p.mass)*(r_ji)/(norm(r_ji)**3))
                # Use gravitational force law and sum components * -G
                next_accel = (-self.G)*np.array([sum(component) for component in zip(*accel_components)])
        return next_accel
    
    def compute_total_energy(self):
        total_energy = self.compute_total_KE() + self.compute_total_PE()
        return total_energy
    
    def compute_total_KE(self):
        total_KE = 0
        for body in self.body_list:
            total_KE += body.KE
        return total_KE
    
    def compute_total_PE(self):
        PE = 0
        for i in self.body_list:
            for j in self.body_list:
                if i == j:
                    continue
                else:
                    r_ij = norm(i.position - j.position)        
                    PE += (self.G)*(i.mass)*(j.mass)/r_ij
        PE = (-1/2)*PE
        return PE
    
    def step_forwards(self, i):
        # First iteration, just use initial values
        if i == 0:
            return self.patch_list
        else:
            # For each body, update position
            for body in range(len(self.patch_list)):
                pos_y = self.body_list[body].position[1]
                self.body_list[body].update_position(self.timestep)
                new_pos_y = self.body_list[body].position[1]
                self.patch_list[body].center = tuple(self.body_list[body].position)
                # If planet has completed orbit, record the number of iterations it took for orbital period computation
                if (pos_y < 0) & (new_pos_y >= 0) & (self.body_list[body].name not in self.orbital_periods.keys()) & (i > 2):
                    self.orbital_periods[self.body_list[body].name] = i
                # Add Total Energy to dataframe
                self.total_energy.append(self.compute_total_energy())
                self.KE.append(self.compute_total_KE)
                self.PE.append(self.compute_total_PE)
                self.iteration.append(i)
                if i%20 == 0:
                    f = open('Energy.txt', 'w')
                    f.write("Total Energy of System: " + str(self.compute_total_energy()) + " J")
                    f.write('\n')
                    f.close()
            # For each body, compute next acceleration, update velocity and acceleration
            for body in range(len(self.patch_list)):
                a = self.compute_next_acceleration(self.body_list[body])
                self.body_list[body].update_velocity(self.timestep, a)
            return self.patch_list
        
    def Display(self):
        # Create Figure and set scales so planets are visible (needs to be changed if plotting different orbital radius')
        fig = plt.figure()
        ax = plt.axes()
        ax.set_xlim(-400000000000, 400000000000)
        ax.set_ylim(-400000000000, 400000000000)
        ax.set_title("Inner Planets of the Solar System")
        ax.set_xlabel("Position (m)")
        ax.set_ylabel("Position (m)")
        # Set space background image
        img = plt.imread('space.jpg')
        imgplot = plt.imshow(img, extent=[-400000000000,400000000000,-400000000000,400000000000])
        # Set x/y pixel count equal so planets are circles not ellipses
        ax.axis("equal")
        # Add each planet to the Display and add legend
        labels = []
        for body in range(len(self.patch_list)):
            ax.add_patch(self.patch_list[body])
            labels.append(self.body_list[body].name)
        plt.legend(self.patch_list, labels)
        # Animate display
        numFrames = self.iterations + 1
        anim = FuncAnimation(fig, self.step_forwards, numFrames, repeat = True, interval = 1, blit = True)
        plt.show()
    
    def compute_orbital_periods(self):
        s = Simulation(self.filename, self.n_bodies)
        for i in range(self.iterations):
            s.step_forwards(i)
            
        # Print Orbital Periods (yrs/days) for each planet
        for body in self.body_list:
            if body.name == 'Sun':
                continue
            else:
                print("Simulated " + body.name + " Orbital Period (yrs): " + str(s.orbital_periods[body.name]/365))
                print("Simulated " + body.name + " Orbital Period (days): " + str(s.orbital_periods[body.name]))
    
    def show_energy_graph(self):
        s = Simulation(self.filename, self.n_bodies)
        for i in range(self.iterations):
            s.step_forwards(i)
        
        # Graph Total Energy by iteration
        fig = plt.figure()
        ax = plt.axes()
        ax.set_title("Total Energy of System")
        ax.set_xlabel("Number of Days (Iteration)")
        ax.set_ylabel("Total Energy (J)")
        ax.plot(s.iteration, s.total_energy)
        
    # Test Hypothesis: "The variance in the total energy of the system will decrease as timestep becomes more exact"
    def energy_variance_decreases_with_timestep(self):
        # 200 different timesteps from 1 second to 1 day (86400)
        timesteps = np.linspace(1,86400,200)
        variances = []
        for timestep in timesteps:
            s = Simulation(self.filename, self.n_bodies)
            s.timestep = timestep
            s.iterations = 365
            for i in range(self.iterations):
                s.step_forwards(i)
            variances.append(np.var(s.total_energy))
        fig = plt.figure()
        ax = plt.axes()
        ax.set_title("Variance of Total Energy by Timestep")
        ax.set_xlabel("Timestep (seconds)")
        ax.set_ylabel("Variance in Total Energy (J^2)")
        ax.plot(timesteps, variances)
    
    # Test Hypothesis: "The variance in the total energy of the system will increase as the number of planets increases"
    def energy_variance_superposition_of_planet_frequencies(self):
        variances = []
        n_bodies = []
        # Total Energy plots for each new body
        for i in range(2, self.n_bodies + 1):
            fig = plt.figure()
            ax = plt.axes()
            ax.set_title("Total Energy of System")
            ax.set_xlabel("Number of Days (Iteration)")
            ax.set_ylabel("Total Energy (J)")
            s = Simulation(self.filename, i)
            for j in range(self.iterations):
                s.step_forwards(j)
            variances.append(np.var(s.total_energy))
            n_bodies.append(i)
            ax.plot(s.iteration, s.total_energy)
        # Show table of number of bodies and corresponding variance
        fig = plt.figure()
        ax = plt.axes()
        fig.patch.set_visible(False)
        ax.axis('off')
        ax.axis('tight')
        table_tuples = list(zip(n_bodies, variances))
        df = pd.DataFrame(np.array(table_tuples), columns=['N-Bodies', 'Variance (J^2)'])
        ax.table(cellText=df.values, colLabels=df.columns, loc='center')
        fig.tight_layout()
        plt.show()
        
        
def main():
    
    S = Simulation('SimParameters.txt', 5)
    S.Display()
    S.compute_orbital_periods()
    S.show_energy_graph()
    #S.energy_variance_decreases_with_timestep()
    #S.energy_variance_superposition_of_planet_frequencies()
    #print("Standard Deviation of Total Energy as percentage: " + str(100*np.sqrt(np.var(S.total_energy))/np.mean(S.total_energy)) + "%")

main()