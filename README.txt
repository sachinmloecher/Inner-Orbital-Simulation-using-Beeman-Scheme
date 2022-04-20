README for Three-step Beeman Scheme Simulation of the Inner Solar System

- File 'SimParameters.txt' for Simulation parameters
- Input in the format 'body6 = [name, color, mass, [x_pos,y_pos], [x_vel,y_vel], planet_patch_size]'
- Simulation('SimParameters.txt', numer_of_bodies) to create simulation
- Simulation.Display() to display solar system
- Simulation.compute_orbital_periods() to print simulated orbital periods
- Simulation.show_energy_graph() to display total energy graph
- Simulation.energy_variance_superposition_of_planet_frequencies() to test hypothesis 2 (get table)
- Simulation.energy_variance_decreases_with_timestep() to test hypothesis 1 (get graph 2)
- Total Energy is written to 'Energy.txt'