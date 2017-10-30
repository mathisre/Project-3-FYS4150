# Project-3-FYS4150

Earth-Sun Euler and Verlet include the algorithms for just the Sun-Earth system using Earth initial conditions (x,y,v_x,v_y)=(0, 1, -2pi, 0). Both of them read files and number of points from the command line. Both of them also use armadillo. Respectively they write eartheuler.txt and earthsun_vervlet.txt. These files contain lines of (x, y, v_x, v_y) for Earth.

The object oriented code runs without any command line arguments and does not use armadillo. It initializes planet objects and adds them to the solver. The verlet solver creates verleto.txt that includes the positions of all the planets. The program also runs a mercury sun two-body planet with the GR correction. That creates perihelion.txt that creates the perihelion angle and the time. Finally the code also runs Forward Euler and it creates a euleroo.txt file. 

The non-perihelion text files from the object oriented code contain lines of x-y coordinates for all the planets. 

There are two python programs that create plots. Both of them read .txt files from the commandline. solar_system_plots.py plots takes data from the object oriented code and plots the positions of the planets. It also calculates the radial distance from the center of mass and plots that as a function of time. perihelion.py calculates the number of lines in the .txt file and then plots the perihelion angle as a function of time.

Hey I'm sorry I didn't include a proper conclusion. Had a lot to do but that didn't happen. That's my bad.
