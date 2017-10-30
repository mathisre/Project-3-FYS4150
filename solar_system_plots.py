from sys import argv
import numpy as np
import matplotlib.pyplot as plt

file = open(argv[6], 'r') #File data read from files
n = int(file.readline())

with file as filename:
    lines = [line.split() for line in filename]
numb = len(lines[1]) #Number of planets times 2
x = np.zeros((n,numb))


#Read planet coordinates into x matrix
for k in range(0,n):
    for m in range(0, numb):
        x[k][m] = float(lines[k][m])
file.close()

#Calculate distance from planet to barycenter
r = np.zeros((n,numb/2))
for k in range(0, n):
    for m in range(0, numb/2):
        r[k][m] = (x[k][2*m]**2 + x[k][2*m+1]**2)**0.2

planets = ['Sun','Mercury', 'Venus', 'Earth', 'Mars',  'Jupiter', 'Saturn', 'Uranus', 'Neptune', 'Pluto', 'Moon']

planets = ['Sun', 'Earth', 'Jupiter']

for k in range(0, numb,2):
    plt.plot(x[:,k], x[:,k+1], label = planets[k/2])
plt.xlabel(r'X position [AU]')
plt.ylabel(r'Y position [AU]')

plt.xlabel(r'X velocity [AU/yr]')
plt.ylabel(r'Y velocity [AU/yr]')

plt.legend(loc = 2, fancybox = True, handlelength = 1, handletextpad = 0.1)
plt.title(r"Velocity of Earth using Forward Euler")
#plt.axis([--50, 50, -50, 50])
plt.show()

time = np.linspace(0, 400,n)
for k in range(0, numb/2):
    plt.plot(time,r[:,k], label = planets[k])
plt.xlabel(r'Time [yrs]')
plt.ylabel(r'Distance from center of mass')
plt.legend(loc = 2, fontsize = 13, scatteryoffsets = [1], markerscale = 100,
           fancybox = True, handlelength = 1, handletextpad = 0.1)
plt.title(r"Planet distance from barycenter with step size 0.0002 years")
plt.show()


#plt.axis([-0.5, 0.5, -0.5, 0.5])
"""""
plt.ion()
for i in range(0, n, 10):
    for k in range(0, numb, 2):
        plt.scatter(x[i,k],x[i,k+1])

    plt.pause(0.00000000000001)

while True:
    plt.pause(0.00000000000001)
"""""



#plt.plot(vx,vy,)
#plt.xlabel('x velocity')
#plt.ylabel('y velocity')
#plt.show()

