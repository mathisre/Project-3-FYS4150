from sys import argv
import numpy as np
import matplotlib.pyplot as plt


file = open(argv[1], 'r')
n = int(file.readline())

#Read number of lines in the text file
with file as f:
    numb = sum(1 for _ in f)

file = open(argv[1], 'r')
n = int(file.readline())

with file as filename:
    lines = [line.split() for line in filename]

ang = np.zeros(numb)
time = np.zeros(numb)

for k in range(0,numb):
    ang[k] = float(lines[k][0])
    time[k] = float(lines[k][1])
file.close()
plt.plot(time,ang)
plt.xlabel(r'Time (years)')
plt.ylabel(r'Perihelion angle')
plt.title(r"Mercury Perihelion angle with relativistic correction")
#plt.axis([--50, 50, -50, 50])
plt.show()

