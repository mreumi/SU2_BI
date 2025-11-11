import numpy as np
import matplotlib.pyplot as plt
import csv

file_name = 'TargetCp2.dat'
data = np.array([[]])  # Initialize an empty 2D array
with open(file_name) as csvfile:
    datareader = csv.reader(csvfile, delimiter=',')
    next(datareader, None) # Skip the header row
    for row in datareader:
        if data.size == 0:
            data = np.array(row, dtype=float).reshape(1,-1)
        else:
            data = np.concatenate((data, np.array(row, dtype=float).reshape(1,-1)), axis=0)


x = data[:, 1]
y = data[:, 2]
cp = data[:, 4]

plt.figure(figsize=(5, 3))
plt.scatter(x, y, c=cp, cmap='viridis', marker='o', label='Target Locations')
plt.colorbar(label='Pressure coefficient (Cp)')
plt.title('Target Locations from TargetCp2.dat')
plt.xlabel('X Coordinate')
plt.ylabel('Y Coordinate')
plt.axis('equal')
plt.grid(True)
plt.savefig('target_locations.png',  bbox_inches='tight') 
plt.show()