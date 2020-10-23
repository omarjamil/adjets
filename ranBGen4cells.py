import numpy as np
import sys
import string
import random

x_cells = 3
y_cells = 3
z_cells = 3
filename = "B_field.par"
B = 2.5e-5
a = 0.0
b = 1.0

dataFile=open(filename, 'w')
openString="//Cell Pos = BField(Tesla) B_theta(degrees) B_phi(degrees)"
dataFile.write(openString+"\n")
for i in range(x_cells):
    for j in range(y_cells):
        for k in range(z_cells):
            xPos = str(i)
            yPos = str(j)
            zPos = str(k)
            B_theta=random.uniform(a,b)*(360.)
            B_phi=random.uniform(a,b)*(180.)
            #B_theta=90.
            #B_phi=90.
            cPos = xPos+yPos+zPos
            dataFile.write(cPos+" "+"="+" "+str(B)+" "+str(B_theta)[:5]+" "+str(B_phi)[:5]+"\n")

dataFile.close()

