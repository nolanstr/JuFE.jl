import numpy as np

R_1y = 10
R_5y = 10

F2_6, F4_7 = 0, 0

"""
Method of Sections
"""

#Section 1

theta = np.arctan(2)

F1_8 = -R_1y / np.sin(theta)
F1_2 = -np.cos(theta) * F1_8

#Section 2

alpha = np.arctan(1)

F8_9 = F1_2 - R_1y
F6_9 = (1/ (2*np.cos(alpha))) * (((R_1y - 2*F1_2)/np.tan(alpha)) - R_1y)
F3_6 = ((R_1y - 2*F1_2)/np.sin(alpha)) - F6_9

#Section 3

F5_10 = -R_1y / np.sin(theta)
F4_5 = -np.cos(theta) * F5_10

#Section 4

F7_9 = (1/(2*np.sin(alpha))) * (((1-np.tan(alpha))*R_5y) - 2*F4_5)
F3_7 = F7_9 + (R_5y/np.cos(alpha))
F9_10 = -F4_5 - np.sin(alpha)*(F7_9 + F3_7)


"""
Method of joints
"""

#Node 3
F3_9 = -np.cos(alpha) * (F3_6 + F3_7)

#Node 2
F2_3 = F1_2

#Node 4
F3_4 = F4_5

#Node 8
beta = np.arctan(0.5)
F6_8 = -np.cos(beta) * F1_8

#Node 10
F7_10 = -np.cos(beta) * F5_10

F = [F1_2,
     F2_3,
     F3_4,
     F4_5,
     F1_8,
     F2_6,
     F3_6,
     F3_9,
     F3_7,
     F4_7,
     F5_10,
     F6_8,
     F6_9,
     F7_9,
     F7_10,
     F8_9,
     F9_10]

for i, F_i in enumerate(F):

    print(f"Force in Truss Member/Element {i+1} = {round(F_i, 3)} kN")

import pdb;pdb.set_trace()
