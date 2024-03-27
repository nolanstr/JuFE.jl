import numpy as np

theta1 = np.arctan(2)
theta2 = np.arctan(1)
theta3 = np.arctan(1/2)

F1 = F2 = -(np.cos(theta1)/np.sin(theta1)) * 10
F3 = F4 = -(np.cos(theta1)/np.sin(theta1)) * 10
F5 = 10 / np.sin(theta1)
F6 = 0
F7 = (-np.cos(theta3) / (2*np.sin(theta1)*np.sin(theta2))) * 10
F10 = 0
F11 = 10/np.sin(theta1)
F12 = (-np.cos(theta3)/np.sin(theta1)) * 10
F13 = (np.cos(theta3)/(2*np.sin(theta1)*np.sin(theta2)))*10
F15 = (-np.cos(theta3)/np.sin(theta1)) * 10
F9 = ((F2-F3)/np.sin(theta2)) + F7 
F14 = -F9
F8 = 20 - np.sin(theta2) * (F13 + F14)
F16 = (np.sin(theta3)/np.sin(theta1)) * 10
F17 = (np.sin(theta3)/np.sin(theta1)) * 10
F = [F1,
     F2,
     F3,
     F4,
     F5,
     F6,
     F7,
     F8,
     F9,
     F10,
     F11,
     F12,
     F13,
     F14,
     F15,
     F16,
     F17]

for i, F_i in enumerate(F):
    print(f"F{i+1} = {F_i}")
