# 1 = noball-bbdoor
# 2 = noball-bblab
# 3 = noball-lbdoor
# 4 =noball-lblab

import numpy as np
import matplotlib.pyplot as plt
import csv
from scipy.optimize import curve_fit
import math as m
import os

values = {f'array_{i}': [] for i in range(1, 5) }
folder_path = 'Data_nB/'

#load data
for i in range(1,5):
    file_name = f"{i}.csv"
    file_path = os.path.join(folder_path, file_name)

    try:
        with open(file_path, 'r') as my_file:
            csv_reader = csv.reader(my_file)
            data_ = np.array(list(csv_reader), dtype=float)
            array_key = f'array_{i}'
            values[array_key].append(data_)
    except FileNotFoundError:
        print(f"Error: {file_name} not found.")

#start analysis

a_val = np.array([])
b_val = np.array([])
c_val = np.array([])
d_val = np.array([])
omega_val = np.array([])
phi_val = np.array([])

a_err = np.array([])
b_err = np.array([])
c_err = np.array([])
d_err = np.array([])
omega_err = np.array([])
phi_err = np.array([])

phiM = np.array([])
phiM_err = np.array([])

for i in range(1,5):
    array_key = f'array_{i}'
    time = values[array_key][0][:, 0]
    angle = values[array_key][0][:, 1]
    plt.plot(time, angle)


    # fit
    def fit(t, a, b, c, d, omega, phi):
        return a * np.exp(d * t) * np.sin(omega * t + phi) + b * t + c


    params, covariance = curve_fit(fit, time, angle, p0=[3, 0.000002, 1, 0.0002, 0.02, -1])
    a, b, c, d, omega, phi = params
    errors = np.sqrt(np.diag(covariance))

    x_fit = np.linspace(0, int(max(time)), 2 * int(max(time)))
    y_fit = fit(x_fit, a, b, c, d, omega, phi)
    plt.plot(x_fit, y_fit)

    # residual plot
    residuals = angle - fit(time, a, b, c, d, omega, phi)
    plt.plot(time, residuals, marker=',', markersize=0.0002)

    # add values
    a_val = np.append(a_val, a)
    b_val = np.append(b_val, b)
    c_val = np.append(c_val, c)
    d_val = np.append(d_val, d)
    omega_val = np.append(omega_val, omega)
    phi_val = np.append(phi_val, phi)

    a_err = np.append(a_err, errors[0])
    b_err = np.append(b_err, errors[1])
    c_err = np.append(c_err, errors[2])
    d_err = np.append(d_err, errors[3])
    omega_err = np.append(omega_err, errors[4])
    phi_err = np.append(phi_err, errors[5])

    phiM_calc = m.pi / 180 * (2 * c + b * time[-1]) / 2
    phiM_err_calc = m.pi / 180 * m.sqrt((errors[2]) ** 2 + (1 / 2 * time[-1] * errors[1]) ** 2)
    phiM = np.append(phiM, phiM_calc)
    phiM_err = np.append(phiM_err, phiM_err_calc)

#final phiM calculation
finalphiM = np.mean(phiM) / 4
finalphiM_err = np.std(phiM) / 4

#print values
result = np.array([])
for i in range(1, 5):
    result = np.append(result, f"{i}")
total = [result, a_val, b_val, c_val, d_val, omega_val, phi_val]
#for i, row in enumerate(zip(*total), start=1):
    #print(f"{', '.join(map(str, row))}")

#omega values
omega_final = np.mean(omega_val)
omega_finalerr = np.std(omega_val)

T_without = 2 * m.pi / omega_final


print("Final phi M : ", finalphiM)
print("Final phi M err", finalphiM_err)
print("Period : ", T_without)
print("Period error : ", 2 * m.pi / omega_final**2 * omega_finalerr)