import numpy as np
import matplotlib.pyplot as plt
import csv
from scipy.optimize import curve_fit
import math as m
import os

values = {f'array_{i}{j}{k}': [] for i in range(1, 5) for j in range(1, 3) for k in range(1, 3)}
folder_path = 'Data_wB/'
#import latex
plt.rcParams['text.usetex'] = True

#load data
for i in range(1, 5):
    for j in range(1, 3):
        for k in range(1, 3):
            file_name = f"{i}{j}{k}.csv"
            file_path = os.path.join(folder_path, file_name)

            try:
                with open(file_path, 'r') as my_file:
                    csv_reader = csv.reader(my_file)
                    data_ = np.array(list(csv_reader), dtype=float)
                    array_key = f'array_{i}{j}{k}'
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

for i in range(1, 5):
    for j in range(1, 3):
        for k in range(1, 3):
            array_key = f'array_{i}{j}{k}'
            time = values[array_key][0][:, 0]
            angle = values[array_key][0][:, 1]
            plt.plot(time, angle, label='Measurements', marker=',')

            #fit
            def fit(t, a, b, c, d, omega, phi):
                return a * np.exp(d * t) * np.sin(omega * t + phi) + b * t + c

            params, covariance = curve_fit(fit, time, angle, p0=[3, 0.000002, 1, 0.0002, 0.02, -1])
            a, b, c, d, omega, phi = params
            errors = np.sqrt(np.diag(covariance))

            x_fit = np.linspace(0, int(max(time)), 2 * int(max(time)))
            y_fit = fit(x_fit, a, b, c, d, omega, phi)

            linear_shift = b * x_fit + c
            damping = a * np.exp(d * x_fit)
            plt.plot(x_fit, linear_shift, label='Linear Shift', linestyle='--')
            plt.plot(x_fit, damping, label='Damping', linestyle=':')

            plt.title("Angle evolution of the torsion pendulum over time")
            plt.ylabel(r"Measured angle $\theta$")
            plt.xlabel(r"Time $t$ [s]")
            plt.plot(x_fit, y_fit, label='Fit')
            plt.legend()
            if i == 1 & j == 1 & k == 1:
                plt.show()

            #residual plot
            residuals = angle - fit(time, a, b, c, d, omega, phi)
            plt.plot(time, residuals, marker=',', markersize=0.0002)
            plt.title("Residual plot of the measurements")
            plt.ylabel("Angle")
            plt.xlabel("Data counts")
            if i == 1 & j == 1 & k == 1:
                plt.show()


            #add values
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
            phiM_err_calc = m.pi / 180 * m.sqrt((errors[2])**2 + (1/2 * time[-1] * errors[1])**2)
            phiM = np.append(phiM, phiM_calc)
            phiM_err = np.append(phiM_err, phiM_err_calc)



#final phiM calculation
phiM_pos = np.array([])
phiM_neg = np.array([])
for i in range(len(phiM)):
    if phiM[i] < 0:
        phiM_neg = np.append(phiM_neg, phiM[i])
    else:
        phiM_pos = np.append(phiM_pos, phiM[i])
finalphiM_pos = np.mean(phiM_pos)
finalphiM_neg = np.mean(phiM_neg)
finalphiM = (-finalphiM_neg + finalphiM_pos) / 4

#print values
result = np.array([])
for i in range(1, 5):
    for j in range(1, 3):
        for k in range(1, 3):
            result = np.append(result, f"{i}{j}{k}")
total = [result, a_val, b_val, c_val, d_val, omega_val, phi_val]
#for i, row in enumerate(zip(*total), start=1):
    #print(f"{', '.join(map(str, row))}")


omega_final = np.mean(omega_val)
omega_finalerr = np.std(omega_val)
#print(omega_final)
#print(omega_finalerr)
T_with = 2 * m.pi / omega_final
#print("T with : ", T_with)
print("Final phi M : ", finalphiM)
finalphiM_err = finalphiM * m.sqrt((np.std(phiM_pos)/4)**2 + np.std(phiM_neg)/4**2)
print("Final phi M err", finalphiM_err)
print("Period : ", T_with)
print("Period error : ", 2 * m.pi / omega_final**2 * omega_finalerr)
#G calculation
M = 10
ma = 0.02
r = 0.024
R = 0.1
beta = 53.18 * m.pi/180
h = 0.2
b = 1/2 * (R/r + r/R)
bprime = 1/2 * (R/r + r/R + h**2/(r*R))
theta = 2.304e-5


G = 4 * m.pi**2 * theta * finalphiM / (ma * M * (2 * r * R)**(-1/2) * m.sin(beta) * ((b - m.cos(beta))**(-3/2) - (bprime + m.cos(beta))**(-3/2))* T_with**2)
print("G : ", G)
Gerr = G * 4 * m.pi**2 / (ma * M * (2 * r * R)**(-1/2) * m.sin(beta) * ((b - m.cos(beta))**(-3/2) - (bprime + m.cos(beta))**(-3/2))) * m.sqrt((finalphiM_err/T_with**2)**2 + (2*finalphiM/T_with)**2)
print("Gerr : ", Gerr)
G_theo = 6.674e-11
print("G relative error : ", (G-G_theo)/G_theo *100)

g0 = 9.816
g0_err = 0.019
Radius = 6378137
EarthM = g0 * Radius**2 / G
EarthM_err = EarthM * m.sqrt((g0 * g0_err)**2 + (Gerr/G)**2)
print("Earth mass : ", EarthM)
print("Earth mass error : ", EarthM_err)

EarthM_theo = 5.972e24
print("Earth mass Relative error : ", (EarthM - EarthM_theo)/EarthM_theo * 100)
