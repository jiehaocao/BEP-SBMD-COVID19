"""Author: J.H.Cao
   Computational Biology group
   Biomedical Engineering
   Eindhoven University of Technology"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from estimating_compmod_params import read_data, seird_model

plt.rcParams['figure.facecolor'] = 'white'

# Read the data (cumulative number of infected)
data = read_data('./data/corona_deaths.csv')
N = 17000000.0
nr_days = 187
shift = 30
nr_days_shift = nr_days + shift

# The initial values of the SIR (system)
S0 = N - 1.0
E0 = 1.0
I0 = 0.0
R0 = 0.0
D0 = 0.0
gamma = 1.0 / 8.0
delta = 1.0 / 4.0
alpha = 1.0 / 18.0
p_i_d = 0.04

# Extracting the estimated time-dependent beta parameters
estimated_sol_seird = seird_model(nr_days_shift, 5.00, 0.52, 0.20, 57.2, 0.02)
D_sol = estimated_sol_seird[5][shift:]  # This is the simulation without Kalman filter applied
R0_values = estimated_sol_seird[6]

# Time-dependent beta
beta = [i * gamma for i in R0_values[shift:]]

# The error covariance matrix

P = np.array([[0.01, 0, 0, 0, 0],
              [0, 0.01, 0, 0, 0],
              [0, 0, 0.01, 0, 0],
              [0, 0, 0, 0.01, 0],
              [0, 0, 0, 0, 0.01]])

# The process noise covariance

Q = np.array([[0.1, 0, 0, 0, 0],
              [0, 0.1, 0, 0, 0],
              [0, 0, 0.1, 0, 0],
              [0, 0, 0, 0.1, 0],
              [0, 0, 0, 0, 0.1]])

# The measurement noise covariance

R = np.array([[0.01, 0, 0, 0, 0],
              [0, 0.01, 0, 0, 0],
              [0, 0, 0.01, 0, 0],
              [0, 0, 0, 0.01, 0],
              [0, 0, 0, 0, 0.01]])

# Observation matrix H

H = np.array([[1.0, 0, 0, 0, 0],
              [0, 1.0, 0, 0, 0],
              [0, 0, 1.0, 0, 0],
              [0, 0, 0, 1.0, 0],
              [0, 0, 0, 0, 1.0]])

# The initial values at k = 0
S_kf_est = [S0]
E_kf_est = [E0]
I_kf_est = [I0]
R_kf_est = [R0]
D_kf_est = [D0]

t_length = np.linspace(0, nr_days - 1, nr_days)

# For-loop over the simulation's time-span
# This approach uses the standard Kalman filter
# The discrete-time SIR model was computed
# the time step-size was chosen to be dt = 1.0

for k in range(0, len(t_length) - 1):
    # state-transition matrix A withe time-depent K
    A_11 = (1 - (beta[k] * I_kf_est[k]) / (2 * N))
    A_13 = -(beta[k] * S_kf_est[k]) / (2 * N)
    A_21 = (beta[k] * I_kf_est[k]) / (2 * N)
    A_22 = (1 - delta)
    A_23 = (beta[k] * S_kf_est[k]) / (2 * N)
    A_32 = delta
    A_33 = (1 - gamma * (1 - p_i_d) - alpha * p_i_d)
    A_43 = gamma * (1-p_i_d)
    A_53 = alpha * p_i_d

    A = np.array([[A_11, 0.0, A_13, 0.0, 0.0],
                  [A_21, A_22, A_23, 0.0, 0.0],
                  [0.0, A_32, A_33, 0.0, 0.0],
                  [0.0, 0.0, A_43, 1.0, 0.0],
                  [0.0, 0.0, A_53, 0.0, 1.0]], dtype=float)

    # The past matrix x_(k)
    x_past = np.array([[S_kf_est[k]],
                       [E_kf_est[k]],
                       [I_kf_est[k]],
                       [R_kf_est[k]],
                       [D_kf_est[k]]], dtype=float)

    # Prediction step x_(k+1) = A * x_(k)
    x_priori = np.dot(A, x_past)
    P = np.dot(np.dot(A, P), A.T) + Q

    # Measurement update equations a.k.a the Kalman Gain calculation
    S = np.dot(np.dot(H, P), H.T) + R
    S.astype(float)
    K = np.dot(np.dot(P, H.T), np.linalg.inv(S))

    # Observed measurements matrix

    if k < len(data) - 1:
        Y = np.array([[x_priori[0][0]],
                      [x_priori[1][0]],
                      [x_priori[2][0]],
                      [x_priori[3][0]],
                      [data[k + 1]]])
    else:
        Y = np.array([[x_priori[0][0]],
                      [x_priori[1][0]],
                      [x_priori[2][0]],
                      [x_priori[3][0]],
                      [x_priori[4][0]]])

    x_diff = Y - np.dot(H, x_priori)
    x_post = x_priori + np.dot(K, x_diff)

    S_kf_est.append(x_post[0][0])
    E_kf_est.append(x_post[1][0])
    I_kf_est.append(x_post[2][0])
    R_kf_est.append(x_post[3][0])
    D_kf_est.append(x_post[4][0])

plt.figure(figsize=(10, 8))
plt.plot(data, "o", markersize=8, color='black')
plt.plot(D_sol, lw=3, color='blue')
plt.plot(D_kf_est, '--', lw=3, color='darkorange')
plt.legend(["Observed measurements", "Simulation without Kalman filter", "Kalman filtered simulation"], loc='best')
plt.xlabel("Time (days)" + "\n\n" + "from February 27, 2020", fontsize=14)
plt.ylabel("Number of Deaths", fontsize=14)
plt.title("'Deceased' from SEIRD model: Kalman filter vs without Kalman filter", fontsize=14, weight='bold')
plt.show()
