"""Author: J.H.Cao
   Computational Biology group
   Biomedical Engineering
   Eindhoven University of Technology"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from estimating_compmod_params import read_data, sir_model


plt.rcParams['figure.facecolor'] = 'white'

# Read the data (cumulative number of infected)
data = read_data('./data/rivm_corona_in_nl_daily.txt')  # Data is only till May 25,2020
N = 17000000.0
nr_days = 187  # This simulation goes all the way till September 1, 2020
shift = 7
nr_days_shift = nr_days + shift

# The initial values of the SIR (system)
S0 = N - 1.0
I0 = 1.0
R0 = 0.0
gamma = 1.0 / 9.0  # the parameter gamma that is assumed to be a constant

# Extracting the estimated time-dependent beta parameters
estimated_sol_sir = sir_model(nr_days_shift, 5.00, 1.02, 0.08, 21.7)
I_sol = estimated_sol_sir[2][shift:]  # This is the simulation without Kalman filter applied
R0_values = estimated_sol_sir[4]


# Time-dependent beta assumed to be known, thus not estimated with the Kalman filter
beta = [i * gamma for i in R0_values[shift:]]

# The error covariance matrix
P = np.array([[0.01, 0, 0],
              [0, 0.01, 0],
              [0, 0, 0.01]])

# The process noise covariance
Q = np.array([[0.1, 0, 0],
              [0, 0.1, 0],
              [0, 0, 0.1]])

# The measurement noise covariance
R = np.array([[0.1, 0, 0],
              [0, 0.1, 0],
              [0, 0, 0.1]])

# Observation matrix H
H = np.array([[1.0, 0, 0],
              [0, 1.0, 0],
              [0, 0, 1.0]])

# The initial values at k = 0
S_kf_est = [S0]
I_kf_est = [I0]
R_kf_est = [R0]

t_length = np.linspace(0, nr_days - 1, nr_days)

# For-loop over the simulation's time-span
# This approach uses the standard Kalman filter
# The discrete-time SIR model was computed
# the time step-size was chosen to be dt = 1.0
# However, the parameter [beta(k)] is assumed to be known beforehand and not estimated with the Kalman filter


for k in range(0, len(t_length) - 1):
    # state-transition matrix A withe time-dependent K

    A = np.array([[(1 - (beta[k] * I_kf_est[k]) / (2 * N)), -(beta[k] * S_kf_est[k]) / (2 * N), 0.0],
                  [(beta[k] * I_kf_est[k]) / (2 * N), 1 + (beta[k] * S_kf_est[k]) / (2 * N) - gamma, 0.0],
                  [0.0, gamma, 1.0]], dtype=float)

    # The past matrix x_(k)
    x_past = np.array([[S_kf_est[k]],
                       [I_kf_est[k]],
                       [R_kf_est[k]]], dtype=float)

    # Prediction step x-_(k+1) = A * x_(k)
    x_priori = np.dot(A, x_past)
    P = np.dot(np.dot(A, P), A.T) + Q

    # Measurement update equations a.k.a the Kalman Gain calculation
    S = np.dot(np.dot(H, P), H.T) + R
    S.astype(float)
    K = np.dot(np.dot(P, H.T), np.linalg.inv(S))

    # Observed measurements matrix

    if k < len(data) - 1:
        Y = np.array([[x_priori[0][0]],
                      [data[k + 1]],
                      [x_priori[2][0]]])
    else:
        Y = np.array([[x_priori[0][0]],
                      [x_priori[1][0]],
                      [x_priori[2][0]]])

    x_diff = Y - np.dot(H, x_priori)

    # x_(k+1) , a posteriori
    x_post = x_priori + np.dot(K, x_diff)

    # Append in a list with the estimations
    S_kf_est.append(x_post[0][0])
    I_kf_est.append(x_post[1][0])
    R_kf_est.append(x_post[2][0])

# Plot the Figure of the simulation's output
plt.figure(figsize=(10, 8))
plt.plot(data, "o", markersize=8, color='black')
plt.plot(I_sol, lw=3, color='blue')
plt.plot(I_kf_est, '--', lw=3, color='darkorange')
plt.legend(["Observed measurements", "Simulation without Kalman filter", "Kalman filtered simulation"], loc='best')
plt.xlabel("Time (days)" + "\n\n" + "from February 27, 2020", fontsize=14)
plt.ylabel("Number of Infected", fontsize=14)
plt.title("'Infected' from SIR model: Kalman filter vs without Kalman filter", fontsize=14, weight='bold')
plt.show()
