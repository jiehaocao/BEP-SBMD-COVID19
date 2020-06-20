"""Author: J.H.Cao
   Computational Biology group
   Biomedical Engineering
   Eindhoven University of Technology"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from kf_filter_seird import kf_filtering_seird, plot_seirdmodel, seird_model_test
from estimating_compmod_params import read_data, seird_model

from scipy.integrate import odeint

plt.rcParams['figure.facecolor'] = 'white'

# Generate some random data that is considered as "Observed measurements"
data = read_data('./data/corona_deaths.csv')
N = 17000000.0
nr_days = 217
shift = 30
nr_days_shift = nr_days + shift

# The parameters of the SEIRD model
beta = 0.33
gamma = 1.0 / 8.0
delta = 1.0 / 4.0
alpha = 1.0 / 18.0
p_i_d = 0.04

beta_2 = 0.23
gamma_2 = 1.0 / 9.0
delta_2 = 1.0 / 3.0
alpha_2 = 1.0 / 20.0
p_i_d_2 = 0.02

# The initial values of the SEIRD (system)
S0 = N - 1.0
E0 = 1.0
I0 = 0.0
R0 = 0.0
D0 = 0.0

y_0 = S0, E0, I0, R0, D0
t_length = np.linspace(0, nr_days - 1, nr_days)

# ========================================================================= #

# The solutions when solving the system of differential equations
seird_system_solution_obs = odeint(seird_model_test, y_0, t_length,
                                   args=(N, beta, delta, gamma, alpha, p_i_d))  # Artificial observed values
seird_system_solution_v2 = odeint(seird_model_test, y_0, t_length,
                                  args=(N, beta_2, delta_2, gamma_2, alpha_2, p_i_d_2))  # Model simulation

X_S, X_E, X_I, X_R, X_D = seird_system_solution_obs.T  # Artificial observed values
X_S_2, X_E_2, X_I_2, X_R_2, X_2_D = seird_system_solution_v2.T  # Model's output (simulation)


# When using the EKF, the seird_model or sir_model output, should be used as the h(x_k)
# and the real data would be the observed measurements. However, if there are no
# real observed measurements, then you should use sir_model of seird_model output as Yk

estimated_sol_seird = seird_model(nr_days_shift, 5.00, 0.52, 0.20, 57.2, 0.02)

S_fit, E_fit, I_fit = estimated_sol_seird[1], estimated_sol_seird[2], estimated_sol_seird[3]
R_fit, D_fit, R0_values = estimated_sol_seird[4], estimated_sol_seird[5], estimated_sol_seird[6]


estimated_sol_seird_shifted = np.array([S_fit[shift:],
                                        E_fit[shift:],
                                        I_fit[shift:],
                                        R_fit[shift:],
                                        D_fit[shift:]])

estimated_sol_seird_obs = estimated_sol_seird_shifted.T.copy()
seird_ode_sol = estimated_sol_seird_shifted.T.copy()

beta_values = []
for i in R0_values[30:]:
    beta_values.append(i * gamma)

# ======================================================================= #

# ====================== Kalman Filter application ====================== #

# Initial state-space
X_pred = np.array([[S0],
                   [E0],
                   [I0],
                   [R0],
                   [D0]])

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

# Observation matrix H, which is the Jacobian of the measurements equations

H = np.array([[1.0, 0, 0, 0, 0],
              [0, 1.0, 0, 0, 0],
              [0, 0, 1.0, 0, 0],
              [0, 0, 0, 1.0, 0],
              [0, 0, 0, 0, 1.0]])

S_est, E_est, I_est, R_est, D_est = kf_filtering_seird(t_length, X_pred, estimated_sol_seird_obs, seird_ode_sol,
                                                       beta_values, data, P, H, Q, R, real_data=True)

plt.plot(data, 'o', color='blue')
plt.plot(D_fit[30:], '--', color='darkorange')
plt.plot(D_est, '*', color='darkgreen')
plt.legend(["Observed measurements", "Simulation without filtering", "Kalman filtered"])
plt.show()
