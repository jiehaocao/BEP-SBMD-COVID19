"""Author: J.H.Cao
   Computational Biology group
   Biomedical Engineering
   Eindhoven University of Technology"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from kf_filter_sir import kf_filtering_sir, plot_model, sir_model_test
from estimating_compmod_params import read_data, sir_model

from scipy.integrate import odeint

plt.rcParams['figure.facecolor'] = 'white'

# Generate some random data that is considered as "Observed measurements"
data = read_data('./data/rivm_corona_in_nl_daily.txt')
N = 17000000.0
nr_days = 120
shift = 7
nr_days_shift = nr_days + shift

# The parameters of the SIR model
beta = 0.30
beta_2 = 0.23
gamma = 1.0 / 7.5
gamma_2 = 1.0 / 9.0

# The initial values of the SIR (system)
S0 = N - 1.0
I0 = 1.0
R0 = 0.0

y_0 = S0, I0, R0
t_length = np.linspace(0, nr_days - 1, nr_days)

# ========================================================================= #

# The solutions when solving the system of differential equations
sir_system_solution_obs = odeint(sir_model_test, y_0, t_length,
                                 args=(N, beta, gamma))  # Artificial observed values
sir_system_solution_v2 = odeint(sir_model_test, y_0, t_length,
                                args=(N, beta_2, gamma_2))  # Model simulation

X_S, X_I, X_R = sir_system_solution_obs.T  # Artificial observed values
X_S_2, X_I_2, X_R_2 = sir_system_solution_v2.T  # Model's output (simulation)

# When using the EKF, the seird_model or sir_model output, should be used as the h(x_k)
# and the real data would be the observed measurements. However, if there are no
# real observed measurements, then you should use sir_model of seird_model output as Yk

estimated_sol_sir = sir_model(nr_days_shift, 5.00, 1.02, 0.08, 21.7)

S_fit, I_fit, R_fit = estimated_sol_sir[1], estimated_sol_sir[2], estimated_sol_sir[3]
R0_values = estimated_sol_sir[4]

estimated_sol_sir_shifted = np.array([S_fit[shift:],
                                      I_fit[shift:],
                                      R_fit[shift:]])

estimated_sol_sir_obs = estimated_sol_sir_shifted.T.copy()
sir_ode_sol = estimated_sol_sir_shifted.T.copy()


beta_values = []
for i in R0_values[7:]:
    beta_values.append(i * gamma)

# ======================================================================= #

# ====================== Kalman Filter application ====================== #

# Initial state-space
X_pred = np.array([[S0],
                   [I0],
                   [R0]])

# The error covariance matrix

P = np.array([[0.01, 0, 0],
              [0, 0.01, 0],
              [0, 0, 0.01]])

# The process noise covariance

Q = np.array([[0.1, 0, 0],
              [0, 0.1, 0],
              [0, 0, 0.1]])

# The measurement noise covariance

R = np.array([[0.01, 0, 0],
              [0, 0.0, 0],
              [0, 0, 0.01]])

# Observation matrix H, which is the Jacobian of the measurements equations

H = np.array([[1.0, 0, 0],
              [0, 1.0, 0],
              [0, 0, 1.0]])

S_est, I_est, R_est = kf_filtering_sir(t_length, X_pred, estimated_sol_sir_obs, sir_ode_sol, beta_values,
                                       data, P, H, Q, R, real_data=True)

plt.plot(data, 'o', color='blue')
plt.plot(I_fit[7:], '--', color='darkorange')
plt.plot(I_est, '*', color='darkgreen')
plt.legend(["Observed measurements", "Simulation without filtering", "Kalman filtered"])
plt.show()
