"""Author: J.H.Cao
   Computational Biology group
   Biomedical Engineering
   Eindhoven University of Technology"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from scipy.integrate import odeint

plt.rcParams['figure.facecolor'] = 'white'


def kf_filtering_sir(t_length, X_pred, y_observed, sim_sol, beta_values, data, P, H, Q, R, real_data=False):
    """A function for the (extended) Kalman filter applied on
       the SIR model (Susceptibles-Infected-Recovered)"""

    # ================================================== #
    # Input:
    # t_length = the length of the estimations (in this case is the time in [days])
    # x_pred = state-matrix (X_pred would be the initial)
    # Y-pred = observed measurements matrix
    # P = error covariance matrix with respect to X_pred
    # H = measurement matrix
    # Q = process noise covariance
    # R = measurement noise covariance

    # Output:
    # Estimated (a posteriori) X_pred
    # ================================================= #

    # Estimated parameters
    beta = 0.23
    N = 17000000
    gamma = 1.0 / 9.0

    # Predictions or estimations (a posteriori)
    S_pred = []
    I_pred = []
    R_pred = []

    # The observed measurements
    S_obs, I_obs, R_obs = y_observed.T

    # The ODE solution of the SIR system
    S_sol, I_sol, R_sol = sim_sol.T

    # Argument if there are real observed measurements available
    if real_data:
        I_obs[:len(data)] = data
        beta = beta_values.copy()

    # The initial values of the state estimations
    S_pred.append(X_pred[0][0])
    I_pred.append(X_pred[1][0])
    R_pred.append(X_pred[2][0])

    # The Kalman Filter algorithm:

    for k in range(0, len(t_length) - 1):
        # Generate the transition matrix, which is the Jacobian of the state
        if real_data:
            # A = np.array([[(-beta[k] * I_pred[k]) / N, (-beta[k] * S_pred[k]) / N, 0],
            #               [(beta[k] * I_pred[k]) / N, (beta[k] * S_pred[k]) / N - gamma, 0],
            #               [0, gamma, 0]], dtype=float)

            A = np.array([[(-beta[k + 1] * I_sol[k + 1]) / N, (-beta[k + 1] * S_sol[k + 1]) / N, 0],
                          [(beta[k + 1] * I_sol[k + 1]) / N, (beta[k + 1] * S_sol[k + 1]) / N - gamma, 0],
                          [0, gamma, 0]], dtype=float)
        else:

            A = np.array([[(-beta * I_pred[k]) / N, (-beta * S_pred[k]) / N, 0],
                          [(beta * I_pred[k]) / N, (beta * S_pred[k]) / N - gamma, 0],
                          [0, gamma, 0]], dtype=float)

        X_priori = np.array([[S_sol[k + 1]],
                             [I_sol[k + 1]],
                             [R_sol[k + 1]]], dtype=float)

        P = np.dot(np.dot(A, P), A.T) + Q

        # Measurement update equations a.k.a the Kalman Gain calculation
        S = np.dot(np.dot(H, P), H.T) + R
        S.astype(float)
        K = np.dot(np.dot(P, H.T), np.linalg.inv(S))

        # Update to a posteriori state

        if k > len(data):
            Y_measure = np.array([[S_pred[-1]],
                                  [I_pred[-1]],
                                  [R_pred[-1]]], dtype=float)
        else:

            Y_measure = np.array([[S_obs[k + 1]],
                                  [I_obs[k + 1]],
                                  [R_obs[k + 1]]], dtype=float)

        # Y_measure = np.array([[0.0],
        #                       [I_obs[k+1]],
        #                       [0.0]], dtype=float)

        # X_diff = Y_measure - np.dot(H, X_priori)
        X_diff = Y_measure - X_priori

        X_post = X_priori + np.dot(K, X_diff)
        P = np.dot((np.eye(3) - np.dot(K, H)), P)

        S_pred.append(X_post[0][0])
        I_pred.append(X_post[1][0])
        R_pred.append(X_post[2][0])

    return S_pred, I_pred, R_pred


def plot_model(t, S, I, R):
    f, ax = plt.subplots(1, 1, figsize=(10, 8))
    ax.plot(t, S, 'b', alpha=0.7, linewidth=2, label='Susceptibles')
    ax.plot(t, I, 'r', alpha=0.7, linewidth=2, label='Infected')
    ax.plot(t, R, 'g', alpha=0.7, linewidth=2, label='Recovered')

    ax.set_xlabel('Time [days]', fontsize=15)
    ax.set_ylabel("Number of People", fontsize=15)
    ax.set_title("SIR (Susceptibles - Infected - Recovered)", fontsize=20, weight='bold')

    ax.yaxis.set_tick_params(length=0)
    ax.xaxis.set_tick_params(length=0)
    ax.grid(b=True, which='major', c='w', lw=2, ls='-')
    # ax.spines['top'].set_visible(False)
    # ax.spines['right'].set_visible(False)

    legend = ax.legend()
    legend.get_frame().set_alpha(1.0)

    plt.show()


def sir_model_test(y, t, N, beta, gamma):
    """ A system of differential equations for the
    SIR (Susceptibles - Infected - Recovered) model"""

    S, I, R = y

    dS_dt = (-beta * I * S) / N
    dI_dt = (beta * I * S) / N - gamma * I
    dR_dt = gamma * I

    return dS_dt, dI_dt, dR_dt
