"""Author: J.H.Cao
   Computational Biology group
   Biomedical Engineering
   Eindhoven University of Technology"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from scipy.integrate import odeint

plt.rcParams['figure.facecolor'] = 'white'


def kf_filtering_seird(t_length, X_pred, y_observed, sim_sol, beta_values, data, P, H, Q, R, real_data=False):
    """A function for the (extended) Kalman filter applied on
       the SEIRD model (Susceptibles-Exposed-Infected-Recovered-Deceased)"""

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
    delta = 1.0 / 3.0
    alpha = 1.0 / 20.0
    p_i_d = 0.02

    # Predictions or estimations (a posteriori)
    S_pred = []
    E_pred = []
    I_pred = []
    R_pred = []
    D_pred = []

    # The observed measurements
    S_obs, E_obs, I_obs, R_obs, D_obs = y_observed.T

    # The ODE solution of the SIR system
    S_sol, E_sol, I_sol, R_sol, D_sol = sim_sol.T

    # Argument if there are real observed measurements available
    if real_data:
        D_obs[:len(data)] = data
        beta = beta_values.copy()

    # The initial values of the state estimations
    S_pred.append(X_pred[0][0])
    E_pred.append(X_pred[1][0])
    I_pred.append(X_pred[2][0])
    R_pred.append(X_pred[3][0])
    D_pred.append(X_pred[4][0])

    # The Kalman Filter algorithm:

    for k in range(0, len(t_length) - 1):
        # Generate the transition matrix, which is the Jacobian of the state

        if real_data:
            # A = np.array([[(-beta[k] * I_pred[k]) / N, 0.0, (-beta[k] * S_pred[k]) / N, 0.0, 0.0],
            #               [(beta[k] * I_pred[k]) / N, -delta, (beta[k] * S_pred[k]) / N, 0.0, 0.0],
            #               [0.0, delta, (-gamma * (1-p_i_d)-(alpha*p_i_d)), 0.0, 0.0],
            #               [0.0, 0.0, (gamma * (1-p_i_d)), 0.0, 0.0],
            #               [0.0, 0.0, (alpha * p_i_d), 0.0, 0.0]])

            A = np.array([[(-beta[k + 1] * I_sol[k + 1]) / N, 0.0, (-beta[k + 1] * S_sol[k + 1]) / N, 0.0, 0.0],
                          [(beta[k + 1] * I_sol[k + 1]) / N, -delta, (beta[k + 1] * S_sol[k + 1]) / N, 0.0, 0.0],
                          [0.0, delta, (-gamma * (1 - p_i_d) - (alpha * p_i_d)), 0.0, 0.0],
                          [0.0, 0.0, (gamma * (1 - p_i_d)), 0.0, 0.0],
                          [0.0, 0.0, (alpha * p_i_d), 0.0, 0.0]])

        else:
            A = np.array([[(-beta * I_pred[k]) / N, 0.0, (-beta * S_pred[k]) / N, 0.0, 0.0],
                          [(beta * I_pred[k]) / N, -delta, (beta * S_pred[k]) / N, 0.0, 0.0],
                          [0.0, delta, (-gamma * (1 - p_i_d) - (alpha * p_i_d)), 0.0, 0.0],
                          [0.0, 0.0, (gamma * (1 - p_i_d)), 0.0, 0.0],
                          [0.0, 0.0, (alpha * p_i_d), 0.0, 0.0]])

        X_state = np.array([[S_pred[k]],
                            [E_pred[k]],
                            [I_pred[k]],
                            [R_pred[k]],
                            [D_pred[k]]], dtype=float)

        # Prediction step of the Kalman filter
        # X_priori = np.dot(A, X_state)

        X_priori = np.array([[S_sol[k + 1]],
                             [E_sol[k + 1]],
                             [I_sol[k + 1]],
                             [R_sol[k + 1]],
                             [D_sol[k + 1]]], dtype=float)

        P = np.dot(np.dot(A, P), A.T) + Q

        # Measurement update equations a.k.a the Kalman Gain calculation
        S = np.dot(np.dot(H, P), H.T) + R
        S.astype(float)
        K = np.dot(np.dot(P, H.T), np.linalg.inv(S))

        # Update to a posteriori state

        if k > len(data):
            Y_measure = np.array([[S_pred[-1]],
                                  [E_pred[-1]],
                                  [I_pred[-1]],
                                  [R_pred[-1]],
                                  [D_pred[-1]]], dtype=float)
        else:
            Y_measure = np.array([[S_obs[k+1]],
                                  [E_obs[k+1]],
                                  [I_obs[k+1]],
                                  [R_obs[k+1]],
                                  [D_obs[k+1]]], dtype=float)


        # Y_measure = np.array([[0.0],
        #                       [0.0],
        #                       [0.0],
        #                       [0.0],
        #                       [D_obs[k+1]]], dtype=float)



        # X_diff = Y_measure - np.dot(H, X_priori)
        X_diff = Y_measure - X_priori
        X_post = X_priori + np.dot(K, X_diff)
        P = np.dot((np.eye(5) - np.dot(K, H)), P)

        S_pred.append(X_post[0][0])
        E_pred.append(X_post[1][0])
        I_pred.append(X_post[2][0])
        R_pred.append(X_post[3][0])
        D_pred.append(X_post[4][0])

    return S_pred, E_pred, I_pred, R_pred, D_pred


def plot_seirdmodel(t, S, E, I, R, D):
    f, ax = plt.subplots(1, 1, figsize=(10, 8))
    ax.plot(t, S, 'b', alpha=0.7, linewidth=2, label='Susceptible')
    ax.plot(t, E, 'y', alpha=0.7, linewidth=2, label='Exposed')
    ax.plot(t, I, 'r', alpha=0.7, linewidth=2, label='Infected')
    ax.plot(t, R, 'g', alpha=0.7, linewidth=2, label='Recovered')
    ax.plot(t, D, 'k', alpha=0.7, linewidth=2, label='Deceased')

    ax.set_xlabel('Time [days]', fontsize=15)
    ax.set_ylabel("Number of people", fontsize=15)
    ax.set_title("SEIRD (Susceptibles-Exposed-Infected-Recovered-Deceased", fontsize=20, weight='bold')

    ax.yaxis.set_tick_params(length=0)
    ax.xaxis.set_tick_params(length=0)
    ax.grid(b=True, which='major', c='w', lw=2, ls='-')
    # ax.spines['top'].set_visible(False)
    # ax.spines['right'].set_visible(False)

    legend = ax.legend()
    legend.get_frame().set_alpha(1.0)

    plt.show()


def seird_model_test(y, t, N, beta, delta, gamma, alpha, prob_I_to_D):
    S, E, I, R, D = y

    dS_dt = -beta * I * (S / N)
    dE_dt = beta * I * (S / N) - delta * E
    dI_dt = delta * E - gamma * (1 - prob_I_to_D) * I - alpha * prob_I_to_D * I
    dR_dt = gamma * (1 - prob_I_to_D) * I
    dD_dt = alpha * prob_I_to_D * I

    return dS_dt, dE_dt, dI_dt, dR_dt, dD_dt
