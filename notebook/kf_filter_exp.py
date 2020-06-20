"""Author: J.H.Cao
   Computational Biology group
   Biomedical Engineering
   Eindhoven University of Technology"""
from typing import List, Any

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy.optimize import curve_fit

plt.rcParams['figure.facecolor'] = 'white'


def read_data(filename):
    """Read the csv file and return a dataframe"""

    df_covid = pd.read_csv(filename, parse_dates=['Datum']).set_index("Datum")
    y_data = df_covid["Aantal"].values

    return y_data


def kf_filtering_exp(t_length, X_pred, y_observed, A, P, H, Q, R):
    """A function for the basic Kalman filter applied on
       the simple exponential growth system"""

    x_post = np.zeros(len(y_observed)+1)
    x_post[0] = y_observed[0]
    x_priori = X_pred.copy()


    for k in range(0, len(t_length) - 1):
        # Prediction step
        P = A * P * A.T + Q

        # Calculating the Kalman Gain
        S = H * P * H.T + R
        K = P * H.T * np.linalg.inv(S)

        # update the a priori state
        x_diff = y_observed[k + 1] - H * x_priori[k + 1]
        x_post[k + 1] = x_priori[k + 1] + K * x_diff

        P = (np.eye(1) - K * H) * P

    return x_post


y_data = read_data('./data/rivm_corona_in_nl_daily.txt')

x = np.arange(len(y_data))
y = y_data

def exponential(t, a, b):
    return a * np.exp(b * t)


p0_exp = np.random.exponential(size=2)

bounds_exp = (0, [1000, 3])  # upper bounds for the given parameters (a, b, c)

(a, b), exp_cov = curve_fit(exponential, x, y, bounds=bounds_exp, p0=p0_exp)

a = 999.99
b = 0.05


def exp_pred(t):
    return a * np.exp(b * t)


A = np.array([[0.05]])
H = np.array([[1.0]])
P = np.array([[1.0]])
Q = np.array([[0.1]])
R = np.array([[0.1]])

X_pred = exp_pred(np.arange(len(x)+1))

x_post = kf_filtering_exp(x, X_pred, y_data, A, P, H, Q, R)
