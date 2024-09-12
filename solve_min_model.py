import numpy as np
from functions_rl import *

# solve the minimal model in 0D using rush larsen, it returns [voltage, time, ...]

def solve_min_model_0D(dt, BCL, num_s1, s2_interval, interval, stim_start, params, print_flag):
    # Time setup
    Num_pts = int((stim_start + num_s1 * BCL + 2 * BCL) / dt)
    t = np.linspace(0, Num_pts*dt, Num_pts+1)
    t_stim = stim_start

    # State Variables
    u = np.zeros(Num_pts+1)
    v = np.zeros(Num_pts+1)
    w = np.zeros(Num_pts+1)
    s = np.zeros(Num_pts+1)

    # Initial conditions
    u[0] = 0.0
    v[0] = 1.0
    w[0] = 1.0
    s[0] = 0.0


    # Loop over time
    for i in range(Num_pts):
        u[i+1] = u[i] + dt * (- (J_fi(u[i], v[i], params) + J_so(u[i], params) + J_si(u[i], w[i], s[i], params)) + I_stim(t[i], t_stim, interval))
        if (tau_v_rl(u[i], params) > 10**-10):
            v[i+1] = v_inf_rl(u[i], params) + (v[i] - v_inf_rl(u[i], params)) * np.exp(-dt / tau_v_rl(u[i], params))
        else:
            v[i+1] = v[i] + dt * ((1.0 - H(u[i] - params["theta_v"])) * (v_inf(u[i], params) - v[i]) / tau_v_minus(u[i], params) - H(u[i] - params["theta_v"]) * v[i] / params["tau_v_plus"])
        if (tau_w_rl(u[i], params) > 10**-10):
            w[i+1] = w_inf_rl(u[i], params) + (w[i] - w_inf_rl(u[i], params)) * np.exp(-dt / tau_w_rl(u[i], params))
        else:
            w[i+1] = w[i] + dt * ((1.0 - H(u[i] - params["theta_w"])) * (w_inf(u[i], params) - w[i]) / tau_w_minus(u[i], params) - H(u[i] - params["theta_w"]) * w[i] / params["tau_w_plus"])
        if (tau_s(u[i], params) > 10**-10):
            s[i+1] = s_inf_rl(u[i], params) + (s[i] - s_inf_rl(u[i], params)) * np.exp(-dt / tau_s_rl(u[i], params))
        else:
            s[i+1] = s[i] + dt * (((1.0 + np.tanh(params["k_s"] * (u[i] - params["u_s"])))) / 2 - s[i]) / tau_s(u[i], params)

        # if time instant is PRE BEAT or S1 -> t_stim follows BCL
        # if time instant is S2 -> t_stim follows BCL - s2_dec until s2_min is reached
        if (t[i] % (BCL + t_stim) == 0 and t[i] != 0 and t[i] < (stim_start + (BCL * (num_s1)))):
            if (print_flag):
                print(f't[i]: {t[i]}, t_stim: {t_stim}, BCL: {BCL}')
            t_stim = t[i]
        elif ((t[i] % (stim_start + (BCL) * (num_s1 - 1) + s2_interval) == 0) and t[i] != 0):
            if (print_flag):
                print(f't[i]: {t[i]}, t_stim: {t_stim}, s2_interval: {s2_interval}')
            t_stim = t[i]
    
    return [u, t]