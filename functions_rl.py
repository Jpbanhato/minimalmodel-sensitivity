import numpy as np

# define the minimal model functions, adapted to rush larsen and h-h format

# Stimulus function -> for denormalized, use a magnitude parameter to multiply the return
def I_stim(t, t_stim, interval):
    if ((t_stim < t < t_stim + interval)):
        return 1.0
    else:
        return 0.0

# Standard Heaviside function
def H(x):
  if (x > 0.0):
    return 1.0
  else:
    return 0.0

# Functions for infinity values
def v_inf(u, params):
  if (u < params["theta_v_minus"]):
    return 1.0
  else:
    return 0.0

def w_inf(u, params):
  return (1.0 - H(u - params["theta_o"])) * (1.0 - u/params["tau_w_inf"]) + H(u - params["theta_o"]) * params["w_inf_star"]

# Functions for time constants
def tau_v_minus(u, params):
  return (1.0 - H(u - params["theta_v_minus"])) * params["tau_v1_minus"] + H(u - params["theta_v_minus"]) * params["tau_v2_minus"]

def tau_w_minus(u, params):
  tau_w1_minus = params["tau_w1_minus"]
  tau_w2_minus = params["tau_w2_minus"]
  k_w_minus = params["k_w_minus"]
  u_w_minus = params["u_w_minus"]
  return tau_w1_minus + (tau_w2_minus - tau_w1_minus) * (1.0 + np.tanh(k_w_minus * (u - u_w_minus))) / 2

def tau_so(u, params):
  k_so = params["k_so"]
  u_so = params["u_so"]
  tau_so1 = params["tau_so1"]
  tau_so2 = params["tau_so2"]
  return tau_so1 + (tau_so2 - tau_so1) * (1.0 + np.tanh(k_so * (u - u_so))) / 2

def tau_s(u, params):
  k_s = params["k_s"]
  u_s = params["u_s"]
  tau_s1 = params["tau_s1"]
  tau_s2 = params["tau_s2"]
  theta_w = params["theta_w"]
  return (1.0 - H(u - theta_w)) * tau_s1 + H(u - theta_w) * tau_s2

def tau_o(u, params):
  tau_o1 = params["tau_o1"]
  tau_o2 = params["tau_o2"]
  theta_o = params["theta_o"]
  return (1.0 - H(u - theta_o)) * tau_o1 + H(u - theta_o) * tau_o2

# Functions for currents

# novo metodo
def J_fi(u, v, params):
  theta_v = params["theta_v"]
  u_u = params["u_u"]
  tau_fi = params["tau_fi"]
  return -v * H(u - theta_v) * (u - theta_v) * (u_u - u) / tau_fi

def J_so(u, params):
  u_o = params["u_o"]
  theta_w = params["theta_w"]
  return (u - u_o) * (1.0 - H(u - theta_w)) / tau_o(u, params) + H(u - theta_w) / tau_so(u, params)

def J_si(u, w, s, params):
  theta_w = params["theta_w"]
  tau_si = params["tau_si"]
  return -H(u - theta_w) * w * s / tau_si

# RUSH LARSEN IMPLEMENTATION:

# s:
def tau_s_rl(u, params):
  k_s = params["k_s"]
  u_s = params["u_s"]
  theta_w = params["theta_w"]
  tau_s1 = params["tau_s1"]
  tau_s2 = params["tau_s2"]
  return (1.0 - H(u - theta_w)) * tau_s1 + H(u - theta_w) * tau_s2

def s_inf_rl(u, params):
  theta_w = params["theta_w"]
  u_s = params["u_s"]
  k_s = params["k_s"]
  return (1.0 + np.tanh(k_s * (u - u_s))) / 2

# v:
def tau_v_rl(u, params):
  theta_v = params["theta_v"]
  theta_v_minus = params["theta_v_minus"]
  tau_v1_minus = params["tau_v1_minus"]
  tau_v2_minus = params["tau_v2_minus"]
  tau_v_plus = params["tau_v_plus"]
  return (tau_v_plus * tau_v_minus(u, params)) / (tau_v_plus - tau_v_plus * H(u - theta_v) + tau_v_minus(u, params) * H(u - theta_v))

def v_inf_rl(u, params):
  theta_v = params["theta_v"]
  theta_v_minus = params["theta_v_minus"]
  tau_v_plus = params["tau_v_plus"]
  return (tau_v_plus * v_inf(u, params) * (1 - H(u - theta_v))) / (tau_v_plus - tau_v_plus * H(u - theta_v) + tau_v_minus(u, params) * H(u - theta_v))

# w:
def tau_w_rl(u, params):
  theta_w = params["theta_w"]
  theta_v = params["theta_v"]
  tau_w1_minus = params["tau_w1_minus"]
  tau_w2_minus = params["tau_w2_minus"]
  k_w_minus = params["k_w_minus"]
  u_w_minus = params["u_w_minus"]
  tau_w_plus = params["tau_w_plus"]
  return (tau_w_plus * tau_w_minus(u, params)) / (tau_w_plus - tau_w_plus * H(u - theta_w) + tau_w_minus(u, params) * H(u - theta_w))

def w_inf_rl(u, params):
  theta_w = params["theta_w"]
  theta_v = params["theta_v"]
  tau_w_plus = params["tau_w_plus"]
  return (tau_w_plus * w_inf(u, params) * (1 - H(u - theta_w))) / (tau_w_plus - tau_w_plus * H(u - theta_w) + tau_w_minus(u, params) * H(u - theta_w))