import numpy as np
import matplotlib.pyplot as plt
import ap_features as apf
import os
from solve_min_model import *
from parameters import *
from plot import *
from sampling import *

def stim(t,BCL):    
    S1_interval = BCL
    S2_interval = BCL
    num_S1 = 8
    duration=0.1
    if t < S1_interval * num_S1:
        Jstim = 1.0 if (t % S1_interval) < duration else 0.0
    else:
        Jstim = 1.0 if (t - S1_interval * num_S1) < duration else 0.0
    return Jstim

output_dir = 'plots'
os.makedirs(output_dir, exist_ok=True)
subfolder_beats = os.path.join(output_dir, 'beats')
os.makedirs(subfolder_beats, exist_ok=True)
subfolder_curves = os.path.join(output_dir, 'curves')
os.makedirs(subfolder_curves, exist_ok=True)

# run the simulation for model parameters and other specified variables -> it may take special treatment for different number of dimensions
dt, BCL, num_s1, interval, stim_start, s2_min, s2_dec, print_flag, model_name = 0.1, 400, 8, 1, 50.0, 300, 10, False, 'tnnp-mayra'
s2_interval_array = np.arange(BCL, s2_min, -s2_dec)
print(f'Running for {model_name}')

params_list, nomes_variaveis = latin_hypercube(define_parameters(model_name), 0.9, 1.1)

iteration = 0
for params_to_dict in params_list:
    params = convert_to_dict(params_to_dict, nomes_variaveis)
    print(f'Iteration {iteration}')
    derivada = []
    last_beat_apd = []
    last_beat_di = []
    for s2_interval in s2_interval_array:
        print(f'S2 = {s2_interval}')
        voltage, time = solve_min_model_0D(dt, BCL, num_s1, s2_interval, interval, stim_start, params, print_flag)
        
        # trace = apf.Beats(y=voltage, t=time)

        ST=[stim(T,BCL)for T in time]
        trace = apf.Beats(y=voltage, t=time,pacing=ST)
          

        beats = trace.beats
        plot_beats(beats, s2_interval, subfolder_beats, iteration)

        # calculate the APD90, DI and maxdVdt
        last_beat = beats[-1]
        last_beat_before = beats[-2]
        last_beat_apd.append(last_beat.apd(90))
        last_beat_di.append(s2_interval - last_beat_before.apd(90))
        dV = (abs(np.diff(voltage) / dt))
        derivada.append(np.max(dV))

    plot_curves(last_beat_di, last_beat_apd, derivada, subfolder_curves, print_flag, iteration)
    iteration += 1