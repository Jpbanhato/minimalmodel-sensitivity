import os
import numpy as np
import matplotlib.pyplot as plt
import ap_features as apf


def plot_beats(beats, s2_interval, output_dir, iteration):
    subfolder = os.path.join(output_dir, 'it_' + str(iteration))
    os.makedirs(subfolder, exist_ok=True)
    fig, ax = plt.subplots()
    for beat in beats:
        ax.plot(beat.t, beat.y, label='Beat ' + str(beats.index(beat) + 1))
        # ax.legend()
    ax.set_title("S2 = " + str(s2_interval))
    filename = f'beats_plot_{s2_interval}.png'
    filepath = os.path.join(subfolder, filename)
    fig.savefig(filepath)
    plt.close(fig)

def plot_curves(beat1_di, beat1_apd, derivada, output_dir, print_flag, iteration):
    subfolder = os.path.join(output_dir, 'it_' + str(iteration))
    os.makedirs(subfolder, exist_ok=True)
    fig, ax = plt.subplots()
    filename = f'beats_plot_apd.png'
    ax.plot(beat1_di, beat1_apd, linewidth=1, marker='x', label='APD x DI')
    filepath = os.path.join(subfolder, filename)
    fig.savefig(filepath)
    plt.close(fig)
    
    if (print_flag):
        print(f'APD90: {last_beat_apd}')
        print(f'DI: {last_beat_di}')
        print(f'Max dVdt: {derivada}')

    fig, ax = plt.subplots()
    filename = f'beats_plot_derivada.png'
    ax.plot(beat1_di, derivada, color='red', label='Derivada')
    ax.legend()
    filepath = os.path.join(subfolder, filename)
    fig.savefig(filepath)
    plt.close(fig)