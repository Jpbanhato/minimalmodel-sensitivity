# !pip install ap_features
# !pip install SALib
# !pip install scikit-learn

import numpy as np
import matplotlib.pylab as plt
import pylab
from math import *
import ap_features as apf
import random
import os
from SALib.sample import saltelli
from SALib.analyze import sobol

u_o = 0
u_u = 1.58
theta_v = u_o + 0.3 / u_u * (u_u - u_o)
theta_w =  u_o + 0.015 / u_u * (u_u - u_o)
theta_v_minus =  u_o + 0.015 / u_u * (u_u - u_o)
theta_o =  u_o + 0.006 / u_u * (u_u - u_o)
tau_v1_minus = 60.0
tau_v2_minus = 1150
tau_v_plus = 1.4506
tau_w1_minus = 70
tau_w2_minus = 20
k_w_minus = 65.0
u_w_minus = u_o + 0.03 / u_u * (u_u - u_o)
tau_w_plus = 280.0
tau_fi = 0.11
tau_o1 = 6
tau_o2 = 6.0
tau_so1 = 43
tau_so2 = 0.2
k_so = 2
u_so = 0.65
tau_s1 = 2.7342
tau_s2 = 3.0
k_s = 2.0994
u_s = u_o +  0.9087 / u_u * (u_u - u_o)
tau_si = 2.8723
tau_w_inf = 0.07
w_inf_star = 0.94

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

def min_model (BCL, params):
    S1_interval = BCL # S1 tem sempre BCL de BCL
    S2_interval = BCL
    BCL_fixed = BCL
    dt = 0.1
    num_S1 = 10
    T = (S1_interval * num_S1) + (BCL * 1) # vai sempre fazer num_S1 S1, depois 1 S2
    N = int(T/dt) + 1
    time= np.zeros(N);
    t=0
    t_ini = 0
    i=0
    voltage = np.zeros(N);
    duration = 1.

    u = 0.
    v = 1.
    w = 1.
    s = 0.

    BCL = S1_interval
    S2_not_done = True

    while(t<T):

        if t>=t_ini and t<=t_ini+duration and S2_not_done:
            Jstim = 1.
        else:
            Jstim = 0

        u = u + dt * (- (J_fi(u, v, params) + J_so(u, params) + J_si(u, w, s, params)) + Jstim)
        if (tau_v_rl(u, params) > 10**-10):
            v = v_inf_rl(u, params) + (v - v_inf_rl(u, params)) * np.exp(-dt / tau_v_rl(u, params))
        else:
            v = v + dt * ((1.0 - H(u - params["theta_v"])) * (v_inf(u, params) - v) / tau_v_minus(u, params) - H(u - params["theta_v"]) * v / params["tau_v_plus"])
        if (tau_w_rl(u, params) > 10**-10):
            w = w_inf_rl(u, params) + (w - w_inf_rl(u, params)) * np.exp(-dt / tau_w_rl(u, params))
        else:
            w = w + dt * ((1.0 - H(u - params["theta_w"])) * (w_inf(u, params) - w) / tau_w_minus(u, params) - H(u - params["theta_w"]) * w / params["tau_w_plus"])
        if (tau_s(u, params) > 10**-10):
            s = s_inf_rl(u, params) + (s - s_inf_rl(u, params)) * np.exp(-dt / tau_s_rl(u, params))
        else:
            s = s + dt * (((1.0 + np.tanh(params["k_s"] * (u - params["u_s"])))) / 2 - s) / tau_s(u, params)

        voltage[i] = u
        time[i] = t
        t = t+dt
        i+=1
        if t > S1_interval * (num_S1 - 1):
          BCL = S2_interval
          if t > S1_interval * num_S1 + 1:
            S2_not_done = False
        if t>t_ini+duration:
            t_ini += BCL
        # if BCL == S2_interval and BCL != S1_interval:
        #   print(BCL)
    # print(f"T = {S1_interval * num_S1} + {BCL * 1} = {S1_interval * num_S1 + BCL * 1} - t = {int(t)}")
    return [voltage, time]

def run_simulation(params):
    resultados = []
    i = 0
    while (i < 1):
        i += 1
        random_params = params
        BCL = 400
        BCL_max_limit = 400
        BCL_min_limit = 350
        number_it = int((BCL_max_limit - BCL_min_limit) / 5) + 1
        # print(number_it)
        j = 0
        derivada_cv = []
        beat1_apd = np.zeros(number_it)
        beat1_di = np.zeros(number_it)
        while(BCL > BCL_min_limit and j < number_it): # até BCL de BCL_min_limit
          # print(BCL)
          voltagem, time = min_model(BCL, random_params)
          trace = apf.Beats(y=voltagem, t=time)
          beats = trace.beats
          tam = len(beats)
          fig, ax = plt.subplots()
          for beat in beats:
              ax.plot(beat.t, beat.y, label='Beat ' + str(beats.index(beat) + 1))
              ax.legend()
          ax.set_title("BCL = " + str(BCL))
          plt.show()
          if tam == number_it:
              valid_simulation = True
              #print(f"Valid simulation. tam: {tam}, beats:{beats}")
              #print("Valid simulation")
              # break
          else:
              valid_simulation = False
              #print(f"Not enough beats {tam}. Re-running with new parameters...")
              break  # Sair do loop atual para gerar novos parâmetros
          # print(f'tam={tam}')
          # if (tam == number_it):
          beat1 = beats[-1]
          beat1_ant = beats[-2]
          # else:
              # beat1 = beats[-2]
              # beat1_ant = beats[-3]
          beat1_apd[j] = beat1.apd(90)
          beat1_di[j] = BCL - beat1_ant.apd(90)
          dt = 0.1
          dV = (abs(np.diff(voltagem) / dt))
          derivada_cv.append(np.max(dV))
          # BCL_array.append(BCL)
          # plt.plot(time[:-1], dV, linewidth=1, label='U x t')
          # plt.show()
          # print(f'(np.max(dV): {(np.max(dV)}')
          # print(f"BCL: {BCL} --- APD: {beat1_apd[j]} + DI: {beat1_di[j]} = {beat1_apd[j] + beat1_di[j]}")
          BCL -= 5  # reduzindo de 5 em 5
          j += 1
        if (not valid_simulation):
            #retornar se a simulação não for valida
            #agora os numeros sao fixos, o sorteio nao altera
            resultados.append(-inf)
            #i -= 1
            #continue
        else:
            # Calcular derivada de APD e CV
            derivada = np.diff(beat1_apd) / np.diff(beat1_di)
            beat1_apd = beat1_apd[:-1]
            beat1_di = beat1_di[:-1]
            # derivada = np.diff(beat1_apd) / np.diff(BCL_array)

            resultados.append(np.max(derivada))
            resultados.append(np.max(beat1_apd))
            resultados.append(np.min(beat1_apd))

            resultados.append(np.max(derivada_cv))
            resultados.append(np.max(derivada_cv))
            resultados.append(np.min(derivada_cv))

            # plt.plot(time, voltagem, linewidth=1, label='U x t')
            # plt.show()

            # plt.plot(beat1_di, beat1_apd, linewidth=1, marker='x', label='APD x DI')
            # plt.show()

            if (len(beat1_di) != len(derivada_cv)):
              derivada_cv = derivada_cv[:-1]
            # plt.plot(beat1_di, derivada_cv, color='red', label='Derivada_cv Numérica')
            # plt.legend()
            # plt.show()

    return resultados

from SALib.sample import latin

def get_parameters(params):
    # Nomes das variáveis correspondentes
    # nomes_variaveis = ['u_u', 'theta_v', 'tau_v_plus', 'tau_w_plus', 'tau_fi', 'tau_so1', 'tau_so2', 'tau_si']
    nomes_variaveis = [
        'u_o', 'u_u', 'theta_v', 'theta_w', 'theta_v_minus', 'theta_o',
        'tau_v1_minus', 'tau_v2_minus', 'tau_v_plus', 'tau_w1_minus',
        'tau_w2_minus', 'k_w_minus', 'u_w_minus', 'tau_w_plus', 'tau_fi',
        'tau_o1', 'tau_o2', 'tau_so1', 'tau_so2', 'k_so', 'u_so', 'tau_s1',
        'tau_s2', 'k_s', 'u_s', 'tau_si', 'tau_w_inf', 'w_inf_star'
    ]
    # Criar dicionário para armazenar os parâmetros
    params_dict = {}

    # Atribuir valores aos parâmetros conforme a ordem especificada
    for i, nome in enumerate(nomes_variaveis):
        params_dict[nome] = params[i]

    return params_dict

nomes_variaveis = [
    'u_o', 'u_u', 'theta_v', 'theta_w', 'theta_v_minus', 'theta_o',
    'tau_v1_minus', 'tau_v2_minus', 'tau_v_plus', 'tau_w1_minus',
    'tau_w2_minus', 'k_w_minus', 'u_w_minus', 'tau_w_plus', 'tau_fi',
    'tau_o1', 'tau_o2', 'tau_so1', 'tau_so2', 'k_so', 'u_so', 'tau_s1',
    'tau_s2', 'k_s', 'u_s', 'tau_si', 'tau_w_inf', 'w_inf_star'
]

# nomes_variaveis = ['u_u', 'theta_v', 'tau_v_plus', 'tau_w_plus', 'tau_fi', 'tau_so1', 'tau_so2', 'tau_si']

bounds = [[0.0, 10**-25], [1.1899916369215005, 2.010026851326774], [0.2050796377375319, 0.3878536344463603], [0.0102736424810174, 0.0197166419480789], [0.0096524036833057, 0.0200332551860258], [0.0043296310212687, 0.0077592846824974], [44.28409980432256, 75.40597799853211], [801.728413047797, 1479.18327219918], [1.1087718035535137, 1.9376055040946905], [44.0720250125579, 88.83283912921195], [13.6430744283549, 25.044583440217792], [47.862248223651335, 84.11232947155908], [0.0215534595717027, 0.0383198647137663], [185.4873168369487, 349.8818128774007], [0.0765754468303557, 0.1393223946515484], [4.102112457419924, 7.600797702604699], [3.725262758973494, 8.020986117754756], [29.54065042765125, 52.75861093039384], [0.1353197333899674, 0.2513137705478079], [1.4767119243399986, 2.502485943593606], [0.4387623742167367, 0.7849352956805291], [1.9201174344257208, 3.432395344994949], [2.182727847433602, 3.9255312327142255], [1.3658962389588614, 2.6543405361019783], [0.754247325495616, 1.2010141500237823], [2.220319358822649, 3.752612632082664], [0.0479787723164891, 0.0911068043932656], [0.6608731540611883, 1.155508648393741]]

# bounds = [[1.1899916369215005, 2.010026851326774],
#  [0.2050796377375319, 0.3878536344463603],
#  [1.1087718035535137, 1.9376055040946905],
#  [185.4873168369487, 349.8818128774007],
#  [0.0765754468303557, 0.1393223946515484],
#  [4.102112457419924, 7.600797702604699],
#  [3.725262758973494, 8.020986117754756],
#  [2.220319358822649, 3.752612632082664]]

nomes_variaveis = [
    'u_o', 'u_u', 'theta_v', 'theta_w', 'theta_v_minus', 'theta_o',
    'tau_v1_minus', 'tau_v2_minus', 'tau_v_plus', 'tau_w1_minus',
    'tau_w2_minus', 'k_w_minus', 'u_w_minus', 'tau_w_plus', 'tau_fi',
    'tau_o1', 'tau_o2', 'tau_so1', 'tau_so2', 'k_so', 'u_so', 'tau_s1',
    'tau_s2', 'k_s', 'u_s', 'tau_si', 'tau_w_inf', 'w_inf_star'
]

nomes_variaveis_amostragem = ['u_u', 'theta_v', 'tau_v_plus', 'tau_w_plus', 'tau_fi', 'tau_so1', 'tau_so2', 'tau_si']
bounds_amostragem = [[0.5*u_u, 1.0*u_u],
[1.0*theta_v, 2*theta_v],
[0.3*tau_v_plus, 1.0*tau_v_plus],
[0.3*tau_w_plus, 1.0*tau_w_plus],
[0.5*tau_fi, 1.2*tau_fi],
[0.5*tau_so1, 1.0*tau_so1],
[0.4*tau_so2, 1.2*tau_so2],
[1.0*tau_si, 1.5*tau_si]]

# u_u: 1.181602385531229 / 1.58 = 0.7478496110957145
# theta_v: 1.8401723027328656 / 0.3 = 6.133907675776219
# tau_v_plus: 7.6624492098163275 / 1.4506 = 5.282261967335122
# tau_w_plus: 256.17272807166955 / 280.0 = 0.9149026002559627
# tau_fi: 0.722066178200383 / 0.11 = 6.5642379836398455
# tau_so1: 42.569425619151346 / 43 = 0.9899866423058452
# tau_so2: 0.08618041550963146 / 0.2 = 0.4309020775481573
# tau_si: 2.92127940366824 / 2.8723 = 1.0170523286802353

valores_padrao = [
    u_o,  # u_o
    None, # u_u (a ser amostrado)
    None, # theta_v (a ser amostrado)
    theta_w,  # theta_w
    theta_v_minus,  # theta_v_minus
    theta_o,  # theta_o
    tau_v1_minus,   # tau_v1_minus
    tau_v2_minus,    # tau_v2_minus
    None, # tau_v_plus (a ser amostrado)
    tau_w1_minus,    # tau_w1_minus
    tau_w2_minus,    # tau_w2_minus
    k_w_minus,  # k_w_minus
    u_w_minus,  # u_w_minus
    None, # tau_w_plus (a ser amostrado)
    None, # tau_fi (a ser amostrado)
    tau_o1,   # tau_o1
    tau_o2,   # tau_o2
    None, # tau_so1 (a ser amostrado)
    None, # tau_so2 (a ser amostrado)
    k_so,  # k_so
    u_so,  # u_so
    tau_s1,  # tau_s1
    tau_s2,   # tau_s2
    k_s,  # k_s
    u_s,   # u_s
    None, # tau_si (a ser amostrado)
    tau_w_inf,  # tau_w_inf
    w_inf_star   # w_inf_star
]

# Definir os limites do problema
problem_amostragem = {
    'num_vars': len(nomes_variaveis_amostragem),
    'names': nomes_variaveis_amostragem,
    'bounds': bounds_amostragem
}

# Gerar amostras usando Latin Hypercube Sampling

num_samples = (int)(18*556 * 10 / 7)  # Número de amostras desejadas = 10008
num_samples = 22

# Fixando a semente para gerar amostras idênticas
seed = 42
np.random.seed(seed)
param_values_amostragem = latin.sample(problem_amostragem, num_samples)
np.random.seed(seed)  # A mesma semente deve ser usada para obter amostras idênticas

# Combinar amostras com os valores padrão
param_values_completos = []
for params_amostrados in param_values_amostragem:
    params_completos = valores_padrao.copy()
    for nome, valor in zip(nomes_variaveis_amostragem, params_amostrados):
        index = nomes_variaveis.index(nome)
        params_completos[index] = valor
    param_values_completos.append(params_completos)

param_values_completos = np.array(param_values_completos)

# Função para converter a lista de parâmetros em um dicionário
def convert_to_dict(params_list, names):
    return {name: value for name, value in zip(names, params_list)}

import csv

header = problem_amostragem['names'] + ['maxDerivada', 'maxAPD', 'minAPD', 'maxDerivadaCV', 'maxCV', 'minCV']

# Nome do arquivo CSV
from datetime import datetime

filename = 'result_latin_hyper_cube.csv'
now = datetime.now()
date_string = now.strftime("%Y-%m-%d_%H-%M-%S")
filename = f"{date_string}_{filename}"

# Criar e escrever os cabeçalhos no arquivo CSV
with open(filename, 'w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(header)

total = 0
certo = 0

parametros = {'u_o': 0,'u_u': 1.58,'theta_v': 0.3,'theta_w': 0.015,'theta_v_minus': 0.015,'theta_o': 0.006,'tau_v1_minus': 60.0,'tau_v2_minus': 1150,'tau_v_plus': 1.4506,'tau_w1_minus': 70,'tau_w2_minus': 20,'k_w_minus': 65.0,'u_w_minus': 0.03,'tau_w_plus': 280.0,'tau_fi': 0.11,'tau_o1': 6,'tau_o2': 6.0,'tau_so1': 43,'tau_so2': 0.2,'k_so': 2,'u_so': 0.65,'tau_s1': 2.7342,'tau_s2': 3.0,'k_s': 2.0994,'u_s': 0.9087,'tau_si': 2.8723,'tau_w_inf': 0.07,'w_inf_star': 0.94}

for params in param_values_completos:
    results = run_simulation(convert_to_dict(params, nomes_variaveis))
    inputs = params
    valores = convert_to_dict(params, nomes_variaveis)
    for nome, valor in valores.items():
        if (nome in nomes_variaveis_amostragem):
          print(f'{nome}: {valor} / {parametros[nome]} = {valor / parametros[nome]}')
    # print(results)

    if results[0] != -inf:
      with open(filename, 'a', newline='') as file:
        writer = csv.writer(file)
        new_row = list(inputs) + list(results)
        writer.writerow(new_row)
      certo+=1

    total+=1
    print(f"{total} de {len(param_values_completos)} - Certos: {certo}")

print(f"Dados salvos em {filename}")
