import numpy as np
from SALib.sample import latin

def latin_hypercube(params, lower_mult, upper_mult):   
    nomes_variaveis = [
    'u_o', 'u_u', 'theta_v', 'theta_w', 'theta_v_minus', 'theta_o',
    'tau_v1_minus', 'tau_v2_minus', 'tau_v_plus', 'tau_w1_minus',
    'tau_w2_minus', 'k_w_minus', 'u_w_minus', 'tau_w_plus', 'tau_fi',
    'tau_o1', 'tau_o2', 'tau_so1', 'tau_so2', 'k_so', 'u_so', 'tau_s1',
    'tau_s2', 'k_s', 'u_s', 'tau_si', 'tau_w_inf', 'w_inf_star'
]
    nomes_variaveis_amostragem = ['u_u', 'theta_v', 'tau_v_plus', 'tau_w_plus', 'tau_fi', 'tau_so1', 'tau_so2', 'tau_si']
    
    bounds_amostragem = [[lower_mult*params["u_u"], upper_mult*params["u_u"]],
    [lower_mult*params["theta_v"], upper_mult*params["theta_v"]],
    [lower_mult*params["tau_v_plus"], upper_mult*params["tau_v_plus"]],
    [lower_mult*params["tau_w_plus"], upper_mult*params["tau_w_plus"]],
    [lower_mult*params["tau_fi"], upper_mult*params["tau_fi"]],
    [lower_mult*params["tau_so1"], upper_mult*params["tau_so1"]],
    [lower_mult*params["tau_so2"], upper_mult*params["tau_so2"]],
    [lower_mult*params["tau_si"], upper_mult*params["tau_si"]]]
    problem_amostragem = {
        'num_vars': len(nomes_variaveis_amostragem),
        'names': nomes_variaveis_amostragem,
        'bounds': bounds_amostragem
    }
    num_samples = 18
    seed = 42
    np.random.seed(seed)
    param_values_amostragem = latin.sample(problem_amostragem, num_samples)
    np.random.seed(seed)

    
    valores_padrao = [
        params["u_o"],  # u_o
        None,  # u_u (a ser amostrado)
        None,  # theta_v (a ser amostrado)
        params["theta_w"],  # theta_w
        params["theta_v_minus"],  # theta_v_minus
        params["theta_o"],  # theta_o
        params["tau_v1_minus"],  # tau_v1_minus
        params["tau_v2_minus"],  # tau_v2_minus
        None,  # tau_v_plus (a ser amostrado)
        params["tau_w1_minus"],  # tau_w1_minus
        params["tau_w2_minus"],  # tau_w2_minus
        params["k_w_minus"],  # k_w_minus
        params["u_w_minus"],  # u_w_minus
        None,  # tau_w_plus (a ser amostrado)
        None,  # tau_fi (a ser amostrado)
        params["tau_o1"],  # tau_o1
        params["tau_o2"],  # tau_o2
        None,  # tau_so1 (a ser amostrado)
        None,  # tau_so2 (a ser amostrado)
        params["k_so"],  # k_so
        params["u_so"],  # u_so
        params["tau_s1"],  # tau_s1
        params["tau_s2"],  # tau_s2
        params["k_s"],  # k_s
        params["u_s"],  # u_s
        None,  # tau_si (a ser amostrado)
        params["tau_w_inf"],  # tau_w_inf
        params["w_inf_star"]  # w_inf_star
    ]
    
    param_values_completos = []
    for params_amostrados in param_values_amostragem:
        params_completos = valores_padrao.copy()
        for nome, valor in zip(nomes_variaveis_amostragem, params_amostrados):
            index = nomes_variaveis.index(nome)
            params_completos[index] = valor
        param_values_completos.append(params_completos)

    param_values_completos = np.array(param_values_completos)
    
    return param_values_completos, nomes_variaveis

def convert_to_dict(params_list, names):
    return {name: value for name, value in zip(names, params_list)}