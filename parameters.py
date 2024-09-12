def define_parameters(cell_name):
    # create a dictionary with parameters names and values to be set by the name
    parameters = {}
    if (cell_name == 'endo'):
        parameters = {
        'u_o': 0.0,
        'u_u': 1.56,
        'theta_v': 0.3,
        'theta_w': 0.13,
        'theta_v_minus': 0.2,
        'theta_o': 0.006,
        'tau_v1_minus': 75.0,
        'tau_v2_minus': 10.0,
        'tau_v_plus': 1.4506,
        'tau_w1_minus': 6.0,
        'tau_w2_minus': 140.0,
        'k_w_minus': 200.0,
        'u_w_minus': 0.016,
        'tau_w_plus': 280.0,
        'tau_fi': 0.1,
        'tau_o1': 470.0,
        'tau_o2': 6.0,
        'tau_so1': 40.0,
        'tau_so2': 1.2,
        'k_so': 2.0,
        'u_so': 0.65,
        'tau_s1': 2.7342,
        'tau_s2': 2.0,
        'k_s': 2.0994,
        'u_s': 0.9087,
        'tau_si': 2.9013,
        'tau_w_inf': 0.0273,
        'w_inf_star': 0.78
    }
    elif (cell_name == 'm'):
        parameters = {
            'u_o': 0.0,
            'u_u': 1.61,
            'theta_v': 0.3,
            'theta_w': 0.13,
            'theta_v_minus': 0.1,
            'theta_o': 0.005,
            'tau_v1_minus': 80.0,
            'tau_v2_minus': 1.4506,
            'tau_v_plus': 1.4506,
            'tau_w1_minus': 70.0,
            'tau_w2_minus': 8.0,
            'k_w_minus': 200.0,
            'u_w_minus': 0.016,
            'tau_w_plus': 280.0,
            'tau_fi': 0.078,
            'tau_o1': 410.0,
            'tau_o2': 7.0,
            'tau_so1': 91.0,
            'tau_so2': 0.8,
            'k_so': 2.1,
            'u_so': 0.6,
            'tau_s1': 2.7342,
            'tau_s2': 4.0,
            'k_s': 2.0994,
            'u_s': 0.9087,
            'tau_si': 3.3849,
            'tau_w_inf': 0.01,
            'w_inf_star': 0.5
        }
    elif (cell_name == 'epi'):
        parameters = {
        'u_o': 0.0,
        'u_u': 1.55,
        'theta_v': 0.3,
        'theta_w': 0.13,
        'theta_v_minus': 0.006,
        'theta_o': 0.006,
        'tau_v1_minus': 60.0,
        'tau_v2_minus': 1150.0,
        'tau_v_plus': 1.4506,
        'tau_w1_minus': 60.0,
        'tau_w2_minus': 15.0,
        'k_w_minus': 65.0,
        'u_w_minus': 0.03,
        'tau_w_plus': 200.0,
        'tau_fi': 0.11,
        'tau_o1': 400.0,
        'tau_o2': 6.0,
        'tau_so1': 30.0181,
        'tau_so2': 0.9957,
        'k_so': 2.0458,
        'u_so': 0.65,
        'tau_s1': 2.7342,
        'tau_s2': 16.0,
        'k_s': 2.0994,
        'u_s': 0.9087,
        'tau_si': 1.8875,
        'tau_w_inf': 0.07,
        'w_inf_star': 0.94
}
    # DEFINE THE OTHER PARAMETERS THAT MIN MODEL FITS TO (TNNP, Torord, etc) -> for denormalized models, it will take special treatment
    elif (cell_name == 'tnnp-mayra'):
        parameters = {
        'u_o': 0,
        'u_u': 1.58,
        'theta_v': 0.3,
        'theta_w':  0.015,
        'theta_v_minus':  0.015,
        'theta_o':  0.006,
        'tau_v1_minus': 60.0,
        'tau_v2_minus': 1150,
        'tau_v_plus': 1.4506,
        'tau_w1_minus': 70,
        'tau_w2_minus': 20,
        'k_w_minus': 65.0,
        'u_w_minus': 0.03,
        'tau_w_plus': 280.0,
        'tau_fi': 0.11,
        'tau_o1': 6,
        'tau_o2': 6.0,
        'tau_so1': 43,
        'tau_so2': 0.2,
        'k_so': 2,
        'u_so': 0.65,
        'tau_s1': 2.7342,
        'tau_s2': 3.0,
        'k_s': 2.0994,
        'u_s': 0.9087,
        'tau_si': 2.8723,
        'tau_w_inf': 0.07,
        'w_inf_star': 0.94
        }
    # elif (cell_name == '...'):
    #     ...
    else: # return a standard min model
        parameters = {

        }
    return parameters