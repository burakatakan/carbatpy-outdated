'''
heat transfer calculation for heat exchanger
'''
# import libraries
import fluid_properties_rp as fprop
import numpy as np
import ht as ht

def alpha_condenser(p, h , fluid, comp, m_dot, d, x_position=None, rel_roughness=None, darcy_friction=None):
    """
    calculates alpha for condensation and superheated/subcooled state before or after
    defines state first, for working fluids mixtures
    uses correlations from ht heat transfer database
    :param p: pressure [Pa]
    :param h: enthalpy [J/(kg*K)]
    :param fluid: fluid name [string]
    :param comp: fluid composition in format [x_a, x_b]
    :param m_dot: mass flow [kg/s]
    :param d: inner diameter
    :param x_position: axial position in tube if available [m]
    :param rel_roughness: relative roughness if available [-]
    :param darcy_friction: darcy friction factor if available [-]
    :return: alpha [W/(m^2*K)]
    """
    # superheated, wet-steam or subcooled
    act_state = fprop.hp_v(h, p, fluid, comp, option=0)
    T_AF = act_state[0]
    sat_state = fprop.p_prop_sat(p, fluid, comp, option=0)
    Re =  4 * m_dot / (act_state[7] * np.pi * d)

    if T_AF < sat_state[1,0]:
        # subcooled
        alpha = ht.conv_internal.Nu_conv_internal(Re, act_state[9], eD=rel_roughness, Di=d, x=x_position, fd=darcy_friction)

    elif T_AF > sat_state[0,0]:
        # superheated
        alpha = ht.conv_internal.Nu_conv_internal(Re, act_state[9], Di=d, x=x_position)

    else:
        # wet-steam
        x = (h - sat_state[1,2]) / (sat_state[0,2] - sat_state[1,2])
        alpha = ht.condensation.Cavallini_Smith_Zecchin(m_dot, x, d, 1/sat_state[1,3], 1/sat_state[0,3],
                                                        sat_state[1,7], sat_state[0,7], sat_state[1,8], sat_state[2,6])
    return alpha
alpha_condenser_vec = np.vectorize(alpha_condenser)



def alpha_tube_annulus():


