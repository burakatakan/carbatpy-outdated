import numpy as np
import matplotlib.pyplot as plt
from carbatpy.fl_props_compressor import z_uv, z_ps, z_Tp, z_Tx, z_mm
from scipy.integrate import solve_bvp
from carbatpy.src import compressor_welp as comwel
from carbatpy import compressor_roskosch_orig_rp as comros

def test_getalp():
    size = 100
    step = np.zeros(size)
    a_help = int(0.25 * size)
    b_help = int(0.5 * size)
    c_help = int(0.75 * size)
    step[0: a_help] = 0
    step[a_help: b_help] = 1
    step[b_help: c_help] = 2
    step[c_help: size+1] = 3
    pV = np.zeros(8)
    pV = [34e-3, 34e-3, 3.5, .04, .06071, 48.916, 50., 50. / 2., 2.]  # parameter see above
    Ti = np.linspace(300, 500, size)
    pi = np.linspace(350, 3500, size)
    theta = np.linspace(0, 2 * np.pi, size)
    dxdtheta = -pV[1] / 2 * np.sin(theta) * (1 + 1/pV[2] * np.cos(theta) * (1 - (1/pV[2] * np.sin(theta)) ** 2) ** -0.5)
    dxdt = (2 * np.pi * pV[7]) * dxdtheta
    alpha_welp = np.zeros(size)
    z_it = np.zeros([size, 16])
    z_it[:, 0] = theta
    z_it[:, 1] = -(pV[1] / 2. * (1. - np.cos(z_it[:, 0]) + pV[2] *
                                 (1. - np.sqrt(1. - (1. / pV[2] * np.sin(z_it[:, 0])) ** 2.)))) + \
                 pV[4] * pV[1] + pV[1]  # piston position, x=0 at UT
    z_it[:, 4] = step
    z_it[:, 5] = Ti
    z_it[:, 6] = pi
    v_p = np.zeros(size)
    alpha_roskosch = np.zeros(100)
    for i in range(0, size):
        alpha_welp[i] = comwel.getalp(pV, step[i], dxdt[i], Ti[i], pi[i])
        z_it = comros.getalp(z_it, i, pV)
        v_p[i] = np.abs(z_it[i, 1] - z_it[i - 1, 1]) / ((z_it[i, 0] - z_it[i - 1, 0]) /
                                                     (2. * np.pi * pV[7]))  # dX/dt

    alpha_roskosch = z_it[:, 13]
    plt.plot(theta, alpha_roskosch, '*r')
    plt.plot(theta, alpha_welp, '*y')
    plt.show()
    assert alpha_welp.all == alpha_roskosch.all

def test_main_function():
    "aims to calculate dm, du and dT for 1 iteration"







