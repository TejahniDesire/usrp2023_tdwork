
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.figure
from astropy import constants as const
from astropy import units as u

# All units are in cgs

# TO DO: given I get Tb brightness.
# Fix Jv
# Iv = size of emitting reion * jv. units of gravity radius, swartchchild radius for s, 6 rs.
# Iv = 5 Rsch * jv
# m = m87blk
# x = nu/nu_c. # n = number density


def temp_from_thetaE(thetae):
    return (thetae * (const.m_e.cgs * (const.c.cgs ** 2)) / const.k_B.cgs).to(u.K)


# Calculate j emission coefficient with frequency
def j_emission(n, b_field, thetae, theta_b, gamma, nu):
    x = x_calc(nu, b_field, theta_b, gamma)
    return ((n * (const.e.esu ** 2) * nu / (2*np.sqrt(3) * const.c.cgs * thetae ** 2)) \
            * intensity(x)).to(u.erg / (u.cm ** 3 * u.s * u.Hz))


# Calculate j emission coefficient with dimensionless quantity x = nu/nu_c
# n: Number density|
def j_emission_with_x(n, b_field, thetae, theta_b, gamma, x):
    nu = x * nu_c(b_field, theta_b, gamma)
    return ((n * (const.e.esu ** 2) * nu / (2*np.sqrt(3) * const.c.cgs * thetae ** 2)) \
            * intensity(x)).to(u.erg / (u.cm ** 3 * u.s * u.Hz))


def intensity(x):
    return 2.5651 * (1 + 1.92 * (x ** (-1/3)) + (0.9977*x**(-2/3)))*np.exp(-1.8899*x**(1/3))


def nu_b(b_field):
    return (const.e.esu * b_field / (2*np.pi * const.m_e.cgs * const.c.cgs)).to(u.Hz)


def nu_c(b_field, theta_b, gamma):
    return 3/2*nu_b(b_field)*np.sin(theta_b)*gamma**2


# Parameters: [n, ν(nu), c]
def theta_e(temp):
    return ((const.k_B.cgs * temp) / (const.m_e.cgs * (const.c.cgs ** 2))).to(u.dimensionless_unscaled)


def x_calc(nu, b_field, theta_b, gamma):
    return nu / nu_c(b_field, theta_b, gamma)


