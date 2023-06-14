
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.figure
from astropy import constants as const
from astropy import units as u

# All units are in cgsa

# TO DO: given I get Tb brightness.
# Fix Jv
# Iv = size of emitting reion * jv. units of gravity radius, swartchchild radius for s, 6 rs.
# Iv = 5 Rsch * jv
# m = m87blk
# x = nu/nu_c. # n = number density
#
def j_emission(n,temp, mass, b_field, theta_b, nu, gamma):
    nu_c_value = nu_c(b_field, mass, theta_b, gamma)
    theta_e_value = theta_e(temp, mass)
    x = x_calc(nu, b_field, mass, theta_b, gamma)
    return ((n * (const.e.esu ** 2) * x * nu_c_value / \
             (2*np.sqrt(3) * const.c.cgs * theta_e_value)) * intensity(nu, b_field, mass, theta_b, gamma)).to(u.erg / (u.cm ** 3 * u.s * u.Hz))
#u.g / (u.s ** 2 * u.cm)


def intensity(nu, b_field, mass, theta_b, gamma):
    x = x_calc(nu, b_field, mass, theta_b, gamma)
    return 2.565 * (1 + 1.92 * (x ** (-1/3)) + (0.9977*x**(-2/3)) )*np.exp(-1.8899*x**(1/3))


def nu_b(b_field, mass):
    return (const.e.esu * b_field / (2*np.pi*mass * const.c.cgs)).to(u.Hz)


def nu_c(b_field, mass, theta_b, gamma):
    return 3/2*nu_b(b_field, mass)*np.sin(theta_b)*gamma**2


# Parameters: [n, ν(nu), c]
def theta_e(temp, mass):
    return ((const.k_B.cgs * temp) / (mass * (const.c.cgs ** 2))).to(u.dimensionless_unscaled)


def x_calc(nu, b_field, mass, theta_b, gamma):
    return nu / nu_c(b_field, mass, theta_b, gamma)


