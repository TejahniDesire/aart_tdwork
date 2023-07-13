from aart_func import *
from params import * 
from astropy import constants as const
from astropy import units as u
from lmfit import Parameters,minimize, fit_report

# All r in units of Rg
# Constants

# TODO Units of the images pixel value
# TODO find values to change j coeff peak to occur between 3-20rg

# Constants
G = const.G.cgs
c = const.c.cgs
me = const.m_e.cgs
mp = const.m_p.cgs
e = const.e.esu
kB = const.k_B.cgs

# Units
cm = u.cm
cmcubed = 1 / (cm ** 3)
kelv = 1 * u.K
grams = 1 * u.g
gauss = 1 * u.cm ** (-1 / 2) * u.g ** (1 / 2) * 1 / (1 * u.s)  # Gauss units in cgs
rads = 1 * u.rad
Hz = 1 * u.Hz

'''KWARGS------------------------------'''
'''Physical'''
# kw_nu0 = 230e9 * Hz
kw_nu0 = 86e9 * Hz
#kw_nu0 = 345e9 * Hz
kw_mass = (MMkg * u.kg).to(u.g)  # m87 mass as put in AART
kw_scale_height = .5
kw_theta_b = 60 * (np.pi / 180) * rads
kw_beta = 1

'''Math Coeff'''
kw_rb_0 = 50
kw_n_th0 = 1.23e4 * cmcubed
kw_t_e0 = 8.1e9 * kelv

#kw_n_th0 = 1.23e6 * cmcubed


'''Exponents'''
kw_p_temp = -.84
kw_p_dens = -.7
# kw_p_temp = -.05
# kw_p_dens = -.05

'''B_field '''
b_0 = None
p_b = None


def full_b_func(r,mass=kw_mass,beta=kw_beta,rb_0=kw_rb_0,n_th0=kw_n_th0,p_dens=kw_p_dens):
    nth = nth_func(r,mass,rb_0,n_th0,p_dens)
    if n_th0.unit != u.cm ** -3:
        raise ValueError('n_th0 needs units of number density')
    return np.sqrt(nth * 8 * np.pi * mp * c ** 2 / (6 * beta * r))


def b_simple_func(params, r, y):
    rg = params['rg']
    rb = params['rb']
    b_0 = params['b_0']
    p_b = params['p_b']
    y_fit = b_0 * (r * rg / rb) ** p_b
    return y_fit-y

# TODO Check if least_squares assumes squares
def b_best_fit(mass=kw_mass,beta=kw_beta,rb_0=kw_rb_0,n_th0=kw_n_th0,p_dens=kw_p_dens):
    # Create r and y data
    r = np.array([2, 15, 50])
    y = full_b_func(r,mass,beta,rb_0,n_th0,p_dens).value
    params = Parameters()
    params.add('b_0', value=1)
    params.add('p_b', value=1)
    params.add('rg', value=rg_func(mass).value, vary=False)
    params.add('rb', value=rb_func(mass).value, vary=False)
    fitted_params = minimize(b_simple_func, params, args=(r,y), method='least_squares')
    return fitted_params.params['b_0'].value, fitted_params.params['p_b'].value


def set_b_params(mass=kw_mass,beta=kw_beta,rb_0=kw_rb_0,n_th0=kw_n_th0,p_dens=kw_p_dens):
    global b_0
    global p_b
    b_0, p_b = b_best_fit(mass,beta,rb_0,n_th0,p_dens)
    b_0 = b_0 * gauss


def rg_func(mass=kw_mass):
    return G * mass / c ** 2


def rb_func(mass=kw_mass, rb_0=kw_rb_0):
    return rb_0 * G * mass / (c ** 2)


def te_func(r, mass=kw_mass,rb_0=kw_rb_0,t_e0=kw_t_e0,p_temp=kw_p_temp):
    if t_e0.unit != u.K:
        raise ValueError('n_th0 needs units of Kelvin')   
    rg = rg_func(mass)
    rb = rb_func(mass,rb_0)
    return t_e0 * (r * rg / rb) ** p_temp


def theta_e_func(r,mass=kw_mass,rb_0=kw_rb_0,t_e0=kw_t_e0,p_temp=kw_p_temp):
    temp = te_func(r,mass,rb_0,t_e0,p_temp)
    #return (kB * temp / (me * c ** 2)+1e-10).to(u.dimensionless_unscaled)
    return (kB * temp / (me * c ** 2)).to(u.dimensionless_unscaled)


def nth_func(r, mass=kw_mass,rb_0=kw_rb_0,n_th0=kw_n_th0,p_dens=kw_p_dens):
    if n_th0.unit != u.cm ** -3:
        raise ValueError('n_th0 needs units of number density')
    rg = rg_func(mass)
    rb = rb_func(mass,rb_0)
    return n_th0 * (r * rg / rb) ** p_dens


def b_func(r,mass=kw_mass,rb_0=kw_rb_0):
    rg = rg_func(mass)
    rb = rb_func(mass,rb_0)
    return b_0  * (r * rg / rb) ** p_b
    # return np.sqrt(nth * 8 * np.pi * mp * c ** 2 / (6 * beta * r))


def nu_c_func(r,mass=kw_mass,theta_b=kw_theta_b,rb_0=kw_rb_0,t_e0=kw_t_e0,p_temp=kw_p_temp):
    b_field = b_func(r,mass,rb_0)
    nu_b = (e * b_field / (2 * np.pi * me * c)).to(u.Hz)
    theta_e = theta_e_func(r,mass,rb_0,t_e0,p_temp)
    return 3 / 2 * nu_b * np.sin(theta_b) * theta_e ** 2


def synchrotron_func(x):
    return 2.5651 * (1 + 1.92 * (x ** (-1 / 3)) + (0.9977 * x ** (-2 / 3))) * np.exp(-1.8899 * x ** (1 / 3))


def profile(r, redshift, nu0=kw_nu0,mass=kw_mass, scale_height=kw_scale_height, theta_b=kw_theta_b, 
            beta=kw_beta,rb_0=kw_rb_0,n_th0=kw_n_th0,t_e0=kw_t_e0,p_dens=kw_p_dens,p_temp=kw_p_temp):
    n = nth_func(r,mass,rb_0,n_th0,p_dens)
    theta_e = theta_e_func(r,mass,rb_0,t_e0,p_temp)
    nu_c = nu_c_func(r,mass,theta_b,rb_0,t_e0,p_temp)
    nu = nu0/redshift
    x = nu / (nu_c)
    # Returns units of [u.erg / (u.cm ** 3 * u.s * u.Hz)]
    jcoeff = n * e ** 2 * nu * synchrotron_func(x) / (2 * np.sqrt(3) * c * theta_e ** 2)
    specific_intensity = r * scale_height * rg_func(mass) * jcoeff
    return ((c ** 2 / (2 * nu ** 2 * kB)) * specific_intensity).to(u.K)

# Assume I is ndarray without astropy units, but in kelvin
def total_jy(I, nu, mass):
    I = I * u.K
    nu = nu * u.Hz
    mass = mass * u.g
    one_M = rg_func(mass).to(u.m)
    M2rads = np.arctan(one_M.value / dBH)
    rads2pxls = (M2rads * 2 * limits) / (I.shape)[0]
    return (I * nu ** 2 * (2 * kB / c ** 2)).to(u.Jy).sum() * rads2pxls ** 2

def ring_radius(I0):
    count = np.zeros([(I0.shape)[0] + 1000, 4])


