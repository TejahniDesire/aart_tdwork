import matplotlib.pyplot as plt
from astropy import units as u
from aart_func import *
from params import *
from astropy import constants as const
from astropy import units as u
from lmfit import Parameters, minimize, fit_report


def radii_of_theta(I0, size):
    x = np.arange(I0.shape[0])
    y = x
    interp = RegularGridInterpolator((x, y), I0.T)
    theta = np.matrix(np.linspace(0, 2 * np.pi, size))  # 1 x size

    rmax = I0.shape[0] * .4
    r = np.matrix(np.arange(rmax)).T  # rmax x 1

    onest = np.matrix(np.ones(r.shape[0])).T  # (rmax) x 1
    onesr = np.matrix(np.ones(size))  # 1 x size

    thetarray = onest @ theta  # (rmax) x size
    rarray = r @ onesr  # (rmax) x size

    # Convert to pixel coords from aart plot coords
    xaart = np.multiply(rarray, np.cos(thetarray))
    yaart = np.multiply(rarray, np.sin(thetarray))
    xprime = xaart + I0.shape[0] / 2
    yprime = yaart + I0.shape[0] / 2

    coords = np.array([xprime, yprime]).T
    peak = np.argmax(interp(coords), 1)
    return (peak * ((limits * 2) / I0.shape[0])), np.ravel(theta)  # units of Rg


size = 200
rsize = 10000


def radii_of_thetaV2(I0, dx=None):
    x = np.arange(I0.shape[0])  # number of pixels
    y = x
    rmax = I0.shape[0] * .4
    interp = RegularGridInterpolator((x, y), I0.T)
    theta = np.matrix(np.linspace(0, 2 * np.pi, size))  # 1 x size

    # rmax = I0.shape[0] * .4
    rmax = I0.shape[0] * .4

    r = np.matrix(np.linspace(0, rmax, rsize)).T  # rsize x 1

    onest = np.matrix(np.ones(r.shape[0])).T  # (rsize) x 1
    onesr = np.matrix(np.ones(size))  # 1 x size

    thetarray = onest @ theta  # (rsize) x size
    rarray = r @ onesr  # (rsize) x size

    xaart = np.multiply(rarray, np.cos(thetarray))
    yaart = np.multiply(rarray, np.sin(thetarray))

    # Convert to pixel coords from aart plot coords
    xprime = xaart + I0.shape[0] / 2
    yprime = yaart + I0.shape[0] / 2

    coords = np.array([xprime, yprime]).T
    peak = np.argmax(interp(coords), 1)

    if not np.all((peak == 0) == False):
        pass

        # TODO What's the best way to remove the point
        #  without changing the array size. Averaging around could include 0,
        #  going to left or right can hit another zero
        # zero_point = np.argmax(peak == 0)
        # peak = np.delete(peak, zero_point)  # Remove any peak listed occuring at index 0
        # theta = np.delete(theta,zero_point)

    if dx is None:
        dx = (limits * 2) / I0.shape[0]

    peaks = np.ravel(r[peak])  # value of r at that argument

    return (peaks * dx), np.ravel(theta)  # units of Rg


def curve_params(varphi, rho):
    """calculate Appendix B parameters for a curve rho(varphi)
       assume varphis are equally spaced!!!"""

    # spacing in varphi
    dvarphi = varphi[-1] - varphi[-2]

    # area
    area = np.trapz(0.5 * rho ** 2, dx=dvarphi)

    # centroid
    mux = np.trapz((rho ** 3 * np.cos(varphi)) / (3 * area), dx=dvarphi)
    muy = np.trapz((rho ** 3 * np.sin(varphi)) / (3 * area), dx=dvarphi)

    # second moment
    Sxx = np.trapz((rho ** 4 * np.cos(varphi) ** 2) / (4 * area), dx=dvarphi) - mux ** 2
    Syy = np.trapz((rho ** 4 * np.sin(varphi) ** 2) / (4 * area), dx=dvarphi) - muy ** 2
    Sxy = np.trapz((rho ** 4 * np.sin(varphi) * np.cos(varphi)) / (4 * area), dx=dvarphi) - mux * muy

    # diagonalize 2nd moment matrix
    D = np.sqrt((Sxx - Syy) ** 2 + 4 * Sxy * Sxy)
    a = np.sqrt(2 * (Sxx + Syy + D))
    b = np.sqrt(2 * (Sxx + Syy - D))

    # radius, eccentricity, position angle
    r = np.sqrt(0.5 * (a ** 2 + b ** 2))
    e = np.sqrt(1 - b ** 2 / a ** 2)
    chi = 0.5 * np.arcsin(2 * Sxy / D)

    # return (area, mux, muy, r, e, chi)
    return r
