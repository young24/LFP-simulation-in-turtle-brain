#!/usr/bin/env python
import numpy as np
from math import sqrt, pow, log

def PS_without_elec_mapping(sigma_T, sigma_S, h, steps, elec_x, elec_y, xmid, ymid, zmid):
    """
    Arguments: sigma_T, sigma_S, h, steps, elec_x, elec_y, xmid, ymid, zmid

    Calculates the potential at the electrode plane elec_z=0 of the MEA from point sources at xmid, ymid, zmid, in a
    slice of thickness h.  To get potential use np.dot(mapping, imem) where imem is the
    transmembrane currents of the point sources. The conductivities sigma are expected to be scalars
    """

    if (zmid > h).any():
        raise RuntimeError("zmid too high")
    if (zmid < 0).any():
        raise RuntimeError("zmid too low")

    mapping = np.zeros((len(elec_x), len(xmid)))
    W = (sigma_T - sigma_S)/(sigma_T + sigma_S)

    for comp in xrange(len(xmid)):
        for elec in xrange(len(elec_x)):
            delta = 1/sqrt(pow(elec_y[elec] - ymid[comp], 2) +
                           pow(elec_x[elec] - xmid[comp], 2) +
                           pow(-zmid[comp], 2))
            n = 1
            while n < steps:
                delta += pow(W, n) * (1/sqrt(pow(elec_y[elec] - ymid[comp], 2) +
                                             pow(elec_x[elec] - xmid[comp], 2) +
                                             pow(+2*n*h - zmid[comp], 2)) +
                                      1/sqrt(pow(elec_y[elec] - ymid[comp], 2) +
                                             pow(elec_x[elec] - xmid[comp], 2) +
                                             pow(-2*n*h - zmid[comp], 2)))
                n += 1
            mapping[elec, comp] = delta
    return 2*mapping/(4*np.pi*sigma_T)


def anisotropic_PS_without_elec_mapping(h, steps, sigma_T, sigma_S, elec_x, elec_y, xmid, ymid, zmid):
    """
    Calculate the moi point source potential at the MEA plane with saline conductivity
    sigma_S is scaled to k * sigma_T. Only for point sources without electrodes.
    Arguments: h, steps, sigma_T, sigma_S, elec_x, elec_y, xmid, ymid, zmid
    """

    if len(sigma_T) != 3:
        raise RuntimeError("Tissue conductivity must have length 3!")
    if len(sigma_S) != 3:
        raise RuntimeError("Saline conductivity must have length 3!")
    if (zmid < 0).any():
        raise RuntimeError("zmid has negative values!")
    if (zmid > h).any():
        raise RuntimeError("zmid outside slice!")

    ratios = sigma_S / sigma_T
    if np.abs(ratios[0] - ratios[2]) <= 1e-15:
        scale_factor = ratios[0]
    elif np.abs(ratios[1] - ratios[2]) <= 1e-15:
        scale_factor = ratios[1]

    sigma_S_scaled = scale_factor * sigma_T
    sigma_T_net = np.sqrt(sigma_T[0] * sigma_T[2])
    sigma_S_net = np.sqrt(sigma_S_scaled[0] * sigma_S_scaled[2])
    anis_W = (sigma_T_net - sigma_S_net)/(sigma_T_net + sigma_S_net)
    mapping = np.zeros((len(elec_x), len(xmid)))

    for comp in xrange(len(xmid)):
        for elec in xrange(len(elec_x)):
            delta = 1/sqrt(sigma_T[0]*sigma_T[2]*pow(elec_y[elec] - ymid[comp], 2) +
                           sigma_T[0]*sigma_T[1]*pow(-zmid[comp], 2) +
                           sigma_T[1]*sigma_T[2]*pow(elec_x[elec] - xmid[comp], 2))
            n = 1
            while n < steps:
                delta += pow(anis_W, n) * (1/sqrt(sigma_T[0]*sigma_T[2]*pow(elec_y[elec] - ymid[comp], 2) +
                                           sigma_T[0]*sigma_T[1]*pow(-2*n*h - zmid[comp], 2) +
                                           sigma_T[1]*sigma_T[2]*pow(elec_x[elec] - xmid[comp], 2)) +
                                          1/sqrt(sigma_T[0]*sigma_T[2]*pow(elec_y[elec] - ymid[comp], 2) +
                                            sigma_T[0]*sigma_T[1]*pow(2*n*h - zmid[comp], 2) +
                                            sigma_T[1]*sigma_T[2]*pow(elec_x[elec] - xmid[comp], 2)))
                n += 1
            mapping[elec, comp] = delta
    return 2/(4*np.pi) * mapping



def PS_with_elec_mapping(sigma_T, sigma_S, h, steps, n_avrg_points, elec_r, elec_x, elec_y, xmid, ymid, zmid):
    """
    Arguments:
           sigma_T, sigma_S, h, steps, n_avrg_points, elec_r, elec_x, elec_y, xmid, ymid, zmid
    """
    if (zmid < 0).any():
        raise RuntimeError("zmid has negative values!")
    if (zmid > h).any():
        raise RuntimeError("zmid outside slice!")

    mapping = np.zeros((len(elec_x), len(xmid)))
    W = (sigma_T - sigma_S)/(sigma_T + sigma_S)

    for comp in xrange(len(xmid)):
        for elec in xrange(len(elec_x)):
            delta = 0
            for pt in xrange(n_avrg_points):
                xpos = (np.random.rand() - 0.5)* 2 * elec_r
                ypos = (np.random.rand() - 0.5)* 2 * elec_r
                # If outside electrode
                while (xpos*xpos + ypos*ypos) > elec_r*elec_r:
                    xpos = (np.random.rand() - 0.5)* 2 * elec_r
                    ypos = (np.random.rand() - 0.5)* 2 * elec_r

                pt_delta = 1/sqrt(pow(elec_y[elec] + ypos - ymid[comp], 2) +
                                  pow(elec_x[elec] + xpos - xmid[comp], 2) +
                                  pow(- zmid[comp], 2))
                n = 1
                while n < steps:
                    pt_delta += pow(W, n) * (
                        1/sqrt(pow(elec_y[elec] + ypos - ymid[comp], 2) +
                               pow(elec_x[elec] + xpos - xmid[comp], 2) +
                               pow(-2*n*h - zmid[comp], 2)) +
                        1/sqrt(pow(elec_y[elec] + ypos - ymid[comp], 2) +
                               pow(elec_x[elec] + xpos - xmid[comp], 2) +
                               pow(2*n*h - zmid[comp], 2)))
                    n += 1
                delta += pt_delta
            mapping[elec, comp] = delta / n_avrg_points
    return 2*mapping/(4*np.pi*sigma_T)



def _LS_omega(a_x, a_y, a_z, comp_length, dx, dy, dz):
    #See Rottman integration formula 46) page 137 for explanation

    factor_a = comp_length*comp_length
    factor_b = - a_x*dx - a_y*dy - a_z*dz
    factor_c = a_x*a_x + a_y*a_y + a_z*a_z
    b_2_ac = factor_b*factor_b - factor_a * factor_c
    if np.abs(b_2_ac) <= 1e-16:
        num = factor_a + factor_b
        den = factor_b
    else:
        num = factor_a + factor_b + \
              comp_length*sqrt(factor_a + 2*factor_b + factor_c)
        den = factor_b + comp_length*sqrt(factor_c)
    return log(num/den)


def LS_without_elec_mapping(sigma_T, sigma_S, h, steps, elec_x, elec_y, xstart, ystart, zstart, xend, yend, zend):
    """
    Arguments:
            sigma_T, sigma_S, h, steps, elec_x, elec_y, xstart, ystart, zstart, xend, yend, zend
    """
    if (zstart < 0).any():
        raise RuntimeError("zstart has negative values")
    if (zstart > h).any():
        raise RuntimeError("zstart above slice")
    if (zend < 0).any():
        raise RuntimeError("zend has negative values")
    if (zend > h).any():
        raise RuntimeError("zend above slice")

    dxs = xend - xstart
    dys = yend - ystart
    dzs = zend - zstart
    length = (dxs**2 + dys**2 + dzs**2)**0.5
    mapping = np.zeros((len(elec_x), len(xstart)))
    W = (sigma_T - sigma_S)/(sigma_T + sigma_S)
    for comp in xrange(len(xstart)):
        for elec in xrange(len(elec_x)):
            a_z0 = -zstart[comp]
            a_x = elec_x[elec] - xstart[comp]
            a_y = elec_y[elec] - ystart[comp]
            n = 1
            delta = _LS_omega(a_x, a_y, a_z0, length[comp], dxs[comp], dys[comp], dzs[comp])
            while n < steps:
                a_z1 = 2*n*h - zstart[comp]
                a_z2 =-2*n*h - zstart[comp]
                delta += pow(W, n) * (_LS_omega(a_x, a_y, a_z1, length[comp], dxs[comp], dys[comp], dzs[comp]) +
                                      _LS_omega(a_x, a_y, a_z2, length[comp], dxs[comp], dys[comp], dzs[comp]))
                n += 1
            mapping[elec, comp] = delta / length[comp]
    return 2 * mapping /(4*np.pi*sigma_T)


def LS_with_elec_mapping(sigma_T, sigma_S, h, steps, n_avrg_points, elec_r, elec_x, elec_y, xstart, ystart, zstart, xend, yend, zend):
    """
    Arguments:
            sigma_T, sigma_S, h, steps, n_avrg_points, elec_r, elec_x, elec_y, xstart, ystart, zstart, xend, yend, zend
    """
    if (zstart < 0).any():
        raise RuntimeError("zstart has negative values")
    if (zstart > h).any():
        raise RuntimeError("zstart above slice")
    if (zend < 0).any():
        raise RuntimeError("zend has negative values")
    if (zend > h).any():
        raise RuntimeError("zend above slice")

    dxs = xend - xstart
    dys = yend - ystart
    dzs = zend - zstart
    length = (dxs**2 + dys**2 + dzs**2)**0.5
    mapping = np.zeros((len(elec_x), len(xstart)))
    W = (sigma_T - sigma_S)/(sigma_T + sigma_S)
    for comp in xrange(len(xstart)):
        for elec in xrange(len(elec_x)):
            delta = 0
            for pt in xrange(n_avrg_points):
                xpos = (np.random.rand() - 0.5)* 2 * elec_r
                ypos = (np.random.rand() - 0.5)* 2 * elec_r
                # If outside electrode
                while (xpos*xpos + ypos*ypos) > elec_r*elec_r:
                    xpos = (np.random.rand() - 0.5)* 2 * elec_r
                    ypos = (np.random.rand() - 0.5)* 2 * elec_r
                a_z0 = - zstart[comp]
                a_x = elec_x[elec] + xpos - xstart[comp]
                a_y = elec_y[elec] + ypos - ystart[comp]
                n = 1
                pt_delta = _LS_omega(a_x, a_y, a_z0, length[comp],
                                     dxs[comp], dys[comp], dzs[comp])
                while n < steps:
                    a_z1 = 2*n*h - zstart[comp]
                    a_z2 =-2*n*h - zstart[comp]
                    pt_delta += pow(W, n) * (_LS_omega(a_x, a_y, a_z1, length[comp], dxs[comp],
                                                       dys[comp], dzs[comp]) +
                                             _LS_omega(a_x, a_y, a_z2, length[comp], dxs[comp],
                                                       dys[comp], dzs[comp]))
                    n += 1
                delta += pt_delta / n_avrg_points
            mapping[elec, comp] = delta / length[comp]
    return 2 * mapping /(4*np.pi*sigma_T)

