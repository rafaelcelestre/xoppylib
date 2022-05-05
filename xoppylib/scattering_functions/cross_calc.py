import numpy
import scipy.constants as codata

from xoppylib.xoppy_xraylib_util import interface_reflectivity

import xraylib
from dabax.dabax_xraylib import DabaxXraylib


# def f1f2_calc(descriptor, energy, theta=3.0e-3, F=0, density=None, rough=0.0, verbose=True,
#               material_constants_library = None,):

def cross_calc(descriptor, energy, calculate=0, unit=None, density=None, verbose=True):
    """
    calculate the atomic cross sections and attenuation coefficients.
    :param descriptor: string with the element symbol or integer with Z
    :param energy: array with energies (eV)
    :param calculate:
            0: total cross section
            1: photoelectric cross section
            2: rayleigh cross serction
            3: compton cross section
            4: total minus raileigh cross section
    :param unit:An flag indicating the unit of the output array
            None (default) return all units in multiple columns
            0: barn/atom (Cross Section calculation)
            1: cm^2 (Cross Section calculation)
            2: cm^2/g (Mass Attenuation Coefficient)
            3: cm^-1 (Linear Attenuation Coefficient)
    :param density: the material density in g/cm^3
    :return:  if unit=None an array (5, npoints) with energy and unit=0 to 4, else returns one-column array
    """


    energy = numpy.array(energy,dtype=float).reshape(-1)
    out = numpy.zeros_like(energy)

    if isinstance(descriptor,str):
        Z = xraylib.SymbolToAtomicNumber(descriptor)
        symbol = descriptor
    else:
        Z = descriptor
        symbol = xraylib.AtomicNumberToSymbol(descriptor)

    tmp = numpy.zeros_like(energy)

    if calculate == 0:
        for i,ienergy in enumerate(energy):
            tmp[i] = xraylib.CSb_Total(Z,1e-3*ienergy)
    elif calculate == 1:
        for i,ienergy in enumerate(energy):
            tmp[i] = xraylib.CSb_Photo(Z,1e-3*ienergy)
    elif calculate == 2:
        for i,ienergy in enumerate(energy):
            tmp[i] = xraylib.CSb_Rayl(Z,1e-3*ienergy)
    elif calculate == 3:
        for i,ienergy in enumerate(energy):
            tmp[i] = xraylib.CSb_Compt(Z,1e-3*ienergy)
    elif calculate == 4:
        for i,ienergy in enumerate(energy):
            tmp[i] = xraylib.CSb_Total(Z,1e-3*ienergy) - xraylib.CSb_Rayl(Z,1e-3*ienergy)

    if density is None:
        density = xraylib.ElementDensity(Z)

    out = numpy.zeros((5,energy.size))
    out[0,:] = energy
    out[1,:] = tmp         # barn/atom (Cross Section calculation)
    out[2,:] = tmp * 1e-24 #  cm^2 (Cross Section calculation)
    out[3,:] = tmp * 1e-24 * codata.Avogadro / xraylib.AtomicWeight(Z)           # cm^2/g (Mass Attenuation Coefficient)
    out[4,:] = tmp * 1e-24 * codata.Avogadro / xraylib.AtomicWeight(Z) * density # cm^-1 (Linear Attenuation Coefficient)

    if unit is None:
        return out
    else:
        return out[1+unit,:].copy()


def cross_calc_mix(descriptor, energy, calculate=0, unit=None, parse_or_nist=0, density=None, verbose=True):
    """
    Same as cross_calc, but for a compund formula
    :param descriptor: a compound descriptor (as in xraylib)
    :param energy: photon energy array in eV
    :param calculate:
            0: total cross section
            1: photoelectric cross section
            2: rayleigh cross serction
            3: compton cross section
            4: total minus raileigh cross section
    :param unit:An flag indicating the unit of the output array
            None (default) return all units in multiple columns
            0: barn/atom (Cross Section calculation)
            1: cm^2 (Cross Section calculation)
            2: cm^2/g (Mass Attenuation Coefficient)
            3: cm^-1 (Linear Attenuation Coefficient)
    :param parse_or_nist: useless. Kept for back-compatibility.
    :param density: the material density in g/cm^3
    :return:
    """

    energy = numpy.array(energy,dtype=float).reshape(-1)
    out = numpy.zeros_like(energy)

    if (density is None):
        raise Exception("Please define density")

    if verbose: print("cross_calc_mix: Using density %g g/cm3"%density)
    tmp = numpy.zeros_like(energy)
    tmp2 = numpy.zeros_like(energy)

    if calculate == 0:
        for i,ienergy in enumerate(energy):
            tmp[i] = xraylib.CSb_Total_CP(descriptor,1e-3*ienergy)
            tmp2[i] = xraylib.CS_Total_CP(descriptor,1e-3*ienergy)
    elif calculate == 1:
        for i,ienergy in enumerate(energy):
            tmp[i] = xraylib.CSb_Photo_CP(descriptor,1e-3*ienergy)
            tmp2[i] = xraylib.CS_Photo_CP(descriptor,1e-3*ienergy)
    elif calculate == 2:
        for i,ienergy in enumerate(energy):
            tmp[i] = xraylib.CSb_Rayl_CP(descriptor,1e-3*ienergy)
            tmp2[i] = xraylib.CS_Rayl_CP(descriptor,1e-3*ienergy)
    elif calculate == 3:
        for i,ienergy in enumerate(energy):
            tmp[i] = xraylib.CSb_Compt_CP(descriptor,1e-3*ienergy)
            tmp2[i] = xraylib.CS_Compt_CP(descriptor,1e-3*ienergy)
    elif calculate == 4:
        for i,ienergy in enumerate(energy):
            tmp[i] = xraylib.CSb_Total_CP(descriptor,1e-3*ienergy) - xraylib.CSb_Rayl_CP(descriptor,1e-3*ienergy)
            tmp2[i] = xraylib.CS_Total_CP(descriptor,1e-3*ienergy) - xraylib.CS_Rayl_CP(descriptor,1e-3*ienergy)


    out = numpy.zeros((5,energy.size))
    out[0,:] = energy
    out[1,:] = tmp # barn/atom (Cross Section calculation)
    out[2,:] = tmp * 1e-24 #  cm^2 (Cross Section calculation)
    out[3,:] = tmp2 # cm^2/g (Mass Attenuation Coefficient)
    out[4,:] = tmp2 * density # cm^-1 (Linear Attenuation Coefficient)

    if unit is None:
        return out
    else:
        return out[1+unit,:].copy()

def cross_calc_nist(descriptor0, energy, calculate=0, unit=None, density=None, verbose=True):
    """
    Same as cross_calc, but for a compund from the NIST compound list
    :param descriptor: a compound descriptor (as in xraylib)
    :param energy: photon energy array in eV
    :param calculate:
            0: total cross section
            1: photoelectric cross section
            2: rayleigh cross serction
            3: compton cross section
            4: total minus rayleigh cross section
    :param unit:An flag indicating the unit of the output array
            None (default) return all units in multiple columns
            0: barn/atom (Cross Section calculation)
            1: cm^2 (Cross Section calculation)
            2: cm^2/g (Mass Attenuation Coefficient)
            3: cm^-1 (Linear Attenuation Coefficient)
    :param density: the material density in g/cm^3. If set to None, use the xraylib value.
    :return:
    """

    energy = numpy.array(energy,dtype=float).reshape(-1)
    out = numpy.zeros_like(energy)


    if isinstance(descriptor0, int):
        nist_compound = xraylib.GetCompoundDataNISTByIndex(descriptor0)
        descriptor = nist_compound["name"]
    else:
        nist_compound = xraylib.GetCompoundDataNISTByName(descriptor0)
        descriptor = descriptor0

    if verbose:
        print("nist compound:")
        for key in nist_compound.keys():
            print("   %s: "%key, nist_compound[key])

    if density is None:
        density = nist_compound["density"]

    if verbose: print("cross_calc_mix: Using density %g g/cm3"%density)
    tmp = numpy.zeros_like(energy)
    tmp2 = numpy.zeros_like(energy)

    if calculate == 0:
        for i,ienergy in enumerate(energy):
            tmp[i] = xraylib.CSb_Total_CP(descriptor,1e-3*ienergy)
            tmp2[i] = xraylib.CS_Total_CP(descriptor,1e-3*ienergy)
    elif calculate == 1:
        for i,ienergy in enumerate(energy):
            tmp[i] = xraylib.CSb_Photo_CP(descriptor,1e-3*ienergy)
            tmp2[i] = xraylib.CS_Photo_CP(descriptor,1e-3*ienergy)
    elif calculate == 2:
        for i,ienergy in enumerate(energy):
            tmp[i] = xraylib.CSb_Rayl_CP(descriptor,1e-3*ienergy)
            tmp2[i] = xraylib.CS_Rayl_CP(descriptor,1e-3*ienergy)
    elif calculate == 3:
        for i,ienergy in enumerate(energy):
            tmp[i] = xraylib.CSb_Compt_CP(descriptor,1e-3*ienergy)
            tmp2[i] = xraylib.CS_Compt_CP(descriptor,1e-3*ienergy)
    elif calculate == 4:
        for i,ienergy in enumerate(energy):
            tmp[i] = xraylib.CSb_Total_CP(descriptor,1e-3*ienergy) - xraylib.CSb_Rayl_CP(descriptor,1e-3*ienergy)
            tmp2[i] = xraylib.CS_Total_CP(descriptor,1e-3*ienergy) - xraylib.CS_Rayl_CP(descriptor,1e-3*ienergy)

    out = numpy.zeros((5,energy.size))
    out[0,:] = energy
    out[1,:] = tmp # barn/atom (Cross Section calculation)
    out[2,:] = tmp * 1e-24 #  cm^2 (Cross Section calculation)
    out[3,:] = tmp2 # cm^2/g (Mass Attenuation Coefficient)
    out[4,:] = tmp2 * density # cm^-1 (Linear Attenuation Coefficient)

    if unit is None:
        return out
    else:
        return out[1+unit,:].copy()



