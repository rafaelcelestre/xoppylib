#
# xoppy functions
#
import numpy
from dabax.dabax_xraylib import DabaxXraylib




# needed by bragg_calc
import xraylib
import scipy.constants as codata

from xoppylib.crystals.bragg_preprocessor_file_io import bragg_preprocessor_file_v2_write
from dabax.common_tools import f0_xop
from dabax.common_tools import bragg_metrictensor, lorentz, atomic_symbols


class XoppyDecorator(object):

    # def f0_calc(self,
    #     MAT_FLAG,
    #     DESCRIPTOR,
    #     GRIDSTART,
    #     GRIDEND,
    #     GRIDN,
    #     FILE_NAME="",
    #     charge=0.0,
    #     ):
    #     raise NotImplementedError
    #
    # def f1f2_calc(self, descriptor, energy, theta=3.0e-3, F=0, density=None, rough=0.0, verbose=True):
    #     raise NotImplementedError
    #
    # def f1f2_calc_mix(self, descriptor, energy, theta=3.0e-3, F=0, density=None, rough=0.0, verbose=True):
    #     raise NotImplementedError
    #
    # def f1f2_calc_nist(self, descriptor, energy, theta=3.0e-3, F=0, density=None, rough=0.0, verbose=True):
    #     raise NotImplementedError
    #
    # def cross_calc(self, descriptor, energy, calculate=0, unit=None, density=None, verbose=True):
    #     raise NotImplementedError
    #
    # def cross_calc_mix(self, descriptor, energy, calculate=0, unit=None, parse_or_nist=0, density=None, verbose=True):
    #     raise NotImplementedError
    #
    # def cross_calc_nist(self, descriptor0, energy, calculate=0, unit=None, density=None, verbose=True):
    #     raise NotImplementedError
    #
    # def xpower_calc(self, energies=numpy.linspace(1000.0, 50000.0, 100), source=numpy.ones(100),
    #                 substance=["Be"], flags=[0], dens=["?"], thick=[0.5], angle=[3.0], roughness=0.0,
    #                 output_file=None):
    #     raise NotImplementedError
    #
    # def bragg_calc(self, descriptor="Si", hh=1, kk=1, ll=1, temper=1.0,
    #                emin=5000.0, emax=15000.0, estep=100.0, fileout=None):
    #     raise NotImplementedError
    #
    # def bragg_calc2(descriptor="YB66", hh=1, kk=1, ll=1, temper=1.0, emin=5000.0, emax=15000.0, estep=100.0,
    #                 ANISO_SEL=0,
    #                 fileout=None,
    #                 do_not_prototype=0,  # 0=use site groups (recommended), 1=use all individual sites
    #                 verbose=True,
    #                 material_constants_library=None,
    #                 ):
    #     raise NotImplementedError
    #
    # def crystal_fh(self, input_dictionary, phot_in, theta=None, forceratio=0):
    #     raise NotImplementedError
    #
    # def mare_calc(self, descriptor, H, K, L, HMAX, KMAX, LMAX, FHEDGE, DISPLAY, lambda1, deltalambda, PHI, DELTAPHI,
    #               verbose=0):
    #     raise NotImplementedError

    def bragg_calc2(self, descriptor="YB66", hh=1, kk=1, ll=1, temper=1.0, emin=5000.0, emax=15000.0, estep=100.0,
                    ANISO_SEL=0,
                    fileout=None,
                    # sourceCryst=2, # 0=xraylib, 1=dabax, 2=auto
                    # sourceF0=2,    # 0=xraylib, 1=dabax, 2=auto
                    do_not_prototype=0,  # 0=use site groups (recommended), 1=use all individual sites
                    verbose=True,
                    # material_constants_library=xraylib,
                    ):
        """
        Preprocessor for Structure Factor (FH) calculations. It calculates the basic ingredients of FH.

        :param descriptor: crystal name (as in xraylib)
        :param hh: miller index H
        :param kk: miller index K
        :param ll: miller index L
        :param temper: temperature factor (scalar <=1.0 )
        :param emin:  photon energy minimum
        :param emax: photon energy maximum
        :param estep: photon energy step
        :param ANISO_SEL: source of temperature factor:
                    0: use scalar value defined in temper
                    1: use isotropic value calculated from keyword UNIANISO_COFF in dabax Crystal.dat file
                    2: use anisotropic value calculated from keyword UNIANISO_COFF in dabax Crystal.dat file
        :param fileout: name for the output file (default=None, no output file)
        :return: a dictionary with all ingredients of the structure factor.
        """
        material_constants_library = self
        output_dictionary = {}

        codata_e2_mc2 = codata.e ** 2 / codata.m_e / codata.c ** 2 / (4 * numpy.pi * codata.epsilon_0)  # in m

        # f = open(fileout,'w')

        txt = ""
        txt += "# Bragg version, Data file type\n"
        txt += "2.6 1\n"

        cryst = material_constants_library.Crystal_GetCrystal(descriptor)

        if cryst is None:
            raise Exception("Crystal descriptor %s not found in material constants library" % descriptor)

        volume = cryst['volume']

        # test crystal data - not needed
        itest = 0
        if itest:

            print("  Unit cell dimensions are %f %f %f" % (cryst['a'], cryst['b'], cryst['c']))
            print("  Unit cell angles are %f %f %f" % (cryst['alpha'], cryst['beta'], cryst['gamma']))
            print("  Unit cell volume is %f A^3" % volume)
            print("  Atoms at:")
            print("     Z  fraction    X        Y        Z")
            for i in range(cryst['n_atom']):
                atom = cryst['atom'][i]
                print("    %3i %f %f %f %f" % (atom['Zatom'], atom['fraction'], atom['x'], atom['y'], atom['z']))
            print("  ")

        volume = volume * 1e-8 * 1e-8 * 1e-8  # in cm^3
        rn = (1e0 / volume) * (codata_e2_mc2 * 1e2)

        dspacing = bragg_metrictensor(cryst['a'], cryst['b'], cryst['c'], cryst['alpha'], cryst['beta'], cryst['gamma'],
                                      HKL=[hh, kk, ll])
        dspacing *= 1e-8  # in cm

        txt += "# RN = (e^2/(m c^2))/V) [cm^-2], d spacing [cm]\n"
        txt += "%e %e \n" % (rn, dspacing)

        output_dictionary["rn"] = rn
        output_dictionary["dspacing"] = dspacing

        atom = cryst['atom']
        number_of_atoms = len(atom)
        list_Zatom = [atom[i]['Zatom'] for i in range(len(atom))]

        list_fraction = [atom[i]['fraction'] for i in range(number_of_atoms)]
        try:
            list_charge = [atom[i]['charge'] for i in range(number_of_atoms)]
        except:
            list_charge = [0.0] * number_of_atoms
        list_x = [atom[i]['x'] for i in range(number_of_atoms)]
        list_y = [atom[i]['y'] for i in range(number_of_atoms)]
        list_z = [atom[i]['z'] for i in range(number_of_atoms)]

        # calculate array of temperature factor for all atoms
        #
        # Consider anisotropic temperature factor
        # X.J. Yu, slsyxj@nus.edu.sg
        # A dummy dictionary Aniso with start =0 if no aniso temperature factor input
        # start
        if 'Aniso' in cryst.keys() and cryst['Aniso'][0]['start'] > 0:  # most crystals have no Anisotropic input
            TFac = TemperFactor(1.0 / (2.0 * dspacing * 1e8), cryst['Aniso'], Miller={'h': hh, 'k': kk, 'l': ll}, \
                                cell={'a': cryst['a'], 'b': cryst['b'], 'c': cryst['c']}, n=len(atom))
            B_TFac = 1
        else:
            B_TFac = 0
        #
        #
        #
        list_temper = []
        list_temper_label = []
        if ANISO_SEL == 0:
            for i in range(number_of_atoms):
                list_temper.append(temper)
                list_temper_label.append(-1)
        elif ANISO_SEL == 1:
            if B_TFac:
                for i in range(number_of_atoms):
                    list_temper.append(TFac[0, i])
                    list_temper_label.append(TFac[2, i])
            else:
                raise Exception("No crystal data to calculate isotropic temperature factor for crystal %s" % descriptor)
        elif ANISO_SEL == 2:
            if B_TFac:
                for i in range(number_of_atoms):
                    list_temper.append(TFac[1, i])
                    list_temper_label.append(TFac[2, i])
            else:
                raise Exception(
                    "No crystal data to calculate anisotropic temperature factor for crystal %s" % descriptor)

        list_AtomicName = []
        for i in range(number_of_atoms):
            s = atomic_symbols()[atom[i]['Zatom']]
            # if sourceCryst == 1: # charge is not available in xraylib
            try:  # charge is not available in xraylib
                if atom[i]['charge'] != 0.0:  # if charge is 0, s is symbol only, not B0, etc
                    s = s + f'%+.6g' % atom[i]['charge']
            except:
                pass
            list_AtomicName.append(s)

        # identify the prototypical atoms
        labels_prototypical = []
        for i in range(number_of_atoms):
            labels_prototypical.append(
                "Z=%d C=%g F=%g T=%g" % (list_Zatom[i], list_charge[i], list_fraction[i], list_temper_label[i]))

        if do_not_prototype:
            indices_prototypical = numpy.arange(number_of_atoms)  # different with diff_pat for complex crystal
        else:
            indices_prototypical = numpy.unique(labels_prototypical, return_index=True)[1]

        number_of_prototypical_atoms = len(indices_prototypical)

        # for i in range(number_of_prototypical_atoms):
        #     print("   >>> ", i, indices_prototypical[i], labels_prototypical[indices_prototypical[i]])
        #
        # for i in indices_prototypical:
        #     print("   >>>>> ", i, labels_prototypical[i])
        #
        # print(">>>>  list_labels", len(labels_prototypical), len(indices_prototypical), labels_prototypical)

        #
        # get f0 coefficients
        #

        # f0coeffs = []
        # if sourceF0 == 0:
        #     for i in indices_prototypical:
        #         f0coeffs.append(f0_xop(atom[i]['Zatom']))
        # elif sourceF0 == 1:
        #     for i in indices_prototypical:
        #             f0coeffs.append(material_constants_library.f0_with_fractional_charge(atom[i]['Zatom'], atom[i]['charge']) )
        # elif sourceF0 == 2:
        #     total_charge_flag = numpy.abs(numpy.array(list_charge)).sum() # note the abs(): to be used as flag...
        #
        #     if total_charge_flag != 0: # Use dabax
        #         for i in indices_prototypical:
        #             f0coeffs.append(material_constants_library.f0_with_fractional_charge(atom[i]['Zatom'], atom[i]['charge']))
        #     else: # use xraylib
        #         if 'AtomicName' not in atom[0].keys():
        #             for i in indices_prototypical:  #normal case come in here
        #                 f0coeffs.append(f0_xop(atom[i]['Zatom']))
        #         else:   #for case with like 'Y3+' entries in f0_xop
        #             import re
        #             for i in indices_prototypical:
        #                 x = atom[i]['AtomicName']
        #                 tmp_x = re.search('(^[a-zA-Z]*)',x)
        #                 if tmp_x.group(0) == x:
        #                     f0coeffs.append(f0_xop(atom[i]['Zatom']))  #neutral atom
        #                 else:
        #                     f0coeffs.append(f0_xop(0,AtomicName=x))    #charged atom

        try:
            f0coeffs = []
            for i in indices_prototypical:
                f0coeffs.append(
                    material_constants_library.f0_with_fractional_charge(atom[i]['Zatom'], atom[i]['charge']))
        except:
            f0coeffs = []
            for i in indices_prototypical:
                f0coeffs.append(f0_xop(atom[i]['Zatom']))
        txt += "# Number of different element-sites in unit cell NBATOM:\n%d \n" % number_of_prototypical_atoms
        output_dictionary["nbatom"] = number_of_prototypical_atoms

        txt += "# for each element-site, the number of scattering electrons (Z_i + charge_i)\n"
        atnum_list = []
        for i in indices_prototypical:
            txt += "%f " % (list_Zatom[i] + list_charge[i])
            atnum_list.append(list_Zatom[i] + list_charge[i])
        txt += "\n"
        output_dictionary["atnum"] = atnum_list

        txt += "# for each element-site, the occupation factor\n"
        unique_fraction = [list_fraction[i] for i in indices_prototypical]
        for z in unique_fraction:
            txt += "%g " % (z)
        txt += "\n"
        output_dictionary["fraction"] = unique_fraction

        txt += "# for each element-site, the temperature factor\n"  # temperature parameter
        unique_temper = []
        for i in indices_prototypical:
            txt += "%g " % list_temper[i]
            unique_temper.append(list_temper[i])
        txt += "\n"
        output_dictionary["temper"] = unique_temper

        #
        # Geometrical part of structure factor:  G and G_BAR
        #
        txt += "# for each type of element-site, COOR_NR=G_0\n"
        list_multiplicity = []
        for i in indices_prototypical:
            # zz = list_AtomicName[i]
            # fraction = list_fraction[i]
            # temper = list_temper[i]
            # count = 0
            # for j in range(len(list_Zatom)):
            #     if (list_AtomicName[j] == zz) and (list_fraction[j] == fraction) and (list_temper[j] == temper): count += 1

            if do_not_prototype:
                txt += "%d " % 1
                list_multiplicity.append(1)
            else:
                count = 0
                for j in range(number_of_atoms):
                    if labels_prototypical[j] == labels_prototypical[i]: count += 1
                txt += "%d " % count
                list_multiplicity.append(count)
        txt += "\n"
        output_dictionary["G_0"] = list_multiplicity

        txt += "# for each type of element-site, G and G_BAR (both complex)\n"
        list_g = []
        list_g_bar = []
        for i in indices_prototypical:

            if do_not_prototype:
                # # ga_item = numpy.exp(2j * numpy.pi * (hh * list_x[i] + kk * list_y[i] + ll * list_z[i]))
                # ga += ga_item
                ga = numpy.exp(2j * numpy.pi * (hh * list_x[i] + kk * list_y[i] + ll * list_z[i]))
            else:
                ga = 0.0 + 0j
                for j in range(number_of_atoms):
                    if labels_prototypical[j] == labels_prototypical[i]:
                        # if list_AtomicName[j] == zz and list_fraction[j] == ff and list_temper[j] == tt:
                        ga_item = numpy.exp(2j * numpy.pi * (hh * list_x[j] + kk * list_y[j] + ll * list_z[j]))
                        ga += ga_item

            txt += "(%g,%g) \n" % (ga.real, ga.imag)
            txt += "(%g,%g) \n" % (ga.real, -ga.imag)
            list_g.append(ga)
            list_g_bar.append(ga.conjugate())
        output_dictionary["G"] = list_g
        output_dictionary["G_BAR"] = list_g_bar

        #
        # F0 part
        #
        txt += "# for each type of element-site, the number of f0 coefficients followed by them\n"
        for f0coeffs_item in f0coeffs:
            txt += "%d " % len(f0coeffs_item)
            for cc in f0coeffs_item:
                txt += "%g " % cc
            txt += "\n"
        output_dictionary["f0coeff"] = f0coeffs

        # X.J. Yu, use ceil to round up, otherwise we may get actual max energy less than emax
        npoint = int(numpy.ceil(((emax - emin) / estep + 1)))
        txt += "# The number of energy points NPOINT: \n"
        txt += ("%i \n") % npoint
        output_dictionary["npoint"] = npoint

        txt += "# for each energy point, energy, F1(1),F2(1),...,F1(nbatom),F2(nbatom)\n"
        list_energy = []
        out_f1 = numpy.zeros((len(indices_prototypical), npoint), dtype=float)
        out_f2 = numpy.zeros((len(indices_prototypical), npoint), dtype=float)
        out_fcompton = numpy.zeros((len(indices_prototypical), npoint), dtype=float)  # todo: is complex?

        for i in range(npoint):
            energy = (emin + estep * i)
            txt += ("%20.11e \n") % (energy)
            list_energy.append(energy)

            for j, jj in enumerate(indices_prototypical):
                f1a = xraylib.Fi(list_Zatom[jj], energy * 1e-3)
                f2a = -xraylib.Fii(list_Zatom[jj], energy * 1e-3)  # TODO: check the sign!!
                txt += (" %20.11e %20.11e 1.000 \n") % (f1a, f2a)
                out_f1[j, i] = f1a
                out_f2[j, i] = f2a
                out_fcompton[j, i] = 1.0

        output_dictionary["energy"] = list_energy
        output_dictionary["f1"] = out_f1
        output_dictionary["f2"] = out_f2
        output_dictionary["fcompton"] = out_fcompton

        if fileout != None:
            with open(fileout, "w") as f:
                f.write(txt)
            if verbose: print("File written to disk: %s" % fileout)

        return output_dictionary

class XraylibDecorated(XoppyDecorator):
    def __init__(self):
        pass

class DabaxDecorated(DabaxXraylib, XoppyDecorator):
    def __init__(self,
                 dabax_repository="http://ftp.esrf.eu/pub/scisoft/DabaxFiles/",
                 file_f0="f0_InterTables.dat",
                 file_f1f2="f1f2_Windt.dat",
                 file_CrossSec="CrossSec_EPDL97.dat",
                 ):
        DabaxXraylib.__init__(self,
                           dabax_repository=dabax_repository,
                           file_f0=file_f0,
                           file_f1f2=file_f1f2,
                           file_CrossSec=file_CrossSec)


if __name__ == "__main__":




    xrl = XraylibDecorated()
    dx = DabaxDecorated(dabax_repository="http://ftp.esrf.fr/pub/scisoft/DabaxFiles/")


    # print(xrl.f1f2_calc("Si", 8000, theta=3.0e-3, F=0, density=None, rough=0.0, verbose=True))

    # print(dx.f1f2_calc("Si", 8000, theta=3.0e-3, F=0, density=None, rough=0.0, verbose=True))


    tmp = dx.bragg_calc2(descriptor="YB66", hh=1, kk=1, ll=1, temper=1.0,
                emin=5000.0, emax=15000.0, estep=100.0,
                ANISO_SEL=0,
                fileout="bragg_yb66.dat",
                # sourceCryst=2, # 0=xraylib, 1=dabax, 2=auto
                # sourceF0=2,    # 0
                      # =xraylib, 1=dabax, 2=auto
                do_not_prototype=0, # 0=use site groups (recommended), 1=use all individual sites
                verbose=True)