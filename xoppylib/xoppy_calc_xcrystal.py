
from xraylib import Crystal_GetCrystalsList
import numpy
import os, platform, six

from xoppylib.xoppy_xraylib_util import bragg_calc
from xoppylib.xoppy_util import locations

def xoppy_calc_xcrystal(
    CRYSTAL_MATERIAL = 32,
    MILLER_INDEX_H = 1,
    MILLER_INDEX_K = 1,
    MILLER_INDEX_L = 1,
    TEMPER = 1,
    MOSAIC = 0,
    GEOMETRY = 0,
    SCAN = 2,
    UNIT = 1,
    SCANFROM = -100,
    SCANTO = 150,
    SCANPOINTS = 200,
    ENERGY = 5000,
    ASYMMETRY_ANGLE = 0.0,
    THICKNESS = 1.0,
    MOSAIC_FWHM = 0,
    RSAG = 0,
    RMER = 0,
    ANISOTROPY = 0,
    POISSON = 0,
    CUT = "",
    FILECOMPLIANCE = "",
    ):

    for file in ["diff_pat.dat", "diff_pat.gle", "diff_pat.par", "diff_pat.xop", "xcrystal.bra"]:
        try:
            os.remove(os.path.join(locations.home_bin_run(), file))
        except:
            pass

    if (GEOMETRY == 1) or (GEOMETRY == 3):
        if ASYMMETRY_ANGLE == 0.0:
            print(
                "xoppy_calc_xcrystal: WARNING: In xcrystal the asymmetry angle is the angle between Bragg planes and crystal surface," +
                "in BOTH Bragg and Laue geometries.")

    descriptor = Crystal_GetCrystalsList()[CRYSTAL_MATERIAL]

    if SCAN == 3:  # energy scan
        emin = SCANFROM - 1
        emax = SCANTO + 1
    else:
        emin = ENERGY - 100.0
        emax = ENERGY + 100.0

    print("Using crystal descriptor: ", descriptor)

    bragg_dictionary = bragg_calc(descriptor=descriptor,
                                  hh=MILLER_INDEX_H, kk=MILLER_INDEX_K, ll=MILLER_INDEX_L,
                                  temper=float(TEMPER),
                                  emin=emin, emax=emax, estep=5.0, fileout="xcrystal.bra")

    with open("xoppy.inp", "wt") as f:
        f.write("xcrystal.bra\n")
        f.write("%d\n" % MOSAIC)
        f.write("%d\n" % GEOMETRY)

        if MOSAIC == 1:
            f.write("%g\n" % MOSAIC_FWHM)
            f.write("%g\n" % THICKNESS)
        else:
            f.write("%g\n" % THICKNESS)
            f.write("%g\n" % ASYMMETRY_ANGLE)

        scan_flag = 1 + SCAN

        f.write("%d\n" % scan_flag)

        f.write("%19.9f\n" % ENERGY)

        if scan_flag <= 3:
            f.write("%d\n" % UNIT)

        f.write("%g\n" % SCANFROM)
        f.write("%g\n" % SCANTO)
        f.write("%d\n" % SCANPOINTS)

        if MOSAIC > 1:  # bent
            f.write("%g\n" % RSAG)
            f.write("%g\n" % RMER)
            f.write("0\n")

            if ((descriptor == "Si") or (descriptor == "Si2") or (descriptor == "Si_NIST") or (
                    descriptor == "Ge") or descriptor == "Diamond"):
                pass
            else:  # not Si,Ge,Diamond
                if ((ANISOTROPY == 1) or (ANISOTROPY == 2)):
                    raise Exception(
                        "Anisotropy data not available for this crystal. Either use isotropic or use external compliance file. Please change and run again'")

            f.write("%d\n" % ANISOTROPY)

            if ANISOTROPY == 0:
                f.write("%g\n" % POISSON)
            elif ANISOTROPY == 1:
                # elas%crystalindex =  irint('CrystalIndex: 0,1,2=Si,3=Ge,4=Diamond: ')
                if descriptor == "Si":
                    f.write("0\n")
                elif descriptor == "Si2":
                    f.write("1\n")
                elif descriptor == "Si_NIST":
                    f.write("2\n")
                elif descriptor == "Ge":
                    f.write("3\n")
                elif descriptor == "Diamond":
                    f.write("4\n")
                else:
                    raise Exception(
                        "Cannot calculate anisotropy data for %s this crystal. Use Si, Ge, Diamond." % descriptor)

                f.write("%g\n" % ASYMMETRY_ANGLE)
                f.write("%d\n" % MILLER_INDEX_H)
                f.write("%d\n" % MILLER_INDEX_K)
                f.write("%d\n" % MILLER_INDEX_L)
            elif ANISOTROPY == 2:
                # elas%crystalindex =  irint('CrystalIndex: 0,1,2=Si,3=Ge,4=Diamond: ')
                if descriptor == "Si":
                    f.write("0\n")
                elif descriptor == "Si2":
                    f.write("1\n")
                elif descriptor == "Si_NIST":
                    f.write("2\n")
                elif descriptor == "Ge":
                    f.write("3\n")
                elif descriptor == "Diamond":
                    f.write("4\n")
                else:
                    raise Exception(
                        "Cannot calculate anisotropy data for %s this crystal. Use Si, Ge, Diamond." % descriptor)

                f.write("%g\n" % ASYMMETRY_ANGLE)
                # TODO: check syntax for CUT: Cut syntax is: valong_X valong_Y valong_Z ; vnorm_X vnorm_Y vnorm_Z ; vperp_x vperp_Y vperp_Z
                f.write("%s\n" % CUT.split(";")[0])
                f.write("%s\n" % CUT.split(";")[1])
                f.write("%s\n" % CUT.split(";")[2])
            elif ANISOTROPY == 3:
                f.write("%s\n" % FILECOMPLIANCE)

    if platform.system() == "Windows":
        command = "\"" + os.path.join(locations.home_bin(), 'diff_pat.exe\" < xoppy.inp')
    else:
        command = "'" + os.path.join(locations.home_bin(), 'diff_pat') + "' < xoppy.inp"
    print("Running command '%s' in directory: %s " % (command, locations.home_bin_run()))
    print("\n--------------------------------------------------------\n")
    os.system(command)
    print("\n--------------------------------------------------------\n")

    # show calculated parameters in standard output
    txt_info = open("diff_pat.par").read()
    for line in txt_info:
        print(line, end="")


if __name__ == "__main__":

    import numpy
    from srxraylib.plot.gol import plot


    tmp = xoppy_calc_xcrystal(
        CRYSTAL_MATERIAL = 32,
        MILLER_INDEX_H = 1,
        MILLER_INDEX_K = 1,
        MILLER_INDEX_L = 1,
        TEMPER = 1,
        MOSAIC = 0,
        GEOMETRY = 0,
        SCAN = 2,
        UNIT = 1,
        SCANFROM = -100,
        SCANTO = 150,
        SCANPOINTS = 200,
        ENERGY = 5000,
        ASYMMETRY_ANGLE = 0.0,
        THICKNESS = 1.0,
        MOSAIC_FWHM = 0,
        RSAG = 0,
        RMER = 0,
        ANISOTROPY = 0,
        POISSON = 0,
        CUT = "",
        FILECOMPLIANCE = "",
    )

    a = numpy.loadtxt("diff_pat.dat", skiprows=5)
    # L  Th-ThB{in} [microrad]  Th-ThB{out} [microrad]  phase_p[rad]  phase_s[rad]  Circ Polariz  p-polarized  s-polarized
    print(a.shape)
    # plot(a[:,0],a[:,-1],title="Energy = %f eV "%energy)

    step = a[1, 0] - a[0, 0]
    h = a[:, -1]
    peak = a[:, -1].max()
    integrated_intensity = a[:, -1].sum() * step
    tt = numpy.where(h >= max(h) * 0.5)
    fwhm = step * (tt[0][-1] - tt[0][0])

    print("Peak value = ", peak)
    print("Integrated intensity [urad] = ", integrated_intensity)
    print("FWHM=", fwhm)

    plot(a[:,0], a[:,-1], xtitle="anggle [urad]", ytitle="reflectivity")


