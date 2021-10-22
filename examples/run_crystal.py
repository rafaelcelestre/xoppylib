from xoppylib.xoppy_calc_xcrystal import xoppy_calc_xcrystal

if __name__ == "__main__":
    import numpy
    from srxraylib.plot.gol import plot

    tmp = xoppy_calc_xcrystal(
        CRYSTAL_MATERIAL=32, # index of xraylib.Crystal_GetCrystalsList()
        MILLER_INDEX_H=1,
        MILLER_INDEX_K=1,
        MILLER_INDEX_L=1,
        TEMPER=1,
        MOSAIC=0,
        GEOMETRY=0,
        SCAN=2,
        UNIT=1,
        SCANFROM=-100,
        SCANTO=150,
        SCANPOINTS=200,
        ENERGY=5000,
        ASYMMETRY_ANGLE=0.0,
        THICKNESS=1.0,
        MOSAIC_FWHM=0,
        RSAG=0,
        RMER=0,
        ANISOTROPY=0,
        POISSON=0,
        CUT="",
        FILECOMPLIANCE="",
    )


    # L  Th-ThB{in} [microrad]  Th-ThB{out} [microrad]  phase_p[rad]  phase_s[rad]  Circ Polariz  p-polarized  s-polarized
    a = numpy.loadtxt("diff_pat.dat", skiprows=5)

    step = a[1, 0] - a[0, 0]
    h = a[:, -1]
    peak = a[:, -1].max()
    integrated_intensity = a[:, -1].sum() * step
    tt = numpy.where(h >= max(h) * 0.5)
    fwhm = step * (tt[0][-1] - tt[0][0])

    print("Peak value = ", peak)
    print("Integrated intensity [urad] = ", integrated_intensity)
    print("FWHM=", fwhm)

    plot(a[:, 0], a[:, -1], xtitle="angle [urad]", ytitle="reflectivity s-pol")
