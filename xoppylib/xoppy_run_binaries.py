import sys
import os
import platform
from xoppylib.xoppy_util import locations
import numpy
import scipy.constants as codata



def xoppy_calc_ws(ENERGY=7.0, CUR=100.0, PERIOD=8.5, N=28.0, KX=0.0, KY=8.739999771118164, \
                  EMIN=1000.0, EMAX=100000.0, NEE=2000, D=30.0, XPC=0.0, YPC=0.0, XPS=2.0, YPS=2.0, NXP=10, NYP=10):
    print("Inside xoppy_calc_ws. ")

    try:
        with open("ws.inp", "wt") as f:
            f.write("inputs from xoppy \n")
            f.write("%f     %f\n" % (ENERGY, CUR))
            f.write("%f  %d  %f  %f\n" % (PERIOD, N, KX, KY))
            f.write("%f  %f   %d\n" % (EMIN, EMAX, NEE))
            f.write("%f  %f  %f  %f  %f  %d  %d\n" % (D, XPC, YPC, XPS, YPS, NXP, NYP))
            f.write("%d  \n" % (4))

        if platform.system() == "Windows":
            command = "\"" + os.path.join(locations.home_bin(), 'ws.exe') + "\""
        else:
            command = "'" + os.path.join(locations.home_bin(), 'ws') + "'"
        print("Running command '%s' in directory: %s \n" % (command, locations.home_bin_run()))
        # TODO try to capture the text output of the external code
        os.system(command)

        # write spec file
        txt = open("ws.out").readlines()
        outFile = os.path.join(locations.home_bin_run(), "ws.spec")
        f = open(outFile, "w")

        f.write("#F ws.spec\n")
        f.write("\n")
        f.write("#S 1 ws results\n")
        f.write("#N 4\n")
        f.write("#L  Energy(eV)  Flux(photons/s/0.1%bw)  Spectral power(W/eV)  Cumulated power(W)\n")
        cum = 0.0
        estep = (EMAX - EMIN) / (NEE - 1)
        for i in txt:
            tmp = i.strip(" ")
            if tmp[0].isdigit():
                tmp1 = numpy.fromstring(tmp, dtype=float, sep=' ')
                cum += tmp1[1] * codata.e * 1e3
                f.write("%f  %g  %f  %f \n" % (tmp1[0], tmp1[1], tmp1[1] * codata.e * 1e3, cum * estep))
            else:
                f.write("#UD " + tmp)
        f.close()
        print("File written to disk: ws.spec")

        # print output file
        # for line in txt:
        #     print(line, end="")
        print("Results written to file: %s" % os.path.join(locations.home_bin_run(), 'ws.out'))

        return outFile
    except Exception as e:
        raise e

def xoppy_calc_xtc(
    ENERGY         = 7.0,
    CURRENT        = 100.0,
    ENERGY_SPREAD  = 0.00096,
    SIGX           = 0.274,
    SIGY           = 0.011,
    SIGX1          = 0.0113,
    SIGY1          = 0.0036,
    PERIOD         = 3.23,
    NP             = 70,
    EMIN           = 2950.0,
    EMAX           = 13500.0,
    N              = 40,
    HARMONIC_FROM  = 1,
    HARMONIC_TO    = 15,
    HARMONIC_STEP  = 2,
    HELICAL        = 0,
    METHOD         = 1,
    NEKS           = 100, ):


    for file in ["tc.inp","tc.out"]:
        try:
            os.remove(os.path.join(locations.home_bin_run(),file))
        except:
            pass

    with open("tc.inp", "wt") as f:
        f.write("TS called from xoppy\n")
        f.write("%10.3f %10.2f %10.6f %s\n"%(ENERGY,CURRENT,ENERGY_SPREAD,"Ring-Energy(GeV) Current(mA) Beam-Energy-Spread"))
        f.write("%10.4f %10.4f %10.4f %10.4f %s\n"%(SIGX,SIGY,SIGX1,SIGY1,"Sx(mm) Sy(mm) Sx1(mrad) Sy1(mrad)"))
        f.write("%10.3f %d %s\n"%(PERIOD,NP,"Period(cm) N"))
        f.write("%10.1f %10.1f %d %s\n"%(EMIN,EMAX,N,"Emin Emax Ne"))
        f.write("%d %d %d %s\n"%(HARMONIC_FROM,HARMONIC_TO,HARMONIC_STEP,"Hmin Hmax Hstep"))
        f.write("%d %d %d %d %s\n"%(HELICAL,METHOD,1,NEKS,"Helical Method Print_K Neks"))
        f.write("foreground\n")


    if platform.system() == "Windows":
        command = "\"" + os.path.join(locations.home_bin(),'tc.exe') + "\""
    else:
        command = "'" + os.path.join(locations.home_bin(), 'tc') + "'"


    print("Running command '%s' in directory: %s "%(command, locations.home_bin_run()))
    print("\n--------------------------------------------------------\n")
    os.system(command)
    print("Output file: %s"%("tc.out"))
    print("\n--------------------------------------------------------\n")

    #
    # parse result files to exchange object
    #


    with open("tc.out","r") as f:
        lines = f.readlines()

    # print output file
    # for line in lines:
    #     print(line, end="")


    # remove returns
    lines = [line[:-1] for line in lines]
    harmonics_data = []

    # separate numerical data from text
    floatlist = []
    harmoniclist = []
    txtlist = []
    for line in lines:
        try:
            tmp = line.strip()

            if tmp.startswith("Harmonic"):
                harmonic_number = int(tmp.split("Harmonic")[1].strip())

                if harmonic_number != HARMONIC_FROM:
                    harmonics_data[-1][1] = harmoniclist
                    harmoniclist = []

                harmonics_data.append([harmonic_number, None])

            tmp = float(line.strip()[0])

            floatlist.append(line)
            harmoniclist.append(line)
        except:
            txtlist.append(line)

    harmonics_data[-1][1] = harmoniclist

    data = numpy.loadtxt(floatlist)

    for index in range(0, len(harmonics_data)):
        # print (harmonics_data[index][0], harmonics_data[index][1])
        harmonics_data[index][1] = numpy.loadtxt(harmonics_data[index][1])

    return data, harmonics_data


def xoppy_calc_xtcap(
        ENERGY        = 6.0,
        CURRENT       = 200.0,
        ENERGY_SPREAD = 0.00138,
        SIGX          = 0.0148,
        SIGY          = 0.0037,
        SIGX1         = 0.0029,
        SIGY1         = 0.0015,
        PERIOD        = 2.8,
        NP            = 84,
        EMIN          = 3217.0,
        EMAX          = 11975.0,
        N             = 50,
        DISTANCE      = 30.0,
        XPS           = 1.0,
        YPS           = 1.0,
        XPC           = 0.0,
        YPC           = 0.0,
        HARMONIC_FROM = 1,
        HARMONIC_TO   = 7,
        HARMONIC_STEP = 2,
        HRED          = 0,
        HELICAL       = 0,
        NEKS          = 100,
        METHOD        = 0,
        BSL           = 0,
    ):

    for file in ["tcap.inp","tcap.out","tcap.log"]:
        try:
            os.remove(os.path.join(locations.home_bin_run(),file))
        except:
            pass

    with open("tcap.inp", "wt") as f:
        f.write("TCAP called from xoppy\n")
        f.write("%10.3f %10.2f %10.6f %s\n"%(ENERGY,CURRENT,ENERGY_SPREAD,"Ring-Energy(GeV) Current(mA) Beam-Energy-Spread"))
        f.write("%10.4f %10.4f %10.4f %10.4f %s\n"%(SIGX,SIGY,SIGX1,SIGY1,"Sx(mm) Sy(mm) Sx1(mrad) Sy1(mrad)"))
        f.write("%10.3f %d %s\n"%(PERIOD,NP,"Period(cm) N"))
        f.write("%10.1f %10.1f %d %s\n"%(EMIN,EMAX,N,"Emin Emax Ne"))
        f.write("%10.3f %10.3f %10.3f %10.3f %10.3f %d %d %s\n"%(DISTANCE,XPC,YPC,XPS,YPS,10,10,"d xpc ypc xps yps nxp nyp"))
        f.write("%d %d %d %d %s\n"%(HARMONIC_FROM,HARMONIC_TO,HARMONIC_STEP,HRED,"Hmin Hmax Hstep Hreduction"))
        f.write("%d %d %d %d %d %s\n"%(HELICAL,METHOD,1,NEKS,BSL,"Helical Method Print_K Neks Bsl-Subtr "))
        f.write("foreground\n")


    if platform.system() == "Windows":
        command = "\"" + os.path.join(locations.home_bin(),'tcap.exe') + "\""
    else:
        command = "'" + os.path.join(locations.home_bin(), 'tcap') + "'"


    print("Running command '%s' in directory: %s "%(command, locations.home_bin_run()))
    print("\n--------------------------------------------------------\n")
    # os.system(command)

    #
    #   catch the optut and write the output to a log file as well as print it.
    #
    retvalue = os.popen(command).read()
    print(retvalue)

    with open("tcap.log", "wt") as f:
        f.write(retvalue)

    print("Output file: '%s/tcap.out'"%(os.getcwd()) )
    print("\n--------------------------------------------------------\n")

    #
    # parse result files to exchange object
    #
    with open("tcap.out","r") as f:
        lines = f.readlines()

    # remove returns
    lines = [line[:-1] for line in lines]
    harmonics_data = []

    # separate numerical data from text
    floatlist = []
    harmoniclist = []
    txtlist = []
    for line in lines:
        try:
            tmp = line.strip()

            if tmp.startswith("Harmonic"):
#                   harmonic_number = int(tmp.split("Harmonic")[1].strip())
                harmonic_number = int(tmp.split()[1])

                if harmonic_number != HARMONIC_FROM:
                    harmonics_data[-1][1] = harmoniclist
                    harmoniclist = []

                harmonics_data.append([harmonic_number, None])

            tmp = float(line.strip()[0])

            floatlist.append(line)
            harmoniclist.append(line)
        except:
            txtlist.append(line)

    harmonics_data[-1][1] = harmoniclist

    data = numpy.loadtxt(floatlist)

    for index in range(0, len(harmonics_data)):
        # print (harmonics_data[index][0], harmonics_data[index][1])
        harmonics_data[index][1] = numpy.loadtxt(harmonics_data[index][1])


    return data, harmonics_data