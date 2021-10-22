from xoppylib.mlayer import MLayer
from srxraylib.plot.gol import plot


if __name__ == "__main__":

    a = MLayer.pre_mlayer(
        interactive=False,
        FILE="pre_mlayer.dat",
        E_MIN=100.0, E_MAX=500.0,
        O_DENSITY=7.19, O_MATERIAL="Cr", #"Water, Liquid", # odd: closer to vacuum
        E_DENSITY=3.00, E_MATERIAL="Sc",  # even: closer to substrate
        S_DENSITY=2.33, S_MATERIAL="Si",  # substrate
        GRADE_DEPTH=0,
        N_PAIRS=50,
        THICKNESS=22.0,
        GAMMA=10.0/22.0,  #  gamma ratio  =  t(even) / (t(odd) + t(even))")
        ROUGHNESS_EVEN=0.0,
        ROUGHNESS_ODD=0.0,
        FILE_DEPTH="myfile_depth.dat",
        GRADE_SURFACE=0,
        FILE_SHADOW="mlayer1.sha",
        FILE_THICKNESS="mythick.dat",
        FILE_GAMMA="mygamma.dat",
        AA0=1.0,AA1=0.0,AA2=0.0,AA3=0.0)

    b = MLayer()
    b.read_preprocessor_file("pre_mlayer.dat")


    #
    # energy scan
    #
    rs, rp, e, t = a.scan(h5file="",
            energyN=100,energy1=300.0,energy2=500.0,
            thetaN=1,theta1=45.0,theta2=45.0)

    print(rs.shape,rp.shape,e.shape,t.shape)

    plot(e,rs[:,0],xtitle="Photon energy [eV]",ytitle="Reflectivity")

    #
    # theta scan
    #
    rs, rp, e, t = a.scan(h5file="",
            energyN=1,energy1=400.0,energy2=401.0,
            thetaN=1000,theta1=40.0,theta2=50.0)

    print(rs.shape,rp.shape,e.shape,t.shape)

    plot(t,rs[0],xtitle="angle [deg]",ytitle="Reflectivity",ylog=False)

    #
    # single point
    #
    a.scan(h5file="",
            energyN=1,energy1=398.0,thetaN=1,theta1=45.0)
