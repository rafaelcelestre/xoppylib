import numpy
from xoppylib.power.power1d_calc import power1d_calc

def xoppy_calc_power(
                energies                   = numpy.linspace(1,100,100),
                source                     = numpy.ones(100),
                substance                  = ['Si']*5,
                thick                      = [1.0]*5,
                angle                      = [0.0]*5,
                dens                       = ['?']*5,
                roughness                  = [0.0]*5,
                flags                      = [1]*5,
                NELEMENTS                  = 1,
                FILE_DUMP                  = 0,
                material_constants_library = None,
):

    substance = substance[0:NELEMENTS+1]
    thick     = thick[0:NELEMENTS+1]
    angle     = angle[0:NELEMENTS+1]
    dens      = dens[0:NELEMENTS+1]
    roughness = roughness[0:NELEMENTS+1]
    flags     = flags[0:NELEMENTS+1]


    if FILE_DUMP:
        output_file = "power.spec"
    else:
        output_file = None

    out_dictionary = power1d_calc(energies=energies, source=source, substance=substance,
                flags=flags, dens=dens, thick=thick, angle=angle, roughness=roughness,
                output_file=output_file, material_constants_library=material_constants_library)

    return out_dictionary

if __name__ == "__main__":
    import xraylib
    print(xoppy_calc_power(material_constants_library=xraylib))