import numpy as np
from xoppylib.fit_gaussian2d import twoD_Gaussian, fit_gaussian2d, info_params
from srxraylib.plot.gol import plot_image

if __name__ == "__main__":
    nx = 200
    ny = 300
    x0 = np.linspace(-20,20,nx)
    y0 = np.linspace(-10,10,ny)
    x, y = np.meshgrid(x0, y0)
    data = twoD_Gaussian((x, y), 100, 0, 0, 10, 5, 0, 0)

    plot_image(data.reshape((nx,ny)),x0,y0,title="DATA")


    popt = fit_gaussian2d(data,x0,y0,p0=None)

    data_fitted = twoD_Gaussian((x, y), *popt)
        
    print(info_params(popt))

    plot_image(data_fitted.reshape((nx,ny)),x0,y0,title="FIT")



    





