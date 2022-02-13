import xml_parser as parser
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np


import os


class surface_gen:
    def __init__(self,xml_file):
        self._par = parser.xml_parser(xml_file)
        return

    def set_x_label(self,x_label):
        self._x_label = x_label
        return

    def set_y_label(self,y_label):
        self._y_label = y_label
        return

    def set_z_label(self,z_label):
        self._z_label = z_label
        return

    def set_title(self,title):
        self._title = title
        return

    def plot(self,clip_last_rows=None):
        xy = self._par.abscissa()
        if self._par.type() == parser.xml_type.SURFACE_ST:
            z = self._par.ordinate()
            hf = plt.figure()
            ha = hf.add_subplot(111, projection='3d')
            # `plot_surface` expects `x` and `y` data to be 2D
            X, Y = np.meshgrid(xy[0], xy[1])
            ha.plot_surface(X, Y, z, rstride=1, cstride=1,
                            cmap='viridis', edgecolor='none')
            ha.set_title(self._title)
            ha.set_xlabel(self._x_label)
            ha.set_ylabel(self._y_label)
            ha.set_zlabel(self._z_label)
            # plt.savefig("")
            plt.show()
        elif self._par.type() == parser.xml_type.SURFACE_SS:
            z = self._par.ordinate(clip_last_rows)
            if not clip_last_rows is None:
                xy[1] = xy[1][:-clip_last_rows]
            hf = plt.figure()
            ha = hf.add_subplot(111, projection='3d')
            # `plot_surface` expects `x` and `y` data to be 2D
            Y, X = np.meshgrid(xy[1], xy[0])
            ha.plot_surface(X, Y, z, rstride=1, cstride=1,
                            cmap='viridis', edgecolor='none')
            ha.set_title(self._title)
            ha.set_xlabel(self._x_label)
            ha.set_ylabel(self._y_label)
            ha.set_zlabel(self._z_label)
            # plt.savefig("")
            plt.show()
            return
        else:
            return

    def plot_2d(self, xy_plane, z_plane):
        hf = plt.figure()
        ha = hf.add_subplot(111, projection='3d')
        # `plot_surface` expects `x` and `y` data to be 2D
        Y, X = np.meshgrid(xy_plane[1], xy_plane[0])
        ha.plot_surface(X, Y, z_plane, rstride=1, cstride=1,
                        cmap='viridis', edgecolor='none')
        ha.set_title(self._title)
        ha.set_xlabel(self._x_label)
        ha.set_ylabel(self._y_label)
        ha.set_zlabel(self._z_label)
        plt.show()
        return


path = os.getcwd()
files_path = "\\".join(path.split("\\")[:-1]) + "\\xmls"
file_crv = files_path+"\\hhw_dsssolver_dr_0_66_cn_srf_numerical.xml"

srf = surface_gen(file_crv)
srf.set_x_label('Spot')
srf.set_y_label('Volatility')
srf.set_z_label('Option Price')
srf.set_title('Heston Vanilla Put PDE (implicit, non-uniform scheme)')
srf.plot(2)


crv = parser.xml_parser(file_crv)
crv_type = crv.type()
crv_y = crv.ordinate()
crv_x = crv.abscissa()
xy_plane = np.asarray([crv_x[0],crv_x[1]])

srf = surface_gen(file_crv)
srf.set_x_label('Spot')
srf.set_y_label('Volatility')
srf.set_z_label('Option Price')
srf.set_title('Heston-Hull-White Call PDE (implicit, uniform scheme)')

srf.plot_2d(xy_plane,crv_y[:,:,1])

crv_y[:,1,:]



