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

    def plot_2d(self, xy_plane, z_plane,clip_last_rows=None):
        hf = plt.figure()
        ha = hf.add_subplot(111, projection='3d')
        if not clip_last_rows is None:
            xy_plane[1] = xy_plane[1][:-clip_last_rows]
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
file_crv = files_path+"\\heston_upoutcall_barrier_thomas_lu_cn_dr_srf_stepping_numerical.xml"

srf = surface_gen(file_crv)
srf.set_x_label('Spot')
srf.set_y_label('Volatility')
srf.set_z_label('Call value')
srf.set_title('SABR PDE (Double Sweeep, nonuniform grid)')
srf.plot(3)


crv = parser.xml_parser(file_crv)
crv_type = crv.type()
crv_y = crv.ordinate()
crv_x = crv.abscissa()
xy_plane = np.asarray([crv_x[0],crv_x[1]])

srf = surface_gen(file_crv)
srf.set_x_label('Spot')
srf.set_y_label('Volatility')
srf.set_z_label('Call option Price')
srf.set_title('SABR PDE (Double Sweeep, nonuniform grid)')

srf.plot_2d(xy_plane,crv_y[180,:,:-3],3)

np.shape(crv_y)

# ================= MULTIPLE SURFACES in one plot ==================

fig = plt.figure()
ax = plt.axes(projection='3d')

plt.rcParams['legend.fontsize'] = 10


# First constraint
spot = crv_x[0]
vol = crv_x[1][:-3]
VOL,SPOT = np.meshgrid(vol,spot)
call_0 = crv_y[0,:,:-3]

ax = fig.gca(projection='3d')
c1 = ax.plot_surface(SPOT,VOL, call_0, label = "Call Price Surface")
c1._facecolors2d=c1._facecolors3d
c1._edgecolors2d=c1._edgecolors3d

# Second
call_100 = crv_y[100,:,:-3]
c2 = ax.plot_surface(SPOT, VOL, call_100, label = "Call Price at Maturity minus 100 time points")
c2._facecolors2d=c2._facecolors3d
c2._edgecolors2d=c2._edgecolors3d

# Third
call_199 = crv_y[130,:,:-3]
c3 = ax.plot_surface(SPOT, VOL, call_199, label="Call Price at Maturity")
c3._facecolors2d=c3._facecolors3d
c3._edgecolors2d=c3._edgecolors3d

'''
# And forth
call_199 = crv_y[199,:,:-3]
c4 = ax.plot_surface(SPOT, VOL, call_199, label="Call Price at Maturity minus 199 time points")
c4._facecolors2d=c4._facecolors3d
c4._edgecolors2d=c4._edgecolors3d
'''

ax.legend() # -> error : 'AttributeError: 'Poly3DCollection' object has no attribute '_edgecolors2d''


# labeling the figure
fig.suptitle("SABR PDE dynamics (Double Sweep, CN, nonuniform grid)")
#plt.xlabel('g2', fontsize=14)
#plt.ylabel('g3', fontsize=14)
ax.set_xlabel(r'$S$', fontsize=15, rotation=60)
ax.set_ylabel('$v$', fontsize=15, rotation=60)
ax.set_zlabel('$Price$', fontsize=15, rotation=60)
#plt.savefig('Constraints.jpg')
plt.show()

