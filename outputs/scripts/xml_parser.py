from bs4 import BeautifulSoup
import os
from enum import Enum
import numpy as np

class xml_type(Enum):
    CURVE = 1
    SURFACE_ST = 2
    SURFACE_SS = 3
    SURFACE_SST = 4
    SURFACE_SSS = 5


class xml_parser:
    def __init__(self,xml_file):
        self._xml = xml_file
        data = None
        with open(self._xml, 'r') as f:
            data = f.read()
        self.bs_obj = BeautifulSoup(data,'xml')
        return

    def type(self):
        if not self.bs_obj.find('CURVE') is None:
            return xml_type.CURVE
        elif not self.bs_obj.find('TYPE') is None:
            type = self.bs_obj.find('TYPE')
            if type.string == 'SPACE_TIME':
                return xml_type.SURFACE_ST
            elif type.string == 'SPACE_SPACE':
                return xml_type.SURFACE_SS
            elif type.string == 'SPACE_SPACE_TIME':
                return xml_type.SURFACE_SST
            elif type.string == 'SPACE_SPACE_SPACE':
                return xml_type.SURFACE_SSS
            else:
                return None
        return None

    def abscissa(self):
        if self.type() == xml_type.CURVE:
            pts = self.bs_obj.find('POINTS')
            x_points = np.asarray([float(s) for s in pts.string.split(',')])
            return x_points
        elif self.type() == xml_type.SURFACE_ST:
            # space_points
            pts = self.bs_obj.find('SPACE_POINTS')
            x = np.asarray([float(s) for s in pts.string.split(',')])
            # time_points
            pts = self.bs_obj.find('TIME_POINTS')
            t = np.asarray([float(s) for s in pts.string.split(',')])
            return np.asarray([x,t])
        elif self.type() == xml_type.SURFACE_SS:
            # space_points_0
            pts = self.bs_obj.find_all('POINTS')
            x_0 = np.asarray([float(s) for s in pts[0].string.split(',')])
            # space_points_1
            x_1 = np.asarray([float(s) for s in pts[1].string.split(',')])
            return np.asarray([x_0,x_1])
        elif self.type() == xml_type.SURFACE_SST:
            # space_points_0
            pts = self.bs_obj.find_all('POINTS')
            x_0 = np.asarray([float(s) for s in pts[0].string.split(',')])
            # space_points_1
            x_1 = np.asarray([float(s) for s in pts[1].string.split(',')])
            # time_points
            pts = self.bs_obj.find('TIME_POINTS')
            t = np.asarray([float(s) for s in pts.string.split(',')])
            return np.asarray([x_0,x_1,t])
        elif self.type() == xml_type.SURFACE_SSS:
            # space_points_0
            pts = self.bs_obj.find_all('POINTS')
            x_0 = np.asarray([float(s) for s in pts[0].string.split(',')])
            # space_points_1
            x_1 = np.asarray([float(s) for s in pts[1].string.split(',')])
            # space_points_2
            x_2 = np.asarray([float(s) for s in pts[2].string.split(',')])
            return np.asarray([x_0,x_1,x_2])
        else:
            return None

    def ordinate(self,clip_last_rows = None):
        if self.type() == xml_type.CURVE:
            pts = self.bs_obj.find('VALUES')
            values = np.asarray([float(s) for s in pts.string.split(',')])
            return values
        elif (self.type() == xml_type.SURFACE_ST
              or self.type() == xml_type.SURFACE_SS) :
            rows = int(self.bs_obj.find('ROW_SIZE').string)
            cols = int(self.bs_obj.find('COLUMN_SIZE').string)
            pts = self.bs_obj.find('VALUES')
            values = np.asarray([float(s) for s in pts.string.split(',')])
            mat = values.reshape((rows, cols))
            if not clip_last_rows is None:
                mat = mat[:, :-clip_last_rows]
            return mat
        elif (self.type() == xml_type.SURFACE_SST \
                or self.type() == xml_type.SURFACE_SSS) :
            rows = int(self.bs_obj.find('ROW_SIZE').string)
            cols = int(self.bs_obj.find('COLUMN_SIZE').string)
            lays = int(self.bs_obj.find('LAYER_SIZE').string)
            pts = self.bs_obj.find('VALUES')
            values = np.asarray([float(s) for s in pts.string.split(',')])
            mat = values.reshape((lays,rows,cols))
            if not clip_last_rows is None:
                mat = mat[:, :,:-clip_last_rows]
            return mat
        else:
            return None



# EXAMPLES:
'''
path = os.getcwd()
files_path = "\\".join(path.split("\\")[:-1]) + "\\xmls"
file_crv = files_path+"\\sabr_double_sweep_cn_srf_numerical_stepping.xml"

crv = xml_parser(file_crv)
crv_type = crv.type()
crv_y = crv.ordinate()
crv_x = crv.abscissa()

np.shape(crv_y[199,:,:])

'''
