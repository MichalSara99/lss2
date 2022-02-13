import xml_parser as parser
import matplotlib.pyplot as plt
import os


class curve_gen:
    def __init__(self,xml_file):
        self._par = parser.xml_parser(xml_file)
        return

    def set_x_label(self,x_label):
        self._x_label = x_label
        return

    def set_title(self,title):
        self._title = title
        return

    def set_legend(self,title,pos):
        self._legend = title
        self._pos = pos
        return

    def plot(self):
        plt.figure()
        plt.plot(self._par.abscissa(),self._par.ordinate())
        plt.title(self._title)
        plt.legend(self._legend,loc=self._pos)
        plt.show()


#path = os.getcwd()
#files_path = "\\".join(path.split("\\")[:-1]) + "\\xmls"
#file_crv = files_path+"\\ode_bvp_neumann_robin.xml"


#crv = curve_gen(file_crv)
#crv.set_x_label('Location')
#crv.set_title('Neumann BVP')
#crv.set_legend('Neumann BVP','lower_right')
#crv.plot()