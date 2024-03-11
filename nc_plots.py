import os
import netCDF4 as nC
import matplotlib.pyplot as plt
import numpy as np

# todo
# complete pcolor mesh plots to use z max/min and other formatting
# add saftey checks to functions which create classes so we dont make plots with dead class vars
# add loops to allow batch generation of plots


class NcVariable:
    def __init__(self, file_ref, var_name, transpose=False, expand=0, offset=False):
        self.name = ''
        self.units = ''
        self.data, self.dataT = np.zeros(0), np.zeros(0)
        self.valid = False

        self.initialise_var(file_ref, var_name)
        if self.valid:
            self.load_data(file_ref, var_name, transpose)
            if expand > 0:
                self.expand_1d_np_array(expand)
            if offset:
                self.offset_1d_values()

    def initialise_var(self, file_id, var):
        try:
            self.name = file_id.variables[var].name
            self.units = file_id.variables[var].units
        except (NameError, KeyError, AttributeError):
            self.valid = False
        else:
            self.valid = True

    def load_data(self, file_id, var, tran=False):
        try:
            self.data = file_id.variables[var][:]
            if self.name == 'time':
                t0 = np.mod(self.data[0], 86400)
                self.data = np.round((self.data - self.data[0]) + t0) / 3600
                self.units = 'hour of day'
        except (NameError, KeyError, AttributeError):
            self.valid = False
        else:
            self.dataT = self.data.T if tran else np.zeros(0)
            self.valid = True

    def expand_1d_np_array(self, num_points):
        last_index = len(self.data) - 1
        increment = self.data[last_index] - self.data[last_index-1]
        new_data = np.copy(self.data[last_index - (num_points - 1):last_index+1]) + (num_points*increment)
        self.data = np.append(self.data, new_data)
        del new_data

    def offset_1d_values(self):
        offset = 0.5*(self.data[1] - self.data[0])
        self.data -= offset


def plot_2d_data(file, var, z_range):
    nc_id = nC.Dataset(file, 'r')
    # add checks to see if x_data.valid, y_data.valid etc are all true
    x_data = NcVariable(nc_id, nc_id.variables[var].dimensions[0], expand=1, offset=False)
    y_data = NcVariable(nc_id, nc_id.variables[var].dimensions[1], expand=1, offset=True)
    z_data = NcVariable(nc_id, var, transpose=True)
    nc_id.close()

    display_2d(x_data, y_data, z_data, z_range, file)


def display_2d(x, y, z, z_min_max, file):
    plt.title(os.path.basename(file))
    plt.xlabel(x.name + ' (' + x.units + ')')
    plt.ylabel(y.name + ' (' + y.units + ')')

    im = plt.pcolormesh(x.data, y.data, z.dataT)
                        #, cmap=cmap, norm=norm)

    #im = plt.imshow(z.dataT, interpolation='None', aspect='auto', origin='lower', vmin=z_min_max[0], vmax=z_min_max[1],
    #                extent=(x.min, x.max, y.min, y.max))
    cbar = plt.colorbar(im)
    cbar.set_label(z.name + ' (' + z.units + ')')
    plt.show(block=False)
    plt.savefig(os.path.splitext(file)[0] + '.png')
    plt.close()

# JC laptop test
path = 'C:/FirsData/MRR/MRR/MRRPro91/PROC_new/MRRPro_20230612_0.nc'

# JC home PC test
#path = 'C:/Users/jonny/Desktop/Data/MRR/202303/20230305/20230305_222000.nc'
plot_2d_data(path, 'Ze', [-10, 30])
