import os
import netCDF4 as nC
import matplotlib.pyplot as plt
import numpy as np


class NcVariable:
    def __init__(self, file_ref, var_name, transpose=False):
        try:
            self.valid = False
            self.name = file_ref.variables[var_name].name
            self.units = file_ref.variables[var_name].units
            self.data = file_ref.variables[var_name][:]
            self.min = 0
            self.max = 0

            if transpose:
                self.data = self.data.T

            if self.name == 'time':
                self.time_override(file_ref)
            else:
                self.min = np.nanmin(self.data)
                self.max = np.nanmax(self.data)
        finally:
            self.valid = True

    def time_override(self, file_ref):
        self.data = nC.num2date(file_ref.variables[self.name][:], file_ref.variables[self.name].units)
        self.min = hours_from_time(self.data[0])
        self.max = hours_from_time(self.data[len(self.data) - 1])
        self.units = 'hour of day'


def hours_from_time(date_time):
    return date_time.hour + (date_time.minute / 60) + (date_time.second / 3600)


def plot_2d_data(file, var, z_range):
    nc_id = nC.Dataset(file, 'r')
    x_data = NcVariable(nc_id, nc_id.variables[var].dimensions[0])
    y_data = NcVariable(nc_id, nc_id.variables[var].dimensions[1])
    z_data = NcVariable(nc_id, var, transpose=True)
    nc_id.close()

    display_2d(x_data, y_data, z_data, z_range, file)


def display_2d(x, y, z, z_min_max, file):
    plt.title(os.path.basename(file))
    plt.xlabel(x.name + ' (' + x.units + ')')
    plt.ylabel(y.name + ' (' + y.units + ')')
    im = plt.imshow(z.data, interpolation='None', aspect='auto', origin='lower', vmin=z_min_max[0], vmax=z_min_max[1],
                    extent=(x.min, x.max, y.min, y.max))
    cbar = plt.colorbar(im)
    cbar.set_label(z.name + ' (' + z.units + ')')
    plt.show(block=False)
    plt.savefig(os.path.splitext(file)[0] + '.png')
    plt.close()


path = '/data/2023/202303'
file_path = [path+'/'+file for file in os.listdir(path) if file.endswith('.nc')][0]
plot_2d_data(file_path, 'Ze', [-10, 30])
