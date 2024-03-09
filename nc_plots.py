import os
import netCDF4 as nC
import matplotlib.pyplot as plt
import numpy as np

# todo
# implement pcolormesh in place of imshow
# alter class definitions and class implementation to allow re-use of class instance for improved performance


class NcVariable:
    def __init__(self, file_ref, var_name, transpose=False):
        self.name = ''
        self.units = ''
        self.data, self.dataT = np.zeros(0), np.zeros(0)
        self.min, self.max = np.nan, np.nan
        self.valid = False

        self.initialise_var(file_ref, var_name)
        if self.valid:
            self.load_data(file_ref, var_name, transpose)

    def initialise_var(self, file_id, var):
        try:
            self.name = file_id.variables[var].name
            self.units = file_id.variables[var].units
        except (NameError, KeyError, AttributeError):
            self.valid = False
        else:
            self.valid = True

    def load_data(self, file_id, var, tran=False):
        if self.name == 'time':
            self.load_time(file_id, var)
        else:
            self.load_numeric(file_id, var, tran)

    def load_numeric(self, file_id, var, tran=False):
        try:
            self.data = file_id.variables[var][:]
        except (NameError, KeyError, AttributeError):
            self.valid = False
        else:
            self.min = np.nanmin(self.data)
            self.max = np.nanmax(self.data)
            self.dataT = self.data.T if tran else np.zeros(0)
            self.valid = True

    def load_time(self, file_id, var):
        try:
            self.data = nC.num2date(file_id.variables[var][:], file_id.variables[var].units)
        except (NameError, KeyError, AttributeError):
            self.valid = False
        else:
            self.min = self.hours_from_time(self.data[0])
            self.max = self.hours_from_time(self.data[len(self.data) - 1])
            self.units = 'hour of day'
            self.dataT = np.zeros(0)
            self.valid = True

    @staticmethod
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
    im = plt.imshow(z.dataT, interpolation='None', aspect='auto', origin='lower', vmin=z_min_max[0], vmax=z_min_max[1],
                    extent=(x.min, x.max, y.min, y.max))
    cbar = plt.colorbar(im)
    cbar.set_label(z.name + ' (' + z.units + ')')
    plt.show(block=False)
    plt.savefig(os.path.splitext(file)[0] + '.png')
    plt.close()


path = 'C:/Users/jonny/Desktop/Data/MRR/202303/20230305/20230305_222000.nc'
#file_path = [path+'/'+file for file in os.listdir(path) if file.endswith('.nc')][0]
plot_2d_data(path, 'Ze', [-10, 30])
