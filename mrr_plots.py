import os
import netCDF4 as nC
import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm


class NcVariable:
    def __init__(self, file_ref, var_name, expand=0, offset=False, t_format="hour"):
        self.valid = False
        try:
            self.name = file_ref.variables[var_name].name
            self.units = file_ref.variables[var_name].units
            self.data = file_ref.variables[var_name][:]
            if self.name == 'time':
                self.format_time(t_format)
            if expand > 0:
                self.expand_1d_np_array(expand)
            if offset:
                self.offset_1d_values()
        except NameError:
            print(f'Error! File ref found.')
        except KeyError:
            print(f'Error! Variable not found: {var_name}')
        except AttributeError:
            print(f'Error! Attribute name or units not found for Variable: {var_name}')
        else:
            self.valid = True

    def format_time(self, time_format):
        if time_format == 'hour':
            self.data = np.round((self.data - self.data[0]) + np.mod(self.data[0], 86400)) / 3600
            self.units = 'hour of day'
        elif time_format == 'seconds':
            self.data = np.round((self.data - self.data[0]) + np.mod(self.data[0], 86400))
            self.units = 'seconds since midnight'
        # flexibility to add other elif's if required to allow other time formats to be shown

    def expand_1d_np_array(self, num_points):
        last_index = len(self.data) - 1
        increment = self.data[last_index] - self.data[last_index-1]
        new_data = np.copy(self.data[last_index - (num_points - 1):last_index+1]) + (num_points*increment)
        self.data = np.append(self.data, new_data)
        del new_data

    def offset_1d_values(self):
        offset = 0.5*(self.data[1] - self.data[0])
        self.data -= offset


def display_2d(x, y, z, file, z_min_max, show=False):
    plt.title(os.path.basename(file))
    plt.subplots(figsize=(12, 4))
    plt.xlabel(x.name + ' (' + x.units + ')')
    plt.ylabel(y.name + ' (' + y.units + ')')
    if z_min_max[0]:
        im = plt.pcolormesh(x.data, y.data, z.data.T, vmin=z_min_max[1], vmax=z_min_max[2], cmap='turbo')
    else:
        im = plt.pcolormesh(x.data, y.data, z.data.T, cmap='turbo')
    cbar = plt.colorbar(im)
    cbar.set_label(z.name + ' (' + z.units + ')')
    if show:
        plt.show(block=False)
    plt.savefig(os.path.splitext(file)[0] + '_' + z.name + '.jpg', dpi=250)
    plt.close()


def plot_2d_data(file_list, var, z_range, display=False, offset_x=False, offset_y=False):
    for j in tqdm(range(len(file_list))):
        try:
            nc_id = nC.Dataset(file_list[j], 'r')
            if len(nc_id.variables[var].dimensions) != 2:
                print(f'Error! Variable {var} does not have correct dims for 2D plot in file: {file_list[j]}')
                nc_id.close()
                break
            x_data = NcVariable(nc_id, nc_id.variables[var].dimensions[0], expand=1, offset=offset_x)
            y_data = NcVariable(nc_id, nc_id.variables[var].dimensions[1], expand=1, offset=offset_y)
            z_data = NcVariable(nc_id, var)
        except FileNotFoundError:
            print(f'Error! File not found: {file_list[j]}')
        except OSError:
            print(f'Error! Not a valid netCDF file: {file_list[j]}')
        except KeyError:
            print(f'Error! Invalid key for File: {file_list[j]}')
        else:
            nc_id.close()
            if x_data.valid and y_data.valid and z_data.valid:
                display_2d(x_data, y_data, z_data, file_list[j], z_range, display)
            else:
                print(f'Errors encountered! See prior message(s) for info. File: {file_list[j]}')
