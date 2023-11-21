import os
import sys
import numpy as np
import netCDF4 as nC
from main import is_valid_data_file

N_DATA_VALS = 7
DATA_Z_UNITS = 0
DATA_SPEC_R = 1
DATA_SPEC_Z = 2
DATA_SPEC_N = 3
DATA_RANG_N = 4
DATA_TIME_DT = 5
DATA_RANG_DR = 6


class NcFileProperties:
    def __init__(self, name, nc_id):
        self.file_name = name
        self.time = nc_id.variables['time'][0]
        self.var_list = [var for var in nc_id.variables]
        self.data = np.empty(7)
        self.data[DATA_Z_UNITS] = 1 if nc_id.variables['Z'].getncattr('units') == "dBZ" else 0
        self.data[DATA_SPEC_R] = 1 if 'spectrum_raw' in self.var_list else 0
        self.data[DATA_SPEC_Z] = 1 if 'spectrum_reflectivity' in self.var_list else 0
        self.data[DATA_SPEC_N] = nc_id.dimensions['spectrum_n_samples'].size
        self.data[DATA_RANG_N] = nc_id.dimensions['range'].size
        self.data[DATA_TIME_DT] = nc_id.variables['time'][1] - nc_id.variables['time'][0]
        self.data[DATA_RANG_DR] = nc_id.variables['range'][1] - nc_id.variables['range'][0]


def read_nc_data(path, filename):
    nc_id = nC.Dataset(path+'/'+filename, "r")
    file_info = NcFileProperties(filename, nc_id)
    nc_id.close()
    return file_info


def dir_check(path, folder_name, length):

    if os.path.isdir(path+'/'+folder_name) is False:
        print(folder_name,"1")
        return False

    if len(folder_name) != length:
        return False

    try:
        year_value = int(folder_name[0:4])
        month_value = int(folder_name[4:6])
        if year_value < 2020 or year_value > 2050:
            return False
        if month_value < 1 or month_value > 12:
            return False
        if length == 8:
            day_value = int(folder_name[6:8])
            if day_value < 1 or day_value > 31:
                return False
        return True
    except ValueError:
        return False


def path_from_argv(arg_val):
    try:
        path = str(arg_val)
        if os.path.exists(path) is False:
            return ""
        if os.path.isdir(path) is False:
            return ""
    except ValueError:
        return ""
    return path


if __name__ == "__main__":

    # if len(sys.argv) < 2:
    #     sys.exit()
    # data_path = path_from_argv([1])
    # if len(data_path) == 0:
        #sys.exit()
    data_path = 'C:/Users/jonny/Desktop/2022'

    mrr_folders = []
    destination_array = np.empty((N_DATA_VALS, 0))

    folders = sorted([folder for folder in os.listdir(data_path) if dir_check(data_path, folder, 6)])
    for folder in folders:
        current_path = os.path.join(data_path, folder)
        sub_folders = [sub_folder for sub_folder in os.listdir(current_path) if dir_check(current_path, sub_folder, 8)]
        sub_folders.sort()
        mrr_folders.extend(['/' + folder + '/' + sub_folder for sub_folder in sub_folders])

    for sub_path in mrr_folders:
        current_path = data_path + sub_path
        file_list = sorted([file for file in os.listdir(current_path) if is_valid_data_file(current_path+'/'+file)])
        data_block = np.empty((N_DATA_VALS, len(file_list)))
        for i, file in enumerate(file_list):
            file_data = read_nc_data(current_path, file)
            data_block[:, i] = file_data.data[:]

        # Copy buffer into destination array
        # Todo: This does not properly copt the shape of the dest data into a continuous blaock
        # Todo: Time stamp data from each file are passed into the file_data class, but not sued in any way. It is needed for all the diagnostics plots!
        destination_array = np.hstack((destination_array, data_block))

        print(destination_array)

    # todo Need to make some basic plots to show when settings are changing
    # todo Output plots need some time resolution/interval so things can be seen
    # Export results/plots