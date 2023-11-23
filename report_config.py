import os
import sys
import numpy as np
import netCDF4 as nC
import matplotlib.pyplot as plt
from main import is_valid_data_file
from tqdm import tqdm

DEBUG_OFF = 0
DEBUG_HOME = 1
DEBUG_LAPTOP = 2

DEBUG = DEBUG_LAPTOP

FIGURE_LABELS = ['Z-dBZ', 'R-Spec', 'Z-Spec', 'N-Spec', 'N-Rng', 'dt', 'dR']
N_DATA_VALS = len(FIGURE_LABELS)
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
        self.time = int(nc_id.variables['time'][0])
        self.var_list = [var for var in nc_id.variables]
        self.data = np.empty(7, dtype=np.int16)
        self.data[DATA_Z_UNITS] = 1 if nc_id.variables['Z'].getncattr('units') == "dBZ" else 0
        self.data[DATA_SPEC_R] = 1 if 'spectrum_raw' in self.var_list else 0
        self.data[DATA_SPEC_Z] = 1 if 'spectrum_reflectivity' in self.var_list else 0
        self.data[DATA_SPEC_N] = nc_id.dimensions['spectrum_n_samples'].size
        self.data[DATA_RANG_N] = nc_id.dimensions['range'].size
        self.data[DATA_TIME_DT] = round(nc_id.variables['time'][1] - nc_id.variables['time'][0])
        self.data[DATA_RANG_DR] = round(nc_id.variables['range'][1] - nc_id.variables['range'][0])


def read_nc_data(path, filename):
    nc_id = nC.Dataset(path+'/'+filename, "r")
    file_info = NcFileProperties(filename, nc_id)
    nc_id.close()
    return file_info


def dir_check(path, folder_name, length):

    if os.path.isdir(path+'/'+folder_name) is False:
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


def plot_data(time_data, diagnostic_data, figure_path, sub_path):

    n_rows = diagnostic_data.shape[0]
    plt.ion()
    fig, axs = plt.subplots(n_rows, sharex=True)
    fig.suptitle('MRR-Pro config from netCDF Archive: '+sub_path)
    for row in range(n_rows):
        axs[row].plot(time_data, diagnostic_data[row], 'o')
        axs[row].set_ylabel(FIGURE_LABELS[row], fontsize=8, labelpad=10)
        if row < (n_rows-1):
            axs[row].tick_params('x', labelbottom=False)
        else:
            axs[row].tick_params('x', labelsize=8)
    plt.show(block=False)
    plt.savefig(figure_path + '/' + sub_path + '.png')
    plt.close()


def daily_data(path, out_path):
    mrr_folders = []
    # data_array = np.empty((N_DATA_VALS, 0), dtype=np.uint16)
    # time_array = np.empty(0, dtype='datetime64[s]')

    folders = sorted([folder for folder in os.listdir(path) if dir_check(path, folder, 6)])
    for folder in folders:
        try:
            os.mkdir(out_path+'/'+folder)
        except FileExistsError:
            pass
        current_path = os.path.join(path, folder)
        sub_folders = [sub_folder for sub_folder in os.listdir(current_path) if dir_check(current_path, sub_folder, 8)]
        sub_folders.sort()
        mrr_folders.extend(['/' + folder + '/' + sub_folder for sub_folder in sub_folders])

    num_folders = len(mrr_folders)
    for j in tqdm(range(num_folders)):
        sub_path = mrr_folders[j]
        current_path = path + sub_path
        file_list = sorted([file for file in os.listdir(current_path) if is_valid_data_file(current_path + '/' + file)])
        data_block = np.empty((N_DATA_VALS, len(file_list)), dtype=np.uint16)
        time_block = np.empty(len(file_list), dtype='datetime64[s]')
        for i, file in enumerate(file_list):
            file_data = read_nc_data(current_path, file)
            data_block[:, i] = file_data.data[:]
            time_block[i] = np.datetime64('1970-01-01') + np.timedelta64(file_data.time, 's')
        plot_data(time_block, data_block, out_path, sub_path)
        # data_array = np.concatenate((data_array, data_block), axis=1)
        # time_array = np.concatenate((time_array, time_block), axis=0)
    # return data_array, time_array


if __name__ == "__main__":

    if DEBUG == DEBUG_LAPTOP:
        data_path = 'C:/FirsData/MRR/MRR/MRRPro91'
    elif DEBUG == DEBUG_HOME:
        data_path = 'C:/Users/jonny/Desktop/2023'
    else:
        if len(sys.argv) < 2:
            sys.exit()
        data_path = path_from_argv([1])
        if len(data_path) == 0:
            sys.exit()

    save_path = data_path + '/diagnostic_plots'
    try:
        os.mkdir(save_path)
    except FileExistsError:
        pass
    else:
        sys.exit()

    daily_data(data_path, save_path)
