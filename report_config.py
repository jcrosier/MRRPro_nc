import os
import sys
import numpy as np
import configparser
import netCDF4 as nC
import matplotlib.pyplot as plt
from mrr_finder import is_valid_data_file
from tqdm import tqdm

CONFIG_FILENAME = './config.ini'
CONFIG_PATHS = 'MRR_PATHS'
CONFIG_INPUT_PATH = 'INPUT_PATH'
CONFIG_OUTPUT_PATH = 'DIAG_PLOT_PATH'

PLOT_BASE_TITLE = 'MRR-Pro config from netCDF Archive: '

DATA_UNIT_CHECK = 0
DATA_VAR_LIST_CHECK = 1
DATA_DIM_SIZE = 2
DATA_DELTA_VAL = 3

DATA_FORMAT = {
    'Z_UNITS': {'var': 'Z', 'type': DATA_UNIT_CHECK, 'label': 'Z-dBZ', 'row': 0, 'units': 'dBZ'},
    'SPEC_R': {'var': 'spectrum_raw', 'type': DATA_VAR_LIST_CHECK, 'label': 'R-SPEC', 'row': 1, 'units': 'dBZ'},
    'SPEC_Z': {'var': 'spectrum_reflectivity', 'type': DATA_VAR_LIST_CHECK, 'label': 'Z-SPEC', 'row': 2, 'units': 'dBZ'},
    'SPEC_N': {'var': 'spectrum_n_samples', 'type': DATA_DIM_SIZE, 'label': 'N-SPEC', 'row': 3, 'units': 'dBZ'},
    'RANG_N': {'var': 'range', 'type': DATA_DIM_SIZE, 'label': 'N-Rng', 'row': 4, 'units': 'dBZ'},
    'TIME_DT': {'var': 'time', 'type': DATA_DELTA_VAL, 'label': 'dt', 'row': 5, 'units': 'dBZ'},
    'RANG_DR': {'var': 'range', 'type': DATA_DELTA_VAL, 'label': 'dR', 'row': 6, 'units': 'dBZ'}
}

N_DATA_ROWS = len(DATA_FORMAT)


class NcFileProperties:
    def __init__(self, name, nc_id):
        self.file_name = name
        self.time_start = int(nc_id.variables['time'][0])
        self.time_end = int(nc_id.variables['time'][nc_id.dimensions['time'].size - 1])
        self.var_list = [var for var in nc_id.variables]
        self.data = np.empty(N_DATA_ROWS, dtype=np.int16)

        for item in DATA_FORMAT:
            var = DATA_FORMAT[item]['var']
            units = DATA_FORMAT[item]['units']
            row = DATA_FORMAT[item]['row']
            method = DATA_FORMAT[item]['type']
            self.data[row] = get_nc_diagnostic_vals(nc_id, var, units, self.var_list, method)


def get_nc_diagnostic_vals(file_id, var_name, unit_str, var_list, data_type):
    if data_type == DATA_UNIT_CHECK:
        if file_id.variables[var_name].getncattr('units') == unit_str:
            return 1
        else:
            return 0
    elif data_type == DATA_VAR_LIST_CHECK:
        if var_name in var_list:
            return 1
        else:
            return 0
    elif data_type == DATA_DIM_SIZE:
        return file_id.dimensions[var_name].size
    elif data_type == DATA_DELTA_VAL:
        try:
            delta_val = round(file_id.variables[var_name][1] - file_id.variables[var_name][0])
        except IndexError:
            delta_val = 0
        return delta_val


def read_nc_diagnostics(path, filename):
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


def plot_data(time_start, time_end, diagnostic_data, figure_path, sub_path):

    plt.ion()
    fig, axs = plt.subplots(N_DATA_ROWS, sharex=True)
    fig.suptitle(PLOT_BASE_TITLE+sub_path)

    for item in DATA_FORMAT:
        row = DATA_FORMAT[item]['row']
        label = DATA_FORMAT[item]['label']
        axs[row].plot(time_start, diagnostic_data[row], 'rx')
        axs[row].plot(time_end, diagnostic_data[row], 'b+')
        axs[row].set_ylabel(label, fontsize=8, labelpad=10)
        if row < (N_DATA_ROWS-1):
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
        data_block = np.empty((N_DATA_ROWS, len(file_list)), dtype=np.uint16)
        start_times = np.empty(len(file_list), dtype='datetime64[s]')
        end_times = np.empty(len(file_list), dtype='datetime64[s]')
        for i, file in enumerate(file_list):
            file_data = read_nc_diagnostics(current_path, file)
            data_block[:, i] = file_data.data[:]
            start_times[i] = np.datetime64('1970-01-01') + np.timedelta64(file_data.time_start, 's')
            end_times[i] = np.datetime64('1970-01-01') + np.timedelta64(file_data.time_end, 's')
        plot_data(start_times, end_times, data_block, out_path, sub_path)
        # data_array = np.concatenate((data_array, data_block), axis=1)
        # time_array = np.concatenate((time_array, time_block), axis=0)
    # return data_array, time_array


def config_okay(config_data):
    try:
        if os.path.isdir(config_data[CONFIG_PATHS][CONFIG_INPUT_PATH]) is False:
            return False
        if os.path.isdir(config_data[CONFIG_PATHS][CONFIG_INPUT_PATH]) is False:
            os.mkdir(config_data[CONFIG_PATHS][CONFIG_OUTPUT_PATH])
    except KeyError:
        return False
    except FileNotFoundError:
        return False
    return True


if __name__ == "__main__":

    config = configparser.ConfigParser()
    config.read(CONFIG_FILENAME)

    if config_okay(config) is False:
        sys.exit()

    data_path = config[CONFIG_PATHS][CONFIG_INPUT_PATH]
    save_path = config[CONFIG_PATHS][CONFIG_OUTPUT_PATH]

    daily_data(data_path, save_path)
