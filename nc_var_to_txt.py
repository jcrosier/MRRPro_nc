import netCDF4 as nC
import numpy
import datetime
import os
import sys
import configparser

CONFIG_FILENAME = './config.ini'
CONFIG_GROUP = 'REALTIME_NC'
CONFIG_INPUT = 'INPUT_PATH'
CONFIG_OUTPUT = 'OUTPUT_PATH'
CONFIG_DT = "AVERAGE_INTERVAL_SECS"
CONFIG_FIELDS = 'FIELD_NAMES'
CONFIG_SCAN = 'SCAN_ARCHIVE'
CONFIG_SCAN_N = 'SCAN_MAX'


def get_latest_file(path_str):
    if os.path.exists(path_str) is False:
        return ""
    file_list = [os.path.join(path_str, filename) for filename in os.listdir(path_str) if filename.endswith(".nc")]
    if len(file_list) == 0:
        return ""
    return max(file_list, key=os.path.getctime)


def get_latest_input(input_path_str, scan, scan_num):
    latest_file = ""
    if scan:
        ref_date = datetime.date.today()
        for i in range(scan_num):
            latest_file = get_latest_file(get_folder_from_date(ref_date, input_path_str))
            if len(latest_file) > 0:
                break
            ref_date -= datetime.timedelta(days=1)
    else:
        latest_file = get_latest_file(input_path_str)
    return latest_file


def get_folder_from_date(date_input, path):
    folder_date_string = str(date_input.year) + str(date_input.month).zfill(2)
    sub_folder_date_string = folder_date_string + str(date_input.day).zfill(2)
    data_path = path + "\\" + folder_date_string + "\\" + sub_folder_date_string + "\\"
    return data_path


def export_data(f_input, f_output, variable, date_str, average_int):
    data_shape = get_dim_tuple(f_input, variable)

    time_delta = int(f_input.variables["time"][1]) - int(f_input.variables["time"][0])
    time_inc = int(average_int / time_delta)
    data_array = numpy.empty(data_shape, dtype=f_input.variables[variable].dtype)
    data_array[:] = f_input.variables[variable][:]

    start_i = 0
    more_output = True

    while more_output:
        time_value = int(f_input.variables["time"][start_i]) % 86400
        f_output.write(date_str + "," + str(time_value) + ",")
        mean_values = numpy.nanmean(data_array[start_i:start_i + time_inc, :], axis=0)

        num_vals = calc_number_of_nums(mean_values)

        for j in range(num_vals + 1):
            val = float(mean_values[j])
            if numpy.isnan(val):
                f_output.write("nan")
            else:
                f_output.write("{:.2f}".format(val))
            if j < num_vals:
                f_output.write(",")

        f_output.write("\r")
        start_i += time_inc
        if start_i >= (data_shape[0]):
            more_output = False


def calc_number_of_nums(values):
    dim_len = len(values)
    last_number = -1
    for j in range(dim_len-1, -1, -1):
        if numpy.isnan(values[j]):
            continue
        else:
            last_number = j
            break
    return last_number


def get_dim_tuple(file_id, var):
    len_list = []
    dim_names_tuple = file_id.variables[var].dimensions
    for dim in dim_names_tuple:
        len_list.append(len(file_id.dimensions[dim]))
    return len_list


def config_okay(config_data):
    try:
        if os.path.isdir(config_data[CONFIG_INPUT]) is False:
            return False
        if os.path.isdir(config_data[CONFIG_OUTPUT]) is False:
            os.mkdir(config_data[CONFIG_OUTPUT])
        if len(config_data[CONFIG_FIELDS]) == 0:
            return False
        mean_val = int(config_data[CONFIG_DT])
        if mean_val < 0 or mean_val > 1000:
            return False
        scan_val = int(config_data[CONFIG_SCAN_N])
        if scan_val < 0 or scan_val > 1000:
            return False
        bool_val = config_data[CONFIG_SCAN]
        if bool_val != "True" and bool_val != "False":
            return False
    except ValueError:
        return False
    except KeyError:
        return False
    except FileNotFoundError:
        return False
    return True


if __name__ == '__main__':

    config = configparser.ConfigParser()
    config.read(CONFIG_FILENAME)
    settings = config[CONFIG_GROUP]

    if config_okay(settings) is False:
        sys.exit()

    field_names = settings[CONFIG_FIELDS].split(',')
    file = get_latest_input(settings[CONFIG_INPUT], settings[CONFIG_SCAN], int(settings[CONFIG_SCAN_N]))

    if len(file) == 0:
        sys.exit()

    date_string = os.path.basename(file)[2:8]
    with nC.Dataset(file, "r", format="NETCDF4_CLASSIC") as input_id:
        input_vars = [var for var in input_id.variables]
        for field in field_names:
            if field not in input_vars:
                continue
            with open(settings[CONFIG_OUTPUT]+'_' + field + '.txt', 'w') as output_id:
                export_data(input_id, output_id, field, date_string, int(settings[CONFIG_DT]))
