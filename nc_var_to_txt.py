import netCDF4 as nC
import numpy
import datetime
import os
import sys


def get_latest_input(input_path_str):
    path_str = get_folder_from_date(datetime.date.today(), input_path_str)
    if os.path.exists(path_str):
        file_list = [os.path.join(path_str, filename) for filename in os.listdir(path_str) if filename.endswith(".nc")]
        if len(file_list) > 0:
            latest_file = max(file_list, key=os.path.getctime)
            return latest_file

    path_str = get_folder_from_date(datetime.date.today() - datetime.timedelta(days=1), input_path_str)
    if os.path.exists(path_str):
        file_list = [os.path.join(path_str, filename) for filename in os.listdir(path_str) if filename.endswith(".nc")]
        if len(file_list) > 0:
            latest_file = max(file_list, key=os.path.getctime)
            return latest_file

    return ""


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


if __name__ == '__main__':

    if len(sys.argv) < 5:
        sys.exit()

    try:
        field_name = str(sys.argv[1])
        average_int_seconds = int(sys.argv[2])
        input_path = str(sys.argv[3])
        output_path = str(sys.argv[4])

    except ValueError:
        sys.exit()

    file = get_latest_input(input_path)
    if len(file) > 0:
        date_string = os.path.basename(file)[2:8]
        with nC.Dataset(file, "r", format="NETCDF4_CLASSIC") as input_id:
            input_vars = [var for var in input_id.variables]
            if field_name in input_vars:
                with open(output_path, 'w') as output_id:
                    export_data(input_id, output_id, field_name, date_string, average_int_seconds)
