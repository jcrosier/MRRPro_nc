import netCDF4 as nC
import numpy
import datetime
import os
import sys
from mrr_finder import files_in_path


# This specifies the required order of commandline arguements
CMD_ARG_FIELDS = 1
CMD_ARG_AVERAGE = 2
CMD_ARG_INPUT_PATH = 3
CMD_ARG_OUTPUT_PATH = 4


def export_data(f_input, f_output, variable, date_str, average_int):
    data_shape = list(f_input.variables[variable].shape)

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

        num_vals = number_of_nums(mean_values)

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


def get_latest_mrr_netcdf(input_path: str, max_days_prior: int = 0) -> str:

    for day in range(max_days_prior + 1):
        newest_file = get_newest_file_for_date(datetime.date.today() - datetime.timedelta(days=day), input_path)
        if len(newest_file) == 0: return newest_file

    print("Error: No source netCDF files found")
    return ""


def get_newest_file_for_date(date: datetime.date, root: str) -> str:

    path_str = mrrpro_folder_from_date(datetime.date.today(), root)
    file_list = files_in_path(path_str)
    if len(file_list) > 0:
        return max(file_list, key=os.path.getctime)
    else:
        return ""
    

def mrrpro_folder_from_date(date: datetime.date, root: str) -> str:

    folder = str(date.year) + str(date.month).zfill(2)
    subfolder = folder + str(date.day).zfill(2)
    return root + "\\" + folder + "\\" + subfolder + "\\"


def number_of_nums(values: numpy.ndarray) -> int:

    dim_len = len(values)
    last_number = 0
    for j in range(dim_len-1, -1, -1):
        if numpy.isnan(values[j]) is False:
            last_number = j
            break

    return last_number


def check_nc_for_var(file_id: nC.Dataset, var: str) -> bool:

    if var_name not in [var for var in file_id.variables]:
        print('Error: Variable not found in source netCDF file')
        return False
    
    return True


def parameter_check(params: list) -> bool:

    if len(params) < 5:
        print('Error: Not enough input params')
        return False
    if list_contains_strings([sys.argv[CMD_ARG_FIELDS], sys.argv[CMD_ARG_INPUT_PATH], sys.argv[CMD_ARG_OUTPUT_PATH]]) is False:
        print('Error: Expected string inputs are not strings')
        return False
    if list_contains_paths(sys.argv[CMD_ARG_INPUT_PATH], sys.argv[CMD_ARG_OUTPUT_PATH]) is False:
        print('Error: Expected path inputs do not point to valid destinations')
        return False
    if list_contains_ints([sys.argv[CMD_ARG_AVERAGE]]) is False:
        print('Error: Expected int inputs are not ints')
        return False
    
    return True


def list_contains_strings(list: list) -> bool:
    return all([isinstance(item, str) for item in list])


def list_contains_ints(list: list) -> bool:
    return all([isinstance(item, int) for item in list])


def list_contains_paths(list: list) -> bool:
    return all([os.path.exists(path) for path in list])


if __name__ == '__main__':

    if parameter_check(sys.argv) is False: sys.exit()

    latest_nc_file = get_latest_mrr_netcdf(sys.argv[CMD_ARG_INPUT_PATH])
    if len(latest_nc_file) == 0: sys.exit()

    with nC.Dataset(latest_nc_file, "r", format="NETCDF4_CLASSIC") as nc_id:

        for var_name in sys.argv[CMD_ARG_FIELDS]:
            if check_nc_for_var(nc_id, var_name) is False: continue

            with open(sys.argv[CMD_ARG_OUTPUT_PATH], 'w') as output_id:
                export_data(nc_id, output_id, var_name, os.path.basename(latest_nc_file)[2:8], sys.argv[CMD_ARG_AVERAGE])
