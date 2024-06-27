import netCDF4 as nC
import datetime
import numpy
import os
import sys

from mrr_finder import files_in_path

# These globals specify the required order of commandline arguements
ARG_FIELDS = 1                  # String containing the names of fields/vars to extract: e.g. "Z,ML,"
ARG_INPUT = 2                   # String specifying the root directory containing source nc data files: e.g. "C:/Users/user1/MRRPro/"
ARG_OUTPUT = 3                  # String specifying the destination path for the exported txt data: e.g. "C:/Users/user1/Desktop/"
ARG_AVERAGE = 4                 # Int specifying any averaging time (seconds) for the extracted data@ e.g. 60

REQUIRED_DIMS_LIST = ["time"]   # We only consider vars with 'time' as the first dimension
REQUIRED_NUM_DIMS = 2           # We only consider vars with 2 dims (2D data)


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


def get_latest_mrr_netcdf(input_path: str, max_days_prior: int = 1) -> str:
    """ Find the newest data-containing source data folder, starting from the current data, and optionally scanning back in time.
    Default is to scan back 1 day to ensure we can get latest data if it is pass midnight, but no new file has been created yet. """
    for day in range(max_days_prior + 1):
        newest_file = get_newest_file_for_date(datetime.date.today() - datetime.timedelta(days=day), input_path)
        if len(newest_file) == 0: return newest_file

    print("Error: No source netCDF files found")
    return ""


def get_newest_file_for_date(date: datetime.date, root: str) -> str:
    """ Find the newest file in a source folder for a given date. """
    path_str = mrrpro_folder_from_date(datetime.date.today(), root)
    file_list = files_in_path(path_str)
    if len(file_list) > 0:
        return max(file_list, key=os.path.getctime)
    else:
        return ""
    

def mrrpro_folder_from_date(date: datetime.date, root: str) -> str:
    """ Construct expected full path to MRR-Pro source folder for a given date. """
    folder = str(date.year) + str(date.month).zfill(2)
    subfolder = folder + str(date.day).zfill(2)
    return root + "\\" + folder + "\\" + subfolder + "\\"


def number_of_nums(values: numpy.ndarray) -> int:
    """ Works out the range of valid data points in a 1-d np array - ignores trailing nans. """
    dim_len = len(values)
    last_number = 0
    for j in range(dim_len-1, -1, -1):
        if numpy.isnan(values[j]) is False:
            last_number = j
            break
    return last_number


def check_var_requirements(nc_id: nC.Dataset, vars: list[str], dims: list[str], n_dims_required: int) -> bool:
    """ Checks that all requested variables can be found and have required dimensional requirements. """
    if nc_vars_present(nc_id, vars + dims) is False:
        print('Error: Variable not found in source netCDF file')
        return False
    if nc_var_dim_dependancy(nc_id, vars, dims) is False:
        print(f'Error: Variable does not have required dependency on dims: {dims}')
        return False
    if nc_var_num_dim_criteria(nc_id, vars, n_dims_required) is False:
        print(f'Error: Variable does not have required number of dims. Must have: {n_dims_required} dims')
        return False
    
    return True


def nc_vars_present(file_id: nC.Dataset, requested_vars: list[str]) -> bool:
    """ Checks an nC file to see if it contains the specificed requested_vars. """
    nc_vars = [nc_var for nc_var in file_id.variables]
    return all([var in nc_vars for var in requested_vars])


def nc_var_dim_dependancy(nc_id: nC.Dataset, requested_vars: list[str], required_dim_list: list[str]) -> bool:
    """ Checks all variables are dependant on specified dimensions. """
    for var in requested_vars:
        var_dims = nc_id.variables[var].dimensions
        if all([required_dim in var_dims for required_dim in required_dim_list]) is False:
            return False
        
    return True


def nc_var_num_dim_criteria(nc_id: nC.Dataset, requested_vars: list[str], num_dims: int) -> bool:
    """ Checks all variables have the required number of dimensions. """
    return [len(nc_id.variables[var].dimensions)==num_dims for var in requested_vars]


def parameter_check(params: list, strings: list[int] = None, paths: list[int] = None, ints: list[int] = None) -> bool:
    """ Checks validity of command line input parameters. """
    string_params = [par for idx, par in enumerate(params) if idx in strings]
    if (strings is not None) and (list_contains_strings(string_params) is False):
        print('Error: Expected string inputs are not strings')
        return False
    path_params = [par for idx, par in enumerate(params) if idx in paths]
    if (paths is not None) and (list_contains_paths(path_params) is False):
        print('Error: Expected path inputs do not point to valid destinations')
        return False
    int_params = [par for idx, par in enumerate(params) if idx in ints]
    if (ints is not None) and (list_contains_ints(int_params) is False):
        print('Error: Expected int inputs are not ints')
        return False
    
    return True


def list_contains_strings(list: list) -> bool:
    """ Checks if all items in a list are of type str. """
    return all([isinstance(item, str) for item in list])


def list_contains_ints(list: list) -> bool:
    """ Check if all items in a list are of type int. """
    return all([isinstance(item, int) for item in list])


def list_contains_paths(list: list) -> bool:
    """ Checks if all items in a list are valid paths. """
    if list_contains_strings(list) is False:
        return False
    return all([os.path.exists(path) for path in list])


if __name__ == '__main__':

    if parameter_check(sys.argv, [ARG_FIELDS], [ARG_INPUT, ARG_OUTPUT], [ARG_AVERAGE]) is False: sys.exit()

    latest_nc_file = get_latest_mrr_netcdf(sys.argv[ARG_INPUT])
    if len(latest_nc_file) == 0: sys.exit()

    nc_id = nC.Dataset(latest_nc_file, "r", format="NETCDF4_CLASSIC")
    variable_list = sys.argv[ARG_FIELDS].split(',')

    if check_var_requirements(nc_id, variable_list, REQUIRED_DIMS_LIST, REQUIRED_NUM_DIMS) is False: nc_id.close(); sys.exit()

    for variable in variable_list:

        with open(sys.argv[ARG_OUTPUT], 'w') as output_id:
            export_data(nc_id, output_id, variable, os.path.basename(latest_nc_file)[2:8], sys.argv[ARG_AVERAGE])

    nc_id.close()
    