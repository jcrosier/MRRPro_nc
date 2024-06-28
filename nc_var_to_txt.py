import netCDF4 as nC
import datetime
import numpy
import os
import sys

from typing import TextIO
from mrr_finder import is_valid_data_folder, files_in_path

# Generic usage:
# C:/Users/User1/code/nc_var_to_txt.py field_list source_path dest_path ave_time_seconds output_fname
# Example usage
# C:/Users/User1/code/nc_var_to_txt.py Ze C:\Users\jonny\Desktop\tempMRR C:\Users\jonny\Desktop\ 60 output


# These globals specify the required order of commandline arguements
ARG_FIELDS = 1                  # String containing the names of fields/vars to extract: e.g. "Z,ML,"
ARG_INPUT = 2                   # String specifying the root directory containing source nc data files: e.g. "C:/Users/user1/MRRPro/"
ARG_OUTPUT = 3                  # String specifying the destination path for the exported txt data: e.g. "C:/Users/user1/Desktop/"
ARG_AVERAGE = 4                 # Int specifying any averaging time (seconds) for the extracted data@ e.g. 60
ARG_NAME = 5                    # basename for the output file
ARG_NUM = 6                     # total number of expected command line inputs

REQUIRED_DIMS_LIST = ["time"]   # We only consider vars with 'time' as the first dimension
REQUIRED_NUM_DIMS = 2           # We only consider vars with 2 dims (2D data)
AVERAGING_DIM = "time"          # This specifies the name of the var to average along

def export_data_text(f_input: nC.Dataset, f_output: TextIO, variable: str, date_str: str, average_val: int):
    """ Load a var and export the data to the output file. """
    ave_inc = get_increment_idx(nc_id, AVERAGING_DIM, average_val)
    data = load_var(f_input, variable)
    start_i = 0  
    while start_i < (data.shape[0]):
        time_value = int(f_input.variables["time"][start_i]) % 86400
        f_output.write(date_str + "," + str(time_value) + ",")
        slice = get_numpy_data_slice(data, start_i, ave_inc)
        output_numpy_slice_to_text(f_output, slice)
        start_i += ave_inc


def get_numpy_data_slice(data: numpy.ndarray, start_idx: int, average: int) -> numpy.ndarray:
    """ Get a 1D numpy array from a 2D numpy array, either by taking a single slice, or by averaging. """
    if average == 1:
        slice = data[start_idx, :]
    else:
        slice = numpy.nanmean(data[start_idx:start_idx + average, :], axis=0)
    return slice


def output_numpy_slice_to_text(f_out: TextIO, data: numpy.ndarray):
    """ Output a 1D numpy float array to a file (output with 2 decimal places). """
    num_vals = trimmed_number_of_nums(data)

    for j in range(num_vals):
        val = float(data[j])
        if numpy.isnan(val):
            f_out.write("nan")
        else:
            f_out.write("{:.2f}".format(val))
        if j < num_vals:
            f_out.write(",")
    f_out.write("\r")
        


        

def load_var(nc_id: nC.Dataset, var: str) -> numpy.ndarray:
    """ Load an entire var into a numpy array. """
    data_shape = list(nc_id.variables[var].shape)
    data_array = numpy.empty(data_shape, dtype=nc_id.variables[var].dtype)
    data_array[:] = nc_id.variables[var][:]
    return data_array


def get_increment_idx(nc_id: nC.Dataset, var: str, average_period: int = None) -> int:
    """ Calculate the number of elements which span an averaging period. """
    if average_period ==0: return 1
    increment = round(nc_id.variables[var][1] - nc_id.variables[var][0])
    if average_period < increment: return 1
    return int(average_period/increment)


def get_latest_mrr_netcdf(input_path: str, max_days_prior: int = 1) -> str:
    """ Find the newest data-containing source data folder, starting from the current data, and optionally scanning back in time.
    Default is to scan back 1 day to ensure we can get latest data if it is pass midnight, but no new file has been created yet. """
    for day in range(max_days_prior + 1):
        newest_file = get_newest_file_for_date(datetime.date.today() - datetime.timedelta(days=day), input_path)
        if len(newest_file) > 0: return newest_file

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


def trimmed_number_of_nums(values: numpy.ndarray) -> int:
    """ Works out the range of valid data points in a 1-d np array - ignores trailing nans. """
    for j in range(len(values)-1, -1, -1):
        if numpy.isnan(values[j]):
            continue
        else:
            return j+1
    return 0


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


def parameter_check(params: list[str], number_expected: int, paths: list[int] = None, ints: list[int] = None) -> bool:
    """ Checks validity of command line input parameters. """
    if len(params) < number_expected:
        print('Error: Not enough Command Line parameters')
        return False
    path_params = [par for idx, par in enumerate(params) if idx in paths]
    if list_contains_paths(path_params) is False:
        print('Error: Expected path inputs do not point to valid destinations')
        return False
    int_params = [par for idx, par in enumerate(params) if idx in ints]
    if list_contains_ints(int_params) is False:
        print('Error: Expected int inputs are not ints')
        return False
    
    return True


def list_contains_strings(list: list) -> bool:
    """ Checks if all items in a list are of type str. """
    return all([isinstance(item, str) for item in list])


def list_contains_ints(list: list) -> bool:
    """ Check if all items in a list can be converted to type int. """
    try:
        [int(item) for item in list]
    except ValueError:
        return False
    return True


def list_contains_paths(list: list) -> bool:
    """ Checks if all items in a list are valid paths. """
    if list_contains_strings(list) is False:
        return False
    return all([os.path.exists(path) for path in list])


if __name__ == '__main__':

    if parameter_check(sys.argv, ARG_NUM, [ARG_INPUT, ARG_OUTPUT], [ARG_AVERAGE]) is False: sys.exit()

    latest_nc_file = get_latest_mrr_netcdf(sys.argv[ARG_INPUT])
    if len(latest_nc_file) == 0: sys.exit()

    nc_id = nC.Dataset(latest_nc_file, "r", format="NETCDF4_CLASSIC")
    variable_list = sys.argv[ARG_FIELDS].split(',')

    if check_var_requirements(nc_id, variable_list, REQUIRED_DIMS_LIST, REQUIRED_NUM_DIMS) is False: nc_id.close(); sys.exit()

    for variable in variable_list:

        with open(sys.argv[ARG_OUTPUT] + sys.argv[ARG_NAME] + '_' + variable + '.txt', 'w') as output_id:
            export_data_text(nc_id, output_id, variable, os.path.basename(latest_nc_file)[2:8], int(sys.argv[ARG_AVERAGE]))

    nc_id.close()
    