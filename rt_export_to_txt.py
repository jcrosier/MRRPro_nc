import netCDF4 as nC
import numpy
import os
import sys

from typing import TextIO
from file_utils import get_latest_netcdf
from nc_util import load_var, get_numpy_data_slice, check_var_requirements
from cli_utils import parameter_check

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


def get_increment_idx(nc_id: nC.Dataset, var: str, average_period: int = None) -> int:
    """ Calculate the number of elements which span an averaging period. """
    if average_period ==0: return 1
    increment = round(nc_id.variables[var][1] - nc_id.variables[var][0])
    if average_period < increment: return 1
    return int(average_period/increment)


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


def trimmed_number_of_nums(values: numpy.ndarray) -> int:
    """ Works out the range of valid data points in a 1-d np array - ignores trailing nans. """
    for j in range(len(values)-1, -1, -1):
        if numpy.isnan(values[j]):
            continue
        else:
            return j+1
    return 0


if __name__ == '__main__':

    if parameter_check(sys.argv, ARG_NUM, [ARG_INPUT, ARG_OUTPUT], [ARG_AVERAGE]) is False: sys.exit()

    latest_nc_file = get_latest_netcdf(sys.argv[ARG_INPUT])
    if len(latest_nc_file) == 0: sys.exit()

    nc_id = nC.Dataset(latest_nc_file, "r", format="NETCDF4_CLASSIC")
    variable_list = sys.argv[ARG_FIELDS].split(',')

    if check_var_requirements(nc_id, variable_list, REQUIRED_DIMS_LIST, REQUIRED_NUM_DIMS) is False: nc_id.close(); sys.exit()

    for variable in variable_list:

        with open(sys.argv[ARG_OUTPUT] + sys.argv[ARG_NAME] + '_' + variable + '.txt', 'w') as output_id:
            export_data_text(nc_id, output_id, variable, os.path.basename(latest_nc_file)[2:8], int(sys.argv[ARG_AVERAGE]))

    nc_id.close()
    