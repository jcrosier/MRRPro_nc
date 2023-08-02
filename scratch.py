import netCDF4 as nC
import numpy
import datetime
import os

INPUT_PATH = "C:\\FirsData\\MRR\\MRR\\MRRPro91\\"
OUTPUT_PATH = "C:\\FirsData\\MRR\\MRR\\realtime\\current.txt"
EXPORT_NAME = "Z"
AVERAGE_TIME = 60

# todo 1: add averging capability


def get_latest_input():
    path_str = get_folder_from_date(datetime.date.today())
    file_list = [filename for filename in os.listdir(path_str) if filename.endswith(".nc")]
    if len(file_list) == 0:
        path_str = get_folder_from_date(datetime.date.today() - datetime.timedelta(days=1))
        file_list = [filename for filename in os.listdir(path_str) if filename.endswith(".nc")]
    if len(file_list) == 0:
        return ""
    else:
        latest_file = max(file_list, key=os.path.getctime)
        return path_str + latest_file


def get_folder_from_date(date_input):
    folder_date_string = str(date_input.year) + str(date_input.month).zfill(2)
    sub_folder_date_string = folder_date_string + str(date_input.day).zfill(2)
    data_path = INPUT_PATH + folder_date_string + "\\" + sub_folder_date_string + "\\"
    return data_path


def export_data(f_input, f_output, variable, date_str):
    data_shape = get_dim_tuple(f_input, variable)

    time_increment = int(f_input.variables["time"][1]) - int(f_input.variables["time"][0])
    time_average_pnts = AVERAGE_TIME / time_increment

    for i in range(data_shape[0]):
        time_value = int(f_input.variables["time"][i]) % 86400
        f_output.write(date_str + "," + str(time_value) + ",")
        for j in range(data_shape[1]):
            val = float(f_input.variables[variable][i][j])
            if numpy.isnan(val):
                f_output.write("nan")
            else:
                f_output.write("{:.2f}".format(val))
            if j < (data_shape[1] - 1):
                f_output.write(",")
        f_output.write("\r")


def get_dim_tuple(file_id, var):
    len_list = []
    dim_names_tuple = file_id.variables[var].dimensions
    for dim in dim_names_tuple:
        len_list.append(len(file_id.dimensions[dim]))
    return len_list


if __name__ == '__main__':
    file = get_latest_input()
    if len(file) > 0:
        date_string = os.path.basename(file)[2:8]
        with nC.Dataset(file, "r", format="NETCDF4_CLASSIC") as input_id:
            with open(OUTPUT_PATH, 'w') as output_id:
                export_data(input_id, output_id, EXPORT_NAME, date_string)
