import os.path
import netCDF4 as nC
import time
import datetime

root_input_folder = "C:\\FirsData\\MRR\\MRR\\MRRPro91"
root_output_folder = "C:\\FirsData\\MRR\\MRR\\realtime"
export_field_list = ["time", "range", "Z"]

nc_file_ext = ".nc"
title_att_name = "title"
title_att_value = "METEK MRR Pro"
export_temp_name = "temp_"


def list_mrrpro_nc(path):
    """
    List at all the valid nc files in a given path.

    Parameters
    ----------
    path : string
        Path to a folder to check for MRR-Pro nc data

    Returns
    -------
    nc_file_list    : list
        a list of mrrpro nc files
    """
    nc_file_list = []
    file_list = os.listdir(path)
    for file in file_list:
        full_path = os.path.join(path, file)
        if valid_data_file(full_path):
            nc_file_list.append(full_path)
    return nc_file_list


def valid_data_file(file_path):
    """
    check the input file is a valid MRR-Pro netcdf:
    1) checks the file exists
    2) checks the suffix is .nc
    3) checks the contents of the "title" global attribute in the file to see presence of "MRR-Pro"

    Parameters
    ----------
    file_path : string
        full path to an individual file

    Returns
    -------
        True    : file PASSES checks: file IS a valid MRR-Pro file
        False   : file FAILS checks: file IS NOT a valid MRR-Pro file
    """
    if os.path.exists(file_path) is False:
        return False
    if file_path.endswith(nc_file_ext) is False:
        return False
    nc_file = nC.Dataset(file_path, "r")
    attribute_title = getattr(nc_file, title_att_name)
    nc_file.close()
    if title_att_value in attribute_title:
        return True
    else:
        return False


def export_fields(variable_list, file_path, output_fold):
    """
    Export subsets of data from nc files
    Look at all the valid nc files in a given path....
    Check if some corresponding file exists in the output path
    If no corresponding output file exists, create it
    Check which files can be merged based on dataset info (variable names, dtypes, units and dimensions)
    Merge the files where possible

    Parameters
    ----------
    variable_list : list
        string containing names of the fields to export
    file_path : string
        full path of source data file
    output_fold   : string
        path to an output directory to store the extracted subset of data

    Returns
    -------
    None

    """

    input_id = nC.Dataset(file_path, "r")

    output_id = nC.Dataset(output_fold, "w", format="NETCDF4_CLASSIC")

    t_dim_input = input_id.dimensions["time"]
    r_dim_input = input_id.dimensions["range"]

    output_id.createDimension("time", len(t_dim_input))
    output_id.createDimension("range", len(r_dim_input))

    for var in variable_list:
        try:
            var_input = input_id.variables[var]
            var_data = var_input[:]
            dimension = input_id.variables[var].dimensions
            data_type = input_id.variables[var].dtype
            var_output = output_id.createVariable(str(var), data_type, dimension, compression='zlib', complevel=9)
            var_output[:] = var_data
        except:
            continue

    input_id.close()
    output_id.close()


# Press the green button in the gutter to run the script.
if __name__ == '__main__':

    now = datetime.datetime.now()
    
    if now.hour < 2:
        number_of_days = 2
    else:
        number_of_days = 1
    
    for day in range(number_of_days):

        current_date = str(datetime.date.today() - datetime.timedelta(days=day))
        year = current_date[0:4]
        month = current_date[5:7]
        day = current_date[8:10]

        input_folder = root_input_folder + "\\" + year + month + "\\" + year + month + day

        if os.path.exists(input_folder) is False:
            continue

        nc_list = list_mrrpro_nc(input_folder)

        for file in nc_list:
            output_path = root_output_folder + "\\" + export_temp_name + os.path.basename(file)
            if os.path.exists(output_path) is False:
                export_fields(export_field_list, file, output_path)

    time_now = time.time()
    for file in os.listdir(root_output_folder):
        file_path = os.path.join(root_output_folder, file)
        time_file = os.stat(file_path).st_mtime
        if time_file < (time_now - number_of_days * 86400):
            if os.path.isfile(file_path):
                os.remove(file_path)
