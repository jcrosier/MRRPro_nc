import os
import netCDF4 as nC
import sys
import datetime
import numpy

NC_FILE_EXT = ".nc"
NC_OUTPUT_PREFIX = "MRRPro_"
NC_MERGE_DIMENSION = "time"
NC_ATT_TITLE_NAME = "title"
NC_ATT_TITLE_VALUE = "METEK MRR Pro"
NC_ATT_FIELDNAMES = "field_names"
NC_MAX_CHUNK_TDIM = 300
NC_COMPRESSION_SETTING = 2

# todo 1: add docstrings where missing


def is_valid_data_folder(path):
    """
    check the input path contains valid MRR-Pro netcdf files:
    1) checks the path exists
    2) get list of all files in the path, and check they are valid files by calling valid_data_file()

    Parameters
    ----------
    path : string
        Full path to a folder which we want to search for valid MRR-Pro files

    Returns
    -------
    Boolean :
        True    : folder PASSES checks: folder DOES contain valid MRR-Pro file(s)
        False   : folder FAILS checks: folder DOES NOT contain valid MRR-Pro file(s)
    """
    if os.path.exists(path) is False:
        return False
    else:
        file_list = os.listdir(path)
        for file in file_list:
            if is_valid_data_file(os.path.join(path, file)):
                return True
        return False


def is_valid_data_file(file):
    """
    check the input file is a valid MRR-Pro netcdf:
    1) checks the file exists
    2) checks the suffix is .nc
    3) checks the contents of the "title" global attribute in the file to see presence of "MRR-Pro"

    Parameters
    ----------
    file : string
        full path to an individual file

    Returns
    -------
        True    : file PASSES checks: file IS a valid MRR-Pro file
        False   : file FAILS checks: file IS NOT a valid MRR-Pro file
    """
    if os.path.exists(file) is False:
        return False

    if file.endswith(NC_FILE_EXT) is False:
        return False

    try:
        nc_file = nC.Dataset(file, "r")
    except OSError:
        return False

    try:
        attribute_title = getattr(nc_file, NC_ATT_TITLE_NAME)
        nc_file.close()
    except AttributeError:
        nc_file.close()
        return False

    if NC_ATT_TITLE_VALUE in attribute_title:
        return True
    else:
        return False


def create_merge_list(path):
    """
    Look at all the valid nc files in a given path....
    Check which files can be merged based on dataset info (variable names, dtypes, units and dimensions)
    Return a nested list which groups together files which can be merged

    Parameters
    ----------
    path : string
        string containing a path pointing to MRR-Pro data

    Returns
    -------
    master_merge_list    : list
        a list of lists, each sublist contains a list of files which can be merged
    """
    nc_file_list = []

    file_list = os.listdir(path)
    for file in file_list:
        full_path = os.path.join(path, file)
        if is_valid_data_file(full_path):
            nc_file_list.append(full_path)

    number_of_nc_files = len(nc_file_list)
    master_merge_list = []
    current_merge_list = []

    for idx, file in enumerate(nc_file_list):
        if len(current_merge_list) == 0:
            current_merge_list.append(file)
        else:
            merge_test_result = merge_test(current_merge_list[0], file)
            if merge_test_result is True:
                current_merge_list.append(file)
            else:
                master_merge_list.append(current_merge_list)
                current_merge_list = [file]
        if idx == (number_of_nc_files - 1):
            master_merge_list.append(current_merge_list)

    return master_merge_list


def merge_test(base_file, ref_file):
    """
    Check to see if the files are compatible for merging
    The following are all checked and must pass check to be suitable for a merge:
        1) identical variable names listed in the fieldnames attributes
        2) identical dimension names and lengths (excluding the global NC_MERGE_DIMENSION)
        3) identical variables and variable attributes

    Parameters
    ----------
    base_file   : string
        full path to a base nc file which is used as the reference in the comparison
    ref_file    : string
        full path to another nc file for the comparison

    Returns
    -------
    result  : Boolean
        True    : files ARE compatible
        False   : files ARE NOT compatible
    """
    result = True

    # open the nc files
    base_file_id = nC.Dataset(base_file, "r")
    ref_file_id = nC.Dataset(ref_file, "r")

    # check the relevant global attribute to check the list of variables is the same
    try:
        base_field_names = base_file_id.getncattr(NC_ATT_FIELDNAMES)
        ref_field_names = ref_file_id.getncattr(NC_ATT_FIELDNAMES)
        if base_field_names != ref_field_names:
            result = False
    except AttributeError:
        result = False

    # check the dimensions are consistent
    if result:
        nc_dims_base = [dim for dim in base_file_id.dimensions]
        nc_dims_ref = [dim for dim in base_file_id.dimensions]
        if nc_dims_base != nc_dims_ref:
            result = False
        else:
            for dim in nc_dims_base:
                if dim == NC_MERGE_DIMENSION:
                    continue
                if base_file_id.dimensions[dim].size != ref_file_id.dimensions[dim].size:
                    result = False

    # get a list of all the variables in the two files and check the attributes are the same
    if result:
        nc_vars_base = [var for var in base_file_id.variables]
        nc_vars_ref = [var for var in ref_file_id.variables]
        for var_base in nc_vars_base:
            if var_base in nc_vars_ref:
                result = compare_var_attributes(base_file_id, ref_file_id, var_base, [])
                if result is False:
                    break
            else:
                result = False
                break

    # close the nc files
    base_file_id.close()
    ref_file_id.close()

    return result


def compare_var_attributes(data1_id, data2_id, var_name, att_list):
    """
    Check to see if specified attributes in the input att_list, for a given variable, are identical or not

    Parameters
    ----------
    data1_id   : netcdf file object
        netcdf file object created using e.g. nC.Dataset("my_netcdf.nc", "r")
    data2_id   : netcdf file object
        netcdf file object created using e.g. nC.Dataset("my_netcdf.nc", "r")
    var_name    : string
        name of a variable found in the nc files
    att_list    : python list
        a list of the attributes to compare
        if att_list is empty (len==0), then automatically check all attributes

    Returns
    -------
    result  : Boolean
        True    : all attributes present and identical values
        False   : var attribute values ARE NOT identical or attributes are missing
    """
    if len(att_list) == 0:
        att_list = data1_id.variables[var_name].ncattrs()

    for att in att_list:
        try:
            att1 = data1_id.variables[var_name].getncattr(att)
            att2 = data2_id.variables[var_name].getncattr(att)
            try:
                if numpy.isnan(att1) and numpy.isnan(att2):
                    continue
            except TypeError:
                pass

            if att1 != att2:
                return False

        except AttributeError:
            return False

    return True


def merge_nc_files(group_list, out_path, base_name):
    """
    Merge all the files listed in the file list, including:
        Attributes: global and variables
        Dimensions:
        Datasets:

    Parameters
    ----------
    group_list   : list
        contains full paths of all files to be merged
    out_path    :   string
        contains path to a new output file to store the results
    base_name    : string
        basename for the output file

    Returns
    -------
    True    : file merge SUCCESSFUL
    False   : file merge FAILED
    """
    for iteration, group in enumerate(group_list):

        full_output_path = out_path + "\\" + NC_OUTPUT_PREFIX + base_name + ("_" + str(iteration) + NC_FILE_EXT)
        output_id = nC.Dataset(full_output_path, "w", format="NETCDF4_CLASSIC")

        # calculate length of the merged time dimension
        new_dim_len = total_dimension_length(group, NC_MERGE_DIMENSION)

        # open the first file to get list of vars and dims
        # create the destination dimensions, global atts, variables, and variable atts
        input_id = nC.Dataset(group[0], "r")
        copy_dimensions(input_id, output_id, NC_MERGE_DIMENSION, new_dim_len)
        copy_attributes(input_id, output_id)
        input_var_list = [var for var in input_id.variables]
        for variable in input_var_list:
            var_dims = input_id.variables[variable].dimensions
            var_type = input_id.variables[variable].dtype
            shape = get_size_from_dims(input_id, var_dims, "", 0)
            if len(var_dims) >= 3 and NC_MERGE_DIMENSION in var_dims:
                shape = get_size_from_dims(input_id, var_dims, NC_MERGE_DIMENSION, min(new_dim_len, NC_MAX_CHUNK_TDIM))
            output_data = output_id.createVariable(str(variable), var_type, var_dims,
                                                   compression='zlib', complevel=NC_COMPRESSION_SETTING,
                                                   chunksizes=shape)
            if NC_MERGE_DIMENSION not in var_dims:
                if len(var_dims) > 0:
                    input_array = numpy.empty(shape, dtype=var_type)
                    input_array[:] = input_id.variables[variable][:]
                    output_data[:] = input_array
            copy_attributes(input_id.variables[variable], output_id.variables[variable])
        input_id.close()

        start_pos, end_pos = 0, 0
        for idx, file in enumerate(group):
            input_id = nC.Dataset(file, "r")
            merge_dim_points = input_id.variables[NC_MERGE_DIMENSION].size
            end_pos += merge_dim_points
            for variable in input_var_list:
                var_type = input_id.variables[variable].dtype
                var_dims = input_id.variables[variable].dimensions
                shape = get_size_from_dims(input_id, var_dims, "", 0)
                output_data = output_id.variables[variable]
                if "time_coverage_end" in variable:
                    if idx == (len(group) - 1):
                        input_array = numpy.empty(shape, dtype=var_type)
                        input_array[:] = input_id.variables[variable][:]
                        output_data[:] = input_array
                if NC_MERGE_DIMENSION in var_dims:
                    input_array = numpy.empty(shape, dtype=var_type)
                    input_array[:] = input_id.variables[variable][:]
                    if len(var_dims) == 1:
                        output_data[start_pos:end_pos] = input_array
                    elif len(var_dims) == 2:
                        output_data[start_pos:end_pos, :] = input_array
                    elif len(var_dims) == 3:
                        output_data[start_pos:end_pos, :, :] = input_array
            input_id.close()
            start_pos += merge_dim_points

        output_id.close()

    return True


def total_dimension_length(file_list, dimension):
    """
    Calculate the total number of data points in a dimension, integrating across all files in a list

    Parameters
    ----------
    file_list   : list
        contains full paths of all files to integrate over
    dimension    : string
        name of a valid netCDF dimension

    Returns
    -------
    total_len : int
        total number of datapoints integrated across the files
    """
    total_len = 0
    for file in file_list:
        input_id = nC.Dataset(file, "r")
        total_len += len(input_id.dimensions[dimension])
        input_id.close()
    return total_len


def copy_dimensions(input_nc, output_nc, override_dim, override_value):
    """
    Copies dimension data from one input netCDF to an output netCDF
    The length of one dimension can be overriden, which facilitates merging of files (time dimension)

    Parameters
    ----------
    input_nc   : netCDF id
        input which contains the original dimension data
    output_nc   : netCDF id
        destination to copy the dimension data into
    override_dim   : string
        name of a dimension to override the default size
    override_value  : int
        new value to use instead of the value found in the input

    Returns
    -------
    None
    """
    input_dims = [dim for dim in input_nc.dimensions]
    for dim in input_dims:
        if str(dim) == override_dim:
            output_nc.createDimension(str(dim), override_value)
        else:
            output_nc.createDimension(str(dim), len(input_nc.dimensions[dim]))


def copy_attributes(input_id, output_id):
    input_attributes = input_id.ncattrs()
    for attr in input_attributes:
        if "_FillValue" in str(attr):
            continue
        output_id.setncattr(str(attr), input_id.getncattr(attr))


def get_size_from_dims(data_id, dim_tuple, override_dim, override_val):
    chunk_list = []
    for item in dim_tuple:
        if item in override_dim:
            chunk_list.append(override_val)
        else:
            chunk_list.append(len(data_id.dimensions[item]))
    return chunk_list


if __name__ == '__main__':

    if len(sys.argv) < 4:
        sys.exit()

    try:
        data_path = str(sys.argv[1])
        start_day = int(sys.argv[2])
        stop_day = int(sys.argv[3])
    except ValueError:
        sys.exit()

    if not os.path.exists(data_path):
        sys.exit()

    if start_day < 0 or stop_day < 0:
        sys.exit()

    for i_day in range(start_day, stop_day + 1):
        date = datetime.date.today() - datetime.timedelta(days=i_day)
        folder_date_string = str(date.year) + str(date.month).zfill(2)
        sub_folder_date_string = folder_date_string + str(date.day).zfill(2)

        input_folder = data_path + "\\" + folder_date_string + "\\" + sub_folder_date_string
        output_folder = data_path + "\\" + folder_date_string
        if is_valid_data_folder(input_folder) is False:
            continue

        process_list = create_merge_list(input_folder)
        merge_result = merge_nc_files(process_list, output_folder, sub_folder_date_string)
