import os
import netCDF4 as nC
import sys
import datetime
import numpy
from tqdm import tqdm
from multiprocessing import Pool

NC_FILE_EXT = ".nc"
NC_OUTPUT_PREFIX = "MRRPro_"
NC_MERGE_DIMENSION = "time"
NC_ATT_TITLE_NAME = "title"
NC_ATT_TITLE_VALUE = "METEK MRR Pro"
NC_ATT_FIELDNAMES = "field_names"
NC_MAX_CHUNK_TDIM = 300
NC_COMP_VAL = 2
N_CPU = 0


def is_valid_data_folder(path):
    """
    check if the input path string contains valid MRR-Pro netcdf files.

    Parameters
    ----------
    path : string
        Full path to a folder which we want to search for valid MRR-Pro files.

    Returns
    -------
    Boolean :
        True    : path contains valid MRR-Pro file(s)
        False   : path does not contain valid MRR-Pro file(s)
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
    check the input file string is a valid MRR-Pro netcdf file:
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
    Look at all the valid nc files in a given path string.
    Check which files can be merged based on dataset info (variable names, dtypes, units and dimensions)
    Return a nested list which groups together files which can be merged

    Parameters
    ----------
    path : string
        string containing a path pointing to a folder contain MRR-Pro data file(s)

    Returns
    -------
    master_merge_list    : list
        a list of lists, each sublist contains a list of files which can be merged
    """

    nc_file_list = [file for file in os.listdir(path) if is_valid_data_file(os.path.join(path, file))]
    number_of_nc_files = len(nc_file_list)

    master_merge_list, current_merge_list = [], []

    for idx, file in enumerate(nc_file_list):
        file_full_path = os.path.join(path, file)
        if len(current_merge_list) == 0:
            current_merge_list.append(file_full_path)
        else:
            if merge_test(current_merge_list[0], file_full_path):
                current_merge_list.append(file_full_path)
            else:
                master_merge_list.append(current_merge_list)
                current_merge_list = [file_full_path]
        if idx == (number_of_nc_files - 1):
            master_merge_list.append(current_merge_list)

    return master_merge_list


def merge_test(base_file, ref_file):
    """
    Check to see if two files are compatible for merging
    The following are all checked and must pass to be suitable for a merge:
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
    result = False

    # open the nc files
    base_file_id = nC.Dataset(base_file, "r")
    ref_file_id = nC.Dataset(ref_file, "r")

    # check the relevant global attribute to check the list of variables is the same
    if merge_test_global_attrs(base_file_id, ref_file_id):
        # check the dimensions are consistent
        if merge_test_dims(base_file_id, ref_file_id):
            # get a list of all the variables in the two files and check the attributes are the same
            if merge_test_var_attrs(base_file_id, ref_file_id):
                result = True

    # close the nc files
    base_file_id.close()
    ref_file_id.close()

    return result


def merge_test_global_attrs(base_nc_id, ref_nc_id):
    """
    Check to see if the files have identical variable names listed in the fieldnames global attributes

    Parameters
    ----------
    base_nc_id   : netcdf id
        base nc file id which is used as the reference in the comparison
    ref_nc_id    : netcdf id
        another nc file id for the comparison

    Returns
    -------
    result  : Boolean
        True    : files ARE compatible
        False   : files ARE NOT compatible
    """
    try:
        base_field_names = base_nc_id.getncattr(NC_ATT_FIELDNAMES)
        ref_field_names = ref_nc_id.getncattr(NC_ATT_FIELDNAMES)
        if base_field_names != ref_field_names:
            return False
    except AttributeError:
        return False
    return True


def merge_test_dims(base_nc_id, ref_nc_id):
    """
    Compare 2 netCDF files: check they have identical dimension names and lengths, excluding the merging (time) dim.

    Parameters
    ----------
    base_nc_id   : netcdf id
        base nc file id which is used as the reference in the comparison
    ref_nc_id    : netcdf id
        another nc file id for the comparison

    Returns
    -------
    result  : Boolean
        True    : dims ARE identical
        False   : dims ARE NOT identical
    """
    nc_dims_base = [dim for dim in base_nc_id.dimensions]
    nc_dims_ref = [dim for dim in ref_nc_id.dimensions]

    if nc_dims_base != nc_dims_ref:
        return False

    for dim in nc_dims_base:
        if dim == NC_MERGE_DIMENSION:
            continue
        if base_nc_id.dimensions[dim].size != ref_nc_id.dimensions[dim].size:
            return False

    return True


def merge_test_var_attrs(base_nc_id, ref_nc_id):
    """
    Check to see if 2 netCDF files have identical attributes for all data variables

    Parameters
    ----------
    base_nc_id   : netcdf id
        base nc file id which is used as the reference in the comparison
    ref_nc_id    : netcdf id
        another nc file id for the comparison

    Returns
    -------
    result  : Boolean
        True    : vars attributes for all vars ARE identical
        False   : vars attributes for all vars ARE NOT identical
    """
    nc_vars_base = [var for var in base_nc_id.variables]
    nc_vars_ref = [var for var in ref_nc_id.variables]

    for var_base in nc_vars_base:
        if var_base in nc_vars_ref:
            if compare_var_attributes(base_nc_id, ref_nc_id, var_base, []) is False:
                return False
        else:
            return False

    return True


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
        True    : attributes ARE identical
        False   : attributes ARE NOT identical or attributes are missing
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
        Variables:

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

        # calculate length of the merged time dimension
        new_dim_len = total_dimension_length(group, NC_MERGE_DIMENSION)

        # specify name of output file
        full_output_path = out_path + "\\" + NC_OUTPUT_PREFIX + base_name + ("_" + str(iteration) + NC_FILE_EXT)

        # create file, create dims, atts and vars. Initialise vars where appropriate (i.e. time invariant datasets)
        output_id = create_merged_nc_using_ref(full_output_path, group[0], new_dim_len)

        # copy var data, file by file, for the time dependant datasets
        start_pos, end_pos = 0, 0
        for idx, file in enumerate(group):
            input_id = nC.Dataset(file, "r")
            merge_dim_points = input_id.variables[NC_MERGE_DIMENSION].size
            end_pos += merge_dim_points

            for variable in input_id.variables:
                var_dims = input_id.variables[variable].dimensions
                if NC_MERGE_DIMENSION in var_dims:
                    copy_data(output_id, input_id, variable, start_pos, end_pos)

            if idx == (len(group) - 1):
                copy_data(output_id, input_id, "time_coverage_end")

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


def create_merged_nc_using_ref(output_path, ref_path, merge_dim_len):
    """
    Create a new netCDF to contain merged data, create and initialise dims, atts and vars.
    Can only initialise "static" vars, any vars with time varying data is copied elsewhere.

    Parameters
    ----------
    output_path   : string
        full path including filename for the merged output file
    ref_path   : string
        full path including filename for one of the input filenames which is used to model the output
    merge_dim_len    : int
        the new length of the merged (time) dimension, which is used to correctly define dimensions of merged vars

    Returns
    -------
    output_id : netcdf if
        a netcdf object reference to the output file
    """
    # create the output file
    output_id = nC.Dataset(output_path, "w", format="NETCDF4_CLASSIC")

    # create the destination dimensions, global atts, variables, and variable atts
    ref_id = nC.Dataset(ref_path, "r")
    copy_dimensions(ref_id, output_id, NC_MERGE_DIMENSION, merge_dim_len)
    copy_attributes(ref_id, output_id)
    copy_variables(ref_id, output_id, NC_MERGE_DIMENSION, merge_dim_len)
    ref_id.close()

    return output_id


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
    """
    Copies attributes from one input netCDF to an output netCDF

    Parameters
    ----------
    input_id   : netCDF id
        input which contains the original data. Can point to either root, for global atts, or a var for var atts
    output_id   : netCDF id
        destination to copy the new data into. Can point to either root, for global atts, or a var for var atts

    Returns
    -------
    None
    """
    input_attributes = input_id.ncattrs()
    for attr in input_attributes:
        if "_FillValue" in str(attr):
            continue
        output_id.setncattr(str(attr), input_id.getncattr(attr))


def copy_variables(f_in, f_out, override_dim, override_val):
    """
    Copies variables, (inc atts and also data where possible) from one input netCDF to an output netCDF

    Parameters
    ----------
    f_in   : netCDF id
        input which contains the original model data.
    f_out   : netCDF id
        destination to copy the new data into.
    override_dim    :   string
        specifies the name of the dimension which is merged/overridden
    override_val    :   int
        specifies the new length of the overriden/merged dimension

    Returns
    -------
    None
    """
    input_var_list = [var for var in f_in.variables]

    for var in input_var_list:

        var_dims = f_in.variables[var].dimensions
        # when we have large 3-d vars, make assertions about the chunking
        if len(var_dims) >= 3 and override_dim in var_dims:
            chunks = get_size_from_dims(f_in, var, NC_MERGE_DIMENSION, min(override_val, NC_MAX_CHUNK_TDIM))
        # for smaller arrays, use the default chunking
        else:
            chunks = get_size_from_dims(f_in, var, "", 0)

        f_out.createVariable(str(var), f_in.variables[var].dtype, var_dims,
                             compression='zlib', complevel=NC_COMP_VAL, chunksizes=chunks)
        copy_attributes(f_in.variables[var], f_out.variables[var])
        # copy the data for the vars which are time invariant
        if NC_MERGE_DIMENSION not in var_dims:
            copy_data(f_out, f_in, var)


def copy_data(f_out, f_in, var, start_pos=0, end_pos=0):
    """
    Copies data from one var in input file into another var in an output file

    Parameters
    ----------
    f_out   : netCDF id
        destination to copy the new data into.
    f_in   : netCDF id
        input which contains the original data.
    var    :   string
        name of the var of interest
    start_pos    :   int
        the starting index for specifying where merged arrays are to be stored in the output.
    end_pos    :   int
        the ending index for specifying where merged arrays are to be stored in the output.

    Returns
    -------
    None
    """
    input_shape = get_size_from_dims(f_in, var, "", 0)
    input_array = numpy.empty(input_shape, dtype=f_in.variables[var].dtype)
    input_array[:] = f_in.variables[var][:]
    input_dims = f_in.variables[var].dimensions
    destination = f_out.variables[var]

    if NC_MERGE_DIMENSION not in input_dims:
        destination[:] = input_array
    elif NC_MERGE_DIMENSION in input_dims:
        if len(input_dims) == 1:
            destination[start_pos:end_pos] = input_array
        elif len(input_dims) == 2:
            destination[start_pos:end_pos, :] = input_array
        elif len(input_dims) == 3:
            destination[start_pos:end_pos, :, :] = input_array


def get_size_from_dims(f_in, var, override_dim, override_val):
    """
    Calculates the "shape" of a var dataset
    Named dimension can be overridden

    Parameters
    ----------
    f_in   : netCDF id
        input netcdf which contains the var dataset.
    var : string
        name of a netcdf var dataset
    override_dim    :   string
        specifies the name of the dimension which is merged/overridden
    override_val    :   int
        specifies the new length of the overriden/merged dimension

    Returns
    -------
    chunk_list  :   list
        list specifying shape of var dataset, after applying any dim overrides
    """
    chunk_list = []
    dim_tuple = f_in.variables[var].dimensions
    if len(dim_tuple) == 0:
        chunk_list.append(f_in.variables[var].size)
    else:
        for item in dim_tuple:
            if item in override_dim:
                chunk_list.append(override_val)
            else:
                chunk_list.append(len(f_in.dimensions[item]))
    return chunk_list


def file_worker(folder):
    """
    Starts the merge process for a given data folder

    Parameters
    ----------
    folder   : string
        full path to a valid MRR data folder

    Returns
    -------
    True    : folder merge SUCCESSFUL
    False   : folder merge FAILED
    """
    output_folder, basename = os.path.split(folder)
    merge_list = create_merge_list(folder)
    merge_result = merge_nc_files(merge_list, output_folder, basename)
    return merge_result


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

    if isinstance(N_CPU, int):
        if N_CPU <= 0:
            n_proc = max(1, os.cpu_count()-2)
        else:
            n_proc = min(N_CPU, os.cpu_count()-2)
    else:
        sys.exit()

    process_list = []

    for i_day in range(start_day, stop_day + 1):
        date = datetime.date.today() - datetime.timedelta(days=i_day)
        folder_date_string = str(date.year) + str(date.month).zfill(2)
        input_folder = data_path + "\\" + folder_date_string + "\\" + folder_date_string + str(date.day).zfill(2)
        if is_valid_data_folder(input_folder):
            process_list.append(input_folder)

    with Pool(processes=n_proc) as pool:
        results = list(tqdm(pool.imap(file_worker, process_list), total=len(process_list)))
