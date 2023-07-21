import os.path
import netCDF4 as nC
from cftime import num2date
import matplotlib.pyplot as plt
import sys
import numpy as np

# import datetime as dt

run_mode = True

nc_file_ext = ".nc"
time_dimension = "time"
range_dimension = "range"
title_att_name = "title"
vars_att_name = "field_names"
title_att_value = "METEK MRR Pro"
output_extension = "_merged"
export_field_list = ["time", "range", "Z"]
export_folder = ".\\dest"
export_temp_name = "temp_"
folder_length_y = 4
folder_length_ym = 6
folder_length_ymd = 8
folder_name_lengths = [folder_length_y, folder_length_ym, folder_length_ymd]


def print_ncattr(nc_fid, key, f_out, opt):
    """
    Prints the NetCDF file attributes for a given key

    Parameters
    ----------
    nc_fid : netCDF4.Dataset
        A netCDF4 dateset object
    key : unicode
        a valid netCDF4.Dataset.variables key
    f_out : file object
        a destination file for outputting attribute information
    opt : variable
        a variable controlling if output should be sent to file (1) or not (0)
    """
    try:
        if opt == 1:
            f_out.write("\t\ttype:" + repr(nc_fid.variables[key].dtype) + "\n")
            for nc_att in nc_fid.variables[key].ncattrs():
                f_out.write('\t\t%s:' % nc_att + repr(nc_fid.variables[key].getncattr(nc_att)) + "\n")
        else:
            print("\t\ttype:", repr(nc_fid.variables[key].dtype))
            for nc_att in nc_fid.variables[key].ncattrs():
                print('\t\t%s:' % nc_att, repr(nc_fid.variables[key].getncattr(nc_att)))

    except KeyError:
        if opt == 1:
            f_out.write("\t\tWARNING: %s does not contain variable attributes" % key + "\n")
        else:
            print("\t\tWARNING: %s does not contain variable attributes" % key)


def ncdump(nc_fid, f_out, opt, verb=True):
    """
    ncdump outputs dimensions, variables and their attribute information.
    The information is similar to that of NCAR's ncdump utility.
    ncdump requires a valid instance of Dataset.

    Parameters
    ----------
    nc_fid : netCDF4.Dataset
        A netCDF4 dateset object
    f_out : python file object
        A destination file created with e.g. open('myfile.txt',w)
    opt :   variable
        A variable with value 1 (output data to file) or 0 (no output to file)
    verb : Boolean
        whether nc_attrs, nc_dims, and nc_vars are printed

    Returns
    -------
    nc_attrs : list
        A Python list of the NetCDF file global attributes
    nc_dims : list
        A Python list of the NetCDF file dimensions
    nc_vars : list
        A Python list of the NetCDF file variables
    """

    # NetCDF global attributes
    nc_attrs = nc_fid.ncattrs()
    if verb:
        if opt == 1:
            f_out.write("NetCDF Global Attributes:" + "\n")
            for nc_attr in nc_attrs:
                f_out.write('\t%s:' % nc_attr + repr(nc_fid.getncattr(nc_attr)) + "\n")
        else:
            print("NetCDF Global Attributes:")
            for nc_attr in nc_attrs:
                print('\t%s:' % nc_attr, repr(nc_fid.getncattr(nc_attr)))
    nc_dims = [dim for dim in nc_fid.dimensions]  # list of nc dimensions
    # Dimension shape information.
    if verb:
        if opt == 1:
            f_out.write("NetCDF dimension information:")
            for dim in nc_dims:
                f_out.write("\tName:" + dim + "\n")
                f_out.write("\t\tsize:" + str(len(nc_fid.dimensions[dim])) + "\n")
                print_ncattr(nc_fid, dim, f_out, opt)
        else:
            print("NetCDF dimension information:")
            for dim in nc_dims:
                print("\tName:", dim)
                print("\t\tsize:", len(nc_fid.dimensions[dim]))
                print_ncattr(nc_fid, dim, f_out, opt)
    # Variable information.
    nc_vars = [var for var in nc_fid.variables]  # list of nc variables
    if verb:
        if opt == 1:
            f_out.write("NetCDF variable information:")
            for var in nc_vars:
                if var not in nc_dims:
                    f_out.write('\tName:' + var + "\n")
                    f_out.write("\t\tdimensions:" + str(nc_fid.variables[var].dimensions) + "\n")
                    f_out.write("\t\tsize:" + str(nc_fid.variables[var].size) + "\n")
                    print_ncattr(nc_fid, var, f_out, opt)
        else:
            print("NetCDF variable information:")
            for var in nc_vars:
                if var not in nc_dims:
                    print('\tName:', var)
                    print("\t\tdimensions:", nc_fid.variables[var].dimensions)
                    print("\t\tsize:", nc_fid.variables[var].size)
                    print_ncattr(nc_fid, var, f_out, opt)
    return nc_attrs, nc_dims, nc_vars


def test_file_manipulation(path):
    nc_file_list = []

    file_list = os.listdir(path)
    for file in file_list:
        if file.endswith(nc_file_ext):
            full_path = os.path.join(path, file)
            if valid_data_file(full_path):
                nc_file_list.append(full_path)

    for file in nc_file_list:
        time_data = netcdf_load_var(file, time_dimension)
        range_data = netcdf_load_var(file, range_dimension)
        time_start = netcdf_load_var(file, "time_coverage_start")
        time_end = netcdf_load_var(file, "time_coverage_end")
        time_ref = netcdf_load_var(file, "time_reference")
        units = 'seconds since 1970-01-01 00:00:00'
        time_start_string = nC.chartostring(time_start)
        time_end_string = nC.chartostring(time_end)
        time_ref_string = nC.chartostring(time_ref)
        print(time_start_string, time_end_string, time_ref_string)
        print(time_data.shape, min(time_data), max(time_data), len(time_data))
        print(range_data.shape, min(range_data), max(range_data), len(range_data))
        print(num2date(time_data[0:3], units))
        print(num2date(time_data[-3:], units))
        netcdf_print_var_attributes(file, "time")
        netcdf_print_dims(file)
        break


def netcdf_lib_version():
    print(nC.getlibversion())


def netcdf_print_groups(file):
    root_group = nC.Dataset(file, "r")
    print("grp list:")
    for group in root_group.groups:
        print(group)
    root_group.close()


def netcdf_print_dims(file):
    # print all the dimensions in a given file / group
    root_group = nC.Dataset(file, "r")
    print("dim list:")
    for dim in root_group.dimensions:
        print(root_group.dimensions[dim].name, len(root_group.dimensions[dim]))
    root_group.close()


def netcdf_print_var_attributes(file, variable):
    root_group = nC.Dataset(file, "r")
    var = root_group.variables[variable]
    for att in var.ncattrs():
        try:
            att_val = var.getncattr(att)
            print(att, att_val)
        except AttributeError:
            print("Attribute not found", att)
            continue
    root_group.close()


def netcdf_print_attributes(file):
    # print all the dimensions in a given file / group
    root_group = nC.Dataset(file, "r")
    print("atr list:")
    for attribute in root_group.ncattrs():
        print(attribute, getattr(root_group, attribute))
    root_group.close()


def netcdf_print_var_names(file):
    root_group = nC.Dataset(file, "r")
    print("var list:")
    for variable in root_group.variables:
        print(variable)
    root_group.close()


def netcdf_plot_2d(file, variable):
    root_group = nC.Dataset(file, "r")
    var = root_group.variables[variable]
    var_data = var[:]
    root_group.close()
    fig, ax = plt.subplots()
    ax.imshow(var_data)
    plt.show()


def netcdf_load_var(file, variable):
    root_group = nC.Dataset(file, "r")
    var = root_group.variables[variable]
    var_data = var[:]
    root_group.close()
    return var_data


def valid_data_file(file):
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
    if file.endswith(nc_file_ext) is False:
        return False
    nc_file = nC.Dataset(file, "r")
    attribute_title = getattr(nc_file, title_att_name)
    nc_file.close()
    if title_att_value in attribute_title:
        return True
    return False


def valid_data_folder(path):
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
    # check the data path is a real input
    if os.path.exists(path) is False:
        return False
    else:
        file_list = os.listdir(path)
        for file in file_list:
            if valid_data_file(os.path.join(path, file)):
                return True
        return False


def search_directories(path):
    """
    Search a given path to find MRR-Pro datafiles
    the search method is based on the known structure of the MRR-Pro data archive:
    */yyyy/yyyymm/yyyymmdd/*.nc
    the function examines the last part of the input path and determines how deep to search for data

    Parameters
    ----------
    path : string
        string containing a path to a MRR-Pro raw data archive

    Returns
    -------
    folders : list
        A list of the full paths to all NetCDFs found in the root input path
    """

    folders = []

    if os.path.exists(path) is False:
        return folders

    base_folder_length = len(os.path.basename(path))
    if base_folder_length not in folder_name_lengths:
        return folders

    if base_folder_length == folder_length_ymd:
        if valid_data_folder(path) is True:
            folders.append(path)
            return folders

    if base_folder_length == folder_length_ym:
        for sub_folder in os.listdir(path):
            full_path = os.path.join(path, sub_folder)
            if valid_data_folder(full_path) is True:
                folders.append(full_path)
        return folders

    if base_folder_length == folder_length_y:
        for sub_folder in os.listdir(path):
            sub_path = os.path.join(path, sub_folder)
            for sub_sub_folder in os.listdir(sub_path):
                full_path = os.path.join(sub_path, sub_sub_folder)
                if valid_data_folder(full_path) is True:
                    folders.append(full_path)
        return folders

    return folders


def compare_attributes(data1_id, data2_id, var_name, att_list):
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

    Returns
    -------
    result  : Boolean
        True    : attribute ARE identical
        False   : files ARE NOT identical
    """
    result = True
    for att in data1_id.variables[var_name].ncattrs():
        if att in att_list:
            att1 = data1_id.variables[var_name].getncattr(att)
            att2 = data2_id.variables[var_name].getncattr(att)
            if att1 != att2:
                result = False
                break
    return result


def merge_test(base_file, ref_file):
    """
    Check to see if the files are compatible for merging
    The following are all checked and must pass check to be suitable for a merge:
        1) identical variable names, including identical dtypes and identical units
        2) identical range dimension
        3) identical spectrum_n_samples dimension
        5) any more checks?.....

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
    base_field_names = getattr(base_file_id, vars_att_name)
    ref_field_names = getattr(ref_file_id, vars_att_name)
    if base_field_names != ref_field_names:
        result = False

    # check the dimensions are consistent
    if result:
        range_comp_list = ["size", "units", "type", "meters_to_center_of_first_gate", "meters_between_gates"]
        result = compare_attributes(base_file_id, ref_file_id, "range", range_comp_list)
        if result:
            n_samples_base = base_file_id.dimensions["spectrum_n_samples"]
            n_samples_ref = ref_file_id.dimensions["spectrum_n_samples"]
            if n_samples_base.size != n_samples_ref.size:
                result = False

    # get a list of all the variables in the two files and check the units and dtype are the same
    if result:
        nc_vars_base = [var for var in base_file_id.variables]
        nc_vars_ref = [var for var in ref_file_id.variables]
        for var_base in nc_vars_base:
            if var_base in nc_vars_ref:

                result = compare_attributes(base_file_id, ref_file_id, var_base, ["units", "type"])
                if not result:
                    break
            else:
                result = False
                break

    # close the nc files
    base_file_id.close()
    ref_file_id.close()

    return result


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

        full_output_path = out_path + "\\" + base_name + ("_" + str(iteration) + ".nc")
        output_id = nC.Dataset(full_output_path, "w", format="NETCDF4_CLASSIC")

        # declare all dimensions and calculate length of the merged time dimension
        time_len = 0
        for idx, file in enumerate(group):
            input_id = nC.Dataset(file, "r")
            time_len += len(input_id.dimensions["time"])
            if idx == 0:
                nc_vars = [var for var in input_id.variables]
                nc_dims = [dim for dim in input_id.dimensions]  # list of nc dimensions
                for dim in nc_dims:
                    if str(dim) == "time":
                        continue
                    new_dim = output_id.createDimension(str(dim), len(dim))
            input_id.close()
        pass
        time_destination = np.empty(time_len, dtype=np.float64)
        t_dim = output_id.createDimension("time", time_len)

        # FINISH ME
        for var in nc_vars:
            start_pos, end_pos = 0, 0
            for idx, file in enumerate(group):
                input_id = nC.Dataset(file, "r")
                if idx == 0:
                    try:
                        var_dims = input_id.variables[var].dimensions
                        # var_type = input_id.variables[var].getncattr("type")
                        print(str(var), var_dims)
                    except:
                        input_id.close()
                        break
                    # if "time is not in the dimensions
                    # create a copy of the var and put it in the nc
                    # copy the attributes
                    # break out of the loop
                    # else:
                    # loop the dimensions to get the dimensions of the arrays
                    # override the size of the time dimension
                    # create a numpy array with the correct dimensions
                # create the nc variable in the destination
                # copy the attributes
                # reference the source variable
                # copy data from source array to the destination array

                input_id.close()

        # this copies and exports the time dataset
        # need to make this more general for all the datasets
        # need to also grab all attributes (globals as well as for the variables themselves)
        # start_pos, end_pos = 0, 0

        # for idx, file in enumerate(group):
        #     input_id = nC.Dataset(file, "r")
        #     time_points = input_id.variables["time"].size
        #     end_pos += time_points
        #
        #     time_destination[start_pos:end_pos] = input_id.variables["time"][:]
        #
        #     input_id.close()
        #     start_pos += time_points
        #
        # t1 = output_id.createVariable("time", "f8", ("time",))
        # t1.units = "days since 2000-01-01"
        # t1.calendar = "standard"
        # t1[:] = time_destination

        output_id.close()

    return False


def export_fields(variable_list, group_list, output_fold, output_file):
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
    group_list : python list
        names of folders to search for source datafiles
    output_fold   : string
        path to an output directory to store the extracted subset of data
    output_file   : string
        prefix to add to filename to output files

    Returns
    -------
    num_files   : variable
        number of files which have been extracted
    """
    num_files = 0

    for iteration, group in enumerate(group_list):

        for idx, file in enumerate(group):
            input_id = nC.Dataset(file, "r")

            full_output_path = output_fold + "\\" + output_file + os.path.basename(file)
            output_id = nC.Dataset(full_output_path, "w", format="NETCDF4_CLASSIC")

            t_dim_input = input_id.dimensions["time"]
            r_dim_input = input_id.dimensions["range"]

            t_dim_output = output_id.createDimension("time", len(t_dim_input))
            r_dim_output = output_id.createDimension("range", len(r_dim_input))

            for var in variable_list:
                var_input = input_id.variables[var]
                var_data = var_input[:]
                dimension = input_id.variables[var].dimensions
                data_type = input_id.variables[var].dtype
                var_output = output_id.createVariable(str(var), data_type, dimension, compression='zlib', complevel=9)
                var_output[:] = var_data

            input_id.close()
            output_id.close()
            num_files += 1

    return num_files


def check_merge(path):
    """
    Look at all the valid nc files in a given path....
    Check which files can be merged based on dataset info (variable names, dtypes, units and dimensions)
    Merge the files where possible

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
        if valid_data_file(full_path):
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


# Press the green button in the gutter to run the script.
if __name__ == '__main__':

    if run_mode:
        data_path = ".\\data\\2023\\202303\\20230322"
    else:
        data_path = (sys.argv[0])

    folder_list = search_directories(data_path)

    for folder in folder_list:
        process_list = check_merge(folder)
        # destination_fold, sub_fold = os.path.split(folder)
        # basename = "MRRPro" + "_" + sub_fold
        # merge_result = merge_nc_files(process_list, destination_fold, basename)
        export_fields(export_field_list, process_list, export_folder, export_temp_name)
