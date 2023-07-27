import os
import netCDF4 as nC
import sys
import datetime
import numpy as np

nc_file_ext = ".nc"
title_att_name = "title"
title_att_value = "METEK MRR Pro"
vars_att_name = "field_names"
output_prefix = "MRRPro_"


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
        time_data = netcdf_load_var(file, "time")
        range_data = netcdf_load_var(file, "range")
        time_start = netcdf_load_var(file, "time_coverage_start")
        time_end = netcdf_load_var(file, "time_coverage_end")
        time_ref = netcdf_load_var(file, "time_reference")
        time_start_string = nC.chartostring(time_start)
        time_end_string = nC.chartostring(time_end)
        time_ref_string = nC.chartostring(time_ref)
        print(time_start_string, time_end_string, time_ref_string)
        print(time_data.shape, min(time_data), max(time_data), len(time_data))
        print(range_data.shape, min(range_data), max(range_data), len(range_data))
        # units = 'seconds since 1970-01-01 00:00:00'
        # print(cftime.num2date(time_data[0:3], units))
        # print(cftime.num2date(time_data[-3:], units))
        break


def netcdf_load_var(file, variable):
    root_group = nC.Dataset(file, "r")
    var = root_group.variables[variable]
    var_data = var[:]
    root_group.close()
    return var_data


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

        full_output_path = out_path + "\\" + output_prefix + base_name + ("_" + str(iteration) + nc_file_ext)
        output_id = nC.Dataset(full_output_path, "w", format="NETCDF4_CLASSIC")

        # calculate length of the merged time dimension
        time_len = 0
        for idx, file in enumerate(group):
            input_id = nC.Dataset(file, "r")
            time_len += len(input_id.dimensions["time"])
            input_id.close()

        # open the first file to get list of vars and dims
        input_id = nC.Dataset(group[0], "r")
        f1_vars = [var for var in input_id.variables]
        f1_dims = [dim for dim in input_id.dimensions]
        for dim in f1_dims:
            if str(dim) == "time":
                output_id.createDimension(str(dim), time_len)
                continue
            output_id.createDimension(str(dim), len(input_id.dimensions[dim]))
        input_id.close()

        destVal = 0

        for var in f1_vars:

            destVal += 1

            input_id = nC.Dataset(group[0], "r")

            var_dims = input_id.variables[var].dimensions
            var_type = input_id.variables[var].dtype
            var_shape_str = get_shape_for_variable(input_id, var, time_len)

            att_names = []
            att_vals = []
            for nc_att in input_id.variables[var].ncattrs():
                att_names.append(str(nc_att))
                att_vals.append(input_id.variables[var].getncattr(nc_att))

            input_id.close()

            nc_att_dict = dict(zip(att_names, att_vals))

            if "_FillValue" in att_names:
                output_data = output_id.createVariable(str(var), var_type, var_dims,
                                                       fill_value=nc_att_dict["_FillValue"],
                                                       compression='zlib', complevel=9)
            else:
                output_data = output_id.createVariable(str(var), var_type, var_dims, compression='zlib', complevel=9)

            print(nc_att_dict)
            for item in nc_att_dict:
                print(item, nc_att_dict[item])
                output_data.setncattr(item, nc_att_dict[item])
            #output_data.setncatts(nc_att_dict)

            try:
                output_data[:] = destVal
            except:
                output_data[:] = str(destVal)

            # t1.units = "days since 2000-01-01"
            # t1.calendar = "standard"
            # t1[:] = time_destination
            print(str(var), var_shape_str, str(var_type), str(var_dims))

            continue

            # start_pos, end_pos = 0, 0
            for idx, file in enumerate(group):
                input_id = nC.Dataset(file, "r")
                var_shape_str = get_shape_for_variable(input_id, var, time_len)
                var_type = input_id.variables[var].dtype
                print(str(var), var_shape_str, str(var_type))
                input_id.close()
                break

                dims_for_variable = input_id.variables[var].dimensions

                if "time" not in dims_for_variable:

                    input_id.close()
                    break

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


def get_shape_for_variable(file_id, var, time_new_len):
    dim_list = file_id.variables[var].dimensions
    shape_string = "("
    num_dims = len(dim_list)

    for idx, dim in enumerate(dim_list):
        if str(dim) == "time":
            len_str = str(time_new_len)
        else:
            len_str = str(len(file_id.dimensions[dim]))
        shape_string += len_str
        if idx == 0 or idx < (num_dims - 1):
            shape_string += ","
    shape_string += ")"

    return shape_string


    # Press the green button in the gutter to run the script.
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

    for i_day in range(start_day, stop_day+1):
        date = datetime.date.today() - datetime.timedelta(days=i_day)
        folder_date_string = str(date.year) + str(date.month).zfill(2)
        sub_folder_date_string = folder_date_string + str(date.day).zfill(2)

        input_folder = data_path + "\\" + folder_date_string + "\\" + sub_folder_date_string
        output_folder = data_path + "\\" + folder_date_string
        if valid_data_folder(input_folder) is False:
            continue

        process_list = create_merge_list(input_folder)
        merge_result = merge_nc_files(process_list, output_folder, sub_folder_date_string)
