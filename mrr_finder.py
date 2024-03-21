import os
import netCDF4 as nC

NC_FILE_EXT = ".nc"
NC_ATT_TITLE_NAME = "title"
NC_ATT_TITLE_VALUE = "METEK MRR Pro"


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


def cpu_checker(n_cpu=0):
    if isinstance(n_cpu, int):
        if n_cpu == 0:
            return os.cpu_count()
        if n_cpu < 0:
            return max(os.cpu_count() + n_cpu, 1)
        if n_cpu > 0:
            return min(os.cpu_count(), n_cpu)
    return os.cpu_count()


def files_in_path(path):
    files = [path + file for file in os.listdir(path) if is_valid_data_file(path + file)]
    return files


def folders_in_path(path):
    folders = [path + folder for folder in os.listdir(path) if is_valid_data_folder(path + folder)]
    return folders
