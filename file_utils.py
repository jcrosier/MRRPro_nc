import netCDF4 as nC
from datetime import date, timedelta
import os


NC_FILE_EXT = ".nc"
NC_ATT_TITLE_NAME = "title"
NC_ATT_TITLE_VALUE = "METEK MRR Pro"


def is_valid_data_folder(path: str) -> bool:
    """ Check if the input path string contains valid MRR-Pro netcdf files. """
    if os.path.isdir(path) is False:
        return False

    file_list = [file for file in os.listdir(path) if file.endswith(NC_FILE_EXT)]

    for file in file_list:
        if is_valid_data_file(os.path.join(path, file)):
            return True

    return False


def is_valid_data_file(file: str) -> bool:
    """ Check the input file string is a valid MRR-Pro netcdf file. """
    if os.path.exists(file) is False:
        return False

    if file.endswith(NC_FILE_EXT) is False:
        return False

    if is_valid_mrr_nc(file) is False:
        return False

    return True


def is_valid_mrr_nc(file: str) -> bool:
    """ Check if a nc data file can be opened with nc lib and contains valid attribute/metadata. """
    try:
        nc_file = nC.Dataset(file, "r")
    except OSError:
        return False

    try:
        title_data = getattr(nc_file, NC_ATT_TITLE_NAME)
        nc_file.close()
    except AttributeError:
        nc_file.close()
        return False

    return NC_ATT_TITLE_VALUE in title_data


def files_in_path(path: str) -> list[str]:
    """ Create a list of all valid MRR-Pro nc files in a given path. """
    if os.path.exists(path) is False: return []
    return [path + file for file in os.listdir(path) if is_valid_data_file(path + file)]


def folders_in_path(path: str) -> list[str]:
    """ Create a list of all 'MRR-Pro data containing' subdirectories in a given path. """
    if os.path.exists(path) is False: return []
    return [path + folder for folder in os.listdir(path) if is_valid_data_folder(path + folder)]


def get_latest_netcdf(input_path: str, max_days_prior: int = 1) -> str:
    """ Find the newest data-containing source data folder, starting from the current data, and optionally scanning back in time.
    Default is to scan back 1 day to ensure we can get latest data if it is pass midnight, but no new file has been created yet. """
    for day in range(max_days_prior + 1):
        newest_file = get_newest_file_for_date(date.today() - timedelta(days=day), input_path)
        if len(newest_file) > 0: return newest_file

    print("Error: No source netCDF files found")
    return ""


def get_newest_file_for_date(date: date, root: str) -> str:
    """ Find the newest file in a source folder for a given date. """
    path_str = mrrpro_folder_from_date(date.today(), root)
    file_list = files_in_path(path_str)
    if len(file_list) > 0:
        return max(file_list, key=os.path.getctime)
    else:
        return ""
    

def mrrpro_folder_from_date(date: date, root: str) -> str:
    """ Construct expected full path to MRR-Pro source folder for a given date. """
    folder = str(date.year) + str(date.month).zfill(2)
    subfolder = folder + str(date.day).zfill(2)
    return root + "\\" + folder + "\\" + subfolder + "\\"


def cpu_checker(cpu_request: int = 0) -> int:
    """ Obtain number of CPU's available for use.
    cpu_request == 0: obtain max possible cpus from OS
    cpu_request < 0: obtain the max number from the OS, but keep some cpu in reserve
    cpu_request > 0: obtain a specific number, checking upper limits from the OS """
    if isinstance(cpu_request, int):
        if cpu_request == 0:
            return os.cpu_count()
        if cpu_request < 0:
            return max(os.cpu_count() + cpu_request, 1)
        if cpu_request > 0:
            return min(os.cpu_count(), cpu_request)
    return os.cpu_count()
