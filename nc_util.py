import numpy
import netCDF4 as nC

def load_var(nc_id: nC.Dataset, var: str) -> numpy.ndarray:
    """ Load an entire var into a numpy array. """
    data_shape = list(nc_id.variables[var].shape)
    data_array = numpy.empty(data_shape, dtype=nc_id.variables[var].dtype)
    data_array[:] = nc_id.variables[var][:]
    return data_array


def get_numpy_data_slice(data: numpy.ndarray, start_idx: int, average: int) -> numpy.ndarray:
    """ Get a 1D numpy array from a 2D numpy array, either by taking a single slice, or by averaging. """
    if average == 1:
        slice = data[start_idx, :]
    else:
        slice = numpy.nanmean(data[start_idx:start_idx + average, :], axis=0)
    return slice


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


def nc_vars_present(nc_id: nC.Dataset, requested_vars: list[str]) -> bool:
    """ Checks an nC file to see if it contains the specificed requested_vars. """
    nc_vars = [nc_var for nc_var in nc_id.variables]
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
