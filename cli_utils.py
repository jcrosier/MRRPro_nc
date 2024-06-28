from os import path

def parameter_check(params: list[str], number_expected: int, paths: list[int] = None, ints: list[int] = None) -> bool:
    """ Checks validity of command line input parameters. """
    if len(params) < number_expected:
        print('Error: Not enough Command Line parameters')
        return False
    path_params = [par for idx, par in enumerate(params) if idx in paths]
    if list_contains_paths(path_params) is False:
        print('Error: Expected path inputs do not point to valid destinations')
        return False
    int_params = [par for idx, par in enumerate(params) if idx in ints]
    if list_contains_ints(int_params) is False:
        print('Error: Expected int inputs are not ints')
        return False
    
    return True


def list_contains_strings(list: list) -> bool:
    """ Checks if all items in a list are of type str. """
    return all([isinstance(item, str) for item in list])


def list_contains_ints(list: list) -> bool:
    """ Check if all items in a list can be converted to type int. """
    try:
        [int(item) for item in list]
    except ValueError:
        return False
    return True


def list_contains_paths(list: list) -> bool:
    """ Checks if all items in a list are valid paths. """
    if list_contains_strings(list) is False:
        return False
    return all([path.exists(path) for path in list])
