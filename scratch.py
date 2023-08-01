import netCDF4 as nC
import numpy

OUTPUT_PATH = "C:\\Users\\jonny\\Desktop\\output.txt"
EXPORT_NAME = "Z"
AVERAGE_TIME = 60

# todo 1: add averging capability

def driver_function():
    file = get_latest_input()
    date_string = get_date_string()
    with nC.Dataset(file, "r", format="NETCDF4_CLASSIC") as input_id:
        with open(OUTPUT_PATH, 'w') as output_id:
            export_data(input_id, output_id, EXPORT_NAME, date_string)


def get_date_string():
    date_string = "221204"
    return date_string


def get_latest_input():
    filename = "C:\\Users\\jonny\\Desktop\\20221204_030000.nc"
    return filename


def export_data(f_input, f_output, variable, date_str):
    data_shape = get_dim_tuple(f_input, variable)

    for i in range(data_shape[0]):
        time_value = int(f_input.variables["time"][i]) % 86400
        f_output.write(date_str+","+str(time_value)+",")
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


def get_chunk_from_dims(data_id, dim_tuple, override_dim, override_val):
    chunk_list = []
    for item in dim_tuple:
        if item in override_dim:
            chunk_list.append(override_val)
        else:
            chunk_list.append(len(data_id.dimensions[item]))
    return chunk_list


def output_tests():
    output_id = nC.Dataset("C:\\Users\\jonny\\Desktop\\test.nc", "w", format="NETCDF4_CLASSIC")

    output_id.createDimension("time", 8640)
    output_id.createDimension("height", 128)
    output_id.createDimension("bins", 64)

    v1_dim = "time",
    v2_dim = "time", "height",
    v3_dim = "time", "height", "bins"

    v1_chunk = get_chunk_from_dims(output_id, v1_dim, "", 0)
    v2_chunk = get_chunk_from_dims(output_id, v2_dim, "", 0)
    v3_chunk = get_chunk_from_dims(output_id, v3_dim, "time", 560)

    var_1d = output_id.createVariable("var_1d", "f8", v1_dim, compression='zlib', complevel=2, chunksizes=v1_chunk)
    var_2d = output_id.createVariable("var_2d", "f8", v2_dim, compression='zlib', complevel=2, chunksizes=v2_chunk)
    var_3d = output_id.createVariable("var_3d", "f8", v3_dim, compression='zlib', complevel=2, chunksizes=v3_chunk)

    var_1d[:] = 0
    var_2d[:] = 0
    var_3d[:] = 0

    output_id.close()


def input_tests():
    print("test")
    input_id = nC.Dataset("C:\\Users\\jonny\\Desktop\\20230323_090000.nc", "r", format="NETCDF4_CLASSIC")

    nc_z = input_id.variables["Z"]
    z_dims = nc_z.dimensions
    data_chunk = get_chunk_from_dims(input_id, z_dims, "", 0)
    output_data = numpy.empty(data_chunk, dtype=numpy.float64, order='C')
    output_data[:] = input_id.variables["Z"][:]

    output_id = nC.Dataset("C:\\Users\\jonny\\Desktop\\test2.h5", "w", format="NETCDF4_CLASSIC")
    output_id.createDimension("time", data_chunk[0])
    output_id.createDimension("range", data_chunk[1])
    output_var = output_id.createVariable("Z", "f8", z_dims, compression='zlib', complevel=2)
    output_var[:] = output_data
    input_id.close()
    output_id.close()

    return "Hello"
    #return output_data
