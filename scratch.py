import netCDF4 as nC
import numpy

# todo need to understand Masked arrays and Fill Values


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
