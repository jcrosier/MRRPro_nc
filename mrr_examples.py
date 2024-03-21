import mrr_plots
import mrr_finder
import mrr_nc_merger

# Example: Get List of all MRR files it a given folder
file_list = mrr_finder.files_in_path('./data/2023/202303/20230322/')

# Example: Plot image of reflectivity vs time, range
mrr_plots.plot_2d_data(file_list, 'Ze', [True, -10, 30], offset_y=True)

# Example: Merge all MRR nc files in a given folder
result = mrr_nc_merger.mrr_merge_folder('./data/2023/202303/20230322/')

# Example: Get a list of all sub-folders containing MRR files (e.g. all days in March 2023)
#folder_list = mrr_finder.folders_in_path('./data/2023/202303/')

# Example: Multithread merge all MRR nc files for each sub-folder in a list of folders
#list_result = mrr_nc_merger.mrr_merge_list(folder_list, n_cpu=6)
