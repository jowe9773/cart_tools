#swath_dat_to_geotiff.py

#import neccesary packaged and modules
from functions import FileFunctions, SickFunctions

#instantiate FileFunctions and SickFunctions
ff = FileFunctions()
sf = SickFunctions()


def swath_dat_to_geotiff(filename, proj_num, out_dir):

    #load the file into a numpy array
    sick_data = sf.load_swath_sick_file(filename)

    #interpolate and fill nulls
    topo_filled = sf.fill_nulls(sick_data[0])

    #export the filled topo as a geotiff
    sf.export_topo_as_geotiff(filename, proj_num, out_dir, topo_filled, sick_data)


if __name__ == "__main__":

    #parameters to change (if applicable)
    PROJ_NUM = 32615

    #select a file to process
    fn = ff.load_fn("Select a .dat file to process")

    out = ff.load_dn("Select a output directory for geotiff")

    #select a directory for storage location of output

    #process the file and export to geotiff
    swath_dat_to_geotiff(fn, PROJ_NUM, out)
