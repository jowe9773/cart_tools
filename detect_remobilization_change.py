#detect_wood.py

#import neccesary packages and modules
from osgeo import gdal
from functions import FileFunctions, SickFunctions


##instantiate classes
ff = FileFunctions()
sf = SickFunctions()

class DetectRemobilizationChange:
    def __init__(self):
        print("Initialized DetectWood class")

    def detect_remobilization_change(self, before_fn, after_fn, parameters):
        #load in these geotiffs so that they can have their change detected
        before_dataset = gdal.Open(before_fn)
        after_dataset = gdal.Open(after_fn)

        #extract the topo information into a numpy raster
        before_topo = before_dataset.GetRasterBand(1).ReadAsArray()
        after_topo = after_dataset.GetRasterBand(1).ReadAsArray()

        #get geotransform information  
        after_geotransform = after_dataset.GetGeoTransform()
        after_width = after_dataset.RasterXSize
        after_height = after_dataset.RasterYSize

        #get xmin, xmax, ymin, and ymax from after geotransform (for export to geotiff at end)
        after = ["", after_geotransform[0],
                after_geotransform[0] + after_width * after_geotransform[1],
                after_geotransform[3] + after_height * after_geotransform[5],
                after_geotransform[3]]

        mask, sieve, connectedness = parameters

        change = sf.detect_change(before_topo, after_topo, mask, sieve, connectedness)

        return after, change



#now run the functions
if __name__ == "__main__":

    #place to adjust parameters for wood detection
    MASK_THRESHOLD = 5
    SIEVE_THRESHOLD = 2650
    CONNECTEDNESS = 4

    #place to adjust parameters for geotiff export
    EPSG = 32615

    params = (MASK_THRESHOLD, SIEVE_THRESHOLD, CONNECTEDNESS)

    #select before and after geotiff files
    before_fn = ff.load_fn("Select a wood topo geotiff file")
    after_fn = ff.load_fn("Select a remobilization topo geotiff file")

    #Select a directory to store output in
    out = ff.load_dn("Select a directory to store remobilization map in")

    #initialize the DetectWood class:
    drc = DetectRemobilizationChange()

    #send these datasets and parameters into the "detect_wood function"
    after, woodmap = drc.detect_remobilization_change(before_fn, after_fn, params)

    #export the wood map as a geotiff
    sf.export_topo_as_geotiff(after_fn, EPSG, out, woodmap, after, marker = "change")
