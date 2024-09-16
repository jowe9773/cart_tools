#functions.py

#import packages and modules
import xml.etree.ElementTree as ET
import numpy as np
from osgeo import gdal,osr
import pandas as pd
import geopandas as gpd
import shapely.geometry
import rasterio



class SickFunctions:
    """class containing methods for sick processing"""
        
    def __init__(self):
        print("initialized")

    def load_processed_sick_file(self, filename):
        """method for taking a .DAT file and turning it into a usable numpy array"""

        file = filename.split('/')[-1]
        # Open the binary file
        with open(filename, 'rb') as f:
            # Read the binary data into a NumPy array
            data = np.fromfile(f, dtype=np.float32)  # Adjust dtype according to your data type

        data[data==-9999] = np.nan

        # Get Grid start and end from filename:
        grid_inds = [file.find('Grid='),file.find('Grid=')+len('Grid=')]
        xy_inds   = [file.find('XY='),file.find('XY=')+len('XY=')]

        grid_string = file[grid_inds[-1]:xy_inds[0]]
        grid_string = grid_string.strip().replace('(','').replace(')','')
        #grid = np.asarray(grid_string.split('x')).astype('float')

        xy_string = file[xy_inds[-1]:]
        xy_string = xy_string.strip().replace('(','').replace(')','').replace('_TopoData.DAT','')
        print(xy_string)
        xy_start_string  = xy_string.split(' to ')[0]
        xy_end_string    = xy_string.split(' to ')[1]
        xy_start         = np.asarray(xy_start_string.split(',')).astype('float')
        xy_end           = np.asarray(xy_end_string.split(',')).astype('float')

        xmin = xy_start[0]
        xmax = xy_end[0]
        ymin = xy_start[1]
        ymax = xy_end[1]

        # Read the Corresponding XML File:
        xml_filename = filename[0:-4]+'.xml'
        xml_data = ET.parse(xml_filename)
        root = xml_data.getroot()
        width = int(root.find('.//parameter[@name="width"]').text)

        topo = data.reshape([int(len(data)/width),width])

        topo = topo[::-1, :]

        print(f'Number of Datapoints: {topo.flatten().shape}')

        return topo, xmin, xmax, ymin, ymax

    def fill_nulls(self, topo):
        """method for filling areas of nodata within a numpy array"""

        #Define the dimensions and data type
        rows, cols = topo.shape
        data_type = gdal.GDT_Float32

        # Create an in-memory raster dataset
        driver = gdal.GetDriverByName('MEM')
        dataset = driver.Create('', cols, rows, 1, data_type)

        # Write the array to the dataset
        band = dataset.GetRasterBand(1)
        no_data_value = -9999
        band.SetNoDataValue(no_data_value)
        topo[np.isnan(topo)] = no_data_value
        band.WriteArray(topo)

        # FillNoData parameters
        mask_band = None
        max_distance = 100
        smoothing_iterations = 0

        # Interpolate NoData values
        gdal.FillNodata(targetBand=band, maskBand=mask_band, maxSearchDist=max_distance, smoothingIterations=smoothing_iterations)

        # Read the interpolated data back into a NumPy array
        interpolated_topo = band.ReadAsArray()

        # Clean up
        band = None
        dataset = None

        print("Interpolation complete.")

        return interpolated_topo


    def extract_wood(self, topo1, topo2, mask_threshold, sieve_threshold, connectedness):

        """method for finding areas of change (presence or movement of wood)"""
        #difference the before and after
        difference = topo2 - topo1

        # Create an in-memory raster dataset
        rows, cols = difference.shape
        data_type = gdal.GDT_Float32

        driver = gdal.GetDriverByName('MEM')
        dataset = driver.Create('', cols, rows, 1, data_type)

        # Write the array to the dataset
        band = dataset.GetRasterBand(1)
        no_data_value = -9999
        band.SetNoDataValue(no_data_value)
        difference[np.isnan(difference)] = no_data_value
        band.WriteArray(difference)

        mask_pos = difference > mask_threshold
        mask_neg = difference < -mask_threshold

        mask = mask_pos | mask_neg

        # Create an in-memory raster dataset
        driver = gdal.GetDriverByName('MEM')
        masked_dataset = driver.Create('', cols, rows, 1, data_type)
        masked_band = masked_dataset.GetRasterBand(1)
        masked_band.WriteArray(mask.astype(np.uint8))

        # Apply the sieve filter to remove small regions
        gdal.SieveFilter(
            srcBand = masked_band,
            maskBand = None,
            dstBand = masked_band,
            threshold = sieve_threshold,
            connectedness = connectedness)

        # Read the sieved data back into a NumPy array
        sieved_mask = masked_band.ReadAsArray().astype(bool)

        # Apply the sieve mask to the difference array
        wood = np.where(sieved_mask, difference, 0)

        return wood


    def export_topo_as_geotiff(self, filename, projection_num, out_directory, topo, sick, wood = False, remobilization = False):
        """method for exporting a particular scan or map as a geotiff file"""

        xmin = sick[1]
        xmax = sick[2]
        ymin = sick[3]
        ymax = sick[4]

        proj = osr.SpatialReference()
        proj.ImportFromEPSG(projection_num)

        datout = np.squeeze(topo)

        datout[np.isnan(datout)] = -9999
        driver = gdal.GetDriverByName('GTiff')
        cols,rows = np.shape(datout)

        out_names = filename.split("/")[-1].split("_")[:3]

        if wood is True:
            output_filename = out_directory + "/" + out_names[0] + "_" + out_names[1] + "_woodmap.tif"

        elif remobilization is True:
            output_filename = out_directory + "/" + out_names[0] + "_" + out_names[1] + "_remobilizationmap.tif"

        else:
            output_filename = out_directory + "/" + out_names[0] + "_" + out_names[1] + "_" + out_names[2] + '.tif'

        ds = driver.Create(output_filename, rows, cols, 1, gdal.GDT_Float32, [ 'COMPRESS=LZW' ] )
        if proj is not None:
            ds.SetProjection(proj.ExportToWkt())

        xres = (xmax - xmin) / float(rows)
        yres = (ymax - ymin) / float(cols)

        geotransform = (xmin, xres, 0, ymax, 0, -yres)

        ds.SetGeoTransform(geotransform)
        ss_band = ds.GetRasterBand(1)
        ss_band.WriteArray(datout)
        ss_band.SetNoDataValue(-9999)
        ss_band.FlushCache()
        ss_band.ComputeStatistics(False)
        ss_band.SetUnitType('m')

class MassaFunctions:
    """This class contains the tools that will come together to process the massa data"""

    def __init__(self):
        print("initialized")

    def load_massa_file(self, filepath, epsg):
        """take massa file and turn it into a geopandas dataframe"""

        # Load the CSV file into a pandas DataFrame
        # Skip the first 28 rows and the line immediately after the header
        df = pd.read_csv(filepath, skiprows=29, header=0)

        # Skip the line immediately after the header by reloading the file and skipping necessary lines
        df = pd.read_csv(filepath, skiprows=list(range(28)) + [29])

        df["Y"] = df["Y"] * (-1)

        # Assuming your CSV has 'latitude' and 'longitude' columns
        # Create a 'geometry' column using the 'longitude' and 'latitude' columns
        geometry = [shapely.geometry.Point(xy) for xy in zip(df['X'], df['Y'])]

        # Convert the DataFrame to a GeoDataFrame
        gdf = gpd.GeoDataFrame(df, geometry=geometry)

        # Set the coordinate reference system (CRS) if known, e.g., WGS84
        gdf.set_crs(epsg = epsg, inplace=True)

        # Now you have a GeoDataFrame ready for spatial operations
        return gdf

    def get_water_depth(self, water_elev, fp_elev, offset = 0):
        """take the loaded massa file (geodataframe) and add a column for floodplain elevation by extracting the nearest elevation point from a geotiff"""

        with rasterio.open(fp_elev) as src:
            # Ensure the CRS of the GeoTIFF and the GeoDataFrame match
            assert src.crs.to_string() == water_elev.crs.to_string(), "CRS mismatch between GeoTIFF and GeoDataFrame"

            # Sample the raster values at the point locations
            coords = [(x,y) for x, y in zip(water_elev.geometry.x, water_elev.geometry.y)]
            floodplain_elevation = [val[0] for val in src.sample(coords)]


            # Append the raster values to the GeoDataFrame
            water_elev["fp_elev"] = floodplain_elevation

            water_elev["Massa Target"] = water_elev["Massa Target"].replace("OutOfRange", -9999)
            water_elev["Massa Target"] = water_elev["Massa Target"].astype(float)

        #differene the floodplain and water elevations to find the water depths
        water_elev["water_depth"] = water_elev["Massa Target"] - water_elev["fp_elev"] + offset

        return water_elev

    def extract_aoi(self, gdf, shp, polygon):
        """take a geodataframe and extract points from and area of interest"""

        # Load the shapefile containing the polygons
        polygons = gpd.read_file(shp)

        # Select the polygons based on an attribute value
        selected_polygons = polygons[polygons['id'] == polygon]

        # Combine the geometries of the selected polygons into a single geometry
        combined_geometry = selected_polygons.unary_union

        # Subset the points to those within the combined geometry
        subset_gdf = gdf[gdf.within(combined_geometry)]

        return subset_gdf
