## This code computes Canopy Height Models, based on LiDAR raster data
# of DSM and DTM derived from the Environment Agency, for East Budleigh
# *note code applied to other study sites, by changing file names* 

# Author: jg857
# Version 3.2
# Date: 17/07/2025
# Enhanced by: ChatGPT

# Import relevant packages
import os
import numpy as np
import pandas as pd
import rasterio
from rasterio.warp import reproject, Resampling
from rasterio.mask import mask
import geopandas as gpd

# Set my working directory
os.chdir('/Users/joegartell/Desktop/python_data/11_06')

# Define a function to fill DSM no data values to corresponding DTM value
def fill_nodata(dsm_path, dtm_path, output_path):
    with rasterio.open(dsm_path) as dsm_src, rasterio.open(dtm_path) as dtm_src:
        dsm_data = dsm_src.read(1).astype('float32')
        dtm_data = dtm_src.read(1).astype('float32')

        # Identify DSM gaps
        if dsm_src.nodata is not None:
            gap_mask = (dsm_data == dsm_src.nodata) | np.isnan(dsm_data)
        else:
            gap_mask = np.isnan(dsm_data)

        filled_dsm = dsm_data.copy()
        filled_dsm[gap_mask] = dtm_data[gap_mask]

        profile = dsm_src.profile.copy()
        profile.update(dtype='float32', nodata=dsm_src.nodata)

        with rasterio.open(output_path, 'w', **profile) as dst:
            dst.write(filled_dsm, 1)

# Define a raster to align raster to template
def template(src_path, template_path, dst_path, resampling_method=Resampling.bilinear):
    with rasterio.open(template_path) as template:
        dst_crs = template.crs
        dst_transform = template.transform
        dst_width = template.width
        dst_height = template.height
        dst_profile = template.profile

    with rasterio.open(src_path) as src:
        dst_array = np.empty((dst_height, dst_width), dtype='float32')
        src_nodata = src.nodata

        reproject(
            source=rasterio.band(src, 1),
            destination=dst_array,
            src_transform=src.transform,
            src_crs=src.crs,
            dst_transform=dst_transform,
            dst_crs=dst_crs,
            resampling=resampling_method,
            src_nodata=src_nodata,
            dst_nodata=src_nodata
        )

        dst_profile.update({
            "height": dst_height,
            "width": dst_width,
            "transform": dst_transform,
            "dtype": 'float32',
            "nodata": src_nodata
        })

        with rasterio.open(dst_path, "w", **dst_profile) as dst:
            dst.write(dst_array, 1)

## Function to clip raster by outline of study site shapefile
# and assign crs 
def clip_shp(raster_path, shapefile_path, output_path):
    gdf = gpd.read_file(shapefile_path)
    gdf = gdf.to_crs("EPSG:27700")  

    with rasterio.open(raster_path) as src:
        out_image, out_transform = mask(src, gdf.geometry, crop=True)
        out_meta = src.meta.copy()
        out_meta.update({
            "driver": "GTiff",
            "height": out_image.shape[1],
            "width": out_image.shape[2],
            "transform": out_transform
        })

        with rasterio.open(output_path, "w", **out_meta) as dest:
            dest.write(out_image)

## Function to calculate Canopy Height model from aligned and clipped DSM and DTM products and on valid pixels
# set a threshold to exclude erronous values
def CHM_Calculation(dsm_path, dtm_path, output_path):
    with rasterio.open(dsm_path) as dsm, rasterio.open(dtm_path) as dtm:
        dsm_data = dsm.read(1).astype('float32')
        dtm_data = dtm.read(1).astype('float32')
        nodata = dsm.nodata

        if nodata is not None:
            valid_mask = (dsm_data != nodata) & (~np.isnan(dsm_data))
        else:
            valid_mask = ~np.isnan(dsm_data)

        chm_data = np.full(dsm_data.shape, nodata, dtype='float32')
        chm_data[valid_mask] = dsm_data[valid_mask] - dtm_data[valid_mask]

        chm_data[chm_data < -2] = 0.0

        profile = dsm.profile
        profile.update(dtype='float32', nodata=nodata)

        with rasterio.open(output_path, 'w', **profile) as dst:
            dst.write(chm_data, 1)

## Run functions 
# Align rasters to original 2010 DSM 
template("trial_join_dsm_2006.tif", "trial_join_dsm_2006.tif", "aligned_2006_dsm.tif", Resampling.bilinear)
template("trial_join_dtm_2006.tif", "trial_join_dsm_2006.tif", "aligned_2006_dtm.tif", Resampling.nearest)
template("trial_join_dsm_2019.tif", "trial_join_dsm_2006.tif", "aligned_2019_dsm.tif", Resampling.bilinear)
template("trial_join_dtm_2019.tif", "trial_join_dsm_2006.tif", "aligned_2019_dtm.tif", Resampling.nearest)

## Fill NoData where relevant 
# based off visual inspection of data quality in ArcGIS Pro
fill_nodata("aligned_2019_dsm.tif", "aligned_2019_dtm.tif", "aligned_2019_dsm_filled.tif")

## Clip the cleaned rasters to the shapefile of East Budleigh
# for other study sites, change name where applicable
clip_shp("aligned_2006_dsm.tif", "EB_STUDY_AREA_.shp", "EB_dsm_2006_clipped.tif")
clip_shp("aligned_2006_dtm.tif", "EB_STUDY_AREA_.shp", "EB_dtm_2006_clipped.tif")
clip_shp("aligned_2019_dsm_filled.tif", "EB_STUDY_AREA_.shp", "EB_dsm_2019_clipped.tif")
clip_shp("aligned_2019_dtm.tif", "EB_STUDY_AREA_.shp", "EB_dtm_2019_clipped.tif")

# Calculate Canopy Height Models
CHM_Calculation("EB_dsm_2006_clipped.tif", "EB_dtm_2006_clipped.tif", "EB_chm_2006_webster.tif")
CHM_Calculation("EB_dsm_2019_clipped.tif", "EB_dtm_2019_clipped.tif", "EB_chm_2019_webster.tif")