## This code computes structural tree metrics of height and coverage from produced Canopy Height Models
# Code also applies to other study sites, by changing file names 
# Input buffer zone to either impact/non impact beaver zones depending on area of interest

# Author: jg857
# Version 2.7
# Date: 05/08/2025
# Enhanced by: ChatGPT

# Import relevant packages
import os
import numpy as np
import pandas as pd
import rasterio
from rasterio.mask import mask
import geopandas as gpd

# Set working directory
os.chdir('/Users/joegartell/Desktop/python_data/11_06')

# Define function to assess metrics of interest within beaver foraging/non foraging zone as the area of interest
def area(chm_path, buffer_shapefile, output_csv, canopy_threshold=2.0):

    # Read buffer shp
    buffer_gdf = gpd.read_file(buffer_shapefile)
    buffer_gdf = buffer_gdf.to_crs("EPSG:27700")

    # Open CHM raster
    with rasterio.open(chm_path) as src:
        pixel_area = abs(src.transform.a * src.transform.e)

        results = []

        # Loop each buffer polygon
        for idx, row in buffer_gdf.iterrows():
            try:
                geom = [row.geometry]
                masked_data, masked_transform = mask(src, geom, crop=True, filled=False)
                masked_array = masked_data[0]

                # Get CHM values
                valid_values = masked_array.compressed()

                if len(valid_values) > 0:
                    # Mean height of vegetation
                    non_zero_values = valid_values[valid_values > 0]
                    mean_height = np.mean(non_zero_values) if len(non_zero_values) > 0 else 0

                    # Canopy coverage (%)
                    canopy_pixels = np.sum(valid_values >= canopy_threshold)
                    total_pixels = len(valid_values)
                    canopy_coverage_percent = (canopy_pixels / total_pixels) * 100

                    # Canopy area (hectares)
                    canopy_area_m2 = canopy_pixels * pixel_area
                    canopy_area_ha = canopy_area_m2 / 10000

                    # Buffer area in hectares
                    buffer_area_m2 = row.geometry.area
                    buffer_area_ha = buffer_area_m2 / 10000

                    # Collect results
                    stats = {
                        'buffer_id': idx,
                        'buffer_area_ha': buffer_area_ha,
                        'mean_height': mean_height,
                        'canopy_pixels': canopy_pixels,
                        'canopy_coverage_percent': canopy_coverage_percent,
                        'canopy_area_ha': canopy_area_ha,
                        'canopy_threshold_m': canopy_threshold
                    }

                    # Include info from buffer shapefile
                    for col in buffer_gdf.columns:
                        if col != 'geometry':
                            stats[f'buffer_{col}'] = row[col]

                    results.append(stats)

                    print(f"Buffer {idx}: Mean height = {mean_height:.2f}m, "
                          f"Canopy coverage = {canopy_coverage_percent:.1f}%, "
                          f"Canopy area = {canopy_area_ha:.2f} ha")

                else:
                    print(f"Buffer {idx}: No valid data found")

            except Exception as e:
                print(f"Error processing buffer {idx}: {str(e)}")

    # Convert to DataFrame
    results_df = pd.DataFrame(results)

    # Raise error in the case of no data
    if results_df.empty:
        raise ValueError("ERROR")

    # Save results to CSV
    results_df.to_csv(output_csv, index=False)

    # Print summary statistics
    print(f"\nSummary Statistics:")
    print(f"Mean vegetation height : {results_df['mean_height'].mean():.2f}m")
    print(f"Average canopy coverage: {results_df['canopy_coverage_percent'].mean():.1f}%")
    print(f"Total canopy area: {results_df['canopy_area_ha'].sum():.2f} hectares")

    return results_df

# Define function which compares before and after beaver impact
def beaver_activity_comparison(chm_2010_path, chm_2020_path, buffer_shapefile, output_csv, canopy_threshold=2.0):
    print("Before Beaver Impact")
    stats_2010 = area(chm_2010_path, buffer_shapefile, 
                                           output_csv=f"before_{output_csv}",
                                           canopy_threshold=canopy_threshold)
    
    print("\nAfter Beaver Impact")
    stats_2020 = area(chm_2020_path, buffer_shapefile, 
                                           output_csv=f"after_{output_csv}",
                                           canopy_threshold=canopy_threshold)
    
    if stats_2010.empty or stats_2020.empty:
        raise ValueError("Data missing")
    
    # Create comparison DataFrame
    comparison_df = pd.DataFrame()
    comparison_df['buffer_id'] = stats_2010['buffer_id']
    comparison_df['buffer_area_ha'] = stats_2010['buffer_area_ha']
    
    # Before beaver impact
    comparison_df['mean_height_2006_m'] = stats_2010['mean_height']
    comparison_df['canopy_coverage_2006_percent'] = stats_2010['canopy_coverage_percent']
    comparison_df['canopy_area_2006_ha'] = stats_2010['canopy_area_ha']
    
    # After beaver impact
    comparison_df['mean_height_2019_m'] = stats_2020['mean_height']
    comparison_df['canopy_coverage_2019_percent'] = stats_2020['canopy_coverage_percent']
    comparison_df['canopy_area_2019_ha'] = stats_2020['canopy_area_ha']
    
    # Changes
    comparison_df['height_change_m'] = stats_2020['mean_height'] - stats_2010['mean_height']
    comparison_df['coverage_change_percent'] = stats_2020['canopy_coverage_percent'] - stats_2010['canopy_coverage_percent']
    comparison_df['area_change_ha'] = stats_2020['canopy_area_ha'] - stats_2010['canopy_area_ha']
    
    # Add buffer attributes
    for col in stats_2010.columns:
        if col.startswith('buffer_'):
            comparison_df[col] = stats_2010[col]
    
    # Save results
    comparison_df.to_csv(output_csv, index=False)
    
    # Print summary
    print(f"\nComparison Summary:")
    print(f"{'='*50}")
    print(f"Average height change: {comparison_df['height_change_m'].mean():.2f}m")
    print(f"Average coverage change: {comparison_df['coverage_change_percent'].mean():.1f} percentage points")
    print(f"Total area change: {comparison_df['area_change_ha'].sum():.2f} hectares")
    
    return comparison_df

## Run the Functions above 
# Define the single buffer zone (either beaver impact or non impact zone)
buffer_shapefile = "hole_eb_nonimpact.shp"
zone_name = "Non Impact Zone"

# Before Beaver activity
results_2006 = area(
    "EB_chm_2006_webster_clipped.tif", 
    buffer_shapefile, 
    "Canopy Pre Beaver_EastBudleigh.csv", 
    canopy_threshold=2.0
)

# After Beaver activity 
results_2019 = area(
    "EB_chm_2019_webster_clipped.tif", 
    buffer_shapefile, 
    "Canopy Post Beaver_EastBudleigh.csv", 
    canopy_threshold=2.0
)

# Comparison between pre and post beaver activity
comparison = beaver_activity_comparison(
    "EB_chm_2006_webster_clipped.tif", 
    "EB_chm_2019_webster_clipped.tif", 
    buffer_shapefile, 
    "Canopy Change_EastBudleigh.csv", 
    canopy_threshold=2.0
)