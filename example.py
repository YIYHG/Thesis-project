# -*- coding: utf-8 -*-
"""
Created on Fri May 19 20:24:57 2023

@author: yiyang huang
"""
import os
import geopandas as gpd
import pandas as pd
import numpy as np
import rasterio as rio
from rasterio.features import geometry_mask
from rasterio.sample import sample_gen
from shapely.geometry import Point,  LineString,shape, box
import matplotlib.pyplot as plt
import rioxarray as rxr
from scipy.spatial import KDTree
import xarray
from rasterio.warp import calculate_default_transform, reproject, Resampling
import pyproj
import math
from rasterio.transform import Affine
from rasterio.crs import CRS
from tqdm import tqdm
import time



# # # Load the transection lines shapefile
# transect_line = gpd.read_file("D:\Thesis\map\example_hn_island/transection_line_hn.shp")

profile_1 = gpd.read_file(r"E:\eu\nl\profile_1.shp")
profile_2 = gpd.read_file(r"E:\eu\nl\profile_2.shp")

# # Set the path to the folder containing the raster files
landcover_folder = r'E:\eu\nl\lc'
# landcover_folder = 'D:\Thesis\map\example_hn_island\lc_two_tile'

# points_gdf = gpd.read_file(r"D:\Thesis\map\example_hn_island/profile_point_10.shp")
# landcover_folder = r'E:\lc_two_tile'

dem = rio.open(r"E:\eu\nl/dem.tif")
dem_array = dem.read(1).astype('float64')

lc = rio.open(r"E:\eu\nl\ESA_WorldCover_10m_2021_V200_N51E003_Map.tif")
lc_array = lc.read(1).astype('int32')

# lc2 = rio.open(r"E:\eu\nl\ESA_WorldCover_10m_2021_V200_N51E003_Map.tif")
# lc2_array = lc.read(1).astype('int32')

# merge geodataframe
points_gdf = pd.concat([profile_1, profile_2], ignore_index=True)

# # Create an empty geodataframe to store the points
# columns = ['LINE_ID', 'ID', 'DIST', 'DIST_SURF', 'X', 'Y', 'Z', 'geometry']
# points_gdf = gpd.GeoDataFrame(columns=columns)

# # Iterate over each line in the transect geodataframe
# for index, row in tqdm(transect_line.iterrows()):
#     line = row['geometry']
#     line_length = line.length
    
#     dist = 0
    
#     while dist <= line_length:
#         point = line.interpolate(dist)
#         x = point.x
#         y = point.y
        
#         # Calculate distance from the beginning of the line
#         dist_to_start = line.project(point)
        
#         # Create a new row for the point and add it to the points geodataframe
#         new_row = {
#             'LINE_ID': row['TR_ID'],
#             'ID': dist / 10,  # Assuming interval of 10 meters
#             'DIST': dist_to_start,
#             'DIST_SURF': dist,
#             'X': x,
#             'Y': y,
#             'Z': 80.0,  # Adjust the Z value as needed
#             'geometry': Point(x, y)
#         }
#         points_gdf = points_gdf.append(new_row, ignore_index=True)
        
#         dist += 10
        



# # Iterate over each point in the geodataframe
# for index, point in tqdm(points_gdf.iterrows()):
#     x = point['X']
#     y = point['Y']
    
#     # Create a point geometry
#     point_geom = Point(x, y)
    
#     # Iterate over each raster file in the folder
#     for filename in os.listdir(landcover_folder):
#         if filename.endswith('.tif'):
#             raster_path = os.path.join(landcover_folder, filename)
            
#             # Open the raster file
#             with rio.open(raster_path) as src:
#                 # Sample the raster value at the given point coordinate
#                 values = list(src.sample([(x, y)]))
#                 pixel_value = values[0][0]  
                
#                 # Create a new column name based on the raster file name
#                 column_name = f'{filename[:-4]}'
                
#                 # Add the pixel value to the geodataframe
#                 points_gdf.at[index, column_name] = pixel_value

# ##option 3

# coords = [(pt.x, pt.y) for pt in points_gdf.geometry]

# # Iterate over each file in the folder
# for filename in os.listdir(landcover_folder):
#     if filename.endswith('.tif'):
#         raster_path = os.path.join(landcover_folder, filename)
        
#         # Open the raster file
#         with rio.open(raster_path) as src:
#             # Sample the raster values at the point coordinates
#             values = list(src.sample(coords))
#             pixel_values = [sample[0] for sample in values]  
            
            # # Add the pixel values to the GeoDataFrame
            # points_gdf['Z'] = pixel_values
# ### option 4
coords = [(x,y) for x, y in zip(points_gdf.X, points_gdf.Y)]
points_gdf['elevation'] = [x for x in dem.sample(coords)]
points_gdf['land_cover'] = [x for x in lc.sample(coords)]
# points_gdf['land_cover'] = [x for x in lc2.sample(coords)]

###########################################################################################
# # Save the updated geodataframe with the extracted values
# points_gdf.to_file('points_with_values.shp')

extract_gdf = gpd.read_file(r"E:\result\extracted_ele_den_2.shp")

sorted_df = extract_gdf.sort_values(['LINE_ID', 'DIST'], ascending=[True, False])

grouped_df = sorted_df.groupby('LINE_ID')

# create empty gdf to store selected_lines
selected_lines = gpd.GeoDataFrame()

# iterate over each group
for name, group in tqdm(grouped_df):
    # reorder the index (starts with 0)
    group = group.reset_index(drop=True)
    
    # Check the last value in the Z column
    last_z_value = group['Z'].iloc[-1]
    
    # Remove the group if the last 'Z' value is 0
    if last_z_value == 80:
        continue
    
    # append the group meets the require
    selected_lines = pd.concat([selected_lines, group])
    
grouped_df = selected_lines.groupby('LINE_ID')

# create the empty gdf
new_lines = gpd.GeoDataFrame()

# iterate over each group
temp =[]
for name, group in tqdm(grouped_df):
    # iterate each row in group
    group = group.reset_index(drop=True) #
    index = group.loc[(group['elevation'] >= 10) | (group['Z'] == 50)].index
    if len(index) == 0:
        new_lines = pd.concat([new_lines,group])
    else:
        first_appear_index = index[0]
        new_lines = pd.concat([new_lines, group.iloc[[first_appear_index]]])
# only add the last row of each group       
last_rows = new_lines.groupby('LINE_ID').last().reset_index()

gdf = gpd.GeoDataFrame(
      last_rows, geometry=gpd.points_from_xy(last_rows.X, last_rows.Y))
output_file = r'E:\results\extracted_den_2.shp'
gdf.to_file(output_file, driver='ESRI Shapefile')

gdf_2 = gpd.GeoDataFrame(
      selected_lines, geometry=gpd.points_from_xy(selected_lines.X, selected_lines.Y))
output_file = r'E:\sensitivity\extracted_den_2.shp'
gdf_2.to_file(output_file, driver='ESRI Shapefile')
