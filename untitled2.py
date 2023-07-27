# -*- coding: utf-8 -*-
"""
Created on Mon Jun 12 14:27:29 2023

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

extract_gdf = gpd.read_file(r"E:\sensitivity_test\extracted_nl.shp")

sorted_df = extract_gdf.sort_values(['LINE_ID', 'DIST'], ascending=[True, False])

# grouped_df = sorted_df.groupby('LINE_ID')

# # Create an empty GeoDataFrame to store selected lines
# selected_lines = gpd.GeoDataFrame()

# # Iterate over each group
# for name, group in tqdm(grouped_df):
#     # Reorder the index (starts with 0)
#     group = group.reset_index(drop=True)
    
#     # Check the last value in the 'Z' column
#     last_z_value = group['Z'].iloc[-1]
    
#     # Remove the group if the last 'Z' value is 0
#     if last_z_value == 80:
#         continue
    
#     # Append the group to the selected lines
#     selected_lines = pd.concat([selected_lines, group])
    
grouped_df = sorted_df .groupby('LINE_ID')

# Create an empty GeoDataFrame to store selected lines
new_lines = gpd.GeoDataFrame()

# iterate over each group
temp =[]
for name, group in tqdm(grouped_df):
    # reorder the index (starts with 0)
    # iterate over each row in the group
    group = group.reset_index(drop=True) #
    index = group.loc[(group['elevation'] >= 6) | (group['Z'] == 50)].index
    if len(index) == 0:
        new_lines = pd.concat([new_lines,group])
    else:
        first_appear_index = index[0]
        new_lines = pd.concat([new_lines, group.iloc[[first_appear_index]]])
        
last_rows = new_lines.groupby('LINE_ID').last().reset_index()

gdf = gpd.GeoDataFrame(
      last_rows, geometry=gpd.points_from_xy(last_rows.X, last_rows.Y))
output_file = r'E:\extracted_nl_6.shp'
gdf.to_file(output_file, driver='ESRI Shapefile')

# gdf_2 = gpd.GeoDataFrame(
#       selected_lines, geometry=gpd.points_from_xy(selected_lines.X, selected_lines.Y))
# output_file = r'E:\sensitivity_test\gree_2.shp'
# gdf_2.to_file(output_file, driver='ESRI Shapefile')

        

        