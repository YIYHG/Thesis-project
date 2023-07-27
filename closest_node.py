# -*- coding: utf-8 -*-
"""
Created on Mon Jun 19 00:11:41 2023

@author: yiyang huang
"""

import numpy as np
from shapely.geometry import Point
import os
import pandas as pd
import geopandas as gpd
import rasterio
from osgeo import gdal, ogr, osr
# from honeybees.library.raster import coords_to_pixels
from scipy.spatial import cKDTree
from shapely.geometry import Point, LineString
import itertools
from operator import itemgetter

# def closest_node(
#     orig: np.ndarray,
#     dest: np.ndarray,
#     verboise=True
# ) -> np.ndarray:
#     '''This function approximates the closest point for an array of coordinates. Used here to determine the coastal section closest to an agent's location
#     Args:
#         orig: array containing latitude and longitude of point_0.
#         dest: array containing latitude and longitude of point_1
#     Returns:
#         closest_node: the index of the closes point_1 for each point_0.'''

#     closest_node = np.full(orig.shape[0], -1)
#     if verboise:
#         print(f'Finding closest point for {orig.shape[0]} points...')
#     # infoprint

#     nr_nodes = np.arange(0, orig.shape[0], 1E4)

#     for i, node in enumerate(orig):
#         # start with 'small' frame
#         j = 0
#         search_rad = 0.5
#         dest_clipped = np.array([])

#         while dest_clipped.size == 0:
#             indices_closest = np.where(np.logical_and(
#                 dest[:, 0] > node[0] - search_rad, dest[:, 0] < node[0] + search_rad))[0]
#             dest_clipped = dest[indices_closest, :]
#             search_rad += 0.5
#             j += 1
#             if j > 100:
#                 print('stuck in while loop')
#                 dest_clipped = dest

#         dist_2 = np.sum((dest_clipped - node)**2, axis=1)
#         closest_node[i] = indices_closest[np.argmin(dist_2)]
#         if i in nr_nodes:
#             print(f'currently at node {i} of {orig.shape[0]}')

#     return closest_node
#######################################ADJUST##########################################


orig_gdf = gpd.read_file(r"E:\RCP45_2100\pro_sew.shp")
dest_gdf = gpd.read_file(r"E:\results\results_10m\extracted_sew.shp")



def ckdnearest(gdfA, gdfB, gdfB_cols=['DIST', 'LINE_ID']):
    A = np.concatenate(
        [np.array(geom.coords) for geom in gdfA.geometry.to_list()])
    B = [np.array(geom.coords) for geom in gdfB.geometry.to_list()]
    B_ix = tuple(itertools.chain.from_iterable(
        [itertools.repeat(i, x) for i, x in enumerate(list(map(len, B)))]))
    B = np.concatenate(B)
    ckd_tree = cKDTree(B)
    dist, idx = ckd_tree.query(A, k=1)
    idx = itemgetter(*idx)(B_ix)
    gdf = pd.concat(
        [gdfA, gdfB.loc[idx, gdfB_cols].reset_index(drop=True),
         pd.Series(dist, name='dist')], axis=1)
    return gdf

c = ckdnearest(orig_gdf, dest_gdf)

c.to_csv(r'E:\closest_node\elevation10_RCP4.5\sew.csv', index=False)
#####################################geodataframe#####################################
# def closest_node(
#     orig: gpd.GeoDataFrame,
#     dest: gpd.GeoDataFrame,
#     verbose=True
# ) -> np.ndarray:
#     '''This function approximates the closest point for each point in the original GeoDataFrame. It returns an array with the closest point indices.
#     Args:
#         orig: GeoDataFrame containing the original points with latitude and longitude.
#         dest: GeoDataFrame containing the destination points with latitude and longitude.
#     Returns:
#         closest_node: the index of the closest point in `dest` for each point in `orig`.'''

#     orig_coords = np.array(orig.geometry.apply(lambda geom: (geom.x, geom.y))).reshape(-1, 2)
#     dest_coords = np.array(dest.geometry.apply(lambda geom: (geom.x, geom.y)))

#     closest_node = np.full(orig_coords.shape[0], -1)
#     if verbose:
#         print(f'Finding closest point for {orig_coords.shape[0]} points...')

#     nr_nodes = np.arange(0, orig_coords.shape[0], 1E4)

#     for i, node in enumerate(orig_coords):
#         # start with 'small' frame
#         j = 0
#         search_rad = 0.5
#         dest_clipped = np.array([])

#         while dest_clipped.size == 0:
#             indices_closest = np.where(np.logical_and(
#                 dest_coords[:, 0] > node[0] - search_rad, dest_coords[:, 0] < node[0] + search_rad))[0]
#             dest_clipped = dest_coords[indices_closest, :]
#             search_rad += 0.5
#             j += 1
#             if j > 100:
#                 print('stuck in while loop')
#                 dest_clipped = dest_coords

#         dist_2 = np.sum((dest_clipped - node)**2, axis=1)
#         closest_node[i] = indices_closest[np.argmin(dist_2)]
#         if i in nr_nodes:
#             print(f'currently at node {i} of {orig_coords.shape[0]}')

#     return closest_node

# Load the GeoDataFrames
# orig_gdf = gpd.read_file(r"E:\projection_eu\pro_bg.shp")
# orig = orig_gdf.to_numpy
# dest_gdf = gpd.read_file(r"E:\results\extracted_bg.shp")
# dest = dest_gdf.to_numpy

# # Call the closest_node function
# closest_node_indices = closest_node(orig, dest)


# df = c

# negative_values = df[df['1'] < 0] #select the erosion value
# negative_values['1'] = negative_values['1'].abs() #convert to absolute value
# negative_values['1'] = negative_values[['1', 'DIST']].max(axis=1) #if the erosion exceeds the limitation, and replace


# df.update(negative_values)
# df.to_csv(r'E:\closest_node\elevation10_RCP8.5\sew.csv', index=False)
