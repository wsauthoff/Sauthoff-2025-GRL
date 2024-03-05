# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.14.1
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# + tags=[]
# TODO
# Add upper slessor lakes (get outlines from matt)
# Add PIG and AP lakes (find citations) if not too close to grounding line
# Should Carter and others, 2013 citation be changed to Fricker and Scambos, 2009? no, redefined as a lake complex? read carter paper or listen to voice memo

# + tags=[]
# This notebook collates Antarctic active subglacial lakes from past outline inventories 
# (Smith and others, 2009; Siegfried & Fricker, 2018) as well as point data of lakes in the latest inventory 
# (Livingstone and others, 2022) as well as individual studies since.
#
# Written 2023-01-17 by W. Sauthoff (sauthoff@mines.edu)

# + tags=[]
# Install package
# %pip install openpyxl

# + tags=[]
# Import packages
import fiona
import geopandas as gpd
import h5py
import math
import os
import numpy as np
import pandas as pd
from shapely.geometry import Point, Polygon, MultiPolygon
from pyproj import CRS, Geod, Transformer

# Define data and script directories dependent on home environment
if os.getenv('HOME') == '/home/jovyan':
    DATA_DIR = '/home/jovyan/data'
    SCRIPT_DIR = '/home/jovyan/repos_my/script_dir'
elif os.getenv('HOME') == '/Users/Wilson': 
    DATA_DIR = '/Volumes/ExtremeSSD/data'
    SCRIPT_DIR = '/Users/Wilson/Documents/0-code/repos_my/script_dir'

# Define utility functions
def ll2ps(lon, lat):
    """
    Transform coordinates from geodetic coordinates (lon, lat)
    to Antarctic Polar Stereograph coordinates (x, y)
    x, y = ll2ps(lon, lat)
    """
    crs_ll = CRS("EPSG:4326")
    crs_xy = CRS("EPSG:3031")
    ll_to_xy = Transformer.from_crs(crs_ll, crs_xy, always_xy = True)
    x, y = ll_to_xy.transform(lon, lat)
    return x, y

def ps2ll(x, y):
    """
    Transform coordinates from Antarctic Polar Stereograph
    to geodetic (lon, lat) coordinates
    
    lon, lat = ps2ll(x, y)
    """
    crs_ll = CRS("EPSG:4326")
    crs_xy = CRS("EPSG:3031")
    xy_to_ll = Transformer.from_crs(crs_xy, crs_ll, always_xy = True)
    lon, lat = xy_to_ll.transform(x, y)
    return lon, lat

def find_intersections(gdf1, gdf2):
    # Create an empty list to store the results
    intersections = []

    # Iterate over each geometry in gdf1
    for index1, geom1 in gdf1.geometry.items():
        # Compare with each geometry in gdf2
        for index2, geom2 in gdf2.geometry.items():
            if geom1.intersects(geom2):
                # If they intersect, add the indices to the list
                intersections.append((index1, index2))

    return intersections


# + [markdown] user_expressions=[]
# # Import datasets

# + tags=[]
# Import active subglacial lake outlines from Smith and others (2009) (S09)
# As released in Smith and others, 2012 dataset (https://doi.org/10.15784/601439)
fiona.drvsupport.supported_drivers['KML'] = 'rw'
S09_outlines_lonlat = gpd.read_file(DATA_DIR + '/boundaries/Smith2009_subglacial_lakes/Antarctic_lakes.kml', driver='KML')

# Ensure GeoDataFrame is in EPSG:4326 for geodesic area calculation
if S09_outlines_lonlat.crs != 'EPSG:4326':
    S09_outlines_lonlat.to_crs('EPSG:4326')

# Create a Geod object for calculating area on the WGS84 ellipsoid
geod = Geod(ellps="WGS84")

# Create empty list to store areas
areas = []

# Iterate through each row in the GeoDataFrame
for index, row in S09_outlines_lonlat.iterrows():
    geom = row.geometry

    # Check if the geometry is a MultiPolygon and handle accordingly
    if isinstance(geom, MultiPolygon):
        # Calculate area for each polygon in the MultiPolygon
        area = sum(abs(geod.geometry_area_perimeter(polygon)[0]) for polygon in geom)
    elif isinstance(geom, Polygon):
        # Calculate the area for a single Polygon
        area = abs(geod.geometry_area_perimeter(geom)[0])
    else:
        # Handle other geometry types (if any) or set to NaN or 0
        area = float('nan')  # or 0, depending on how you want to handle errors

    areas.append(area)

# Convert to crs to epsg:3031
S09_outlines = S09_outlines_lonlat.to_crs(3031)

# Append list of areas to geodataframe
S09_outlines['area (m^2)'] = areas

# Delete original geodataframe
del S09_outlines_lonlat

# + tags=[]
S09_outlines

# + tags=[]
# Import active subglacial lake outlines from Siegfried & Fricker (2018) (SF18)
# Original pub: https://doi.org/10.1017/aog.2017.36 
# Code for loading lake outlines available in code bank associated with Siegfried & Fricker (2021), https://doi.org/10.1029/2020GL091089: 
# https://github.com/mrsiegfried/Siegfried2021-GRL/blob/main/data/outlines/load_lakes.ipynb

# import subglacial lake outlines (Siegfried & Fricker, 2018)
h5f = h5py.File(DATA_DIR + '/boundaries/SiegfriedFricker2018_subglacial_lakes/SiegfriedFricker2018-outlines.h5', 'r')
outline_geometries = [] # store polygons
citations = [] # store citation information

# we're going to calculate geodesic lake area because that is often screwed up
# and occasionally incorrect in the literature
areas = []

# we're going to need to do some coordinate transforms for the geodesic area
# define CRS for Antarcica and make a converter from xy to ll
CRS_LL = "EPSG:4326" # wgs84 in lon,lat
CRS_XY = h5f.attrs.get('proj_crs') # get projection from hdf5 file
XY_TO_LL = Transformer.from_crs(CRS_XY, CRS_LL, always_xy = True) # make coord transformer
geod = CRS(CRS_LL).get_geod() # geod object for calculating geodesic area on defined ellipsoid

# look through each lake and load all of it's info
for lake in h5f.keys():
    outline_x = h5f[lake]['x'][:]
    outline_y = h5f[lake]['y'][:]
    outlines_xy = np.stack((outline_x, outline_y),axis=2).reshape(outline_x.shape[1], 2)

    # A single lake with multiple polygons is NaN broken---need to identify and
    # load as a MultiPolygon. Otherwise it's easy (just load as polygon)
    if np.isnan(outlines_xy)[:,0].sum() == 0:
        geometry = Polygon(outlines_xy)
        lon, lat = XY_TO_LL.transform(outlines_xy[:,0], outlines_xy[:,1])
        this_area = abs(geod.polygon_area_perimeter(lon,lat)[0])
    else:
        this_area = 0
        # break at NaN values and load each as separate polygons
        idx = np.where(np.isnan(outlines_xy[:,0]))[0]

        # grab outline of first lake before getting into the loop
        this_outline = outlines_xy[0:idx[0],:]
        pgons = [Polygon(this_outline)] # put the first polygon in a list
        lon,lat = XY_TO_LL.transform(this_outline[:,0], this_outline[:,1])
        this_area += abs(geod.polygon_area_perimeter(lon,lat)[0])/1e6 # add its area
        for i in np.arange(0,len(idx)):
            if i == len(idx)-1:
                this_outline = outlines_xy[idx[i]+1:,:]
            else:
                this_outline = outlines_xy[idx[i]+1:idx[i+1]]

            pgons.append(Polygon(this_outline))
            lon,lat = XY_TO_LL.transform(this_outline[:,0], this_outline[:,1])
            this_area += abs(geod.polygon_area_perimeter(lon,lat)[0])/1e6
        geometry = MultiPolygon(pgons)

    # append all the results in the right place
    outline_geometries.append(geometry)
    citations.append(h5f[lake].attrs.get('citation')[0].decode('UTF-8'))
    areas.append(this_area)

# make a pandas dataframe with all the necessary info
df = pd.DataFrame(zip(h5f.keys(), outline_geometries, areas, citations),
                  columns=['name', 'geometry', 'area (m^2)', 'cite'])
# convert to geopands geodataframe
SF18_outlines = gpd.GeoDataFrame(df, crs=CRS_XY, geometry=outline_geometries)
# close HDF5 file
h5f.close()

# + tags=[]
SF18_outlines

# + tags=[]
# Read subglacial lake point data from Livingstone and others (2022) (L22), https://doi.org/10.1038/s43017-021-00246-9
url = 'https://static-content.springer.com/esm/art%3A10.1038%2Fs43017-021-00246-9/MediaObjects/43017_2021_246_MOESM1_ESM.xlsx'
use_cols = ['Name', 'Lat.  oN', 'Lon. oE', 'Lake Type', 'References', 'Prior Inventory']
import_rows = np.arange(0,676)
L22_points = pd.read_excel(url, sheet_name='Antarctica', usecols=use_cols, skiprows = lambda x: x not in import_rows)

# View just the active lakes of the pandas dataset
L22_activelake_points = L22_points[L22_points['Lake Type'] == 'Active']

# Reset the index, dropping the old one
L22_activelake_points = L22_activelake_points.reset_index(drop=True)

# Display geodataframe
L22_activelake_points

# + tags=[]
# Isolate lakes not previously included in Smith and others (2009; 2017) or Siegfried and Fricker (2018)
L22_activelake_points_new = L22_activelake_points.copy(deep=True)
L22_activelake_points_new.drop(L22_activelake_points_new.loc[L22_activelake_points['References'].str.contains('Smith et al.|Siegfried & Fricker')].index, inplace=True)
L22_activelake_points_new = L22_activelake_points_new.reset_index(drop=True)
# print(len(L22_activelake_points_new))
L22_activelake_points_new

# + tags=[]
# Copy Siegfried & Fricker (2018) outline inventory to add entries from more recent publications
lake_locations = SF18_outlines.copy(deep=True)

# Create conditional list of lakes within CryoSat-2 SARIn coverage
CS2SARIn_lakes = ['ConwaySubglacialLake', 'Cook_E1', 'Cook_E2', 'David_1', 'David_s1', 
    'EngelhardtSubglacialLake', 'Foundation_1', 'Institute_E1', 'Institute_W1', 'KT1', 'KT2', 'KT3',
    'Lake10', 'Lake12', 'Lake78', 'Lambert_1', 'LennoxKing_1', 'Mac1', 'Mac2', 'Mac3', 
    'MercerSubglacialLake', 'Ninnis_1', 'Ninnis_2', 'Rec1', 'Rutford_1', 'Slessor_1', 
    'Slessor_23', 'Totten_1', 'Totten_2', 'Thw_70', 'Thw_124', 'Thw_142', 'Thw_170', 
    'UpperSubglacialLakeConway', 'Whillans_6', 'WhillansSubglacialLake', 'Wilkes_1']

# Add new column and populate with whether lake is within CryoSat-2 (CS2) InSAR coverage
lake_locations['CS2_SARIn'] = np.where(lake_locations['name'].isin(CS2SARIn_lakes),'True', 'False')

# Display dataframe preview
lake_locations

# + tags=[]
# Change journal abbreviations to ISO4 standard
# First look at citations
lake_locations['cite'].unique()

# + tags=[]
# Replace TC with Cryosphere
lake_locations = lake_locations.replace('Kim and others, 2016, TC, doi:10.5194/tc-10-2971-2016', 'Kim and others, 2016, Cryosphere, doi:10.5194/tc-10-2971-2016')
lake_locations = lake_locations.replace('Smith and others, 2017, TC, doi:10.5194/tc-11-451-2017', 'Smith and others, 2017, Cryosphere, doi:10.5194/tc-11-451-2017')
# Replace GRL with Geophys. Res. Lett.
lake_locations = lake_locations.replace('McMillan and others, 2013, GRL, doi:10.1002/grl.50689', 'McMillan and others, 2013, Geophys. Res. Lett., doi:10.1002/grl.50689')
# Ensure replacements worked as expected
lake_locations['cite'].unique()

# + tags=[]
# Find area of previously identified lakes inventoried in SF18
# to use as guesstimate for lakes without this information or figures 
# allowing for closer guesstimate
np.mean(SF18_outlines['area (m^2)'])/1e6

# + tags=[]
# Add entries for newer lakes from publications not included in Siegfried & Fricker (2018) inventory 
# (using approximated centroid point when outline is unavailable)

# Smith 2009 Recovery_8 (dropped from Siegfried & Fricker (2018) outlines) to see if there's activity
lake_gdf = S09_outlines[S09_outlines['Name'] == 'Recovery_8']
name = lake_gdf['Name'].values[0]
# S09 outline inventory uses 3D polygons with z dimension vs. 2D polygons in SF18 inventory
# Extract the point values that define the perimeter of the polygon to make polygon without third z dimension
# Extract 2D coordinates (X, Y) from the 3D polygon (X, Y, Z)
xy_coords = [(x, y) for x, y, z in lake_gdf.geometry.values[0].exterior.coords]
# Create a new 2D polygon from these coordinates
lake_poly_2d = Polygon(xy_coords)
geometry = lake_poly_2d
area = lake_gdf['area (m^2)'].values[0]
# Store citation info from another lake from the same S09 study in the SF18 citation format
cite = SF18_outlines[SF18_outlines['name'] == 'Bindschadler_1'].cite.values[0]
# Make entry into pandas dataframe to concatenate to inventory
df = pd.DataFrame([[name, geometry, area, cite, 'False']], columns=lake_locations.columns)
# Convert the DataFrame to a GeoDataFrame
# gdf = gpd.GeoDataFrame(df, crs='EPSG:3031')
# Ensure that new entry isn't already in inventory before adding
df_diff = df[~df['name'].isin(lake_locations['name'])]
# gdf_diff = gdf[~gdf['name'].isin(lake_locations['name'])]
# Add entry to inventory
lake_locations = pd.concat([lake_locations, df_diff], ignore_index=True)

# + tags=[]
# Livingstone and others, 2022 active lakes points (no outlines available):
# Wingham and others, 2006 (within Livingstone and others, 2022 inventory)
# Lakes L1, U1, U2, U3
areas = [600e6, 200e6, 200e6, 400e6]  # Lake L1 area estimated in paper, but U1-3 were guestimated from Fig. 1
for i in range(0, 4):
    name = L22_activelake_points_new.iloc[i]['Name']
    lon = L22_activelake_points_new.iloc[i]['Lon. oE']
    lat = L22_activelake_points_new.iloc[i]['Lat.  oN']
    geometry = Point(ll2ps(lon, lat)).buffer(math.sqrt(areas[i] / math.pi))
    area = areas[i]
    cite = 'Wingham and others, 2006, Nature, doi:10.1038/nature04660'
    CS2_SARIn = 'False'
    # Make entry into pandas dataframe to concatenate to inventory
    df = pd.DataFrame([[name, geometry, area, cite, CS2_SARIn]], columns=lake_locations.columns)
    # Ensure that new entry isn't already in inventory before adding
    df_diff = df[~df['name'].isin(lake_locations['name'])]
    # Add entry to inventory
    lake_locations = pd.concat([lake_locations, df_diff], ignore_index=True)
# Delete list
del areas

# N. Young personal comm. (within Wright & Siegert, 2012 and Livingstone and others, 2022 inventories)
# Lakes Site A, Site B, Site C
area = 200e6
for i in range(6, 9):
    name = L22_activelake_points_new.iloc[i]['Name']
    lon = L22_activelake_points_new.iloc[i]['Lon. oE']
    lat = L22_activelake_points_new.iloc[i]['Lat.  oN']
    geometry = Point(ll2ps(lon, lat)).buffer(math.sqrt(area / math.pi))
    cite = 'Wright & Siegert, 2012, Antarct. Sci., doi:10.1017/S095410201200048X'
    CS2_SARIn = 'True'
    # Make entry into pandas dataframe to concatenate to inventory
    df = pd.DataFrame([[name, geometry, area, cite, CS2_SARIn]], columns=lake_locations.columns)
    # Ensure that new entry isn't already in inventory before adding
    df_diff = df[~df['name'].isin(lake_locations['name'])]
    # Add entry to inventory
    lake_locations = pd.concat([lake_locations, df_diff], ignore_index=True)

# Hoffman and others, 2020 (within Livingstone and others, 2022 inventory)
# Western Thwaites (WT) lake
area = SF18_outlines[SF18_outlines['name'] == 'Thw_142']['area (m^2)'].values[0]  # Using nearby Thwaites lake of similar area from Fig. 1
for i in range(9, 10):
    name = L22_activelake_points_new.iloc[i]['Name']
    lon = L22_activelake_points_new.iloc[i]['Lon. oE']
    lat = L22_activelake_points_new.iloc[i]['Lat.  oN']
    geometry = Point(ll2ps(lon, lat)).buffer(math.sqrt(area / math.pi))
    cite = 'Hoffman and others, 2020, Cryosphere, doi:10.5194/tc-14-4603-2020'
    CS2_SARIn = 'True'
    # Make entry into pandas dataframe to concatenate to inventory
    df = pd.DataFrame([[name, geometry, area, cite, CS2_SARIn]], columns=lake_locations.columns)
    # Ensure that new entry isn't already in inventory before adding
    df_diff = df[~df['name'].isin(lake_locations['name'])]
    # Add entry to inventory
    lake_locations = pd.concat([lake_locations, df_diff], ignore_index=True)

# Livingstone and others, 2022 inventory contains a point for Haynes Glacier system of four active lakes detailed in 
# Hoffman and others, 2020; however, only one point for four lakes and point is adjacent to four lakes, not one of four lakes, so adding manually
# Haynes Glacier (HG) lakes point location and area estimated from Hoffman and others, 2020 supplement Fig. S1 (https://doi.org/10.5194/tc-14-4603-2020-supplement)
# TL96
name = 'TL96'
area = 25e6
geometry = Point(-1439000, -544000).buffer(math.sqrt(area / math.pi))
cite = 'Hoffman and others, 2020, Cryosphere, doi:10.5194/tc-14-4603-2020'
CS2_SARIn = 'True'
# Make entry into pandas dataframe to concatenate to inventory
df = pd.DataFrame([[name, geometry, area, cite, CS2_SARIn]], columns=lake_locations.columns)
# Ensure that new entry isn't already in inventory before adding
df_diff = df[~df['name'].isin(lake_locations['name'])]
# Add entry to inventory
lake_locations = pd.concat([lake_locations, df_diff], ignore_index=True)

# TL108
name = 'TL108'
area = 15e6
geometry = Point(-1427500, -542000).buffer(math.sqrt(area / math.pi))
# Make entry into pandas dataframe to concatenate to inventory
df = pd.DataFrame([[name, geometry, area, cite, CS2_SARIn]], columns=lake_locations.columns)
# Ensure that new entry isn't already in inventory before adding
df_diff = df[~df['name'].isin(lake_locations['name'])]
# Add entry to inventory
lake_locations = pd.concat([lake_locations, df_diff], ignore_index=True)

# TL115
name = 'TL115'
area = 5e6
geometry = Point(-1422500, -537500).buffer(math.sqrt(area / math.pi))
# Make entry into pandas dataframe to concatenate to inventory
df = pd.DataFrame([[name, geometry, area, cite, CS2_SARIn]], columns=lake_locations.columns)
# Ensure that new entry isn't already in inventory before adding
df_diff = df[~df['name'].isin(lake_locations['name'])]
# Add entry to inventory
lake_locations = pd.concat([lake_locations, df_diff], ignore_index=True)

# TL122
name = 'TL122'
area = 10e6
geometry = Point(-1419000, -532500).buffer(math.sqrt(area / math.pi))
# Make entry into pandas dataframe to concatenate to inventory
df = pd.DataFrame([[name, geometry, area, cite, CS2_SARIn]], columns=lake_locations.columns)
# Ensure that new entry isn't already in inventory before adding
df_diff = df[~df['name'].isin(lake_locations['name'])]
# Add entry to inventory
lake_locations = pd.concat([lake_locations, df_diff], ignore_index=True)

# Siegfried and Fricker, 2021 (publication not included in Livingstone and others, 2022 inventory)
# Lower Subglacial Lake Mercer (LSLM)
name = 'LowerMercerSubglacialLake'
area = 10e6
geometry = Point(-308000, -509000).buffer(math.sqrt(area / math.pi))
cite = 'Siegfried and Fricker, 2021, Geophys. Res. Lett., doi:10.1029/2020GL091089'
CS2_SARIn = 'True'
# Make entry into pandas dataframe to concatenate to inventory
df = pd.DataFrame([[name, geometry, area, cite, CS2_SARIn]], columns=lake_locations.columns)
# Ensure that new entry isn't already in inventory before adding
df_diff = df[~df['name'].isin(lake_locations['name'])]
# Add entry to inventory
lake_locations = pd.concat([lake_locations, df_diff], ignore_index=True)

# Lower Subglacial Lake Conway (LSLC)
name = 'LowerSubglacialLakeConway'
area = 10e6
geometry = Point(-282000, -502000).buffer(math.sqrt(area / math.pi))
cite = 'Siegfried and Fricker, 2021, Geophys. Res. Lett., doi:10.1029/2020GL091089'
CS2_SARIn = 'True'
# Make entry into pandas dataframe to concatenate to inventory
df = pd.DataFrame([[name, geometry, area, cite, CS2_SARIn]], columns=lake_locations.columns)
# Ensure that new entry isn't already in inventory before adding
df_diff = df[~df['name'].isin(lake_locations['name'])]
# Add entry to inventory
lake_locations = pd.concat([lake_locations, df_diff], ignore_index=True)

# Drop outlines from lake_locations that are duplicative of SiegfiedFricker2018_outlines
lake_locations_nonSF18 = lake_locations.copy(deep=True)
lake_locations_nonSF18.drop(lake_locations_nonSF18.iloc[0:131].index, inplace=True)

# Sort geopandas geodataframe to alphabetize by lake name and reset index
# lake_locations_nonSF18 = lake_locations_nonSF18.sort_values(by=['name'], ignore_index=True)

# Reset the index, dropping the old one
lake_locations_notSF18 = lake_locations_nonSF18.reset_index(drop=True)

# + tags=[]
# Display geodataframe
lake_locations

# + tags=[]
# Display geodataframe
lake_locations_notSF18
# -





# +
# Code to to explore and reconcile outline inventories

# + tags=[]
S09_outlines

# + tags=[]
for idx in range(len(SF18_outlines['cite'].unique())): 
    length = len(SF18_outlines[SF18_outlines['cite'].str.contains(SF18_outlines['cite'].unique()[idx])])
    print(SF18_outlines['cite'].unique()[idx])
    print(length)

# + tags=[]
SF18_outlines[SF18_outlines['cite'] == 'Fricker & Scambos, 2009, J. Glac., doi:10.3189/002214309788608813']

# + tags=[]
SF18_outlines[SF18_outlines['cite'] == 'McMillan and others, 2013, GRL, doi:10.1002/grl.50689']

# + tags=[]
SF18_outlines[SF18_outlines['cite'] == 'Kim and others, 2016, TC, doi:10.5194/tc-10-2971-2016']

# + tags=[]
SF18_outlines[SF18_outlines['cite'] == 'Carter and others, 2013, J. Glac., doi:10.3189/2013JoG13J085']

# + tags=[]
SF18_outlines[SF18_outlines['cite'] == 'Fricker and others, 2010, J. Glac., doi:10.3189/002214310791968557']

# + tags=[]
SF18_outlines[SF18_outlines['cite'] == 'Fricker and others, 2014, J. Glac., doi:10.3189/2014JoG14J063']

# + tags=[]
SF18_outlines[SF18_outlines['cite'] == 'Siegfried & Fricker, 2018, Ann. Glac., doi:10.1017/aog.2017.36']

# +
# Siegfried and Fricker (2018) only contains 97 Smith and others (2009) citations due to: 
# 124 lakes in Smith and others, 2009, J. Glac., doi:10.3189/002214309789470879
# -7 Fricker & Scambos, 2009, J. Glac., doi:10.3189/002214309788608813 redelineates Mercer/Whillans lakes
# -8 Fricker and others, 2010, J. Glac., doi:10.3189/002214310791968557 redelineates MacAyeal Ice Stream lakes
# -1 Carter and others, 2013, J. Glac., doi:10.3189/2013JoG13J085 redelineates Mercer_1 as Lake78
# -9 Fricker and others, 2014, J. Glac., doi:10.3189/2014JoG14J063 redelineates Recovery Glacier lakes
# -2 Siegfried & Fricker, 2018, Ann. Glac., doi:10.1017/aog.2017.36 redelineates Slessor2 and Slessor3 as Slessor23
# 97 citations to Smith and others, 2009, J. Glac., doi:10.3189/002214309789470879 in Siegfried and Fricker, 2018 inventory
# -



# + tags=[]
for idx in range(len(L22_activelake_points['References'].unique())): 
    length = len(L22_activelake_points[L22_activelake_points['References'].str.contains(L22_activelake_points['References'].unique()[idx], regex=False)])
    print(L22_activelake_points['References'].unique()[idx])
    print(length)

# + tags=[]
L22_activelake_points.loc[L22_activelake_points['References'].str.contains('Smith et al. (2009)', regex=False)]
# -

SF18_outlines['name'].tolist()



# + tags=[]
# View head of Smith and others, 2009 static outline geodataframe
S09_outlines.head()

# + tags=[]
S09_outlines[S09_outlines.Name == 'Recovery_8']

# + tags=[]


# + tags=[]
# Find outlines in SiegfriedFricker2018_outlines static outline geodataframe with Smith and others, 2009 citation
# These lakes have not been redefined with a new outline and can be analysed as is
SF18_outlines_S09only = SF18_outlines[SF18_outlines['cite'].isin(['Smith and others, 2009, J. Glac., doi:10.3189/002214309789470879'])]
SF18_outlines_S09only

# + tags=[]
# Find outlines in SiegfriedFricker2018_outlines geodataframe with citations not matching Smith and others, 2009
# These are eitherly newly identified lakes or lakes that have been redefined with new static outlines 
SF18_outlines_SF18only = SF18_outlines[~SF18_outlines['cite'].isin(['Smith and others, 2009, J. Glac., doi:10.3189/002214309789470879'])]
SF18_outlines_SF18only

# + tags=[]
len(SF18_outlines) == len(SF18_outlines_S09only) + len(SF18_outlines_SF18only)

# + tags=[]
# Find spatial intersections between two static outline geodataframes to indicate groups of outlines that are redefined outlines of one lake
S09_SF18_intersections = find_intersections(S09_outlines, SF18_outlines)

# + tags=[]
# Extract intersecting indicies from each static outline geodataframe for later use
S09_intersection_indices = [inner_list[0] for inner_list in S09_SF18_intersections]
SF18_intersection_indices = [inner_list[1] for inner_list in S09_SF18_intersections]

# + tags=[]
# Store the lake names from Smith and others, 2009 and Siegfried & Fricker, 2018 that intersect
# Store 'name' column for rows in gdf1 corresponding to indices in indices1
S09_intersections_names = S09_outlines.iloc[S09_intersection_indices]['Name']

# Store 'name' column for rows in gdf2 corresponding to indices in indices2
SF18_intersections_names = SF18_outlines.iloc[SF18_intersection_indices]['name']

# Concatenate into dataframe
intersecting_lake_names = pd.concat([S09_intersections_names.reset_index(drop=True), SF18_intersections_names.reset_index(drop=True)], axis=1)
intersecting_lake_names.columns = ['S09', 'SF18']

# Filtering rows where column1 and column2 differ indicating lake was redelineated
intersecting_lake_names_diff = intersecting_lake_names[intersecting_lake_names['S09'] != intersecting_lake_names['SF18']]

# + tags=[]
print(intersecting_lake_names.to_string())

# + tags=[]
print(intersecting_lake_names_diff)
# -

# Extract rows with intersection between the Smith and others, 2009 static outlines and the static outlines in Siegfried & Fricker, 2018
# These are lakes that were defined in Smith and others, 2009 and subsequently renamed or redelineated by publications after Smith and others, 2009
# We will analyze as a unified polygon of all past static outlines for each lake
S09_intersections = S09_outlines[S09_outlines.index.isin(S09_intersection_indices)] 
S09_intersections

# + tags=[]
# Do the same for the Siegfried & Fricker, 2018 static outlines
SF18_intersections = SF18_outlines[SF18_outlines.index.isin(SF18_intersection_indices)] 
SF18_intersections

# + tags=[]
# Extract rows with no intersection between the Smith and others, 2009 and the renamed or redefined static outlines in Siegfried & Fricker, 2018
# These are new lakes that were delineated by Smith and others, 2009 and not subsequently rename or redelineated
# These lakes only have one static outline and can be analysed as is
S09_non_intersections = S09_outlines[~S09_outlines.index.isin(S09_intersection_indices)] 
S09_non_intersections

# + tags=[]
# Extract rows with no intersection between the Smith and others, 2009 and the renamed or redefined static outlines in Siegfried & Fricker, 2018
# These are new lakes that were added by publications subsequent to Smith and others, 2009
# These lakes only have one static outline and can be analysed as is
SF18_non_intersections = SF18_outlines[~SF18_outlines.index.isin(SF18_intersection_indices)] 
SF18_non_intersections

# + tags=[]
len(S09_outlines) == len(S09_intersections) + len(S09_non_intersections)

# + tags=[]
len(S09_intersections)

# + tags=[]
len(S09_non_intersections)  # This is Recovery_8 that was dropped in SF18 inventory

# + tags=[]
len(SF18_outlines) == len(SF18_intersections) + len(SF18_non_intersections)

# + tags=[]
len(SF18_intersections)  # Lakes from S09 inventory that were either keep as is or redelinated; there is one more than S09_intersections because 

# + tags=[]
len(SF18_non_intersections)  # These are new lakes added in the SF18 inventory from other publications since S09 inventory
# -



# Create list of S09 lakes redelineated in SF18
S09_SF18_redelineated_lake_names = [   
    ['Cook_E2', 'Cook_E2'],
    ['KambTrunk_1', 'KT1'],
    ['Macayeal_1', 'Mac1'],
    ['Macayeal_2', 'Mac2'],
    ['Mercer_1', 'Lake78'],
    ['Mercer_2', 'MercerSubglacialLake'],
    ['Recovery_3', 'Rec2'],
    ['Recovery_4', 'Rec3'],
    ['Recovery_5', 'Rec4'], 
    ['Recovery_6', 'Rec5'],
    ['Recovery_7', 'Rec6'],
    ['Recovery_9', 'Rec8'],
    ['Recovery_10', 'Rec9'],
    ['Recovery_11', 'Rec10'],     
    ['Whillans_1', 'EngelhardtSubglacialLake'],            
    ['Whillans_2a', 'Lake12'],
    ['Whillans_2b', 'Lake10'],
    ['Whillans_3', 'WhillansSubglacialLake'],
    ['Whillans_4', 'ConwaySubglacialLake'],
    ['Whillans_5', 'UpperSubglacialLakeConway'],
    # Two S09 lakes converted to one SF18 lake
    [['Slessor_2', 'Slessor_3'], 'Slessor_23'],
    [['Recovery_1', 'Recovery_2'], 'Rec1'],
    # One S09 lake converted to two SF18 lakes
    ['Macayeal_3', ['Mac4', 'Mac5']]
]



# +
# Add rows for lake candidates I've observed in ICESat-2 data (using approx centroid point until outline is established)
# New observation at ~–450 km, –540 km (polar stereographic x, y)
# Sauthoff2024_outlines.loc[len(Sauthoff2024_outlines)] = ['lower_Whillans6', 'POINT (-450000, -540000)', 'nan', 'nan', 'Sauthoff and others, in prep', 'True'] 

# display dataframe preview
# Sauthoff2024_outlines
