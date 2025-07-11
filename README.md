# Sauthoff-2025-GRL
Code for the data analysis and figure plotting in Sauthoff et al., submitted to _Geophysical Research Letters_. 

## Licenses
- **Code**: Licensed under GNU GPL v3 (see LICENSE-CODE)
- **Data**: Licensed under CC BY 4.0 (see LICENSE-DATA)

## Inputs
* `CryoSat2_mode_masks` folder contains mode masks used to determine when SARIn mode expanded.
* `lake_outlines` folder contains active lake outlines from publicly available data sets or via correspondence with authors.

## Analysis and plotting workflow

### 0_lake_locations.ipynb
* Notebook collates Antarctic active subglacial lakes from past outline inventories as well as point data of lakes in the latest inventory as well as individual studies that are not included past inventories to generate the most recent active subglacial lake inventory.

### 0_preprocess_data.ipynb
* Notebook pre-processes altimetry data sets to construct a multi-mission CryoSat-2 to ICESat-2 (2010 to present) time series.

### Fig1_subglacial_lake_distribution.ipynb
* Notebok generates Fig. 1.

### FigS1_lake_reexamination_methods.ipynb
* Notebook does data analysis to re-examine previously identified active subglacial lakes and creates Fig. S1 plotting the lake re-examination methods.

### Figs23_S23_lake_reexamination_results.ipynb
* Notebook generates Figs. 2, 3, S2, and S3.

## Outputs
* `cycle_dates.csv` is dataframe listing satellite cycle start and end datetimes from the multi-mission altimetry data set used for temporal analysis.
* `CryoSat2_SARIn_mode_masks` folder contains polygons of the CryoSat-2 SARIn mode coverage areas used in Fig. 1.
* `geometric_calcs` folder contains csv files of geometric variables (e.g., active area, dh, dV) for each re-examined active subglacial lake and continentally integrated summation files using four analysis approaches stored in subfolders:
    * `evolving_outlines_geom_calc`: evolving outlines, evolving outlines (forward filled)
    * `stationary_outline_geom_calc`: stationary outlines, evolving outlines union.]
* `lake_outlines` folder contains geojson files of lake outlines organized in subfolders:
    * `evolving_outlines`: evolving outlines for each re-examined lakes (including a 'forward_fill' subfolder for that analysis approach).
    * `stationary_outlines`: five files of stationary outlines served in geojson format
        * Smith and others, 2009 inventory
        * Siegfried and Fricker, 2018 inventory
        * `stationary_outlines_gdf`: stationary outlines collated from two prior outline inventories, an inventory with only point locations, and individual studies since inventories 
        * `reexamined_stationary_outlines_gdf`: revised version of the original stationary outlines that only includes lakes re-examined in this study
