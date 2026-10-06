[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15758712.svg)](https://doi.org/10.5281/zenodo.15758712)
[![DOI](https://img.shields.io/badge/DOI-10.1029%2F2025GL117121-blue.svg)](https://doi.org/10.1029/2025GL117121)

# Sauthoff-2025-GRL
Code and Data for Sauthoff et al., (2026) "Dynamic Boundaries of Antarctic Active Subglacial Lakes Reveal Underestimated Water Volume Change and Overestimated Lakebed Active Area" in _Geophysical Research Letters_.

## Versions
* **v1.1** (in preparation): renames the six Neckel and others (2021) Jutulstraumen Glacier lakes. That study labeled its lakes only by figure panel, which v1.0 used as `JG_` + panel label. In v1.1 each is named `{basin}_{distance}`: basin is the MEaSUREs refined drainage basin whose grounding-line segment is nearest the lake's subglacial-water-routing pour point, and distance is the straight-line distance (km, rounded) from the lake outline's centroid to the nearest point on that segment. Lake data and geometric calculations are otherwise unchanged from v1.0.

    | v1.0 | v1.1 |
    |---|---|
    | `JG_D2_a` | `Jutulstraumen_162` |
    | `JG_Combined_D2_b_E1` | `Jutulstraumen_163` |
    | `JG_D1_b` | `Jutulstraumen_166` |
    | `JG_D1_a` | `Jutulstraumen_169` |
    | `JG_Combined_E2_F2` | `Jutulstraumen_174` |
    | `JG_F1` | `Jutulstraumen_187` |

    The same mapping is in `output/lake_outlines/renamed_lakes.csv`.
* **v1.0**: version archived with the paper when it was submitted for peer review.

## Licenses
- **Code**: Licensed under GPL-3.0 (see LICENSE-CODE)
- **Data**: Licensed under CC-BY-SA-4.0 (see LICENSE-DATA)

## /Input
* `CryoSat2_mode_masks` folder contains mode masks used to determine when SARIn mode expanded.
* `lake_outlines` folder contains active lake outlines from publicly available data sets or via correspondence with authors.

## Analysis and plotting workflow

### 0_lake_locations.ipynb
* Notebook collates Antarctic active subglacial lakes from past outline inventories as well as point data of lakes in the latest inventory as well as individual studies that are not included past inventories to generate the most recent active subglacial lake inventory.

### 0_preprocess_data.ipynb
* Notebook pre-processes altimetry data sets to construct a multi-mission CryoSat-2 to ICESat-2 (2010 to present) time series.

### Fig1_subglacial_lake_distribution.ipynb
* Notebook generates Fig. 1.

### FigS1_lake_reexamination_methods.ipynb
* Notebook does data analysis to re-examine previously identified active subglacial lakes and creates Fig. S1 plotting the lake re-examination methods.

### Figs23_S23_lake_reexamination_results.ipynb
* Notebook generates Figs. 2, 3, S2, and S3.

## /Output
* `cycle_dates.csv` is dataframe listing satellite cycle start and end datetimes from the multi-mission altimetry data set used for temporal analysis.
* `CryoSat2_SARIn_mode_masks` folder contains polygons of the CryoSat-2 SARIn mode coverage areas used in Fig. 1.
* `geometric_calcs` folder contains csv files of geometric variables (e.g., active area, dh, dV) for each re-examined active subglacial lake and continentally integrated summation files using four analysis approaches stored in subfolders:
    * `evolving_outlines_geom_calc`: evolving outlines, evolving outlines (forward filled)
    * `stationary_outline_geom_calc`: stationary outlines, evolving outlines union.]
* `lake_outlines` folder contains geojson files of lake outlines organized in subfolders, plus `renamed_lakes.csv`, which maps lake names changed since v1.0 (old name, new name, version, naming rule):
    * `evolving_outlines`: evolving outlines for each re-examined lakes (including a 'forward_fill' subfolder for that analysis approach).
    * `stationary_outlines`: five files of stationary outlines served in geojson format
        * Smith and others, 2009 inventory
        * Siegfried and Fricker, 2018 inventory
        * `stationary_outlines_gdf`: stationary outlines collated from two prior outline inventories, an inventory with only point locations, and individual studies since inventories 
        * `reexamined_stationary_outlines_gdf`: revised version of the original stationary outlines that only includes lakes re-examined in this study

## Misc
* ./clean_commit.sh shell script is used to ensure Jupyter notebooks are saved without cell outputs. In terminal, run `./clean_commit.sh` to execute, which will use a hidden file, .pre-commit-config.yaml, to accomplish this.
* .gitignore lists file types not committed to Git

## Notes
* Throughout notebooks the "updated stationary outline" is referred to using "evolving outlines union outline" or a variant that was based on its methodological construction. During peer review "evolving outlines union outline" was changed to "updated stationary outline" for simplicity and clarity, but was mostly not changed in the codebased to avoid breaking the code.
