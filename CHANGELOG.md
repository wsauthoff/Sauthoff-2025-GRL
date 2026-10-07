# Changelog

Versions of the code and data for Sauthoff and others (2026), *Geophysical Research Letters*
(doi:10.1029/2025GL117121), archived on Zenodo. The concept DOI, doi:10.5281/zenodo.15758711,
always resolves to the latest version; each version also has its own DOI.

## v1.1 (2026-10-06)

Version for the published paper, with the revisions made in response to peer review.

### Geometric calculations
- Reprocessed all geometric calculations.
- Some cycle mid-point datetimes were missed because their timestamps differed from
  `cycle_dates.csv`; they are now included.
- Forward fill carried only the first evolving outline where there were several; it now carries
  all of them.
- Missing output is now NaN rather than zero.

### Figures
- New supporting-information figures:
  - CryoSat-2 compared with ICESat-2 geometric variables;
  - individual-lake dV and dV bias time series;
  - Fig. S5.
- Revised Figs. 3 and S1.
- Renamed `Figs23_S23_lake_reexamination_results.ipynb` to
  `Figs23_S2-6_lake_reexamination_results.ipynb`, since it now generates Figs. S2–S6.

### GRL cover image
- Added `GRL_cover_image.ipynb`, the notebook that generated the cover image submitted with the
  paper (previously kept outside this repository).

### Terminology
- Terminology follows the published paper: "updated stationary outline" (see the README Notes).
  Output files were renamed to match.

### Lake names
`output/lake_outlines/renamed_lakes.csv` maps each renamed lake's earlier label to its name in
v1.1. An earlier label is the name the lake had in v1.0 or, for a lake added in v1.1, its label in
the source study. Labels can repeat across studies (L1 is both a Wingham and others, 2006 lake
and an Arthur and others, 2025 lake), so match on `old_name` together with `cite`.

- **Site_B and Site_C.** The re-examination products (evolving outlines and geometric
  calculations) treat them as one lake, Site_BC. The inventory keeps them as separate rows.
- **Arthur and others (2025).**
  - Six lakes that v1.0 named by that study's figure labels are renamed
    `{ice shelf}_{distance}`, where distance is the lake's distance from the grounding line (km)
    as reported in Arthur and others (2025), Table 1:

    | v1.0 | v1.1 |
    |---|---|
    | `M1` | `Muninisen_5` |
    | `M2` | `Muninisen_15` |
    | `R1` | `Roi_Baudouin_19` |
    | `R2` | `Roi_Baudouin_115` |
    | `R3` | `Roi_Baudouin_136` |
    | `V1` | `Vigridisen_54` |
    | `L1` (not in v1.0) | `Lazarevisen_32` |

  - Added `Lazarevisen_32`, published as L1 and renamed by the same rule. v1.0 left it out of
    the inventory and the re-examined lakes by mistake, because its label duplicates a Wingham
    and others (2006) lake, L1, which keeps its name.
  - `0_lake_locations.ipynb` now stops with an error when a source lake shares a name with a
    lake already in the inventory, instead of skipping it as a duplicate.
- **Neckel and others (2021), Jutulstraumen Glacier.**
  - That study labeled its six lakes only by figure panel, which v1.0 used as `JG_` + panel
    label.
  - Each is now named `{basin}_{distance}`:
    - basin is the MEaSUREs refined drainage basin whose grounding-line segment is nearest the
      lake's subglacial-water-routing pour point (Sauthoff and others, in prep.);
    - distance is the straight-line distance (km, rounded) from the lake outline's centroid to
      the nearest point on that segment.
  - The rename changes no data.
  - The published Supporting Information (Table S1) describes these names as a "JG" prefix
    with the distance from the grounding line. v1.1 uses the basin name, `Jutulstraumen`, as
    the prefix instead, matching the `{ice shelf}_{distance}` names of the Arthur and others
    (2025) lakes.

    | v1.0 | v1.1 |
    |---|---|
    | `JG_D2_a` | `Jutulstraumen_162` |
    | `JG_Combined_D2_b_E1` | `Jutulstraumen_163` |
    | `JG_D1_b` | `Jutulstraumen_166` |
    | `JG_D1_a` | `Jutulstraumen_169` |
    | `JG_Combined_E2_F2` | `Jutulstraumen_174` |
    | `JG_F1` | `Jutulstraumen_187` |

- **lower Conway, lower Mercer and upper Engelhardt subglacial lakes.** Renamed
  `LowerConwaySubglacialLake` → `lowerConwaySubglacialLake`,
  `LowerMercerSubglacialLake` → `lowerMercerSubglacialLake` and
  `UpperEngelhardtSubglacialLake` → `upperEngelhardtSubglacialLake`, in the inventories and in
  every output file named after them. This matches the published Supporting Information,
  Table S1, which writes "lower" and "upper" in lowercase.

## v1.0 (2025-06-27)

Version archived with the paper when it was submitted for peer review
(doi:10.5281/zenodo.15758712).
