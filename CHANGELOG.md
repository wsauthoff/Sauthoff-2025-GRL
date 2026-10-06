# Changelog

Versions of the code and data for Sauthoff and others (2026), *Geophysical Research Letters*
(doi:10.1029/2025GL117121), archived on Zenodo. The concept DOI, doi:10.5281/zenodo.15758711,
always resolves to the latest version; each version also has its own DOI.

## v1.1 (in preparation)

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

### Terminology
- Terminology follows the published paper: "updated stationary outline" (see the README Notes).
  Output files were renamed to match.

### Lake names
`output/lake_outlines/renamed_lakes.csv` maps every v1.0 lake name changed in v1.1 to its new
name.

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

  - Added `Lazarevisen_32`, published as L1. v1.0 left it out of the inventory and the
    re-examined lakes by mistake, because its label duplicates a Wingham and others (2006) lake,
    L1. It is new in v1.1 rather than renamed, so it is not in `renamed_lakes.csv`.
- **Neckel and others (2021), Jutulstraumen Glacier.**
  - That study labeled its six lakes only by figure panel, which v1.0 used as `JG_` + panel
    label.
  - Each is now named `{basin}_{distance}`:
    - basin is the MEaSUREs refined drainage basin whose grounding-line segment is nearest the
      lake's subglacial-water-routing pour point (Sauthoff and others, in prep.);
    - distance is the straight-line distance (km, rounded) from the lake outline's centroid to
      the nearest point on that segment.
  - The rename changes no data.

    | v1.0 | v1.1 |
    |---|---|
    | `JG_D2_a` | `Jutulstraumen_162` |
    | `JG_Combined_D2_b_E1` | `Jutulstraumen_163` |
    | `JG_D1_b` | `Jutulstraumen_166` |
    | `JG_D1_a` | `Jutulstraumen_169` |
    | `JG_Combined_E2_F2` | `Jutulstraumen_174` |
    | `JG_F1` | `Jutulstraumen_187` |

## v1.0 (2025-06-27)

Version archived with the paper when it was submitted for peer review
(doi:10.5281/zenodo.15758712).
