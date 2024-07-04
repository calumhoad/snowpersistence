# Repository for the manuscript: *Snow persistence influences vegetation metrics central to Arctic greening analyses*

## Content:
This repository contains the code and data necessary to replicate data analysis, figures and tables in:

Hoad, C. G., Myers-Smith, I. H., Kerby, J. T., Colesie, C. and Assmann, J. J. *Snow persistence influences vegetation metrics central to Arctic greening analyses*

## Contact:
Calum G. Hoad

Email: calum.hoad [at] ed.ac.uk

Website: https://calumhoad.github.io/

## Acknowledgements

(from the manuscript)

We would like to thank everyone who helped with field data collection in the Canadian Yukon during 2022 and in Greenland during 2023, including Joseph Everest, Erica Zaja, Jiri Subrt, Sian Williams and Mariana García Criado. For help with drones and sensors, particular thanks go to Tom Wade at the University of Edinburgh Airborne Research and Innovation facility, and Alex Merrington, Jack Gillespie, Craig Atkins and Robbie Ramsay at the NERC Field Spectroscopy Facility. Additional thanks to Alan Hobbs, Colin Kay and Graham Mitchell from the NERC Geophysical Equipment Facility.

We thank Tim Gyger for support and consultation on our statistical methods, Gwenn Flowers for the time taken to provide climate data for Kluane and Kirsten Schmidt-Pedersen for sharing her extensive knowledge of the people, plants and animals of Qeqertarsuaq, Greenland.

Funding for this research was provided by NERC through a SENSE CDT studentship (NE/T00939X/1), Tundra Time (NE/W006448/1), a 2023 UK - Greenland Bursary, a Geophysical Equipment Facility loan (1152), and a Field Spectroscopy Facility loan (891.0111). Additional funding was provided by a Scottish Alliance for GeoScience, Environment and Society (SAGES) small grant scheme award.

We thank Kluane First Nation, Champagne and Aishihik First Nations and the people of Qeqertarsuaq for their permission to work on their lands. We thank Outpost Research Station and Arctic Station for logistical support.

## Data

Due to the limitations of GitHub, imagery from drones and satellites are not stored in this repository. We intend to make drone data publicly available upon conclusion of Calum Hoad's PhD, however this may be made available upon request by contacting Calum directly (see contact details above). All satellite imagery can be downloaded from either [Copernicus browser](https://browser.dataspace.copernicus.eu/?zoom=5&lat=50.16282&lng=20.78613&themeId=DEFAULT-THEME&visualizationUrl=https%3A%2F%2Fsh.dataspace.copernicus.eu%2Fogc%2Fwms%2Fa91f72b5-f393-4320-bc0f-990129bd9e63&datasetId=S2_L2A_CDAS&demSource3D=%22MAPZEN%22&cloudCoverage=30&dateMode=SINGLE) or [NASA EarthData Search](https://search.earthdata.nasa.gov/search), as per the manuscript. 

The data we derived from satellite and drone imagery are stored within this repository in tabular (.csv) format and enable the statistical analyses to be run.

### Data locations:
**NDVI time series derived from Sentinel-2 and NASA HLSS30 data:**
[/data/ndvi/](/data/ndvi)
- s2-blaesedalen-ndvi-ts-pt.csv
- s2-kluane-low-ndvi-ts-pt.csv
- s2-kluane-high-ndvi-ts-pt.csv
- s30-blaesedalen-ndvi-ts-pt.csv

**NDVI curves modelled using smoothed-spline and Beck (2006):**
[/data/ndvi/](/data/ndvi)
- s2-bl-beck.csv
- s2-kl-beck.csv
- s2-kh-beck.csv
- s2-bl-smooth.csv
- s2-kl-smooth.csv
- s2-kh-smooth.csv
- s30-bl-smooth.csv

**Snow metrics derived from drone imagery, within Sentinel-2 and HLSS30 pixels at each plot:**
[/data/snow/](/data/snow/)
- snow-cover-10m-blaesedalen.csv
- snow-cover-10m-kluane-low.csv
- snow-cover-10m-kluane-high.csv
- snow-cover-30m-blaesedalen.csv

**Snow and smoothed-spline NDVI curves combined to a single dataframe per plot and satellite data product:**
[/data/combined-ndvi-snow/](/data/combined-ndvi-snow/)
- s2-bl-smooth-joined.csv
- s2-kl-smooth-joined.csv
- s2-kh-smooth-joined.csv
- s30-bl-smooth-joined.csv

**Statistical output tables for all models run as part of this manuscript:**
[/plots/stats-tables/](/plots/stats-tables/)
- linear-models-max.html
- linear-models-kh-ln-max.html
- linear-models-doy.html
- linear-models-kh-ln-doy.html
- inla-doy-summaries.html
- inla-max-summaries.html
- inla-bl-max-break-summaries.html

## Scripts:
Scripts for running all analyses included in the manuscript can be found at the following locations. Numbers indicate the sequential order in which scripts should be run.

1. [1-ndvi-ts-extraction-s2.R](/scripts/r/1-ndvi-ts-extraction-s2.R) Script for extracting an NDVI time series from Sentinel-2 data at all three plots.
2. [2-ndvi-ts-extraction-s30.R](/scripts/r/2-ndvi-ts-extraction-s30.R) Script for extracting an NDVI time series from NASA HLSS30 data at Blaesedalen.
3. [3-ndvi-apply-curve-fit.R](/scripts/r/3-ndvi-apply-curve-fit.R) Script which takes the functions in [a-ndvi-curve-fitting-functions.R](/scripts/r/a-ndvi-curve-fitting-functions.R) and applied them to the data from [1] and [2] to create NDVI time series modelled using smoothed-spline and Beck (2006). The script then extracts the peak NDVI and its timing.
4. [4-snow-calculate-coverage-metrics.R](/scripts/r/4-snow-calculate-coverage-metrics.R) Script which classifies areas of snow cover in drone imagery time series, calculates the percentage of satellite pixels (Sentinel-2, HLSS30) covered by snow at each drone time step, then integrates the extent and duration of snow cover within satellite pixels to a single metric representative of snow persistence.
5. [5-combine-ndvi-snow-data.R](/scripts/r/5-combine-ndvi-snow-data.R) Script which joins the smoothed-spline NDVI time series with the calculated snow metrics for every satellite pixel.
6. [6-fit-variograms.R](/scripts/r/6-fit-variograms.R) Script which fits variograms to each variable (snow, peak NDVI, peak NDVI DoY) in Sentinel-2 data at Kluane plots and both Sentinel-2 and HLSS30 data at Blaesedalen.
7. [7-fit-all-models.R](/scripts/r/7-fit-all-models.R) Script which fits all statistical models (OLS, INLA Matern 2D) to the data, using functions defined in [b-modelling-helper-functions.R](/scripts/r/b-modelling-helper-functions.R) and outputs model tables to [/plots/stats-tables/](/plots/stats-tables/).
