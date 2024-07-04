# Repository for the manuscript: *Snow persistence influences vegetation metrics central to Arctic greening analyses*

## Content:
This repository contains the code and data necessary to replicate data analysis, figures and tables in:

Hoad, C. G., Myers-Smith, I. H., Kerby, J. T., Colesie, C. and Assmann, J. J. *Snow persistence influences vegetation metrics central to Arctic greening analyses*

## Contact:
Calum G. Hoad

Email: calum.hoad [at] ed.ac.uk

Website: https://calumhoad.github.io/


## Background, aims and questions:
The Arctic is shown in a large number of studies to have 'greened' since 1982. This 'greening' is characterised as a generally positive trend in vegetation indices (principly maxNDVI) over the period of observation, with vegetation growth driven by increasing Arctic temperatures frequently cited as the likely cause of the trend. 

However, there has been increasing recognition in the literature of the complexities and heterogeneity of Arctic greening (and browning) trends. In particular, spatial and temporal resolution of the EO data used in greening analyses have been cited as complicating factors - with shifts in either the spatial or temporal resolution of analyses having the potential to reverse the direction of a trend.

Over the same period as Arctic greening trends, snow cover extent and duration have undergone change in the Arctic (AMAP, 2017). There are unanswered questions on the relationship between the analyses of Arctic snow cover change and analyses of Arctic greening, especially regarding sub-pixel scale summer snow patches.

Chapter 1 of the PhD shall therefore focus on:

***The impact of within pixel snow cover on the vegetation indices derived from Earth Obvservation pixles (Landsat, Sentinel-2). Could within pixel snow cover be contributing to Arctic greening/browning trends?***

## Repo structure
- [**data**](data): Containing all data used in the analyses
    - [**lsatTS-output**](data/lsatTS-output): Output CSV from the LandsatTS script over the three field sites
    - [**sensor-info**](data/sensor-info): Sensor characteristics for the MAIA S2, Mavic 3 Multispectral
    - [**uav**](data/uav): Orthomosaics created from field data
- [**plots**](plots): As it says on the tin
- [**scripts**](scripts): Scripts related to analyses and data vis
    - [**r**](scripts/r): Scripts in R
        - [**LandsatTS**](scripts/r/LandsatTS): Modified functions from LandsatTS package
    - [**Python**](scripts/Python): Scripts in Python
