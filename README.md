# Chapter 1 analyses repo
This repository contains the analyses and data for Chapter 1 of Calum Hoad's PhD. 

## Background, aims and questions:
The Arctic is shown in a large number of studies to have 'greened' since 1982. This 'greening' is characterised as a generally positive trend in vegetation indices (principly maxNDVI) over the period of observation, with vegetation growth driven by increasing Arctic temperatures frequently cited as the likely cause of the trend. 

However, there has been increasing recognition in the literature of the complexities and heterogeneity of Arctic greening (and browning) trends. In particular, spatial and temporal resolution of the EO data used in greening analyses have been cited as complicating factors - with shifts in either the spatial or temporal resolution of analyses having the potential to reverse the direction of a trend.

Over the 

Chapter 1 of the PhD shall therefore focus on:

***The impact of within pixel snow cover on the vegetation indices derived from Earth Obvservation pixles (Landsat, Sentinel-2).***

## Repo structure
- [data](data): Containing all data used in the analyses
    - [lsatTS-output](lsatTS-output): Output CSV from the LandsatTS script over the three field sites
    - [sensor-info](sensor-info): Sensor characteristics for the MAIA S2, Mavic 3 Multispectral
    - [uav](uav): Orthomosaics created from field data
- [plots](plots): As it says on the tin
- [scripts](scripts): Scripts related to analyses and data vis
    - [r](r): Scripts in R
        - [LandsatTS](LandsatTS): Modified functions from LandsatTS package
    - [Python](Python): Scripts in Python