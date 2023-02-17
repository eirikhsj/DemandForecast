# DemandForecast

This R packages contains the code for the Master Thesis project in Data Science titled:
Energy Demand Forecasting using dimensionality reduced temperature fields (Name Pending) by Eirik Sjåvik, supervised by Alex Lenkoski (NR). 

The package can be installed by:
1: Cloning this repo to a local folder. 
2: Make sure devtools package is installed, or install it through CRAN by install.packages('devtools').
3: Load devtools package by library(devtools).
3: Then load this package using the command devtools::install('PATH/DemandForecast')in e.g. the RStudio console (PATH is the path to the local folder). 
4: Load package by Library(DemandForecast)

The files can roughly be divided up into 3 categories: 
1) Convenience or helper functions, which loads and prepares data. 
2) Model scripts, that takes this data, and perform 
3) R_markdown scripts which provides an overview of test runs and results. 

The data used in this project has a magnitude of over 100 GB, and is not stored on github. 
It is accesible through UiO Data Science hub storage, (upon request, etc. )

A small sample of the data is attached to provide an example of core functionality. 


# The Project

The thesis covers two related problems: 
First: How can we make improvements on a structural model of energy demand by using temperature?
Our results show that by using a gam pca approach we can achieve subtantial improvements (a skill of 0.17 for a 30 day forecast period) compared to a climatology approach. 
This result holds for a structural model (e.g. under the assumption that we know the future temperature).

The main functions used for this problem is demand_forecast.R, which is a rolling cv function which enables the testing and comparison of a host of different energy demand forecast models. 
A detailed breakdown of tests run are found in the RMarkdown documents: 
-
-

Of course, in real life we don't know the future temperatures, which lead us into the second problem.

Second: Can we utilize different temperature forecasts as inputs to our demand model?
Here we focus on five different sub-seasonal(??) forecasting models, where our period of interest is the range 0-30 days. 

The models are:
1) NWP forecast model
2) Reweighted NWP forecast model 
3) Stratospheric Wind model
4) Sea level pressure model 
5) Different combinations of these models

The goal is to find the most accurate forecast which can be utilized in the Energy Demand model. 

The main findings are:
- 
- 

The example data provided thus includes the PC of observed ERA-temperatures, and the predicted temperatures based on the models mentioned above, 
which provides an overview of the project as a whole. 



