# neighbourhood_forecast


Aims

The aim of this code is twofold:

1. Post-process the deterministic forecasts to get probabilistic forecast using an inexpensive method (Theis et al., 2005)

2. Perform a neighbourood based forecast verification using the fractions skill score



##### post-processing deterministic forecast


'1. define a rainfall accumulation threshold  (physical threshold or percentile) over a given accumulation period (3hours, 6 hours, 12 hours and so on)'

'2. convert the rainfall field to a binary field (ones or zeros) depending whether a grid-point exceed or not the rainfall threshold of point 1'

'3. For each grid-point define a squared neighbourhood (e.g. number of grid points of the side of the square)'

'4. For each neighbourhood size around a given grid-point compute the fractions of grid-points exceeding the threshold within the neighbourhood'

'5. The fraction of grid-point just computed can be interpreted the spatial probability of the event occurrence within that neighbourhood'

'The same methodology can be applied to the observations field repeating the same step'

This process will generate the fractions of the rainfall field


##### fractions skill score

1. to perform forecast verification before applying the steps described above you should REGRID the forecast field into the observed grid'



##### python scripts 


The file nbood_code.py contains all the functions used to perform steps 1-5.

The file neighbourhood_forecasts.py contains the main code and is the one to be run from terminal. You can choose some parameters as specified in the script.

