# pecan
R framework for confidence and prediction interval calculation

lib.pecan:
- Folder with "pecan" functions
- pecan.R: main function for parameter estimation and confidence and prediction interval calculation
- others: functions used by pecan

lib.models:
- Folder with model functions that can be used as inputs to pecan()
- pk models start with "pk."

demo.R
- script running examples for pecan usage and functionality
- examples on PK using ADVAN style implementation can only be run if corresponding function are provided in "lib.advan" folder that is not provided in this repository
