# Take archived Minnesota moose telemetry data and create a binary example
# dataset that will be used for the uhcplots package.
library(data.table)
moose12687 <- fread("http://conservancy.umn.edu/bitstream/handle/11299/181607/moose12687.csv?sequence=66&isAllowed=y")
table(moose12687$year)
# Divide step length by 1000 (mm) to make the coefficients similar in SSF
moose12687$step <- moose12687$step/1000
devtools::use_data(moose12687, overwrite = T)
