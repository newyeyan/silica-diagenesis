# silica-diagenesis
# chertDiagenesis

1D diagenetic model for chert formation in marine sediments, revise from Tatzel and Frings (2022)

## Description
This repository contains the scripts used to define a model that simulates the formation of chert (crystalline quartz via amorphous silica and opal-CT precursors) in marine sediments. 

## Information and instructions
Running the requires the pre-installation of Matlab and an active license (/https://de.mathworks.com/products/matlab.html). These files have been tested in Matlab R2024a.

To reproduce the results in the manuscript, download the whole directory and add to your Matlab search path and run the script 'chertKineticModel.m' via the command line or Matlabs GUI. Change the parameters you're interested in (see below).

## Parameters
In the manuscript, results are presented based on varying one or more of the following parameters contained in the first few lines of the script 'chertKineticModel.m'. These are:

Sed_Si_pc: percent surface sediment drymass that is amorphous silica.
heatFlow: in W/m2,
mineralContent: percentage of the sediment that is detrital aluminosilicates
