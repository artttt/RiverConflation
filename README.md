# RiverConflation
A python based tool for spatially conflating 2 GIS river networks

DISCLAIMER
This tool is under development and is considered a pre alpha-release.

## Description
Takes 2 stream networks and identifies congruent locations between them based on similarity of catchments.
Requires 2 vector stream networks and their associated vector subcatchments with an attribute that allows for a oin of the stream link data to the subcatchment data. Since this data can come in many different forms you will probably need to do some massaging to get this into the form required. some examples are given.

For mor detail at this stage see the code

## Developing
I am developing RiverConflation using using ipython notebook started with the --script command included. This gives a .py copy of each notebook so they can be easily imported.

All editing should be in the ipynb files as the .py files are overwritten automatically.

## Requirements
At this stage this is simpily an Ipython notebook developed an minimally tested with python 2.7 (Anaconda 64bit on windows).

See the imports in the code for a full list of modules required.

Memory and CPU time requirements can be significant for large river networks.
