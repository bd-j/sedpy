sedpy
======

modules for storing and operating on astronomical source spectral energy distributions

installation & setup:
	1) download code to some directory 'hdir'. 
	2) add setenv sedpy = 'hdir' to your .bashrc. 
	5)  cp any desired filters from the kcorrect filter directory to 'hdir'/data/filters/

observate.py has methods for generating synthetic photometry through any filters, and classes for dealing with filters generally. Functionality for spectra is being added slowly. With a huge debt to Mike Blanton's k_correct code .

yanny.py (from Erin Sheldon) is for reading filter curves.  

attenuation.py contains dust extenction and attenuation classs and methods (very preliminary).  

modelgrid.py is a module with classes for the storage and interpolation of model SEDs (using numpy structured arrays, usually, and Delaunay triangulation)
