# ncwrap
A Fortran 90-2003 style wrapper module for reading and writing netcdf data. 
It calles several netcdf library subroutines to read or right the netcdf files with ease. 
It requires netcdf static library distributed by Unidata. It work with netcdf version 3.6 or higher. 

## Compile
Please edit makefile. Note that the required netcdf library are different between netcdf version 3.x (-lnetcdf) and version 4.x (-lnetcdf -lnetcdff).

## Copyright and License
Copyright (C) 2015, Takuto Maeda, All rights reserved. All source codes included in this archive are released under the MIT License. 

The netcdf library is copyrighted by (C) Unidata Program Center/University Corporation for Atmospheric Research. 

