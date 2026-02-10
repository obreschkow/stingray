# stingray path
setwd('~/Dropbox/Code/Fortran/stingray/stingray')

# libraries
library(cooltools)
library(rglplus)
library(sphereplot)
source("R/routines.R")

# set file name of mock sky
file = 'example/mocksky.hdf5'

# read partial data
mock = readhdf5(file, list(galaxies=list(ra='',dec='',dc=''),
                           mapping='*',
                           parameters='*'))

# plot galaxy positions in 3D                
x = sphereplot::sph2car(mock$galaxies$ra,mock$galaxies$dec,mock$galaxies$dc) # [Mpc/h]
dc.max = max(mock$galaxies$dc)
rgl.new(aspect=1, fixed=FALSE)
color = lightness(spectrumcolors(1000,rev=TRUE), 0.4)
col = color[mock$galaxies$dc/dc.max*999+1]
spheres3d(x, col = col, radius=3)
spheres3d(0, 0, 0, col = 'black', radius=6)

# add grid
rgl.hold()
rgl.sphgrid(radius=dc.max, radaxis=FALSE, longtype='D', col.long='grey', col.lat='grey', add=TRUE)
rgl.camera(c(0,0,1e4),fov = 4)
rgl.draw()

# add tiles to plot
draw.tiling(mock$mapping$tiles, mock$mapping$shells, mock$parameters$box_side, col='black')