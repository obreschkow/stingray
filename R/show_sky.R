# user
path = '/Users/do/Data/SURFS/stingray/shark/output/'
file = paste0(path,'mocksky.hdf5')
color = rainbow(100,end=2/3)

# load data
gal = h5read(file,'galaxies',bit64conversion='bit64')
tile = h5read(file,'tiling')
para = h5read(file,'parameters')

# make cartesian coordinate
gal$x = sph2car(gal$ra,gal$dec,gal$dc) # [Mpc/h]

# rotation matrix
rotationmatrix = t(cbind(c(0,0,1),c(1,0,0),c(0,1,0)))%*%para$skyrotation

# coloring
color_property = gal$snapshot
transparency_property = gal$mag

# plot
rgl.closeall()
rgl.tiling(tile,rotationmatrix,para$box_length)
alpha = 1-(transparency_property-min(transparency_property))/(max(transparency_property)-min(transparency_property))
col = color[(1-(color_property-min(color_property))/(max(color_property)-min(color_property)))*99+1]
points3d(gal$x, col = col, alpha = alpha)
