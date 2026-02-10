# 3d plotting routines

draw.tiling = function(tile,shell,zoom=1,col='grey') {
  draw.cube = function(x,y,z,rotationmatrix=diag(c(1,1,1)),zoom=1,col='grey') {
    rot.lines3d = function(x,y,z,rot,zoom,...) {
      x = rot%*%rbind(x,y,z)*zoom
      lines3d(t(x),...)
    }
    rot.lines3d(x+c(0.5,0.5,-0.5,-0.5,0.5),y+c(0.5,-0.5,-0.5,0.5,0.5),z+c(0.5,0.5,0.5,0.5,0.5),rotationmatrix,zoom,col=col)
    rot.lines3d(x+c(0.5,0.5,-0.5,-0.5,0.5),y+c(0.5,-0.5,-0.5,0.5,0.5),z-c(0.5,0.5,0.5,0.5,0.5),rotationmatrix,zoom,col=col)
    rot.lines3d(x+c(0.5,0.5),y+c(0.5,0.5),z+c(-0.5,0.5),rotationmatrix,zoom,col=col)
    rot.lines3d(x+c(0.5,0.5),y-c(0.5,0.5),z+c(-0.5,0.5),rotationmatrix,zoom,col=col)
    rot.lines3d(x-c(0.5,0.5),y+c(0.5,0.5),z+c(-0.5,0.5),rotationmatrix,zoom,col=col)
    rot.lines3d(x-c(0.5,0.5),y-c(0.5,0.5),z+c(-0.5,0.5),rotationmatrix,zoom,col=col)
  }
  
  rgl.hold()
  for (itile in seq(length(tile$tile_id))) {
    ishell = tile$shell[itile]
    sxx = shell$transformation$rotation_xx[ishell]
    sxy = shell$transformation$rotation_xy[ishell]
    sxz = shell$transformation$rotation_xz[ishell]
    syx = shell$transformation$rotation_yx[ishell]
    syy = shell$transformation$rotation_yy[ishell]
    syz = shell$transformation$rotation_yz[ishell]
    szx = shell$transformation$rotation_zx[ishell]
    szy = shell$transformation$rotation_zy[ishell]
    szz = shell$transformation$rotation_zz[ishell]
    rotationmatrix = rbind(c(sxx,sxy,sxz),c(syx,syy,syz),c(szx,szy,szz))
    draw.cube(tile$center_x[itile],tile$center_y[itile],tile$center_z[itile],rotationmatrix,as.numeric(zoom),col=col)
  }
  rgl.draw()
}

draw.shells = function(shell,zoom=1) {
  rgl.hold()
  if (shell$dc_min[1]>0) rgl.ball(0,0,0, radius = shell$dc_min[1]*as.numeric(zoom))
  for (i in seq_along(shell$shell_id)) {
    rgl.ball(0,0,0, radius = as.numeric(shell$dc_max[i])*as.numeric(zoom),col='#000000',alpha=0.05)
  }
  rgl.draw()
}