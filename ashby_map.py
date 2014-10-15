# -*- coding: utf-8 -*-
"""
Plot Ashby Charts.

@author: Nicolas Guarin Zapata
@date: October 15, 2014
"""
import numpy as np
from scipy.spatial import ConvexHull
import matplotlib.pyplot as plt
from matplotlib.path import Path
import matplotlib.patches as patches
from matplotlib import rcParams

rcParams['font.family'] = 'serif'
rcParams['font.size'] = 16


def convex_poly(points, color, inc=1.2, rad=0.3, lw=2):
    """
    Plot the convex hull around a set of points as a 
    shaded polygon.
    """
    hull = ConvexHull(points)
    
    
    cent = np.mean(points, 0)
    pts = []
    for pt in points[hull.simplices]:
        pts.append(pt[0].tolist())
        pts.append(pt[1].tolist())
    
    pts.sort(key=lambda p: np.arctan2(p[1] - cent[1],
                                    p[0] - cent[0]))
    pts = pts[0::2]  # Deleting duplicates
    pts.insert(len(pts), pts[0])
    
    
    verts = inc*(np.array(pts)- cent) + cent
    verts2 = np.zeros((3*verts.shape[0]-2,2))
    verts2[0::3] = verts
    verts2[1::3,:] = (1-rad)*verts[0:-1,:] + rad*verts[1:,:]
    verts2[2::3,:] = rad*verts[0:-1,:] + (1-rad)*verts[1:,:]
    verts2[0:-1] = verts2[1:]
    verts2[-1] = verts2[0]
    
    
    codes = [Path.MOVETO, Path.LINETO, Path.CURVE3,]
    for j in range(len(pts)-2):
        codes.extend([Path.CURVE3, Path.LINETO, Path.CURVE3,])
    codes.append(Path.CURVE3)
    
    
    path = Path(verts2, codes)
    patch = patches.PathPatch(path, facecolor=color, lw=0, alpha=0.2)
    edge = patches.PathPatch(path, edgecolor=color, facecolor='none', lw=lw)
    plt.gca().add_patch(patch)
    plt.gca().add_patch(edge)



inc = 1.2
rad = 0.3
lw = 2
colors = ['blue', 'green', 'red', 'orange']

plt.close('all')

##
plt.figure()
for k in range(4):
    points = 1.5*(np.random.rand(30, 2) - 0.5) + k
    plt.plot(points[:,0], points[:,1], 'o', ms=5, color=colors[k])
    convex_poly(points, colors[k], inc=inc, rad=rad, lw=lw)
    
plt.grid(True)
plt.xlabel(r"$x$", size=18)
plt.ylabel(r"$y$", size=18)
  
##  
E = {}
E["poly"] = np.loadtxt('young_poly.txt')
E["metals"] = np.loadtxt('young_metals.txt')
E["comp"] = np.loadtxt('young_comp.txt')
E["ceramic"] = np.loadtxt('young_ceramic.txt')

rho = {}
rho["poly"] = np.loadtxt('dens_poly.txt')
rho["metals"] = np.loadtxt('dens_metals.txt')
rho["comp"] = np.loadtxt('dens_comp.txt')
rho["ceramic"] = np.loadtxt('dens_ceramic.txt')

plt.figure()
for k, key  in enumerate(E.keys()):
    x = rho[key][:,0] * 1000
    y = E[key][:,0] * 1e9
    points = np.vstack([x,y]).T
    convex_poly(points, colors[k], inc=inc, rad=0.3, lw=lw)
    plt.loglog(x, y, 'o', ms=5, color=colors[k])
    
plt.xlabel(r"Density $\rho$ (kg/m$^3$)", size=18)
plt.ylabel(r"Young Modulus $E$ (GPa)", size=18)

plt.grid(True)
plt.show()




