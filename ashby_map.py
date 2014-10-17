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


def poly_enclose(points, color, inc=1.2, rad=0.3, lw=2):
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


def ellip_enclose(points, color, inc=1.2, lw=2, nst=2):
    """
    Plot the minimum ellipse around a set of points.
    
    https://github.com/joferkington/oost_paper_code/blob/master/error_ellipse.py
    """
    
    def eigsorted(cov):
        vals, vecs = np.linalg.eigh(cov)
        order = vals.argsort()[::-1]
        return vals[order], vecs[:,order]
        
    def rot_ellip(w, h, center, angle, n=26):
        angle = angle*np.pi/180
        t = np.linspace(0, 2*np.pi, n-1)
        verts = np.zeros((n,2))
        x = verts[:,0]
        y = verts[:,1]
        x[:-1] = w*np.cos(t)
        x[-1] = x[-2]
        y[:-1] = h*np.sin(t)
        y[-2] = y[-2]
        rot_mat = np.array([[np.cos(theta), -np.sin(theta)],
                             [np.sin(theta), np.cos(theta)]])     
        verts = np.dot(verts, rot_mat.T)
        
        verts = verts + center
              
        codes = [4 for j in range(n)]
        codes[0] = 1
        codes[-1] = 79
        return Path(verts, codes)
        
    
    x = points[:,0]
    y = points[:,1]
    cov = np.cov(x, y)
    vals, vecs = eigsorted(cov)
    theta = np.degrees(np.arctan2(*vecs[:,0][::-1]))
    print theta
    w, h = 2 * nst * np.sqrt(vals)        
    center = np.mean(points, 0)
    ell = patches.Ellipse(center, width=inc*w, height=inc*h, angle=theta,
                          facecolor=color, alpha=0.2, lw=0)
    edge = patches.Ellipse(center, width=inc*w, height=inc*h, angle=theta,
                          facecolor='none', edgecolor=color, lw=lw)
#    path = ell.get_path()
#    print path, len(path.vertices),len(path.codes)
#    plt.plot(path.vertices[:,0], path.vertices[:,1])
#    path = rot_ellip(w, h, center, theta, 26)
#    patch = patches.PathPatch(path, facecolor=color, lw=2, alpha=0.2)
#    plt.gca().add_patch(patch)
    plt.gca().add_artist(ell)
    plt.gca().add_artist(edge)


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
#    plt.loglog(points[:,0], points[:,1], 'o', ms=5, color=colors[k])
#    poly_enclose(points, colors[k], inc=inc, rad=rad, lw=lw)
    ellip_enclose(points, colors[k], inc=inc, lw=lw)
    
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
    poly_enclose(points, colors[k], inc=inc, rad=0.3, lw=lw)
#    ellip_enclose(points, colors[k], inc=inc, lw=lw)
#    plt.loglog(x, y, 'o', ms=5, color=colors[k])
    plt.plot(x, y, 'o', ms=5, color=colors[k])
    
plt.xlabel(r"Density $\rho$ (kg/m$^3$)", size=18)
plt.ylabel(r"Young Modulus $E$ (GPa)", size=18)

plt.grid(True)
plt.show()




