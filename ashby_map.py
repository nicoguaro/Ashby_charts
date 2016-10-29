# -*- coding: utf-8 -*-
"""
Plot Ashby Charts.

@author: Nicolas Guarin Zapata
@date: October, 2014
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
    points = np.log(points)
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
    patch._path._vertices = np.exp(patch._path._vertices)
    return patch, edge


def ellip_enclose(points, color, inc=1.2, lw=2, nst=2):
    """
    Plot the an ellipse around a set of points using Principal
	Component Analysis [1, 2].
    
    	
	References
	==========
	[1] Jofer Kington https://github.com/joferkington/oost_paper_code/blob/master/error_ellipse.py
	
	[2] Wikipedia contributors. "Principal component analysis." 
	Wikipedia, The Free Encyclopedia. Wikipedia, The Free Encyclopedia,
	29 Oct. 2014. Web. 29 Oct. 2014. 
    """
    
    def eigsorted(cov):
        vals, vecs = np.linalg.eigh(cov)
        order = vals.argsort()[::-1]
        return vals[order], vecs[:, order]
        
    points = np.log(points)
    x = points[:, 0]
    y = points[:, 1]
    cov = np.cov(x, y)
    vals, vecs = eigsorted(cov)
    theta = np.degrees(np.arctan2(*vecs[:,0][::-1]))
    w, h = 2 * nst * np.sqrt(vals)        
    center = np.mean(points, 0)
    patch = patches.Ellipse(center, width=inc*w, height=inc*h, angle=theta,
                          facecolor=color, alpha=0.2, lw=0)
    edge = patches.Ellipse(center, width=inc*w, height=inc*h, angle=theta,
                          facecolor='none', edgecolor=color, lw=lw)
    edge._path._vertices = np.exp(edge._path._vertices)
    return patch, edge


def plot_group(points, color, inc=1.2, lw=2, rad=0.3, nst=2, enclosing="poly",
               scaling="log"):
    plt.plot(points[:,0], points[:,1], 'o', ms=8, color=color,
             mfc="white", mec=color)
    if enclosing=="poly":
        patch, edge = poly_enclose(points, color, inc=inc,
                                   rad=rad, lw=lw)
    elif enclosing=="ellipse":
        patch, edge = ellip_enclose(points, color, inc=1,
                                    lw=lw, nst=nst)
    else:
        raise Exception("Not know enclosing option.")

#    plt.xscale(scaling)
#    plt.yscale(scaling)
    plt.gca().add_patch(patch)
    plt.gca().add_patch(edge)

if __name__ == "__main__":
    inc = 1.5
    rad = 0.3
    lw = 2
    colors = ['blue', 'green', 'red', 'orange']
    

    E = {}
    E["poly"] = np.loadtxt('props/young_poly.txt')
    E["metals"] = np.loadtxt('props/young_metals.txt')
    E["comp"] = np.loadtxt('props/young_comp.txt')
    E["ceramic"] = np.loadtxt('props/young_ceramic.txt')
    
    rho = {}
    rho["poly"] = np.loadtxt('props/dens_poly.txt')
    rho["metals"] = np.loadtxt('props/dens_metals.txt')
    rho["comp"] = np.loadtxt('props/dens_comp.txt')
    rho["ceramic"] = np.loadtxt('props/dens_ceramic.txt')
    
    plt.figure()
    for k, key  in enumerate(E.keys()):
        x = rho[key][:,0] * 1000  # Units to be in kg*m**-3
        y = E[key][:,0] *1e9  # Units to be in Pa
        points = np.vstack([x,y]).T
        plot_group(points, colors[k], enclosing="ellipse", inc=inc, rad=rad)
    
    plt.xlabel(r"Density $\rho$ (kg/m$^3$)", size=18)
    plt.ylabel(r"Young Modulus $E$ (Pa)", size=18)
    plt.grid(True)
    
    
    #plt.savefig("log-example.svg")
    #plt.savefig("log-example.pdf")
    #plt.savefig("log-example.png", dpi=300)
    
    plt.show()




