# -*- coding: utf-8 -*-
"""
Created on Mon Nov  9 17:06:32 2020

@author: Dorian
"""

from mpl_toolkits.mplot3d import Axes3D # Required for 3d projection
import numpy as np 
from numpy import cos, sin, sqrt
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, RadioButtons
from PoleFigure import add_polefigure

# Global parameters
scale=0.2                               # Size of fondamental lattice
axcolor = 'lightgoldenrodyellow'        # Color of sliders and so
color_cycle=('red', 'green', 'blue')    # Colors of directions


ranges={'Triclinic':(360,180,360), 'Cubic':(360,90,90), 'Hexagonal':(360,90,60),}

def euler2mat(phi1,Phi,phi2):
    '''Compute rotation matrix from Euler angles'''
    g1=np.array([[cos(phi1), sin(phi1), 0],
                  [-sin(phi1), cos(phi1), 0],
                  [0,0,1]])
    g2=np.array([[1, 0, 0],
                 [0, cos(Phi), sin(Phi)],
                 [0, -sin(Phi), cos(Phi)]])
    g3=np.array([[cos(phi2), sin(phi2), 0],
                  [-sin(phi2), cos(phi2), 0],
                  [0,0,1]])
    return np.matmul(g3, np.matmul(g2, g1))

def cart2sph(x,y,z):
    '''Convert cartesian coordinates (x,y,z) to spherical coordinates 
    (radius, azimuth, colatitude)'''
    r=np.sqrt(x**2+y**2+z**2)
    phi=np.arctan2(y,x)
    theta=np.arccos(z/r)
    return r, phi, theta

## The most simple way to update a 3D arrow is to clear it...
def clear_axis(ax):
    ax.cla()
    ax.set_xlabel('x',labelpad=-15)
    ax.set_ylabel('y',labelpad=-15)
    ax.set_zlabel('z',labelpad=-15)
    ax.set_xticks(())
    ax.set_yticks(())
    ax.set_zticks(())

## Wireframe for cubic crystal
Pts_cube=np.array([[0,0,0],
              [1,0,0],
              [1,1,0],
              [0,1,0],
              [0,0,0],
              [0,0,1],
              [1,0,1],
              [1,1,1],
              [0,1,1],
              [0,0,1],
              [np.nan,0,0],
              [1,0,0],
              [1,0,1],
              [np.nan,0,0],
              [1,1,0],
              [1,1,1],
              [np.nan,0,0],
              [0,1,0],
              [0,1,1]])
    
## Wireframe for hexagonal crystal
theta=np.linspace(0,2*np.pi, 7)
pl1=np.concatenate([ [np.cos(theta)], [np.sin(theta)],[np.zeros(shape=(7))]]).T
pl2=pl1+np.array([0,0,2*np.sqrt(2/3)])
Pts_hex=np.concatenate((pl1,pl2), axis=0)
sep=np.array([[np.nan, 0,0]])
Pts_hex=np.concatenate((Pts_hex, sep, Pts_hex[[1,8],:], sep, Pts_hex[[2,9],:], sep, Pts_hex[[3,10],:], sep, Pts_hex[[4,11],:], sep, Pts_hex[[5,12],:]), axis=0)

## Wireframe for triclinic crystal
def triclin_mat():
    a,b,c = 1, 1.1, 1.5
    alpha=np.deg2rad(80)
    beta=np.deg2rad(85)
    gamma=np.deg2rad(95)
    n=(cos(alpha)-cos(gamma)*cos(beta))/sin(gamma)
    m33=c*sqrt(sin(beta)**2-n**2)
    M=np.array([[a, b*cos(gamma), c*cos(beta)],[0, b*sin(gamma), c*n],[0,0,m33]])
    return M

fig = plt.figure(figsize = [9, 4.5])
ax = fig.add_subplot(121, projection='3d')
plt.subplots_adjust(left=0.25, bottom=0.25)

# Sliders for Euler angles
ax_phi1 = plt.axes([0.3, 0.15, 0.6, 0.03], facecolor=axcolor)
ax_Phi = plt.axes([0.3, 0.1, 0.6, 0.03], facecolor=axcolor)
ax_phi2 = plt.axes([0.3, 0.05, 0.6, 0.03], facecolor=axcolor)
sphi1 = Slider(ax_phi1, '$\phi_1$', 0, 360, valinit=0, valstep=1, valfmt='%0.0f°', closedmax=False)
sPhi = Slider(ax_Phi, '$\Phi$', 0, 180, valinit=0, valstep=1, valfmt='%0.0f°', closedmax=False)
sphi2 = Slider(ax_phi2, '$\phi_2$', 0, 360, valinit=0, valstep=1, valfmt='%0.0f°', closedmax=False)

## Init the pole figure
ax_pf = add_polefigure(fig, 122, projection='stereographic')
scats=[0,0,0]
for i in range(0,3):
    scats[i] = ax_pf.scatter(0, 0, color=color_cycle[i])
ax_pf.set_rlim(0.0, np.pi / 2)
ax_pf.set_theta_zero_location("N")
ax_pf.set_xticks(np.arange(0,2*np.pi, np.pi/2))
ax_pf.set_xticklabels(('x', 'y', '-x', '-y'))
ax_pf.tick_params(axis='x', which='major', pad=-5)

# Radio buttons for crystal symmetry
rax_cs = plt.axes([0.025, 0.05, 0.2, 0.15], facecolor=axcolor, title='Crystal symmetry')
radio_cs = RadioButtons(rax_cs, ('Triclinic','Cubic', 'Hexagonal'), active=0)

# Radio buttons for PF projection
rax_pf = plt.axes([0.025, 0.3, 0.2, 0.15], facecolor=axcolor, title='PF projection')
radio_pf = RadioButtons(rax_pf, ('Stereographic','Lambert'), active=0)


q=[0,0,0]
def show_crystal(val):
    angles=np.deg2rad([sphi1.val, sPhi.val,sphi2.val])
    mat=euler2mat(*angles).T
    uvw=np.eye(3)
    if radio_cs.value_selected == 'Cubic':
        Pts0 = Pts_cube
    elif radio_cs.value_selected == 'Triclinic':
        Pts0 = np.dot(triclin_mat(), Pts_cube.T).T
        uvw = triclin_mat()
    else:
        Pts0 = Pts_hex
    uvw=np.matmul(mat,uvw)
    Pts=np.dot(mat, Pts0.T).T*scale
    dC=np.array([0.5,0.5,0.5])
    Pts += dC  
    clear_axis(ax)
    ax.plot3D(*Pts_cube.T, color='grey')
    ax.plot3D(*Pts.T, color='black', marker='.')
    for i in range(0,uvw.shape[1]):
        uvw_i=uvw[:,i]
        q[i]=ax.quiver3D(*dC,*uvw_i, length=0.5, color=color_cycle[i])
        r, phi, theta = cart2sph(*uvw_i)
        if theta > np.pi/2:
            phi += np.pi
            theta = np.pi - theta
        scats[i].set_offsets((phi,theta))
    if radio_cs.value_selected == 'Hexagonal':
        legend_entries=[r'$[11\bar{2}0]$', r'$[\bar{1}100]$', r'$[0001]$']
    else:
        legend_entries=['[100]', '[010]', '[001]']
    ax.legend(q, legend_entries, bbox_to_anchor=(0, 1), loc='upper right')
    fig.canvas.draw_idle()

show_crystal(0)
             
def switch_geom(geom):
    range=ranges[geom]
    sPhi.valmax = range[1]
    sphi2.valmax = range[2]
    sPhi.set_val(sPhi.val % range[1])
    sphi2.set_val(sphi2.val % range[2])
    show_crystal(0)
    
def switch_proj(proj):
    ax_pf.set_rscale(proj)
    fig.canvas.draw_idle()
    
radio_cs.on_clicked(switch_geom)
radio_pf.on_clicked(switch_proj)       
sphi1.on_changed(show_crystal)
sPhi.on_changed(show_crystal)
sphi2.on_changed(show_crystal)

plt.show()
