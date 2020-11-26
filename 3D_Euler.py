# -*- coding: utf-8 -*-
"""
Created on Mon Nov  9 17:06:32 2020

@author: Dorian
"""

from mpl_toolkits.mplot3d import Axes3D # Required for 3d projection
import numpy as np 
from numpy import cos, sin, sqrt
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, RadioButtons, CheckButtons
from PoleFigure import add_polefigure

# Global parameters
scale=0.2                               # Size of fondamental lattice
axcolor = 'lightgoldenrodyellow'        # Color of sliders and so
color_cycle=('red', 'green', 'blue')    # Colors of directions


ranges={'Triclinic':(360,180,360), 'Cubic':(360,90,90), 'Hexagonal':(360,90,60),}

def euler2mat(phi1,Phi,phi2):
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
    r=x**2+y**2+z**2
    phi=np.arctan2(y,x)
    theta=np.arccos(z/r)
    return r, phi, theta

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

## The most simple way to update a 3D arrow is to clear it...
def clear_axis():
    ax.cla()
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    ax.set_xticks(())
    ax.set_yticks(())
    ax.set_zticks(())


# Sliders for Euler angles
ax_phi1 = plt.axes([0.3, 0.15, 0.6, 0.03], facecolor=axcolor)
ax_Phi = plt.axes([0.3, 0.1, 0.6, 0.03], facecolor=axcolor)
ax_phi2 = plt.axes([0.3, 0.05, 0.6, 0.03], facecolor=axcolor)
sphi1 = Slider(ax_phi1, '$\phi_1$', 0, 359, valinit=0, valstep=1, valfmt='%0.0f°')
sPhi = Slider(ax_Phi, '$\Phi$', 0, 179, valinit=0, valstep=1, valfmt='%0.0f°')
sphi2 = Slider(ax_phi2, '$\phi_2$', 0, 359, valinit=0, valstep=1, valfmt='%0.0f°')

## Init the pole figure
ax2 = add_polefigure(fig, 122, projection='stereographic')
scats=[0,0,0]
for i in range(0,3):
    scats[i] = ax2.scatter(0, 0, color=color_cycle[i])
ax2.set_rlim(0.0, np.pi / 2)
ax2.set_theta_zero_location("N")
ax2.set_xticks(np.arange(0,2*np.pi, np.pi/2))
ax2.set_xticklabels(('x', 'y', '-x', '-y'))

# Radio buttons for crystal symmetry
rax = plt.axes([0.025, 0.05, 0.2, 0.15], facecolor=axcolor)
radio = RadioButtons(rax, ('Triclinic','Cubic', 'Hexagonal'), active=0)


q=[0,0,0]
def show_crystal(val):
    angles=np.deg2rad([sphi1.val, sPhi.val,sphi2.val])
    mat=euler2mat(*angles).T
    uvw=np.eye(3)
    if radio.value_selected == 'Cubic':
        Pts0 = Pts_cube
    elif radio.value_selected == 'Triclinic':
        Pts0 = np.dot(triclin_mat(), Pts_cube.T).T
        uvw = triclin_mat()
    else:
        Pts0 = Pts_hex
    uvw=np.matmul(mat,uvw)
    Pts=np.dot(mat, Pts0.T).T*scale
    dC=np.array([0.5,0.5,0.5])
    Pts += dC  
    clear_axis()
    ax.plot3D(*Pts_cube.T, color='grey')
    ax.plot3D(*Pts.T, color='black', marker='.')
    for i in range(0,3):
        uvw_i=uvw[:,i]
        q[i]=ax.quiver3D(*dC,*uvw_i, length=0.5, color=color_cycle[i])
        r, phi, theta = cart2sph(*uvw_i)
        if theta > np.pi/2:
            phi += np.pi
            theta = np.pi - theta
        scats[i].set_offsets((phi,theta))
    if radio.value_selected == 'Hexagonal':
        legend_entries=[r'$[11\bar{2}0]$', r'$[1\bar{1}00]$', r'$[0001]$']
    else:
        legend_entries=['[100]', '[010]', '[001]']
    ax.legend(q, legend_entries, bbox_to_anchor=(0, 1), loc='upper right')
    fig.canvas.draw_idle()

show_crystal(0)
             
def switch_geom(geom):
    range=ranges[geom]
    sPhi.valmax = range[1]-1
    sphi2.valmax = range[2]-1  
    sPhi.set_val(sPhi.val % range[1])
    sphi2.set_val(sphi2.val % range[2])
    show_crystal(0)
    
radio.on_clicked(switch_geom)        
sphi1.on_changed(show_crystal)
sPhi.on_changed(show_crystal)
sphi2.on_changed(show_crystal)

plt.show()
