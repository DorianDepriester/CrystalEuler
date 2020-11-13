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
a,b,c = 1, 1.1, 1.5
alpha=np.deg2rad(70)
beta=np.deg2rad(75)
gamma=np.deg2rad(80)
n=(cos(alpha)-cos(gamma)*cos(beta))/sin(gamma)
Omega=a*b*c*sqrt(1-cos(alpha)**2-cos(beta)**2-cos(gamma)**2+2*cos(alpha)*cos(beta)*cos(gamma))
M=np.array([[a, b*cos(gamma), c*cos(beta)],[0, b*sin(gamma), c*n],[0,0,Omega/(a*b*sin(gamma))]])
Pts_triclin = np.dot(M, Pts_cube.T).T

def plotCube(geom='cube', scale=1, centered=False, color='grey', marker='.'):
    if geom=='cube':
        Pts0 = Pts_cube
    else:
        Pts0 = Pts_triclin
    Pts=Pts0*scale
    if centered:
        C= np.mean(Pts[np.all(np.isfinite(Pts), axis=1),:],axis=0)
        Pts += np.array([0.5,0.5,0.5])-C
    ax.plot3D(*Pts.T, color, marker=marker)

fig, ax = plt.subplots()
plt.subplots_adjust(left=0.25, bottom=0.25)
ax = plt.axes(projection='3d')

def clear_axis():
    ax.cla()
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    ax.set_xticks(())
    ax.set_yticks(())
    ax.set_zticks(())
    
clear_axis()
scale=0.2
axcolor = 'lightgoldenrodyellow'
ax_phi1 = plt.axes([0.25, 0.2, 0.65, 0.03], facecolor=axcolor)
ax_Phi = plt.axes([0.25, 0.15, 0.65, 0.03], facecolor=axcolor)
ax_phi2 = plt.axes([0.25, 0.1, 0.65, 0.03], facecolor=axcolor)
sphi1 = Slider(ax_phi1, '$\phi_1$', 0, 359, valinit=0, valstep=1, valfmt='%0.0f°')
sPhi = Slider(ax_Phi, '$\Phi$', 0, 179, valinit=0, valstep=1, valfmt='%0.0f°')
sphi2 = Slider(ax_phi2, '$\phi_2$', 0, 259, valinit=0, valstep=1, valfmt='%0.0f°')

rax = plt.axes([0.025, 0.5, 0.2, 0.15], facecolor=axcolor)
radio = RadioButtons(rax, ('Triclinic','Cubic', 'Hexagonal'), active=0)

def show_crystal(val):
    angles=np.deg2rad([sphi1.val, sPhi.val,sphi2.val])
    mat=euler2mat(*angles).T
    if radio.value_selected == 'Cubic':
        Pts0 = Pts_cube
    elif radio.value_selected == 'Triclinic':
        Pts0 = Pts_triclin
    else:
        Pts0 = Pts_hex
    Pts=np.dot(mat, Pts0.T).T*scale
    dC=np.array([0.5,0.5,0.5])
    Pts += dC  
    clear_axis()
    plotCube()
    ax.plot3D(*Pts.T, color='black', marker='.')
    u=np.dot(mat, np.array([1,0,0]))
    v=np.dot(mat, np.array([0,1,0]))
    w=np.dot(mat, np.array([0,0,1]))    
    q1 =ax.quiver3D(*dC,*u, length=0.5, color='red')
    q2 =ax.quiver3D(*dC,*v, length=0.5, color='green')
    q3 =ax.quiver3D(*dC,*w, length=0.5, color='blue')
    if radio.value_selected == 'Hexagonal':
        legend_entries=[r'$[11\bar{2}0]$', r'$[1\bar{1}00]$', r'$[0001]$']
    else:
        legend_entries=['[100]', '[010]', '[001]']
    rax.legend([q1, q2, q3], legend_entries, bbox_to_anchor=(0, 0), loc='upper left')
    fig.canvas.draw_idle()

show_crystal(0)

def switch_geom(geom):
    if geom == 'Cubic':
        sPhi.valmax = 89
        sphi2.valmax = 89        
    elif geom == 'Hexagonal':
        sPhi.valmax = 89
        sphi2.valmax = 59
    else:
        sPhi.valmax = 179
        sphi2.valmax = 359
    sPhi.set_val(sPhi.val % (sPhi.valmax + 1))
    sphi2.set_val(sphi2.val % (sphi2.valmax +1))
    show_crystal(1)
    fig.canvas.draw_idle()
    
radio.on_clicked(switch_geom)        
sphi1.on_changed(show_crystal)
sPhi.on_changed(show_crystal)
sphi2.on_changed(show_crystal)
plt.show()
