# -*- coding: utf-8 -*-
"""
Created on Thu Sep 19 09:36:34 2019

@author: Dan
"""

import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import io
import base64
import numpy as np


# https://www.lifewire.com/contrasting-foreground-background-colors-4061363
# Good color contrast is:
#     Gold nodes #F1C40F
#     Orange nodes border #F39C12
#     Black text # 000000
#     White background # FFFFFF
#     Darker gold connection lines #D4AC0D
#     "Python blue" is #1F77B4


def graph_stl(mesh, figsize=(5, 5), view=(0, 0, 'z')):
    # Create a new plot
    figure = plt.figure(figsize=figsize)
    axes = mplot3d.Axes3D(figure, rect=[-0.27, -0.35, 1.5, 1.5])


    # Load the STL files and add the vectors to the plot
    axes.add_collection3d(mplot3d.art3d.Poly3DCollection(mesh.vectors))

    # Auto scale to the mesh size
    xx = []
    yy = []
    zz = []
    for f in mesh.points:
        v0 = f[0:3]
        v1 = f[3:6]
        v2 = f[6:9]
        for v in [v0, v1, v2]:
            xx.append(v[0])
            yy.append(v[1])
            zz.append(v[2])

    minlim = min(min(xx), min(yy), min(zz))
    maxlim = max(max(xx), max(yy), max(zz))
    axlim = max(abs(minlim), abs(maxlim))
    # print("limits=", axlim)
    plt.xlim([-axlim, axlim])
    plt.ylim([-axlim, axlim])
    axes.set_zlim([-axlim, axlim])
    plt.subplots_adjust(left=0.5)

    # toggle grid and axis from display
    plt.axis('off')

    # set viewing angle
    axes.view_init(elev=view[0], azim=view[1])

    # Draw guides
    # linex = [0, 0]
    # liney = [0, 0]
    # linez = [-10, 10]
    # axes.plot(linex, liney, linez)  # vertical guide
    #
    # linex = [0, 0]
    # liney = [-10, 10]
    # linez = [0, 0]
    # axes.plot(linex, liney, linez)  # horizontal guide

    # save the data
    img = io.BytesIO()
    plt.savefig(img, format='png')
    img.seek(0)
    graph_url = base64.b64encode(img.getvalue()).decode()
    plt.close()

    return '{}'.format(graph_url)


def build_graph(pts=None, halftrace=None, fulltrace=None, axis=False, figsize=(2, 2)):
    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111)
    img = io.BytesIO()
    if pts is not None:
        ptsx = pts[:, 0]
        ptsy = pts[:, 1]
        plt.scatter(ptsx, ptsy, c='#1f77b4')
    if axis:
        margin = .3
        ymax = max(ptsy) + margin
        ymin = min(ptsy) - margin
        ax.plot([0, 0], [ymin, ymax], 'k:', linewidth=5)
    if fulltrace is not None:
        lsarr = np.array(list(fulltrace.coords))
        lsx = lsarr[:, 0]
        lsy = lsarr[:, 1]
        ax.plot(lsx, lsy, 'r', linewidth=5, color="#00FF00")
    if halftrace is not None:
        lsarr = np.array(list(halftrace.coords))
        lsx = lsarr[:, 0]
        lsy = lsarr[:, 1]
        ax.plot(lsx, lsy, 'r', linewidth=5, color="#FF0000")
    plt.savefig(img, format='png')
    img.seek(0)
    graph_url = base64.b64encode(img.getvalue()).decode()
    plt.close()
    return 'data:image/png;base64,{}'.format(graph_url)
