"""
Created by Daniel Fu (Reif Lab, Duke University) at 3/29/2020

Name        : custom.py
Project     : cadaxisdna
Description : A collection of functions for custom creation of shapes and modulex
Interpreter : Python 3.7.4
"""

import numpy as np

from ..routing import modules


class Shape:
    def __init__(self):
        raise NotImplementedError("This is only a default constructor.")

    def construct(self):
        raise NotImplementedError("This is only a default constructor.")


class CustomInput:
    def __init__(self, rings, filename):
        self.rings = rings
        self.inputs = {filename: True}

    def construct(self, idx):
        all_modules = []
        all_planes = []
        for line in self.rings:
            bps = line[0]
            height = line[1]
            dirbit = line[2]
            p = modules.Plane(np.array([0, 0, height]), pitch=0, yaw=0)
            (dnaobj, idx) = p.add_circle(idx, bps, dirbit)
            all_planes.append(p)
            all_modules.append(p.modules[0])
            _ = modules.JoinedStrand([p.modules[0]])
        return all_planes, all_modules

class Sphere2L:
    def __init__(self):
        self.inputs = {'ReMushroom': "Yes"}

    def construct(self, idx):
        all_modules = []
        all_planes = []
        p = modules.Plane(np.array([0, 0, 30.6]), pitch=0, yaw=0)
        (dnaobj, idx) = p.add_circle(idx, 72, True)
        all_planes.append(p)
        all_modules.append(p.modules[0])
        _ = modules.JoinedStrand([p.modules[0]])
        p = modules.Plane(np.array([0, 0, 29.75]), pitch=0, yaw=0)
        (dnaobj, idx) = p.add_circle(idx, 120, False)
        all_planes.append(p)
        all_modules.append(p.modules[0])
        _ = modules.JoinedStrand([p.modules[0]])
        p = modules.Plane(np.array([0, 0, 28.5]), pitch=0, yaw=0)
        (dnaobj, idx) = p.add_circle(idx, 160, True)
        all_planes.append(p)
        all_modules.append(p.modules[0])
        _ = modules.JoinedStrand([p.modules[0]])
        p = modules.Plane(np.array([0, 0, 26.9]), pitch=0, yaw=0)
        (dnaobj, idx) = p.add_circle(idx, 200, False)
        all_planes.append(p)
        all_modules.append(p.modules[0])
        _ = modules.JoinedStrand([p.modules[0]])
        p = modules.Plane(np.array([0, 0, 24.95]), pitch=0, yaw=0)
        (dnaobj, idx) = p.add_circle(idx, 232, True)
        all_planes.append(p)
        all_modules.append(p.modules[0])
        _ = modules.JoinedStrand([p.modules[0]])
        p = modules.Plane(np.array([0, 0, 22.75]), pitch=0, yaw=0)
        (dnaobj, idx) = p.add_circle(idx, 260, False)
        all_planes.append(p)
        all_modules.append(p.modules[0])
        _ = modules.JoinedStrand([p.modules[0]])
        p = modules.Plane(np.array([0, 0, 20.35]), pitch=0, yaw=0)
        (dnaobj, idx) = p.add_circle(idx, 280, True)
        all_planes.append(p)
        all_modules.append(p.modules[0])
        _ = modules.JoinedStrand([p.modules[0]])
        p = modules.Plane(np.array([0, 0, 17.8]), pitch=0, yaw=0)
        (dnaobj, idx) = p.add_circle(idx, 288, False)
        all_planes.append(p)
        all_modules.append(p.modules[0])
        _ = modules.JoinedStrand([p.modules[0]])
        p = modules.Plane(np.array([0, 0, 15.2]), pitch=0, yaw=0)
        (dnaobj, idx) = p.add_circle(idx, 292, True)
        all_planes.append(p)
        all_modules.append(p.modules[0])
        _ = modules.JoinedStrand([p.modules[0]])
        p = modules.Plane(np.array([0, 0, 12.6]), pitch=0, yaw=0)
        (dnaobj, idx) = p.add_circle(idx, 288, False)
        all_planes.append(p)
        all_modules.append(p.modules[0])
        _ = modules.JoinedStrand([p.modules[0]])
        p = modules.Plane(np.array([0, 0, 10.05]), pitch=0, yaw=0)
        (dnaobj, idx) = p.add_circle(idx, 276, True)
        all_planes.append(p)
        all_modules.append(p.modules[0])
        _ = modules.JoinedStrand([p.modules[0]])
        p = modules.Plane(np.array([0, 0, 7.7]), pitch=0, yaw=0)
        (dnaobj, idx) = p.add_circle(idx, 252, False)
        all_planes.append(p)
        all_modules.append(p.modules[0])
        _ = modules.JoinedStrand([p.modules[0]])
        p = modules.Plane(np.array([0, 0, 5.55]), pitch=0, yaw=0)
        (dnaobj, idx) = p.add_circle(idx, 224, True)
        all_planes.append(p)
        all_modules.append(p.modules[0])
        _ = modules.JoinedStrand([p.modules[0]])
        p = modules.Plane(np.array([0, 0, 3.7]), pitch=0, yaw=0)
        (dnaobj, idx) = p.add_circle(idx, 192, False)
        all_planes.append(p)
        all_modules.append(p.modules[0])
        _ = modules.JoinedStrand([p.modules[0]])
        p = modules.Plane(np.array([0, 0, 2.15]), pitch=0, yaw=0)
        (dnaobj, idx) = p.add_circle(idx, 152, True)
        all_planes.append(p)
        all_modules.append(p.modules[0])
        _ = modules.JoinedStrand([p.modules[0]])
        p = modules.Plane(np.array([0, 0, 1]), pitch=0, yaw=0)
        (dnaobj, idx) = p.add_circle(idx, 108, False)
        all_planes.append(p)
        all_modules.append(p.modules[0])
        _ = modules.JoinedStrand([p.modules[0]])
        p = modules.Plane(np.array([0, 0, 30.6]), pitch=0, yaw=0)
        (dnaobj, idx) = p.add_circle(idx, 120, False)
        all_planes.append(p)
        all_modules.append(p.modules[0])
        _ = modules.JoinedStrand([p.modules[0]])
        p = modules.Plane(np.array([0, 0, 29.75]), pitch=0, yaw=0)
        (dnaobj, idx) = p.add_circle(idx, 168, True)
        all_planes.append(p)
        all_modules.append(p.modules[0])
        _ = modules.JoinedStrand([p.modules[0]])
        p = modules.Plane(np.array([0, 0, 28.5]), pitch=0, yaw=0)
        (dnaobj, idx) = p.add_circle(idx, 208, False)
        all_planes.append(p)
        all_modules.append(p.modules[0])
        _ = modules.JoinedStrand([p.modules[0]])
        p = modules.Plane(np.array([0, 0, 26.9]), pitch=0, yaw=0)
        (dnaobj, idx) = p.add_circle(idx, 248, True)
        all_planes.append(p)
        all_modules.append(p.modules[0])
        _ = modules.JoinedStrand([p.modules[0]])
        p = modules.Plane(np.array([0, 0, 24.95]), pitch=0, yaw=0)
        (dnaobj, idx) = p.add_circle(idx, 280, False)
        all_planes.append(p)
        all_modules.append(p.modules[0])
        _ = modules.JoinedStrand([p.modules[0]])
        p = modules.Plane(np.array([0, 0, 22.75]), pitch=0, yaw=0)
        (dnaobj, idx) = p.add_circle(idx, 308, True)
        all_planes.append(p)
        all_modules.append(p.modules[0])
        _ = modules.JoinedStrand([p.modules[0]])
        p = modules.Plane(np.array([0, 0, 20.35]), pitch=0, yaw=0)
        (dnaobj, idx) = p.add_circle(idx, 328, False)
        all_planes.append(p)
        all_modules.append(p.modules[0])
        _ = modules.JoinedStrand([p.modules[0]])
        p = modules.Plane(np.array([0, 0, 17.8]), pitch=0, yaw=0)
        (dnaobj, idx) = p.add_circle(idx, 336, True)
        all_planes.append(p)
        all_modules.append(p.modules[0])
        _ = modules.JoinedStrand([p.modules[0]])
        p = modules.Plane(np.array([0, 0, 15.2]), pitch=0, yaw=0)
        (dnaobj, idx) = p.add_circle(idx, 340, False)
        all_planes.append(p)
        all_modules.append(p.modules[0])
        _ = modules.JoinedStrand([p.modules[0]])
        p = modules.Plane(np.array([0, 0, 12.6]), pitch=0, yaw=0)
        (dnaobj, idx) = p.add_circle(idx, 336, True)
        all_planes.append(p)
        all_modules.append(p.modules[0])
        _ = modules.JoinedStrand([p.modules[0]])
        p = modules.Plane(np.array([0, 0, 10.05]), pitch=0, yaw=0)
        (dnaobj, idx) = p.add_circle(idx, 324, False)
        all_planes.append(p)
        all_modules.append(p.modules[0])
        _ = modules.JoinedStrand([p.modules[0]])
        p = modules.Plane(np.array([0, 0, 7.7]), pitch=0, yaw=0)
        (dnaobj, idx) = p.add_circle(idx, 300, True)
        all_planes.append(p)
        all_modules.append(p.modules[0])
        _ = modules.JoinedStrand([p.modules[0]])
        p = modules.Plane(np.array([0, 0, 5.55]), pitch=0, yaw=0)
        (dnaobj, idx) = p.add_circle(idx, 272, False)
        all_planes.append(p)
        all_modules.append(p.modules[0])
        _ = modules.JoinedStrand([p.modules[0]])
        p = modules.Plane(np.array([0, 0, 3.7]), pitch=0, yaw=0)
        (dnaobj, idx) = p.add_circle(idx, 240, True)
        all_planes.append(p)
        all_modules.append(p.modules[0])
        _ = modules.JoinedStrand([p.modules[0]])
        p = modules.Plane(np.array([0, 0, 2.15]), pitch=0, yaw=0)
        (dnaobj, idx) = p.add_circle(idx, 200, False)
        all_planes.append(p)
        all_modules.append(p.modules[0])
        _ = modules.JoinedStrand([p.modules[0]])
        p = modules.Plane(np.array([0, 0, 1]), pitch=0, yaw=0)
        (dnaobj, idx) = p.add_circle(idx, 156, True)
        all_planes.append(p)
        all_modules.append(p.modules[0])
        _ = modules.JoinedStrand([p.modules[0]])

        return all_planes, all_modules
