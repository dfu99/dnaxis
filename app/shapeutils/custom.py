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


class Flat:
    def __init__(self, name):
        self.inputs = {'Name': name}

    def construct(self, idx):
        all_modules = []
        all_planes = []

        input = [[0, 0, 315, True, 85.75, 0, 0],
                 [0, 2.6, 315, False, -85.75, 1, 1],
                 [0, 5.2, 315, True, 85.75, 2, 2],
                 [0, 7.8, 315, False, -85.75, 3, 3],
                 [0, 10.4, 315, True, 85.75, 4, 4],
                 [0, 13, 315, False, -85.75, 5, 5],
                 [0, 15.6, 315, True, 85.75, 6, 6],
                 [0, 18.2, 315, False, -85.75, 7, 7],
                 [0, 20.8, 315, True, 85.75, 8, 8],
                 [0, 23.4, 315, False, -85.75, 9, 9],
                 [0, 26, 315, True, 85.75, 10, 10],
                 [0, 28.6, 315, False, -85.75, 11, 11],
                 [0, 31.2, 315, True, 85.75, 12, 12],
                 [0, 33.8, 315, False, -85.75, 13, 13],
                 [0, 36.4, 315, True, 85.75, 14, 14],
                 [0, 39, 315, False, -85.75, 15, 15],
                 [0, 41.6, 315, True, 85.75, 16, 16],
                 [0, 44.2, 315, False, -85.75, 17, 17],
                 [0, 46.8, 315, True, 85.75, 18, 18],
                 [0, 49.4, 315, False, -85.75, 19, 19],
                 [0, 52, 315, True, 85.75, 20, 20],
                 [0, 54.6, 315, False, -85.75, 21, 21]]

        for line in input:
            z = line[1]
            x = line[0]
            bp = line[2]
            dirBit = line[3]
            offset = line[4]
            p = modules.Plane(np.array([0, 0, z]), pitch=0, yaw=0)
            (_, idx) = p.add_line(idx, np.array([x, 0, z]), np.array([0, 1, 0]), bp, dirBit, offset=offset)
            all_planes.append(p)
            for m in p.modules:
                all_modules.append(m)
            _ = modules.JoinedStrand(p.modules)

        return all_planes, all_modules


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


class FromFile:
    def __init__(self, filename):
        from app.routing.filehandler import symfile
        self.inputs = {"Name": filename}
        data = symfile.readfile(filename)
        self.data = symfile.list2DToFloat(data)
        for line in self.data:
            line[0] = int(line[0])

        print(self.data)

    def construct(self, idx):
        all_modules = []
        all_planes = []

        for line in self.data:
            bps = line[0]
            z = line[1]
            dirBit = line[2]
            p = modules.Plane(np.array([0, 0, z]), pitch=0, yaw=0)
            (_, idx) = p.add_circle(idx, bps, dirBit)
            all_planes.append(p)
            for m in p.modules:
                all_modules.append(m)
            _ = modules.JoinedStrand(p.modules)

        return all_planes, all_modules


class ManualInput:
    def __init__(self, name):
        self.inputs = {'Name': name}

    def construct(self, idx):
        all_modules = []
        all_planes = []

        input = [[2.25167, 0, 168, False, 38.5, 1, 0],
                 [6.755, 0, 168, False, 38.5, 3, 1],
                 [11.25833, 0, 168, False, 38.5, 5, 2],
                 [0, 1.3, 168, True, -30, 0, 3],
                 [4.50333, 1.3, 168, True, -30, 2, 4],
                 [9.00666, 1.3, 168, True, -30, 4, 5],
                 [13.50999, 1.3, 168, True, 90, 6, 6],
                 [13.50999, 3.9, 168, False, -81.5, 7, 7],
                 [9.00666, 3.9, 168, False, 158.5, 9, 8],
                 [4.50333, 3.9, 168, False, 158.5, 11, 9],
                 [0, 3.9, 168, False, 158.5, 13, 10],
                 [11.25833, 5.2, 168, True, 210, 8, 11],
                 [6.755, 5.2, 168, True, 210, 10, 12],
                 [2.25167, 5.2, 168, True, 210, 12, 13],
                 [-2.25167, 5.2, 168, True, 90, 14, 14],
                 [-2.25167, 7.8, 168, False, -81.5, 15, 15],
                 [2.25166, 7.8, 168, False, 38.5, 17, 16],
                 [6.75499, 7.8, 168, False, 38.5, 19, 17],
                 [11.25832, 7.8, 168, False, 38.5, 21, 18],
                 [0, 9.1, 168, True, -30, 16, 19],
                 [4.50333, 9.1, 168, True, -30, 18, 20],
                 [9.00666, 9.1, 168, True, -30, 20, 21]]

        for line in input:
            z = line[1]
            x = line[0]
            bp = line[2]
            dirBit = line[3]
            offset = line[4]
            p = modules.Plane(np.array([0, 0, z]), pitch=0, yaw=0)
            (_, idx) = p.add_line(idx, np.array([x, 0, z]), np.array([0, 1, 0]), bp, dirBit, offset=offset)
            all_planes.append(p)
            for m in p.modules:
                all_modules.append(m)
            _ = modules.JoinedStrand(p.modules)

        return all_planes, all_modules


class EllipseStraight:
    def __init__(self):
        self.inputs = {'Name': "EllipseStraight"}

    def construct(self, idx):
        all_modules = []
        all_planes = []

        # (0, 0)
        p = modules.Plane(np.array([0, 0, 0]), pitch=0, yaw=0)
        (_, idx) = p.add_line(idx, np.array([1.3, 0, 0]), np.array([0, 1, 0]), 168, True)
        (_, idx) = p.add_arc(idx, 198, True, 0, 180, np.array([-22.192584 + 2.6, 55.776]))
        (_, idx) = p.add_arc(idx, 198, True, 180, 180, np.array([-22.192584 + 2.6, 0]))
        (_, idx) = p.add_line(idx, np.array([1.3 - 44.385168, 0, 0]), np.array([0, 1, 0]), 168, True)
        all_planes.append(p)
        for m in p.modules:
            all_modules.append(m)
        _ = modules.JoinedStrand(p.modules)

        # (0, 1)
        p = modules.Plane(np.array([0, 0, 0]), pitch=0, yaw=0)
        (_, idx) = p.add_line(idx, np.array([3.9, 0, 0]), np.array([0, 1, 0]), 168, False)
        (_, idx) = p.add_arc(idx, 222, False, 0, 180, np.array([-22.192584 + 2.6, 55.776]))
        (_, idx) = p.add_arc(idx, 222, False, 180, 180, np.array([-22.192584 + 2.6, 0]))
        (_, idx) = p.add_line(idx, np.array([3.9 - 44.385168, 0, 0]), np.array([0, 1, 0]), 168, False)
        all_planes.append(p)
        for m in p.modules:
            all_modules.append(m)
        _ = modules.JoinedStrand(p.modules)

        # (0, 2)
        p = modules.Plane(np.array([0, 0, 2.25167]), pitch=0, yaw=0)
        (_, idx) = p.add_line(idx, np.array([5.2, 0, 2.25167]), np.array([0, 1, 0]), 168, True)
        (_, idx) = p.add_arc(idx, 235, True, 0, 180, np.array([-22.192584 + 2.6, 55.776]))
        (_, idx) = p.add_arc(idx, 235, True, 180, 180, np.array([-22.192584 + 2.6, 0]))
        (_, idx) = p.add_line(idx, np.array([5.2 - 44.385168, 0, 2.25167]), np.array([0, 1, 0]), 168, True)
        all_planes.append(p)
        for m in p.modules:
            all_modules.append(m)
        _ = modules.JoinedStrand(p.modules)

        # (0, 3)
        p = modules.Plane(np.array([0, 0, 4.50333]), pitch=0, yaw=0)
        (_, idx) = p.add_line(idx, np.array([3.9, 0, 4.50333]), np.array([0, 1, 0]), 168, False)
        (_, idx) = p.add_arc(idx, 222, False, 0, 180, np.array([-22.192584 + 2.6, 55.776]))
        (_, idx) = p.add_arc(idx, 222, False, 180, 180, np.array([-22.192584 + 2.6, 0]))
        (_, idx) = p.add_line(idx, np.array([3.9 - 44.385168, 0, 4.50333]), np.array([0, 1, 0]), 168, False)
        all_planes.append(p)
        for m in p.modules:
            all_modules.append(m)
        _ = modules.JoinedStrand(p.modules)

        # (0, 4)
        p = modules.Plane(np.array([0, 0, 4.50333]), pitch=0, yaw=0)
        (_, idx) = p.add_line(idx, np.array([1.3, 0, 4.50333]), np.array([0, 1, 0]), 168, True)
        (_, idx) = p.add_arc(idx, 198, True, 0, 180, np.array([-22.192584 + 2.6, 55.776]))
        (_, idx) = p.add_arc(idx, 198, True, 180, 180, np.array([-22.192584 + 2.6, 0]))
        (_, idx) = p.add_line(idx, np.array([1.3 - 44.385168, 0, 4.50333]), np.array([0, 1, 0]), 168, True)
        all_planes.append(p)
        for m in p.modules:
            all_modules.append(m)
        _ = modules.JoinedStrand(p.modules)

        # (0, 5)
        p = modules.Plane(np.array([0, 0, 2.25167]), pitch=0, yaw=0)
        (_, idx) = p.add_line(idx, np.array([0, 0, 2.25167]), np.array([0, 1, 0]), 168, False)
        (_, idx) = p.add_arc(idx, 185, False, 0, 180, np.array([-22.192584 + 2.6, 55.776]))
        (_, idx) = p.add_arc(idx, 185, False, 180, 180, np.array([-22.192584 + 2.6, 0]))
        (_, idx) = p.add_line(idx, np.array([-44.385168, 0, 2.25167]), np.array([0, 1, 0]), 168, False)
        all_planes.append(p)
        for m in p.modules:
            all_modules.append(m)
        _ = modules.JoinedStrand(p.modules)

        return all_planes, all_modules


class Seat:
    def __init__(self):
        self.inputs = {'Name': "Seat"}

    def construct(self, idx):
        all_modules = []
        all_planes = []

        p = modules.Plane(np.array([0, 0, 0]), pitch=0, yaw=0)
        (_, idx) = p.add_arc(idx, 140, True, 0, 180, np.array([0, 0]))
        all_planes.append(p)
        for m in p.modules:
            all_modules.append(m)
        _ = modules.JoinedStrand(p.modules)

        p = modules.Plane(np.array([0, 0, 0]), pitch=0, yaw=0)
        (_, idx) = p.add_arc(idx, 116, False, 0, 180, np.array([0, 0]))
        all_planes.append(p)
        for m in p.modules:
            all_modules.append(m)
        _ = modules.JoinedStrand(p.modules)

        p = modules.Plane(np.array([0, 0, 2.6]), pitch=0, yaw=0)
        (_, idx) = p.add_arc(idx, 140, False, 0, 180, np.array([0, 0]))
        all_planes.append(p)
        for m in p.modules:
            all_modules.append(m)
        _ = modules.JoinedStrand(p.modules)

        p = modules.Plane(np.array([0, 0, 2.6]), pitch=0, yaw=0)
        (_, idx) = p.add_arc(idx, 116, True, 0, 180, np.array([0, 0]))
        all_planes.append(p)
        for m in p.modules:
            all_modules.append(m)
        _ = modules.JoinedStrand(p.modules)

        p = modules.Plane(np.array([-1.3 + 13.529450, 0, -13.529450]), pitch=-90, yaw=90)
        (_, idx) = p.add_arc(idx, 78, True, 0, 90, np.array([0, 0]))
        all_planes.append(p)
        for m in p.modules:
            all_modules.append(m)
        _ = modules.JoinedStrand(p.modules)

        p = modules.Plane(np.array([-1.3 + 13.529450, 0, -13.529450]), pitch=-90, yaw=90)
        (_, idx) = p.add_arc(idx, 66, False, 0, 90, np.array([0, 0]))
        all_planes.append(p)
        for m in p.modules:
            all_modules.append(m)
        _ = modules.JoinedStrand(p.modules)

        p = modules.Plane(np.array([1.3 + 13.529450, 0, -13.529450]), pitch=-90, yaw=90)
        (_, idx) = p.add_arc(idx, 78, False, 0, 90, np.array([0, 0]))
        all_planes.append(p)
        for m in p.modules:
            all_modules.append(m)
        _ = modules.JoinedStrand(p.modules)

        p = modules.Plane(np.array([1.3 + 13.529450, 0, -13.529450]), pitch=-90, yaw=90)
        (_, idx) = p.add_arc(idx, 66, True, 0, 90, np.array([0, 0]))
        all_planes.append(p)
        for m in p.modules:
            all_modules.append(m)
        _ = modules.JoinedStrand(p.modules)

        p = modules.Plane(np.array([-1.3 - 13.529450, 0, -13.529450]), pitch=-90, yaw=90)
        (_, idx) = p.add_arc(idx, 78, True, 0, 90, np.array([0, 0]))
        all_planes.append(p)
        for m in p.modules:
            all_modules.append(m)
        _ = modules.JoinedStrand(p.modules)

        p = modules.Plane(np.array([-1.3 - 13.529450, 0, -13.529450]), pitch=-90, yaw=90)
        (_, idx) = p.add_arc(idx, 66, False, 0, 90, np.array([0, 0]))
        all_planes.append(p)
        for m in p.modules:
            all_modules.append(m)
        _ = modules.JoinedStrand(p.modules)

        p = modules.Plane(np.array([1.3 - 13.529450, 0, -13.529450]), pitch=-90, yaw=90)
        (_, idx) = p.add_arc(idx, 78, False, 0, 90, np.array([0, 0]))
        all_planes.append(p)
        for m in p.modules:
            all_modules.append(m)
        _ = modules.JoinedStrand(p.modules)

        p = modules.Plane(np.array([1.3 - 13.529450, 0, -13.529450]), pitch=-90, yaw=90)
        (_, idx) = p.add_arc(idx, 66, True, 0, 90, np.array([0, 0]))
        all_planes.append(p)
        for m in p.modules:
            all_modules.append(m)
        _ = modules.JoinedStrand(p.modules)

        p = modules.Plane(np.array([0, -1.3 - 13.529450, -13.529450]), pitch=0, yaw=270)
        (_, idx) = p.add_arc(idx, 140, True, 0, 180, np.array([0, 0]))
        all_planes.append(p)
        for m in p.modules:
            all_modules.append(m)
        _ = modules.JoinedStrand(p.modules)

        p = modules.Plane(np.array([0, -1.3 - 13.529450, -13.529450]), pitch=0, yaw=270)
        (_, idx) = p.add_arc(idx, 116, False, 0, 180, np.array([0, 0]))
        all_planes.append(p)
        for m in p.modules:
            all_modules.append(m)
        _ = modules.JoinedStrand(p.modules)

        p = modules.Plane(np.array([0, 1.3 - 13.529450, -13.529450]), pitch=0, yaw=270)
        (_, idx) = p.add_arc(idx, 140, False, 0, 180, np.array([0, 0]))
        all_planes.append(p)
        for m in p.modules:
            all_modules.append(m)
        _ = modules.JoinedStrand(p.modules)

        p = modules.Plane(np.array([0, 1.3 - 13.529450, -13.529450]), pitch=0, yaw=270)
        (_, idx) = p.add_arc(idx, 116, True, 0, 180, np.array([0, 0]))
        all_planes.append(p)
        for m in p.modules:
            all_modules.append(m)
        _ = modules.JoinedStrand(p.modules)

        return all_planes, all_modules
