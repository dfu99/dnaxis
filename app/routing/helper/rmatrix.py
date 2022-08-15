# -*- coding: utf-8 -*-
"""
Created on Sat Oct 13 22:27:29 2018

@author: Dan

rotater.py

rotates a vector by an angle (converted to a rotation matrix)
"""

import numpy as np
import math


def rotate(v, a):
    """
    Rotate a vector 'v' by angle 'a' (degrees)
    :param v: vector
    :param a: float degrees
    :return: new vector
    """
    rm = rmat(a)
    v_ = rm.dot(v)
    return v_


def rmat(a):
    """
    Converts an angle (degrees) to a rotation matrix
    :param a: float angle (degrees)
    :return: rotation matrix
    """
    a = math.radians(a)
    rm = np.array([
            [math.cos(a), -math.sin(a)],
            [math.sin(a), math.cos(a)]
            ])
    return rm
