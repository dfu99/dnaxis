"""
exceptions.py

A list of custom exceptions for CADAxiSDNA

Author: Daniel Fu
"""


class NoValidXoverThreshold(Exception):
    pass


class NoValidSpacing(Exception):
    pass


class HnickError(Exception):
    pass


class DistanceError(Exception):
    pass

class NuclNotFound(Exception):
    pass
