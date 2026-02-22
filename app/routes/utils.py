from flask import session
import numpy as np
from app.routing.helper.mymath import bp2radius
from app.scaffold_library import ALLOWED_EXTENSIONS


def allowed_file(filename):
    """Check whether file has an allowed extension."""
    return '.' in filename and \
           filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS


def is_connected(graph, num):
    """Check whether graph is fully connected."""
    nums = [i for i in range(num)]
    initial_node = list(graph.keys())[0]
    stack = []
    stack.append(initial_node)
    nums.remove(initial_node)
    while stack:
        next_node = stack.pop()
        for target_node in graph[next_node]:
            if target_node != next_node:
                try:
                    nums.remove(target_node)
                    stack.append(target_node)
                except ValueError:
                    pass
    if nums:
        return False
    else:
        return True


def has_cycle(graph):
    """Check for cycle in graph."""
    initial_node = list(graph.keys())[0]
    stack = []
    stack.append(initial_node)
    visited = []
    while stack:
        next_node = stack.pop()
        visited.append(next_node)
        for target_node in graph[next_node]:
            if target_node in visited:
                return True
            elif target_node in stack:
                pass
            else:
                stack.append(target_node)
                graph[target_node].remove(next_node)
    return False


def is_lteq_degree(graph, deg):
    """Check if every node has at most deg edges."""
    for node in graph:
        if len(graph[node]) > deg:
            return False
    return True


def edges_to_linkedlist(edges):
    """Convert edges to linked list (adjacency list)."""
    linkedlist = {}
    for e in edges:
        try:
            linkedlist[e[0]].append(e[1])
        except KeyError:
            linkedlist[e[0]] = [e[1]]
        try:
            linkedlist[e[1]].append(e[0])
        except KeyError:
            linkedlist[e[1]] = [e[0]]
    return linkedlist


def edge_len_bound(edges, lower, upper):
    """Check if any edge is within a distance bound."""
    nodes = session['ringdata']
    for e in edges:
        n1 = nodes[e[0]]
        n2 = nodes[e[1]]

        n1xy = np.array([bp2radius(n1[0]), n1[1]])
        n2xy = np.array([bp2radius(n2[0]), n2[1]])

        d = np.linalg.norm(n2xy - n1xy)
        if lower <= d <= upper:
            return True
        else:
            return False
