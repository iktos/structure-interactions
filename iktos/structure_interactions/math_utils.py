from __future__ import absolute_import

from itertools import product
from logging import getLogger

import numpy as np

logger = getLogger(__name__)


def get_euclidean_distance_3d(v1, v2):
    """
    Faster implementation of euclidean distance for the 3D case
    """
    if isinstance(v1, list):
        v1 = np.array(v1)
    if isinstance(v2, list):
        v2 = np.array(v2)

    return np.sqrt(np.sum((v1 - v2) ** 2, axis=len(v1.shape) - 1))


def get_vector(p1, p2):
    """
    Calculate vector from p1 to p2
    p1: coordinates of point p1
    p2: coordinates of point p2
    Return numpy array with vector coordinates
    """
    if not (len(p1) == 3 and len(p2) == 3):
        logger.error(f'Invalid vector format: {p1}, {p2}')
        raise ValueError()
    return (
        None
        if len(p1) != len(p2)
        else np.array([p2[i] - p1[i] for i in range(len(p1))])
    )


def get_vector_angle(v1, v2, deg=True):
    """
    Calculate the angle between two vectors
    v1: coordinates of vector v1
    v2: coordinates of vector v2
    Return angle in degree (range 0-180) or rad
    Note: round value to prevent floating point errors
    """
    if np.array_equal(v1, v2):
        return 0.0
    dm = np.dot(v1, v2)
    cm = np.linalg.norm(v1) * np.linalg.norm(v2)
    angle = np.arccos(round(dm / cm, 10))
    return np.degrees([angle])[0] if deg else angle


def normalize_vector(v):
    """
    Take a vector and return the normalized vector
    """
    norm = np.linalg.norm(v)
    return v / norm if not norm == 0 else v


def get_centroid(coo):
    """
    Calculate the centroid from a 3D point cloud
    coo: array of coordinates
    Return centroid coordinates as list
    """
    center = np.array([0.0, 0.0, 0.0])
    for i, j in product(range(3), range(len(coo))):
        center[i] += coo[j][i]
    center = center / float(len(coo))
    return center


def project_on_plane(pnormal1, ppoint, tpoint):
    """
    Project point onto plane
    pnormal1: normal of plane
    ppoint: coordinates of point in the plane
    tpoint: coordinates of point to be projected
    Return coordinates of point orthogonally projected on the plane
    """
    # Choose the plane normal pointing to the point to be projected
    pnormal2 = [coo * (-1) for coo in pnormal1]
    d1 = get_euclidean_distance_3d(tpoint, pnormal1 + ppoint)
    d2 = get_euclidean_distance_3d(tpoint, pnormal2 + ppoint)
    pnormal = pnormal1 if d1 < d2 else pnormal2
    # Calculate the projection of tpoint to the plane
    sn = -np.dot(pnormal, get_vector(ppoint, tpoint))
    sd = np.dot(pnormal, pnormal)
    sb = sn / sd
    return [c1 + c2 for c1, c2 in zip(tpoint, [sb * pn for pn in pnormal])]
