from __future__ import absolute_import

from itertools import product

try:
    from iktos.logger import getLogger
except ImportError:
    from logging import getLogger

import numpy as np
import numpy.typing as npt

LOGGER = getLogger(__name__)


def matrix_multiply(*matrices):
    if len(matrices) == 1:
        return matrices
    else:
        m_other = matrix_multiply(*matrices[1:])
        return np.matmul(matrices[0], m_other)


def rotate(
    point: npt.NDArray, axis_1: npt.NDArray, axis_2: npt.NDArray, theta: float
) -> npt.NDArray:
    """Rotates point around axis 1-2 by angle theta (radian)."""
    p = [[pp] for pp in point] + [[1]]
    x1, y1, z1 = axis_1
    x2, y2, z2 = axis_2

    U = [x2 - x1, y2 - y1, z2 - z1]
    U = np.array(U) / np.sqrt(np.dot(U, U))
    a, b, c = U
    d = np.sqrt(b**2 + c**2)

    T = [[1, 0, 0, -x1], [0, 1, 0, -y1], [0, 0, 1, -z1], [0, 0, 0, 1]]
    T_inv = [[1, 0, 0, x1], [0, 1, 0, y1], [0, 0, 1, z1], [0, 0, 0, 1]]

    R_x = [[1, 0, 0, 0], [0, c / d, -b / d, 0], [0, b / d, c / d, 0], [0, 0, 0, 1]]
    R_x_inv = [[1, 0, 0, 0], [0, c / d, b / d, 0], [0, -b / d, c / d, 0], [0, 0, 0, 1]]

    R_y = [[d, 0, -a, 0], [0, 1, 0, 0], [a, 0, d, 0], [0, 0, 0, 1]]
    R_y_inv = [[d, 0, a, 0], [0, 1, 0, 0], [-a, 0, d, 0], [0, 0, 0, 1]]

    ct = np.cos(theta)
    st = np.sin(theta)
    R_z = [[ct, st, 0, 0], [-st, ct, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]]

    p2 = matrix_multiply(T_inv, R_x_inv, R_y_inv, R_z, R_y, R_x, T, p)
    return p2[0][:3].T[0]


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
    """Calculate vector from p1 to p2.

    Args:
        p1: coordinates of point p1.
        p2: coordinates of point p2.

    Returns:
        Numpy array with vector coordinates.
    """
    if not (len(p1) == 3 and len(p2) == 3):
        LOGGER.error(f'Invalid vector format: {p1}, {p2}')
        raise ValueError()
    return (
        None
        if len(p1) != len(p2)
        else np.array([p2[i] - p1[i] for i in range(len(p1))])
    )


def get_vector_angle(v1, v2, deg=True):
    """Calculate the angle between two vectors.

    Note:
        Round value to prevent floating point errors.

    Args:
        v1: coordinates of vector v1.
        v2: coordinates of vector v2.
        deg: Whether convert angles from radians to degrees (default True).

    Returns:
        Angle in degree (range 0-180) or rad.
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
