"""
ic/io
=====

Holds functions for writing and reading data for initial conditions.
"""
__all__ = [
    "write_cubeNCD",
    "read_cubeNCD"
]

import numpy as np

def write_cube(filename, data, NCD,
               include_header=False):
    """
    Write a data cube to file.

    Arguments
    ---------
    filename : str
        Output file path.

    data : dict[str, ndarray]
        Dictionary of fields to write. All arrays must have shape NCD.

    NCD : tuple of int
        Grid dimensions (nx, ny, nz).
    
    include_header : bool, optional
        If True, include header at top of file.
    """
    nx, ny, nz = NCD

    # Validate shapes
    for key, arr in data.items():
        if arr.shape != (nx, ny, nz):
            raise ValueError(f"{key} has shape {arr.shape}, expected {NCD}")

    columns = []
    formats = []
    header = []

    # Indices
    i, j, k = np.indices((nx, ny, nz))
    columns.extend([i.ravel(), j.ravel(), k.ravel()])
    formats.extend(['%3d', '%3d', '%3d'])
    header.extend(["i", "j", "k"])

    for key in data:
        columns.append(data[key].ravel())
        formats.append('%15.7e')
        header.append(key)

    out = np.column_stack(columns)
    fmt = " ".join(formats)
    header = " ".join(header)

    with open(filename, "w") as f:
        f.write(f"# {nx} {ny} {nz}\n")

        if include_header:
            f.write("# " + header + "\n")

        np.savetxt(f, out, fmt=fmt)

def read_cube(filename):
    """
    Read a data cube from file.

    Arguments
    ---------
    filename : str
        Input file path.

    Returns
    -------
    data : dict[str, ndarray]
        Dictionary of fields with shape (nx, ny, nz).

    NCD : tuple of int
        Grid dimensions (nx, ny, nz).
    """
    with open(filename, "r") as f:
        lines = f.readlines()

    if not lines[0].startswith("#"):
        raise ValueError("First line must contain grid dimensions.")

    nx, ny, nz = map(int, lines[0][1:].split())
    NCD = (nx, ny, nz)

    has_header = len(lines) > 1 and lines[1].startswith("#")
    data_start = 2 if has_header else 1

    if has_header:
        names = lines[1][1:].split()
    else:
        names = None

    raw = np.loadtxt(lines[data_start:])
    ncols = raw.shape[1]

    if names is None:
        names = ["i", "j", "k"] + [f"col{i}" for i in range(ncols - 3)]

    if len(names) != ncols:
        raise ValueError("Header column count does not match data.")

    data = {}
    data["i"] = raw[:, 0].astype(int).reshape(NCD)
    data["j"] = raw[:, 1].astype(int).reshape(NCD)
    data["k"] = raw[:, 2].astype(int).reshape(NCD)

    for idx, name in enumerate(names[3:], start=3):
        data[name] = raw[:, idx].reshape(NCD)

    return data, NCD