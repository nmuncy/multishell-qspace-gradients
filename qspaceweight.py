#!/usr/bin/env python
r"""Convert samples to siemens DVS file.

Write a siemens.dvs file to the same directory as the unitary-schema input.

Example:
    python qspaceweight.py \
        --unitary-schema samples.txt \
        --num-b0 5 \
        --interspersed \
        --bvalues 1000 2000 3000

"""
# %%
import os
import sys
import pandas as pd
import numpy as np
from argparse import ArgumentParser, RawTextHelpFormatter

type PT = str | os.PathLike


# %%
def get_args():
    """Get and parse arguments."""
    parser = ArgumentParser(
        description=__doc__, formatter_class=RawTextHelpFormatter
    )
    parser.add_argument(
        "--unitary-schema",
        type=str,
        help="Sample.txt generated from http://www.emmanuelcaruyer.com/q-space-sampling.php",
    )
    parser.add_argument(
        "--num-b0",
        type=int,
        help="Number of B0 in the begining of the acquisition.",
    )
    parser.add_argument(
        "--interspersed",
        action="store_true",
        help="B0 volumes should be interspersed.",
    )
    parser.add_argument(
        "--bvalues",
        nargs="+",
        type=int,
        help="B-Values vector separated by spaces (i.e: 1000 2000 3000). "
        + "Number must match the shells in [unitary_schema]",
    )

    if len(sys.argv) <= 1:
        parser.print_help(sys.stderr)
        sys.exit(0)

    return parser


# %%
def intersperse_b0(directions: list, num_b0: int) -> list:
    """
    Function to intersperse a given number of b-value=0 volumes in
    between the original directions.
    It assumes the scanner will collect a first b-value=0 volume at
    the beginning of the acquisition.

    Parameters
    ----------
    directions : list
    num_b0 : int
        Number of b-value=0 volumes

    Returns
    -------
    interspersed_directions : list

    """
    # To intersperse num_b0 b=0 volumes more or less equidistant in between the N
    # DW images we should divide the N volumes in (num_b0+1) blocks, with each of
    # the num_b0 first blocks having a b=0 volumes at the end, and the last block
    # not having any.

    N = len(directions)  # number of original directions
    n_blocks = N // num_b0 + 1  # number of blocks (equivalent of "ceil")
    vec_per_block = N // int(
        n_blocks
    )  # number of points (directions) per block (equiv. of "floor")

    interspersed_directions = [None] * (N + num_b0)
    oi = 0  # counter for out_dirs
    ii = 0  # counter for input_dirs
    for _ in range(n_blocks):

        # copy the next block from the input:
        interspersed_directions[oi : oi + vec_per_block] = directions[
            ii : ii + vec_per_block
        ]
        oi += vec_per_block
        ii += vec_per_block

        # introduce a b-value = 0:
        interspersed_directions[oi] = "( 0.000, 0.000, 0.000 )"
        oi += 1

    # last block (it might not have 'vec_per_block' points,
    # but just a few left):
    vec_left = (
        N + num_b0 - oi
    )  # how many points/vectors/directions we have left
    for _ in range(vec_left):
        interspersed_directions[oi] = directions[ii]
        oi += 1
        ii += 1

    return interspersed_directions


# %%
def read_unitary(unitary_schema: PT) -> np.ndarray:
    """Convert unitary schem to np.ndarray."""
    df_sample = pd.read_table(
        unitary_schema,
        sep="\t",
        comment="#",
        header=None,
        names=["shell", "x", "y", "z"],
    )
    df_sample = df_sample.round({"x": 6, "y": 6, "z": 6})
    schema = df_sample.to_numpy()
    return schema


# %%
def organize_schema(bvalues: list, schema: np.ndarray) -> list:
    """Organize schema np.ndarray into list of directions."""
    directions = []
    maxb = max(bvalues)
    for _, dir_list in enumerate(schema):
        bvalue = bvalues[int(dir_list[0] - 1)]
        weight = np.sqrt(float(bvalue) / float(maxb))
        u = dir_list[1:4]
        v = u * weight
        v = v.round(6)
        directions.append(f"( {v[0]}, {v[1]}, {v[2]} )")
    return directions


# %%
def write_siemens(
    directions: list, num_b0: int, interspersed: bool, out_dir: PT
):
    """Write siemens.dvs file from directions."""
    siemens_file = os.path.join(out_dir, "siemens.dvs")
    with open(siemens_file, "w") as fd:
        fd.write(f"[directions={len(directions) + num_b0}]\n")
        fd.write("CoordinateSystem = xyz\n")
        fd.write("Normalisation = none\n")

        if interspersed:
            directions_interspersed = intersperse_b0(directions, num_b0)
            for n, dir in enumerate(directions_interspersed):
                fd.write(f"Vector[{n}] = {dir}\n")
            return

        # Non-interspersed b0 files are first.
        n = 0
        for _ in range(num_b0):
            fd.write(f"Vector[{n}] = ( 0.000, 0.000, 0.000 )\n")
            n = n + 1
        for dir in directions:
            fd.write(f"Vector[{n}] = {dir}\n")
            n = n + 1


# %%
def main():
    """Convert unitary to siemens schema."""

    # Capture args
    args = get_args().parse_args()
    unitary_schema = args.unitary_schema
    interspersed = args.interspersed
    bvalues = np.array(args.bvalues)
    num_b0 = args.num_b0

    # # For testing
    # unitary_schema = "/home/nmuncy2/Desktop/samples.txt"
    # interspersed = True
    # bvalues = np.array([1000, 2000, 3000], dtype=int)
    # num_b0 = 5

    # Read unitary
    out_dir = os.path.dirname(unitary_schema)
    schema = read_unitary(unitary_schema)

    # Check number of shells in input schema:
    input_nshells = int(max([item[0] for item in schema]))
    given_bvalues = bvalues.size
    if input_nshells != given_bvalues:
        raise ValueError(
            "Given B-values number and provided Sample.txt "
            + "shells number doesn't match."
        )

    # Organize schema, write siemens file
    directions = organize_schema(bvalues, schema)
    write_siemens(directions, num_b0, interspersed, out_dir)


if __name__ == "__main__":
    main()
