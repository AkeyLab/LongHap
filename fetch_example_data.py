#!/usr/bin/env python3
"""Download the LongHap example data from GitHub."""

import argparse
import io
import sys
import tarfile
import urllib.request
from pathlib import Path

REPO = "AkeyLab/LongHap"
BRANCH = "main"
TARBALL_URL = f"https://github.com/{REPO}/archive/refs/heads/{BRANCH}.tar.gz"
EXAMPLE_PREFIX = f"LongHap-{BRANCH}/example/"


def fetch_example(output_dir: Path) -> None:
    if output_dir.exists():
        print(f"Error: '{output_dir}' already exists. Remove it first or choose a different directory.", file=sys.stderr)
        sys.exit(1)

    print(f"Downloading example data from {REPO} ({BRANCH} branch)...")
    with urllib.request.urlopen(TARBALL_URL) as response:
        tarball = response.read()

    print("Extracting example/ directory...")
    output_dir.mkdir(parents=True, exist_ok=True)
    with tarfile.open(fileobj=io.BytesIO(tarball), mode="r:gz") as tar:
        for member in tar.getmembers():
            if not member.name.startswith(EXAMPLE_PREFIX) or member.name == EXAMPLE_PREFIX:
                continue
            member.name = member.name[len(EXAMPLE_PREFIX):]
            tar.extract(member, path=output_dir, filter="data")

    print(f"Example data saved to {output_dir}/")


def main():
    parser = argparse.ArgumentParser(
        prog="longhap_fetch_example_data",
        description="Download the LongHap example dataset from GitHub.",
    )
    parser.add_argument(
        "-o", "--output",
        type=Path,
        default=Path("example"),
        help="Output directory (default: ./example)",
    )
    args = parser.parse_args()
    fetch_example(args.output)


if __name__ == "__main__":
    main()
