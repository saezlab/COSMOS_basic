"""Run the NCI60 practical from any working directory."""

import argparse
from dataclasses import fields
import os
from pathlib import Path
import tempfile

# Make cache writes local to a standard temporary directory on managed desktops.
os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "cosmos-grn-matplotlib"))
os.environ.setdefault("NUMBA_CACHE_DIR", str(Path(tempfile.gettempdir()) / "cosmos-grn-numba"))

import matplotlib
matplotlib.use("Agg")
from grn_benchmark import Config, run_benchmark


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for field in fields(Config):
        parser.add_argument("--" + field.name.replace("_", "-"),
                            type=type(field.default), default=field.default)
    parser.add_argument("--output-dir", default=None,
                        help="Defaults to GRN_comparison/results; existing result files are replaced.")
    args = vars(parser.parse_args())
    output_dir = args.pop("output_dir")
    try:
        result = run_benchmark(Config(**args), output_dir)
    except (ValueError, FileNotFoundError) as error:
        parser.error(str(error))
    print(result["network_stats"].to_string(index=False))
    print(result["prediction_summary"].to_string(index=False))
    print(f"Results: {result['output_dir'].resolve()}")
    print(f"Analysis time: {result['run_info']['elapsed_seconds']:.1f} seconds (excluding imports)")


if __name__ == "__main__":
    main()
