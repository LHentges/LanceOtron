import argparse
import sys

from .modules import find_and_score_peaks, call_peaks_with_input, score_bed


def cli():
    parser = argparse.ArgumentParser(
        description="Sort significantly enriched regions of ChIP-seq singnals using a CNN"
    )
    subparsers = parser.add_subparsers(help="sub-command help")

    findandscore = subparsers.add_parser(
        "callPeaks", help="Call peaks from a bigWig file."
    )
    findandscore.add_argument("file", help="bigwig file")
    findandscore.add_argument(
        "-c",
        "--cutoff",
        type=float,
        default=0,
        help="Peak score cutoff to report region in output file.",
    )
    findandscore.add_argument(
        "--format",
        type=str,
        default="web",
        help="Format of the output file. 'Web' will mimic the output from the web app, 'Bed' will create a standard three column bed file.",
    )
    findandscore.add_argument(
        "-t",
        "--threshold",
        type=float,
        default=4,
        help="initial threshold used for selecting candidate peaks; default=4",
    )
    findandscore.add_argument(
        "-w",
        "--window",
        type=int,
        default=400,
        help="window size for rolling mean to use for selecting candidate peaks; default=400",
    )
    findandscore.add_argument(
        "-f",
        "--folder",
        type=str,
        default="./",
        help="folder to write results to; default=current directory",
    )
    findandscore.add_argument(
        "--skipheader", default=False, action="store_true", help="skip writing header"
    )
    findandscore.set_defaults(func=find_and_score_peaks)

    fas_input = subparsers.add_parser(
        "callPeaksInput",
        help="Call peaks from a bigWig file with an input track for control.",
    )
    fas_input.add_argument("file", help="bigwig file")
    fas_input.add_argument(
        "-i",
        "--input",
        type=str,
        help="control input track used to calculate Poisson-based significance of peaks",
    )
    fas_input.add_argument(
        "-t",
        "--threshold",
        type=float,
        default=4,
        help="initial threshold used for selecting candidate peaks; default=4",
    )
    fas_input.add_argument(
        "-w",
        "--window",
        type=int,
        default=400,
        help="window size for rolling mean to use for selecting candidate peaks; default=400",
    )
    fas_input.add_argument(
        "-f",
        "--folder",
        type=str,
        default="./",
        help="folder to write results to; default=current directory",
    )
    fas_input.add_argument(
        "--skipheader", default=False, action="store_true", help="skip writing header"
    )
    fas_input.set_defaults(func=call_peaks_with_input)

    score = subparsers.add_parser(
        "scoreBed",
        help="Score an existing bed file using Lanceotron's model and a coverage track",
    )
    score.add_argument("file", help="bigwig file")
    score.add_argument(
        "-b",
        "--bed",
        type=str,
        help="bed file of regions to be scored using L-tron's neural network",
    )
    score.add_argument(
        "-f",
        "--folder",
        type=str,
        default="./",
        help="folder to write results to; default=current directory",
    )
    score.add_argument(
        "--skipheader", default=False, action="store_true", help="skip writing header"
    )
    score.set_defaults(func=score_bed)

    # Parse the arguments and quit if no file specified.
    # This is mainly to catch lanceotron being called with no additional args.
    args = parser.parse_args()
    if "file" not in args:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args.func(**vars(args))
