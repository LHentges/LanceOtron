# LanceOTron Command Line Interface

[![CircleCI](https://circleci.com/gh/Chris1221/lanceotron/tree/main.svg?style=svg&circle-token=bf3f78a54437e63368f5b9dc1c536d7f32f32393)](https://circleci.com/gh/Chris1221/lanceotron/tree/main)

A bare-bones interface to the trained LanceOTron (LoT) model from the command line. 

LoT is an all-in-one peak caller that identifies peak regions from a coverage track and uses a convolutional neural network to classify them based on their shape. This algorithm has a web client available at https://lanceotron.molbiol.ox.ac.uk/ where users can upload coverage tracks, call peaks, and visualize them using [multi locus view](https://lanceotron.readthedocs.io/en/latest/multi_locus_view/multi_locus_view.html), a powerful visualization engine. This web client will do most of the heavy lifting for most of the people that want to use the tool, but for those who need to call peaks in batch mode we provide a command line interface in this package. This document details the installation of LoT as well as typical use cases. See the left side of this page for tutorials on the three main modules. 


## Installation

1. Clone the repository.
2. Install dependencies with pip.
3. Install the package.
4. Run tests to ensure that everything is working.

```{sh}
git clone git@github.com:Chris1221/lanceotron.git; cd lanceotron # Step 1
pip install -r requirements.txt # Step 2
pip install -e . # Step 3
python -m unittest
```

## Usage

To see available commands, use the `--help` flag.

```
lanceotron --help
```

## Call Peaks

To call peaks from a bigWig track, use the `callPeaks` command.

::: lanceotron.find_and_score_peaks


## Call Peaks with Input

To call peaks from a bigWig track with an input file, use the `callPeaks_Input` command.

::: lanceotron.call_peaks_with_input

## Score a Bed file

To score the peaks in an existing Bed file, use the `scoreBed` command.

::: lanceotron.score_bed

## Examples

There is a basic bigWig file included in the `test` subdirectory. To try out the caller, execute it on this file. 

```sh
lanceotron callPeaks test/chr22.bw -f output_folder
```

## Citation

```{bibtex}
@article {Hentges2021.01.25.428108,
	author = {Hentges, Lance D. and Sergeant, Martin J. and Downes, Damien J. and Hughes, Jim R. and Taylor, Stephen},
	title = {LanceOtron: a deep learning peak caller for ATAC-seq, ChIP-seq, and DNase-seq},
	year = {2021},
	doi = {10.1101/2021.01.25.428108},
	publisher = {Cold Spring Harbor Laboratory},
	URL = {https://www.biorxiv.org/content/early/2021/01/27/2021.01.25.428108},
	journal = {bioRxiv}
}
```

## Bug Reports and Improvement Suggestions

Please [raise an issue](https://github.com/Chris1221/lanceotron/issues/new/choose) if there is anything you wish to ask or contribute. 