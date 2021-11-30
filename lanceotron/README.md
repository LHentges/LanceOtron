# LanceOTron CLI

[![PyPI version](https://badge.fury.io/py/lanceotron.svg)](https://badge.fury.io/py/lanceotron) [![Downloads](https://pepy.tech/badge/lanceotron)](https://pepy.tech/project/lanceotron) [![CircleCI](https://circleci.com/gh/Chris1221/lanceotron/tree/main.svg?style=svg&circle-token=bf3f78a54437e63368f5b9dc1c536d7f32f32393)](https://circleci.com/gh/Chris1221/lanceotron/tree/main) [![codecov](https://codecov.io/gh/Chris1221/lanceotron/branch/main/graph/badge.svg?token=yhL3YI00UP)](https://codecov.io/gh/Chris1221/lanceotron)

A bare-bones interface to the trained LanceOTron (LoT) model from the command line.

## Installation

```{sh}
pip install lanceotron
```

### Local installation 

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

| Option          | Description                                            | Default |
|-----------------|--------------------------------------------------------|---------|
| file            | BigWig Track to analyse                                |         |
| -t, --threshold | Threshold for selecting candidate peaks                | 4       |
| -w, --window    | Window size for rolling mean to select candidate peaks | 400     |
| -f, --folder    | Output folder                                          | "./"    |
| --skipheader    | Skip writing the header                                | False   |


## Call Peaks with Input

To call peaks from a bigWig track with an input file, use the `callPeaks_Input` command.

| Option          | Description                                            | Default |
|-----------------|--------------------------------------------------------|---------|
| file            |  BigWig track to analyse                                |         |
| -i, --input     | Control input track to calculate significance of peaks                               |         |
| -t, --threshold | Threshold for selecting candidate peaks                | 4       |
| -w, --window    | Window size for rolling mean to select candidate peaks | 400     |
| -f, --folder    | Output folder                                          | "./"    |
| --skipheader    | Skip writing the header                                | False   |

## Score a Bed file

To score the peaks in an existing Bed file, use the `scoreBed` command.

| Option          | Description                                            | Default |
|-----------------|--------------------------------------------------------|---------|
| file            | BigWig Track to analyse                                |         |
| -b, --bed | Bed file of regions to be scored                |        |
| -f, --folder    | Output folder                                          | "./"    |
| --skipheader    | Skip writing the header                                | False   |


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

## Building the documentation

To serve the documentation locally, use

```
python -m mkdocs serve
```

## Bug Reports and Improvement Suggestions

Please [raise an issue](https://github.com/Chris1221/lanceotron/issues/new/choose) if there is anything you wish to ask or contribute. 
