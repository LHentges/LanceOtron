<p align="center">
    <img src="LanceOtron_logo_shadow_dark.png" alt="LanceOtron logo">
</p>


**LanceOtron** is a machine learning, genomic data extraction and analysis tool trained for ATAC-seq, ChIP-seq, and DNase-seq peak calling. A freely available and fully-featured webtool version, utilising the graphical user interface [MLV](https://mlv.molbiol.ox.ac.uk) and hosted at the [MRC WIMM Centre of Computational Biology, University of Oxford](https://www.imm.ox.ac.uk/research/units-and-centres/mrc-wimm-centre-for-computational-biology), can be found at [LanceOtron.molbiol.ox.ac.uk](https://lanceotron.molbiol.ox.ac.uk).

## Python Packages

LanceOtron uses Python 3 (3.8.3) and TensorFlow 2. The models have been saved such that a TensorFlow 1 setup could be used making only minor amendments to the scripts. Additional packages were used for benchmarking LanceOtron (N.B. [bedtools](https://github.com/arq5x/bedtools2) needs to be installed to use the python implementation). See [requirements.txt](requirements.txt) for specific versions used. 

**Required Python Packages for LanceOtron:**
> * scipy
> * numpy
> * pyBigWig
> * scikit\_learn==0.23.1
> * tensorflow==2.X.X

**Additional Python Packages for Benchmarking:**
> * pandas
> * matplotlib
> * pybedtools
> * seaborn

## Command Line Installation

We recommend using a fresh virtual environment with Python 3.7+, and installing the required packages with pip.

1. Clone/download repository
1. pip install requirements
> `pip install -r /path/to/repo/LanceOtron-master/requirements.txt`

## Citation

Please see our [bioRxiv article](https://www.biorxiv.org/content/10.1101/2021.01.25.428108v1) for further details on this project.

Hentges, L. D., Sergeant, M. J., Downes, D. J., Hughes, J. R., & Taylor, S. (2021). LanceOtron : a deep learning peak caller for ATAC-seq , ChIP-seq , and DNase-seq. *bioRxiv*. https://www.biorxiv.org/content/10.1101/2021.01.25.428108v1
