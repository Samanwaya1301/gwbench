# gwbench

## Acknowledgment

We request that any academic report, publication, or other academic
disclosure of results derived from the use of this software acknowledge
the use of the software by an appropriate acknowledgment or citation.

The gwbench software can be cited from [arXiv:2010.15202](https://arxiv.org/abs/2010.15202), with INSPIRE BibTeX entry:
```
@article{Borhanian:2020ypi,
    author = "Borhanian, Ssohrab",
    title = "{gwbench: a novel Fisher information package for gravitational-wave benchmarking}",
    eprint = "2010.15202",
    archivePrefix = "arXiv",
    primaryClass = "gr-qc",
    month = "10",
    year = "2020"
}
```

## Installation

These instructions are still developing and written for users with access to the Conda
package manager (e.g. via [Miniconda](https://docs.conda.io/en/latest/miniconda.html)).
The required dependencies are specified in `requirements_conda.txt`, if an installation via
other means is preferred.

Clone this repository and follow the next steps.

#### Source Oasis Conda - *do this first; only on LIGO clusters needed*
```
source /cvmfs/oasis.opensciencegrid.org/ligo/sw/conda/etc/profile.d/conda.sh  
which conda
```

The last line should print `/cvmfs/oasis.opensciencegrid.org/ligo/sw/conda/condabin/conda` or similar.

#### Setup conda virtual environment
```
conda create -y --name gwbench python=3.7  
conda activate gwbench  
conda install -y -c conda-forge --file requirements_conda.txt  
```

#### Install while virtual environment `gwbench` active
```
pip install .
```

#### Uninstall
```
pip uninstall gwbench
```

### Test
```
cd example_scripts  
python quick_start.py
```