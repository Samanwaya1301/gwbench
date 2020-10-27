**Source Oasis Conda - *do this first, on LIGO clusters only***:  
source /cvmfs/oasis.opensciencegrid.org/ligo/sw/conda/etc/profile.d/conda.sh  
which conda (this should print `/cvmfs/oasis.opensciencegrid.org/ligo/sw/conda/condabin/conda` or similar)  

**Setup conda VE:**  
conda create -y --name gwbench python=3.7  
conda activate gwbench  
conda install -y -c conda-forge --file requirements_conda.txt  

**Install while *VE gwbench ACTIVE* with:**  
pip install .

**Test:**  
cd example_scripts  
python quick_start.py

**Uninstall:**  
pip uninstall gwbench
