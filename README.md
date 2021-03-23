# centriole_data
scripts for scraping annotation data from CATMAID and computing some stats

# Usage
## 1) Authentication
First replace the api_token, username, and password fields in config.py

## 2) Create a python environment
These scripts were tested on python3.8. Feel free to use other versions but
there is a chance that they won't work as expected.
```
conda create -n centrioles python=3.8
conda activate centrioles
```

## 3) Install dependencies
```
pip install tqdm
pip install networkx
pip install numpy
pip install requests
```

## 4) Run the script
This script is somewhat slow. I have added a progress bar to help.
```
python get_data.py
```
