# centriole_data
scripts for scraping annotation data from CATMAID and computing some stats

# Usage
## 1) Authentication
First replace `catmaid_url`, `project_id`, `api_token`, `username`, and
`password` fields in `config.py`.

Note that the `catmaid_url` should be a fully qualified URL, including the
protocol and a possible sub-directory, e.g.
"https://spaces.itanna.io/catmaid/itanna/".

The `project_id` is an integer number that can be found by opening the
respective project in CATMAID and create a URL to the current view with the link
buttons in the upper right corner of the window. If you paste this URL, e.g.
into an address bar, you will see the part `pid=<a-number>`. The number after
the equals sign is the project ID.

The `api_token` of a user for a project can be obtained through the user menu
that appears when hovering the mouse cursor over the own name in the upper right
corner, after signing in. The option "Get API token" will show a pop-up dialog
asking for a password confirmation. Once entered, the API token will be visible.

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
