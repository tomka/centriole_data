import requests
import pickle
from pathlib import Path
import collections
from math import atan2

from tqdm import tqdm

from api_requests import AUTH_TOKEN, CATMAID_URL, PROJECT_ID

import numpy as np


def get_vector_coords(object_id):
    centriole_vector = requests.post(
        f"{CATMAID_URL}/{PROJECT_ID}/{object_id}/0/1/compact-skeleton",
        verify=False,
        auth=AUTH_TOKEN,
        data="",
    ).json()
    ends = centriole_vector[0]
    assert len(ends) == 2
    a_loc = np.array(ends[0][3:6])
    b_loc = np.array(ends[1][3:6])
    dx, dy, dz = b_loc - a_loc

    # get magnitude
    r = np.linalg.norm(a_loc - b_loc)
    # get mother coordinates
    x, y, z = a_loc
    # get theta/phi
    theta = atan2((dx**2 + dy**2)**0.5, dz)
    phi = atan2(dx, dy)
    return x, y, z, theta, phi, r


neurons_with_annotations = requests.post(
    f"{CATMAID_URL}/{PROJECT_ID}/annotations/query-targets",
    verify=False,
    auth=AUTH_TOKEN,
    headers={"content-type": "application/x-www-form-urlencoded; charset=UTF-8"},
    data="sort_by=name&sort_dir=asc&types[0]=neuron&with_annotations=true&with_timestamps=true",
)

neurons_with_annotations = neurons_with_annotations.json()["entities"]

vectors = [
    skeleton
    for skeleton in neurons_with_annotations
    if "centriole vector line"
    in list(annotation["name"] for annotation in skeleton["annotations"])
]


cells = []
for vector in tqdm(vectors, desc="Getting cells: "):
    name = vector["name"]
    # get object id and assure we only have one
    object_id = vector["skeleton_ids"]
    assert len(object_id) == 1
    object_id = object_id[0]

    # For robustness we don't care about the case of the annoation and convert
    # everything to lower case.
    annotation_list = [
        annotation["name"].lower() for annotation in vector["annotations"]
    ]
    annotations = set()
    cell_num = None
    for annotation in annotation_list:
        if "cell #" in annotation:
            cell_num = int(annotation[-3:])
        else:
            annotations.add(annotation)

    # ignore objects not assigned to a cell
    if cell_num is None:
        continue

    if len(annotations) == 0:
        print(f"No annotations for: {name}: {annotation_list}")
        continue

    cells.append((cell_num, object_id))

rows = []
for cell_num, object_id in tqdm(cells, desc="Fetching cell stats: "):
    x, y, z, theta, phi, r = get_vector_coords(object_id)
    rows.append((cell_num, x, y, z, theta, phi, r))

with open("centriole_vectors.txt", "w") as f:
    f.write("\n".join([",".join([str(x) for x in row]) for row in rows]))
