import requests
import pickle
from pathlib import Path
import collections

from tqdm import tqdm

from api_requests import AUTH_TOKEN, CATMAID_URL, PROJECT_ID

from cell_stats import (
    get_cell,
    get_basal_body,
    get_distance_a,
    get_depth_a,
    get_distance_b,
    get_depth_b,
    get_distance_ab,
    get_cilia_length,
    get_centriole,
    get_cilium,
    get_cell_type,
    get_location,
    get_cell_cycle,
    get_migrating,
)

columns = {
    "Cell": get_cell,
    "Basal body": get_basal_body,
    "Distance A": get_distance_a,
    "Depth A": get_depth_a,
    "Distance B": get_distance_b,
    "Depth B": get_depth_b,
    "Distance AB": get_distance_ab,
    "Cilia length": get_cilia_length,
    "Centriole": get_centriole,
    "Cilium": get_cilium,
    "Cell type": get_cell_type,
    "Location": get_location,
    "Cell Cycle": get_cell_cycle,
    "migrating": get_migrating,
}
columns = [(k, v) for k, v in columns.items()]

all_annotations = set()


neurons_with_annotations = requests.post(
    f"{CATMAID_URL}/{PROJECT_ID}/annotations/query-targets",
    verify=False,
    auth=AUTH_TOKEN,
    headers={"content-type": "application/x-www-form-urlencoded; charset=UTF-8"},
    data="sort_by=name&sort_dir=asc&types[0]=neuron&with_annotations=true&with_timestamps=true",
)

neurons_with_annotations = neurons_with_annotations.json()["entities"]

# data I need: neuron_id -> (cell number, [annotation_names])

# cell_num -> [(object_id, annotations, name, annotation_list)]
cells = collections.defaultdict(list)

for neuron in tqdm(neurons_with_annotations, desc="Getting cells: "):
    name = neuron["name"]
    object_id = neuron["skeleton_ids"]
    assert len(object_id) == 1
    object_id = object_id[0]
    # For robustness we don't care about the case of the annoation and convert
    # everything to lower case.
    annotation_list = [annotation["name"].lower() for annotation in neuron["annotations"]]
    annotations = set()
    cell_num = None
    for annotation in annotation_list:
        if "cell #" in annotation:
            cell_num = int(annotation[-3:])
        else:
            annotations.add(annotation)
    all_annotations = all_annotations.union(annotations)

    # ignore objects not assigned to a cell
    if cell_num is None:
        continue

    if len(annotations) == 0:
        print(f"No annotations for: {name}: {annotation_list}")
        continue

    cells[cell_num].append((object_id, annotations, name))

rows = []
for cell_num, objects in tqdm(cells.items(), desc="Fetching cell stats: "):
    cell = [None] * len(columns)
    for i, (k, v) in enumerate(columns):
        cell[i] = v(objects, cell_num)
    rows.append(cell)

with open("cell_data.csv", "w") as f:
    f.write(",".join([k for k, v in columns]) + "\n")
    for row in tqdm(rows, "Writing to file: "):
        f.write(",".join([str(x) for x in row]) + "\n")

print("Cell data has been written to 'cell_data.csv'")