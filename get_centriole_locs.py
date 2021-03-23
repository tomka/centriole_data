import requests
import pickle
from pathlib import Path

import collections

from api_requests import AUTH_TOKEN

from cell_stats import (
    get_cell,
    get_centriole_coords,
)

columns = {
    "Cell": get_cell,
    "coords": get_centriole_coords,
}
columns = [(k, v) for k, v in columns.items()]

all_annotations = set()

data_file = Path("centrioles.obj")
if data_file.exists():
    rows = pickle.load(data_file.open("rb"))
else:

    neurons_with_annotations = requests.post(
        "https://jls.janelia.org/catmaid/49/annotations/query-targets",
        verify=False,
        auth=AUTH_TOKEN,
        data="sort_by=name&sort_dir=asc&types[0]=neuron&with_annotations=true&with_timestamps=true",
    )
    neurons_with_annotations = neurons_with_annotations.json()["entities"]

    # data I need: neuron_id -> (cell number, [annotation_names])

    # cell_num -> [(object_id, annotations, name, annotation_list)]
    cells = collections.defaultdict(list)

    for neuron in neurons_with_annotations:
        name = neuron["name"]
        object_id = neuron["skeleton_ids"]
        assert len(object_id) == 1
        object_id = object_id[0]
        annotation_list = [annotation["name"] for annotation in neuron["annotations"]]
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
    for cell_num, objects in cells.items():
        cell = [None] * len(columns)
        for i, (k, v) in enumerate(columns):
            cell[i] = v(objects, cell_num)
        cell = get_cell(objects, cell_num)
        locs = get_centriole_coords(objects, cell_num)
        for loc in locs:
            rows.append([float(x) for x in loc])

    pickle.dump(rows, data_file.open("wb"))

print(f"Found {len(rows)} centrioles!")
with open("centriole_locs.csv", "w") as f:
    for row in rows:
        f.write(" ".join([str(x) for x in row]) + "\n")
