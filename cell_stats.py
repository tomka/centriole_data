import networkx as nx
import numpy as np

import requests
import warnings
import logging

from api_requests import AUTH_TOKEN, CATMAID_URL, PROJECT_ID

logger = logging.getLogger(__name__)

warnings.filterwarnings("ignore")

# Store a *lower case* version of all annotations
all_annotations = list(map(lambda a: a.lower(), [
    # locations
    "ML",
    "PCL",
    "EGL",
    "EGL boundary",
    "IGL",
    # cilium
    "pocket cilium",
    "concealed cilium",
    "incomplete cilium",
    "pre-ciliary structure",
    "surface cilium",
    "cilium",
    # cell type
    "Bergmann glia",
    "unsure cell type",
    "Granule Cell",
    "Purkinge Cell",
    # cell cycle
    "mitotic",
    "metaphase/anaphase",
    "telophase",
    "cytokinesis",
    "prophase",
    "prometaphase",
    "S/G2",
    # migrating
    "migrating",
    # centriole
    "docked centriole",
    "duplicating centriole",
    "Mother Centriole",
    "Daughter Centriole",
    "centriole",
    "tethered centriole",
    # other
    "bipolar",
    "midbody",
]))


def normann(a):
    """Normalize an annotation by converting it to lower case."""
    return a.lower()


def get_cell(cell_objects, cell):
    return cell


def get_basal_body(cell_objects, cell):
    basal_body = False
    for (object_id, annotations, name) in cell_objects:
        basal_body = basal_body or ("basal body" in annotations)
    return basal_body


def get_distance_a(cell_objects, cell):
    mother_id = None
    daughter_id = None
    for (object_id, annotations, name) in cell_objects:
        if (
            "mother centriole" in annotations
            and "b" not in annotations
            and "b" not in name
        ):
            assert mother_id is None
            mother_id = object_id
        elif (
            "daughter centriole" in annotations
            and "b" not in annotations
            and "b" not in name
        ):
            assert (
                daughter_id is None
            ), f"{daughter_id}, {object_id}, {name}, {annotations}"
            daughter_id = object_id

    if mother_id is None or daughter_id is None:
        logger.info(
            f"Cell {cell}, Centriole pair A, has mother_id: "
            f"{mother_id} and daughter_id: {daughter_id}"
        )
        return np.nan

    mother_skeleton = requests.post(
        f"{CATMAID_URL}/{PROJECT_ID}/{mother_id}/0/1/compact-skeleton",
        verify=False,
        auth=AUTH_TOKEN,
        data="",
    ).json()
    daughter_skeleton = requests.post(
        f"{CATMAID_URL}/{PROJECT_ID}/{daughter_id}/0/1/compact-skeleton",
        verify=False,
        auth=AUTH_TOKEN,
        data="",
    ).json()
    a_nodes = mother_skeleton[0]
    b_nodes = daughter_skeleton[0]
    assert len(a_nodes) == 1
    assert len(b_nodes) == 1
    a_loc = np.array(a_nodes[0][3:6])
    b_loc = np.array(b_nodes[0][3:6])
    return np.linalg.norm(a_loc - b_loc)


def get_centriole_coords(cell_objects, cell):
    centriole_locations = []
    for (object_id, annotations, name) in cell_objects:
        if "centriole" in annotations:
            centriole_id = object_id

            centriole_skeleton = requests.post(
                f"{CATMAID_URL}/{PROJECT_ID}/{centriole_id}/0/1/compact-skeleton",
                verify=False,
                auth=AUTH_TOKEN,
                data="",
            ).json()
            nodes = centriole_skeleton[0]
            loc = np.array(nodes[0][3:6])
            centriole_locations.append([float(x) for x in loc[::-1]] + [centriole_id])

    return centriole_locations


def get_depth_a(cell_objects, cell):
    mother_id = None
    for (object_id, annotations, name) in cell_objects:
        if (
            "mother centriole" in annotations
            and "b" not in annotations
            and "b" not in name
        ):
            assert mother_id is None
            mother_id = object_id

    if mother_id is None:
        logger.info(f"Cell {cell} has no mother centriole A")
        return np.nan

    mother_skeleton = requests.post(
        f"{CATMAID_URL}/{PROJECT_ID}/{mother_id}/0/1/compact-skeleton",
        verify=False,
        auth=AUTH_TOKEN,
        data="",
    ).json()
    a_nodes = mother_skeleton[0]
    a_loc = np.array(a_nodes[0][3:6])
    return a_loc[1]


def get_distance_b(cell_objects, cell):
    mother_id = None
    daughter_id = None
    for (object_id, annotations, name) in cell_objects:
        if "mother centriole" in annotations and ("b" in annotations or " b " in name):
            assert (
                mother_id is None
            ), f"Found mother centriole b with id: {object_id}, but previously found with id: {mother_id}.\nObjects in cell: {cell_objects}"
            mother_id = object_id
        elif "daughter centriole" in annotations and (
            "b" in annotations or " b " in name
        ):
            assert daughter_id is None
            daughter_id = object_id

    if mother_id is None or daughter_id is None:
        if mother_id is not None or daughter_id is not None:
            logger.info(
                f"Cell {cell}, centriole pair B, has mother_id: "
                f"{mother_id} and daughter_id: {daughter_id}"
            )
        return np.nan

    mother_skeleton = requests.post(
        f"{CATMAID_URL}/{PROJECT_ID}/{mother_id}/0/1/compact-skeleton",
        verify=False,
        auth=AUTH_TOKEN,
        data="",
    ).json()
    daughter_skeleton = requests.post(
        f"{CATMAID_URL}/{PROJECT_ID}/{daughter_id}/0/1/compact-skeleton",
        verify=False,
        auth=AUTH_TOKEN,
        data="",
    ).json()
    a_nodes = mother_skeleton[0]
    b_nodes = daughter_skeleton[0]
    assert len(a_nodes) == 1
    assert len(b_nodes) == 1
    a_loc = np.array(a_nodes[0][3:6])
    b_loc = np.array(b_nodes[0][3:6])
    return np.linalg.norm(a_loc - b_loc)


def get_depth_b(cell_objects, cell):
    mother_id = None
    for (object_id, annotations, name) in cell_objects:
        if "mother centriole" in annotations and ("b" in annotations or " b " in name):
            assert mother_id is None
            mother_id = object_id

    if mother_id is None:
        return np.nan

    mother_skeleton = requests.post(
        f"{CATMAID_URL}/{PROJECT_ID}/{mother_id}/0/1/compact-skeleton",
        verify=False,
        auth=AUTH_TOKEN,
        data="",
    ).json()
    a_nodes = mother_skeleton[0]
    a_loc = np.array(a_nodes[0][3:6])
    return a_loc[1]


def get_distance_ab(cell_objects, cell):
    mother_id_a = None
    mother_id_b = None
    for (object_id, annotations, name) in cell_objects:
        if "mother centriole" in annotations and ("b" in annotations or " b " in name):
            assert mother_id_b is None
            mother_id_b = object_id
        elif (
            "mother centriole" in annotations
            and "b" not in annotations
            and "b" not in name
        ):
            assert mother_id_a is None
            mother_id_a = object_id

    if mother_id_a is None or mother_id_b is None:
        return np.nan

    mother_skeleton_a = requests.post(
        f"{CATMAID_URL}/{PROJECT_ID}/{mother_id_a}/0/1/compact-skeleton",
        verify=False,
        auth=AUTH_TOKEN,
        data="",
    ).json()
    mother_skeleton_b = requests.post(
        f"{CATMAID_URL}/{PROJECT_ID}/{mother_id_b}/0/1/compact-skeleton",
        verify=False,
        auth=AUTH_TOKEN,
        data="",
    ).json()
    a_nodes = mother_skeleton_a[0]
    b_nodes = mother_skeleton_b[0]
    a_loc = np.array(a_nodes[0][3:6])
    b_loc = np.array(b_nodes[0][3:6])
    return np.linalg.norm(a_loc - b_loc)


def get_cilia_length(cell_objects, cell):
    cilia_length = 0
    for (object_id, annotations, name) in cell_objects:
        if "cilium" in annotations:
            assert cilia_length == 0
            skeleton = requests.post(
                f"{CATMAID_URL}/{PROJECT_ID}/{object_id}/0/1/compact-skeleton",
                verify=False,
                auth=AUTH_TOKEN,
                data="",
            ).json()
            nodes = skeleton[0]
            graph = nx.Graph()
            for node in nodes:
                loc = np.array(node[3:6])
                graph.add_node(node[0], location=loc)
                if node[1] is not None:
                    graph.add_edge(node[0], node[1])
            for u, v in graph.edges:
                cilia_length += np.linalg.norm(
                    graph.nodes[u]["location"] - graph.nodes[v]["location"]
                )
    return cilia_length


def get_centriole(cell_objects, cell):
    centriole_type = False
    for (object_id, annotations, name) in cell_objects:
        if "tethered centriole" in annotations:
            assert (
                centriole_type is False or centriole_type == "tethered"
            ), centriole_type
            centriole_type = "tethered"
        elif "docked centriole" in annotations:
            assert centriole_type is False or centriole_type == "docked", centriole_type
            centriole_type = "docked"
    return centriole_type


def get_cilium(cell_objects, cell):
    cilium_type = False
    for (object_id, annotations, name) in cell_objects:
        if "pre-ciliary structure" in annotations:
            assert cilium_type is False
            cilium_type = "pre-ciliary"
        elif "incomplete cilium" in annotations:
            assert cilium_type is False
            cilium_type = "incomplete"
        elif "concealed cilium" in annotations:
            assert cilium_type is False
            cilium_type = "concealed"
        elif "pocket cilium" in annotations:
            assert cilium_type is False
            cilium_type = "pocket"
        elif "surface cilium" in annotations:
            assert cilium_type is False
            cilium_type = "surface"

    return cilium_type


def get_cell_type(cell_objects, cell):
    cell_type = "unsure"
    for (object_id, annotations, name) in cell_objects:
        if "granule cell" in annotations:
            assert cell_type == "unsure" or cell_type == "granule cell", cell_type
            cell_type = "granule cell"
        elif "bergmann glia" in annotations:
            assert cell_type == "unsure" or cell_type == "bergmann glia", cell_type
            cell_type = "bergmann glia"
        elif "purkinge cell" in annotations:
            assert cell_type == "unsure" or cell_type == "purkinge cell", cell_type
            cell_type = "purkinge cell"

    return cell_type


def get_location(cell_objects, cell):
    location = None
    for (object_id, annotations, name) in cell_objects:
        if "ml" in annotations:
            assert (
                location is None or location == "ml"
            ), f"ml in cell annotations, but already found {location}"
            location = "ml"
        elif "egL" in annotations:
            assert (
                location is None or location == "egl"
            ), f"egl in cell annotations, but already found {location}"
            location = "egl"
        elif "igl" in annotations:
            assert (
                location is None or location == "igl"
            ), f"igl in cell annotations, but already found {location}"
            location = "igl"
        elif "pcl" in annotations:
            assert (
                location is None or location == "pcl"
            ), f"pcl in cell annotations, but already found {location}"
            location = "pcl"

    return location


def get_cell_cycle(cell_objects, cell):
    cell_cycle = False
    cell_stage = False
    for (object_id, annotations, name) in cell_objects:
        if "mitotic" in annotations:
            assert cell_stage is False or cell_stage == "mitotic", cell_stage
            cell_stage = "mitotic"
        elif "s/g2" in annotations:
            assert cell_stage is False or cell_stage == "s/g2", cell_stage
            cell_stage = "s/g2"

        if "prophase" in annotations:
            assert cell_cycle is False or cell_cycle == "prophase"
            cell_cycle = "prophase"
        elif "prometaphase" in annotations:
            assert cell_cycle is False or cell_cycle == "prometaphase"
            cell_cycle = "prometaphase"
        elif "telophase" in annotations:
            assert cell_cycle is False or cell_cycle == "telophase"
            cell_cycle = "telophase"
        elif "metaphase/anaphase" in annotations:
            assert cell_cycle is False or cell_cycle == "metaphase/anaphase"
            cell_cycle = "metaphase/anaphase"
        elif "cytokinesis" in annotations:
            assert cell_cycle is False or cell_cycle == "cytokinesis"
            cell_cycle = "cytokinesis"

    return cell_cycle


def get_migrating(cell_objects, cell):
    migrating = False
    for (object_id, annotations, name) in cell_objects:
        if "migrating" in annotations:
            migrating = True

    return migrating
