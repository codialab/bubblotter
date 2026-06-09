#!/usr/bin/env python3
import logging
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
import argparse
import numpy as np
import matplotlib.colors as mcolors
import re
import json
import random
from scipy.cluster.hierarchy import linkage, leaves_list
from collections import Counter


def get_core_nodes(haplotypes, threshold=0.9):
    total_haps = len(haplotypes)
    min_required = total_haps * threshold

    # 1. Count occurrences
    # We use set(hap) so that if a node appears twice in the *same* haplotype,
    # it only counts as 1 presence for that haplotype.
    node_counts = Counter(node for hap in haplotypes for node in set(hap))

    # 2. Filter based on the threshold
    core_nodes = [node for node, count in node_counts.items() if count > min_required]

    return core_nodes


def calculate_distance(hap1, hap2):
    """
    Custom distance function.
    Calculates how many elements differ between two haplotypes.
    """
    # Count mismatches for overlapping parts
    mismatches = sum(1 for a, b in zip(hap1, hap2) if a != b)

    # Add a penalty for length differences (if your haplotypes vary in length)
    mismatches += abs(len(hap1) - len(hap2))

    return mismatches


def sort_haplotypes_by_similarity(haplotypes):
    n = len(haplotypes)

    # Edge case: 1 or 0 items
    if n <= 1:
        return haplotypes

    # 1. Build a condensed distance matrix
    # SciPy expects a flat 1D array of the upper triangle of the
    # distance matrix
    distances = []
    for i in range(n):
        for j in range(i + 1, n):
            distances.append(calculate_distance(haplotypes[i], haplotypes[j]))

    distances = np.array(distances)
    # If all haplotypes are identical, distances will be all 0s.
    # Linkage fails on this, so just return the original list.
    if np.all(distances == 0):
        return haplotypes

    # 2. Perform Hierarchical Clustering
    # 'average' linkage is highly stable for sequence data.
    # optimal_ordering=True is the magic step: it rotates the branches of the
    # clustering tree to ensure adjacent items in the final list are as
    # similar as possible.
    Z = linkage(distances, method='average', optimal_ordering=True)

    # 3. Extract the sorted order of indices
    sorted_indices = leaves_list(Z)

    # 4. Reorder the original list based on the new indices
    sorted_haplotypes = [haplotypes[i] for i in sorted_indices]

    return sorted_haplotypes


def add_custom_arrow(ax, start, end, y, height=0.4, head_length=1.0,
                     color="steelblue", label=None):
    """
    Draws a gggenes-style arrow on a Matplotlib axis.
    """
    # Determine direction based on start and end coordinates
    direction = 1 if end > start else -1

    # Ensure the head length isn't longer than the arrow itself
    actual_head_length = min(head_length, abs(end - start)) * direction

    tail_end = end - actual_head_length
    half_h = height / 2

    # Define the 5 points of the arrow polygon
    vertices = [
        (start, y - half_h),      # Bottom start
        (tail_end, y - half_h),   # Bottom of head base
        (end, y),                 # Arrow tip
        (tail_end, y + half_h),   # Top of head base
        (start, y + half_h)       # Top start
    ]
    # Create and add the patch
    poly = Polygon(vertices, facecolor=color, edgecolor='black', linewidth=1,
                   zorder=2)
    ax.add_patch(poly)

    # Optional: Add text inside or above the arrow
    if label:
        ax.text((start + tail_end) / 2, y, label, ha='center', va='center',
                color='white', fontweight='bold', zorder=3)


def reverse(sequence):
    return list(map(lambda x: (x[0], not x[1]), sequence[::-1]))


def rev_node(node):
    return (node[0], not node[1])


def apply_complex_replacements(haplotype, replacement_rules):
    operations = []

    # 1. Find the boundaries using conditional matching
    for borders, new_value in replacement_rules:
        try:
            # Find index where the first element of the tuple matches the border
            idx1 = next(i for i, val in enumerate(haplotype) if val[0] == borders[0])
            idx2 = next(i for i, val in enumerate(haplotype) if val[0] == borders[1])

            start_idx = min(idx1, idx2)
            end_idx = max(idx1, idx2)

            if start_idx == idx1:
                operations.append((start_idx, end_idx, new_value))
            else:
                operations.append((start_idx, end_idx, (new_value[0],
                                                        not new_value[1])))

        except StopIteration:
            # next() raises StopIteration if it reaches the end without finding a match.
            # This safely skips rules where the borders aren't found.
            continue

    # 2. Sort Right-to-Left
    operations.sort(key=lambda x: x[0], reverse=True)

    # 3. Apply the replacements
    new_haplotype = haplotype.copy()

    for start_idx, end_idx, new_value in operations:
        new_haplotype[start_idx:end_idx + 1] = [new_value]

    return new_haplotype


def main():
    log = logging.getLogger(__name__)
    parser = argparse.ArgumentParser()
    parser.add_argument("filename")
    parser.add_argument("bubble_chains")

    args = parser.parse_args()
    filename = args.filename
    bubble_chains = args.bubble_chains

    node_lengths = {}
    haplotypes = []
    haplotype_to_idx = {}
    start = set()
    end = set()
    anchors = set()
    first = True

    chains = json.load(open(bubble_chains))
    replacements = [(chain['ends'], (f"C{key}", True)) for (key, chain) in chains.items()]
    # Add a length for all replacements
    for ends, replacement in replacements:
        node_lengths[replacement[0]] = 10

    with open(filename) as f:
        for line in f:
            if not line:
                continue
            marker = line[0]

            if marker == 'S':
                fields = line.strip().split("\t")
                name = fields[1]
                length = len(fields[2])
                node_lengths[name] = length
            elif marker == 'P' or marker == 'W':
                if marker == 'P':
                    fields = line.strip().split("\t")
                    haplotype = list(map(lambda x: (x[:-1], True) if x[-1] == '+' else (x[:-1], False),fields[2].split(",")))
                elif marker == 'W':
                    fields = line.strip().split("\t")
                    nodes = re.split('([><][^><]*)', fields[6])[1::2]
                    haplotype = [(node[1:], True) if node[0] == '>' else (node[1:], False) for node in nodes]

                haplotype = apply_complex_replacements(haplotype, replacements)
                # Sequence is backward
                if rev_node(haplotype[-1]) in start or rev_node(haplotype[0]) in end:
                    haplotype = reverse(haplotype)
                elif haplotype[0] not in start and haplotype[-1] not in end:
                    log.warning(f"Could not orient haplotype: {fields[1]}")

                nodes_of_haplotype = set(haplotype)
                if first:
                    anchors = nodes_of_haplotype
                    first = False
                else:
                    anchors = set.intersection(anchors, nodes_of_haplotype)
                if haplotype[0] not in start:
                    start.add(haplotype[0])
                if haplotype[-1] not in end:
                    end.add(haplotype[-1])
                haplotype_key = tuple(haplotype)
                if haplotype_key not in haplotype_to_idx:
                    haplotypes.append(haplotype)
                    haplotype_to_idx[haplotype_key] = len(haplotypes) - 1

    unique_ids = list(node_lengths.keys())
    random.shuffle(unique_ids)
    num_ids = len(unique_ids)
    cmap = plt.get_cmap('nipy_spectral')
    sampled_colors = cmap(np.linspace(0, 1, num_ids))
    color_dict = {id_name: mcolors.to_hex(sampled_colors[i]) for i, id_name in enumerate(unique_ids)}

    sorted_haps = sort_haplotypes_by_similarity(haplotypes)

    print(get_core_nodes(haplotypes, threshold=0.95))

    _, ax = plt.subplots(figsize=(10, 10))

    max_x = 0
    for h_idx, haplotype in enumerate(sorted_haps):
        pos = 0
        for node, orientation in haplotype:
            next_pos = pos + node_lengths[node]
            if orientation:
                start = pos
                end = next_pos
            else:
                start = next_pos
                end = pos
            if end > max_x:
                max_x = end
            add_custom_arrow(ax, start, end, h_idx + 1, color=color_dict[node], label=None)
            pos = next_pos

    ax.set_xlabel("Genomic Position")
    ax.set_xlim(0, max_x)
    ax.set_ylim(len(haplotypes) + 1, 0)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
