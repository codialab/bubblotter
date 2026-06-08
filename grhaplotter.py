#!/usr/bin/env python3
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
import argparse
import numpy as np
import matplotlib.colors as mcolors

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


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("filename")

    args = parser.parse_args()
    filename = args.filename

    node_lengths = {}
    haplotypes = []
    haplotype_to_idx = {}
    start = set()
    end = set()

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
            elif marker == 'P':
                fields = line.strip().split("\t")
                haplotype = list(map(lambda x: (x[:-1], True) if x[-1] == '+' else (x[:-1], False),fields[2].split(",")))
                # Sequence is backward
                print(start, end, haplotype[0], haplotype[-1])
                if rev_node(haplotype[-1]) in start or rev_node(haplotype[0]) in end:                    
                    haplotype = reverse(haplotype)
                if haplotype[0] not in start:
                    start.add(haplotype[0])
                if haplotype[-1] not in end:
                    end.add(haplotype[-1])
                haplotype_key = tuple(haplotype)
                if haplotype_key not in haplotype_to_idx:        
                    haplotypes.append(haplotype)
                    haplotype_to_idx[haplotype_key] = len(haplotypes) - 1

    unique_ids = list(node_lengths.keys())
    num_ids = len(unique_ids)
    cmap = plt.get_cmap('nipy_spectral')
    sampled_colors = cmap(np.linspace(0, 1, num_ids))
    color_dict = {id_name: mcolors.to_hex(sampled_colors[i]) for i, id_name in enumerate(unique_ids)}

    _, ax = plt.subplots(figsize=(10, 10))

    max_x = 0
    for h_idx, haplotype in enumerate(haplotypes):
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
            add_custom_arrow(ax, start, end, h_idx + 1, color=color_dict[node], label=node)
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
