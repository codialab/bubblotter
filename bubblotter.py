#!/usr/bin/env python3
import logging
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, Rectangle
import matplotlib.image as mpimg
from matplotlib.gridspec import GridSpec
import argparse
import numpy as np
import matplotlib.colors as mcolors
import re
import json
import random
import os
from scipy.cluster.hierarchy import linkage, leaves_list
import subprocess
from pathlib import Path
import colorsys


log = logging.getLogger(__name__)


def hsv_to_hex(h, s, v):
    r, g, b = (int(255 * component) for component in colorsys.hsv_to_rgb(h / 255, s / 255, v / 255))
    return f"#{r:02x}{g:02x}{b:02x}"


def calculate_distance(haplotype_1, haplotype_2):
    mismatches = sum(1 for a, b in zip(haplotype_1, haplotype_2) if a != b)
    mismatches += abs(len(haplotype_1) - len(haplotype_2))
    return mismatches


def sort_haplotypes_by_similarity(haplotypes):
    n = len(haplotypes)
    if n <= 1:
        return haplotypes
    distances = []
    for i in range(n):
        for j in range(i + 1, n):
            distances.append(calculate_distance(haplotypes[i][1],
                                                haplotypes[j][1]))
    distances = np.array(distances)
    if np.all(distances == 0):
        return haplotypes
    linkages = linkage(distances, method='average', optimal_ordering=True)
    sorted_indices = leaves_list(linkages)
    sorted_haplotypes = [haplotypes[i] for i in sorted_indices]
    return sorted_haplotypes


def sort_with_reference(haplotypes, reference_index=0):
    if len(haplotypes) <= 1:
        return haplotypes
    if len(haplotypes) <= 2:
        ref = haplotypes.pop(reference_index)
        return [ref] + haplotypes

    reference = haplotypes[reference_index]

    others = [h for i, h in enumerate(haplotypes) if i != reference_index]

    sorted_others = sort_haplotypes_by_similarity(others)

    # Calculate the distances of the reference group to the ends of
    # the ordering
    dist_to_start = calculate_distance(reference, sorted_others[0])
    dist_to_end = calculate_distance(reference, sorted_others[-1])

    # If the end is closer, flip the ordering
    if dist_to_end < dist_to_start:
        sorted_others.reverse()

    return [reference] + sorted_others


def add_custom_arrow(ax, start, end, y, height=0.4, head_length=1.0,
                     color="black", label=None):
    direction = 1 if end > start else -1

    actual_head_length = min(head_length, abs(end - start)) * direction

    tail_end = end - actual_head_length
    half_h = height / 2

    vertices = [
        (start, y - half_h),
        (tail_end, y - half_h),
        (end, y),
        (tail_end, y + half_h),
        (start, y + half_h)
    ]
    poly = Polygon(vertices, facecolor=color, edgecolor='black', linewidth=1,
                   zorder=2)
    ax.add_patch(poly)

    if label:
        ax.text((start + tail_end) / 2, y, label, ha='center', va='center',
                color='white', fontweight='bold', zorder=3)


def reverse(sequence):
    return list(map(lambda x: (x[0], not x[1]), sequence[::-1]))


def rev_node(node):
    return (node[0], not node[1])


def apply_complex_replacements(haplotype, replacement_rules):
    operations = []

    for borders, new_value in replacement_rules:
        try:
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
            continue

    operations.sort(key=lambda x: x[0], reverse=True)
    new_haplotype = haplotype.copy()
    for start_idx, end_idx, new_value in operations:
        new_haplotype[start_idx:end_idx + 1] = [new_value]
    return new_haplotype


def get_haplotypes(filename, bubble_chains=None):
    node_lengths = {}
    haplotypes = []
    haplotype_to_idx = {}
    start = set()
    end = set()

    if bubble_chains is not None:
        chains = json.load(open(bubble_chains))
        replacements = [(chain['ends'],
                         (f"C{key}",
                          True)) for (key, chain) in chains.items()]
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
                    haplotype = list(map(lambda x: (x[:-1], True) if x[-1] == '+' else (x[:-1], False),
                                         fields[2].split(",")))
                    name = fields[1]
                elif marker == 'W':
                    fields = line.strip().split("\t")
                    nodes = re.split('([><][^><]*)', fields[6])[1::2]
                    haplotype = [(node[1:],
                                  True) if node[0] == '>' else (node[1:],
                                                                False) for node in nodes]
                    name = fields[1]

                if bubble_chains is not None:
                    haplotype = apply_complex_replacements(haplotype,
                                                           replacements)
                # Sequence is backward
                if rev_node(haplotype[-1]) in start or rev_node(haplotype[0]) in end:
                    haplotype = reverse(haplotype)
                elif haplotype[0] not in start and haplotype[-1] not in end:
                    log.warning(f"Could not orient haplotype: {fields[1]}")

                if haplotype[0] not in start:
                    start.add(haplotype[0])
                if haplotype[-1] not in end:
                    end.add(haplotype[-1])
                haplotype_key = tuple(haplotype)
                if haplotype_key not in haplotype_to_idx:
                    haplotypes.append(([name], haplotype))
                    haplotype_to_idx[haplotype_key] = len(haplotypes) - 1
                else:
                    haplotype_idx = haplotype_to_idx[haplotype_key]
                    haplotypes[haplotype_idx][0].append(name)
    return (haplotypes, node_lengths)


def draw_plot(filename, annotations, full_coords,
              start_node, end_node,
              reference="", reference_start=0, reference_end=0,
              output_filename=None, plot_bandage=True):
    if output_filename is None:
        output_filename = f"{filename}.png"
    (haplotypes, node_lengths) = get_haplotypes(filename)

    if len(haplotypes) == 0:
        log.warning(f"No haplotypes for bubble {full_coords}")
        return

    unique_ids = list(node_lengths.keys())
    random.shuffle(unique_ids)
    num_ids = len(unique_ids)
    sampled_colors = []
    for h in np.linspace(0, 1, num_ids):
        hex_color = hsv_to_hex(int(h * 255), 130, 200)
        sampled_colors.append(hex_color)
    color_dict = {id_name: mcolors.to_hex(sampled_colors[i]) for i, id_name in enumerate(unique_ids)}
    color_dict[str(start_node)] = "#00ee00"
    color_dict[str(end_node)] = "#ee0000"

    if plot_bandage:
        colors_file = f"{filename}_colors.csv"
        with open(f"{filename}_colors.csv", "w") as f:
            f.write("node,color\n")
            for key, value in color_dict.items():
                f.write(f"{key},{value}\n")
        bandage_image = f"{filename}_bandage.png"
        subprocess.run(["Bandage", "image", filename, bandage_image,
                        "--colors", colors_file,
                        "--width", "1000", "--height", "1000"])

    reference_idx = next((h_id for h_id, (names, _haplo) in enumerate(haplotypes) if reference is not None and reference != "" and any(s.startswith(reference) for s in names)), 0)
    sorted_haps = sort_with_reference(haplotypes,
                                      reference_index=reference_idx)

    if not plot_bandage:
        fig, ax = plt.subplots(figsize=(10, 10))
    else:
        fig = plt.figure(figsize=(20, 12))

        gs = GridSpec(12, 20, figure=fig)

        ax = fig.add_subplot(gs[:10, :10])
        ax1 = fig.add_subplot(gs[:10, 10:])
        ax2 = fig.add_subplot(gs[-2:, :])

        ax2.plot(full_coords, np.zeros(len(full_coords)), c="black")
        ax2.set_ylim(bottom=-1, top=1)
        for annotation in annotations:
            order = 0
            color = "#000000"
            if annotation[2] == "gene" or annotation[2] == "ncRNA_gene":
                order = 0
                color = "#cccccc"
            elif annotation[2] == "exon":
                order = 1
                color = "#555555"
            else:
                continue
            start = int(annotation[3])
            end = int(annotation[4])
            length = end - start
            ax2.add_patch(Rectangle((start, -0.5),
                                    length,
                                    1, color=color, zorder=order))
        ax2.add_patch(Rectangle((reference_start, -0.5),
                                reference_end - reference_start,
                                1, color='r', zorder=5))
        ax2.spines['top'].set_visible(False)
        ax2.spines['right'].set_visible(False)
        ax2.spines['left'].set_visible(False)
        ax2.yaxis.set_visible(False)
        ax2.set_xlim(full_coords[0], full_coords[1])

        img = mpimg.imread(bandage_image)
        ax1.imshow(img)
        ax1.axis('off')

    max_x = 0
    for h_idx, (_names, haplotype) in enumerate(sorted_haps):
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
            add_custom_arrow(ax, start, end, h_idx + 1, color=color_dict[node],
                             label=node)
            pos = next_pos

    ax.set_xlabel("Genomic Position")
    ax.set_xlim(0, max_x)
    ax.set_ylim(len(sorted_haps) + 1, 0)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    # ax.yaxis.set_visible(False)

    y_positions = list(range(1, len(sorted_haps) + 1))
    ax.set_yticks(y_positions)
    multi_line_labels = [
        f"Group {h_id}*:\n{len(names)} Haplotypes" if reference is not None and reference != "" and any(s.startswith(reference) for s in names) else f"Group {h_id}:\n{len(names)} Haplotypes"
        for h_id, (names, _haplo) in enumerate(sorted_haps)
    ]
    ax.set_yticklabels(multi_line_labels)
    ax.tick_params(axis='y', length=0, pad=15)
    fig.suptitle(f"{reference}:{reference_start}-{reference_end}")

    for tick in ax.get_yticklabels():
        tick.set_horizontalalignment('right')

    fig.text(
        x=0.01,
        y=0.01,
        s="*: this group contains the reference",
        ha="left",
        va="bottom",
        fontsize=9,
        style="italic",
        color="dimgray"
    )

    plt.tight_layout()
    # plt.show()
    plt.savefig(output_filename)
    plt.close(fig)


def get_reference_coords(filename, reference):
    if reference is None:
        return (0, 0)
    reference_start = 0
    reference_end = 0
    with open(filename) as f:
        for line in f:
            marker = line[0]
            if marker == 'P':
                fields = line.strip().split()
                if fields[1].startswith(reference):
                    naming_parts = fields[1].split("#")
                    coords = naming_parts[-1].split(":")
                    if len(coords) != 2:
                        log.error("ERROR: path does not contain coordinates")
                    else:
                        reference_start = int(coords[0])
                        reference_end = int(coords[1])
                        break
            elif marker == 'W':
                fields = line.strip().split()
                name = "#".join(fields[1:4])
                if name.startswith(reference):
                    reference_start = int(fields[4])
                    reference_end = int(fields[5])
                    break
    return (reference_start, reference_end)


def get_maximum_node_id(filename):
    maximum_id = 0
    with open(filename) as f:
        for line in f:
            marker = line[0]
            if marker == 'S':
                fields = line.strip().split()
                id = int(fields[1])
                if id > maximum_id:
                    maximum_id = id
    return maximum_id


def main():
    random.seed(1)

    parser = argparse.ArgumentParser()
    parser.add_argument("filename")
    parser.add_argument("work_dir")
    parser.add_argument("-r", "--reference",
                        help="Reference that is used to assign positions to the output, please use the PanSN spec format if the sample name is not unique")
    parser.add_argument("-a", "--reference_annotation",
                        help="GFF3 annotation of the reference, used to mark positions")
    parser.add_argument("-s", "--reference_seq_id",
                        help="seqid in the GFF3 annotation of the reference")
    parser.add_argument("-i", "--include_ins", action="store_true", help="If this is set simple insertions and deletions (so ins/dels consisting of only a single node) will be plotted as well")

    args = parser.parse_args()
    filename = args.filename
    work_dir = args.work_dir
    reference = args.reference
    reference_annotation = args.reference_annotation
    reference_seq_id = args.reference_seq_id
    include_ins = args.include_ins

    reference_start, reference_end = get_reference_coords(filename, reference)

    annotations = []
    if reference_annotation:
        with open(reference_annotation) as f:
            for line in f:
                if line[0] == "#":
                    continue
                fields = line.strip().split()
                chr = fields[0]
                start_pos = int(fields[3])
                end_pos = int(fields[4])
                if (reference_seq_id is None or chr == reference_seq_id) and start_pos <= reference_end and end_pos >= reference_start:
                    annotations.append(fields)

    Path(work_dir).mkdir(parents=True, exist_ok=True)

    # Sort ids
    sorted_vg = os.path.join(work_dir, f"{filename}_sorted.vg")
    with open(sorted_vg, "wb") as outfile:
        subprocess.run(["vg", "ids", "--sort", filename], stdout=outfile)
    sorted_gfa = os.path.join(work_dir, f"{filename}_sorted.gfa")
    with open(sorted_gfa, "w") as outfile:
        subprocess.run(["vg", "view", sorted_vg], stdout=outfile)

    # Call bubble_gun

    bubbles_pb = os.path.join(work_dir, f"{filename}_sorted_bubbles.pb")
    with open(bubbles_pb, "wb") as outfile:
        subprocess.run(["vg", "snarls", "-a", sorted_vg], stdout=outfile)
    bubbles_json = os.path.join(work_dir, f"{filename}_sorted_bubbles.json")
    with open(bubbles_json, "w") as outfile:
        subprocess.run(["vg", "view", "-j", "-R", bubbles_pb], stdout=outfile)

    node_lengths = {}
    with open(sorted_gfa) as f:
        for line in f:
            if not line:
                continue
            marker = line[0]
            if marker == 'S':
                fields = line.strip().split()
                node_lengths[fields[1]] = len(fields[2])

    ends = []
    with open(bubbles_json) as f:
        for line in f:
            if not line:
                continue
            bubble = json.loads(line)
            bubble_end = int(bubble["end"]["node_id"])
            bubble_start = int(bubble["start"]["node_id"])
            # Check if we are in a subbubble
            if "parent" in bubble:
                continue
            # Check if this is an insertion or deletion
            if not include_ins and bubble_end - bubble_start <= 2:
                continue
            # Check if this is a small bubble consisting of only 1bp nodes
            if bubble_end - bubble_start <= 5:
                all_are_1bp = all(map(lambda x: node_lengths[str(x)] < 2,
                                      range(bubble_start + 1,
                                            bubble_end)))
                if all_are_1bp:
                    continue
            ends.append([bubble_start, bubble_end])

    sorted_ends = sorted(ends, key=lambda x: x[0])
    maximum_length = len(str(sorted_ends[-1][-1]))
    chunk_files = []
    for i in range(0, len(sorted_ends)):
        start_in_between = sorted_ends[i][0]
        end_in_between = sorted_ends[i][1]
        log.info(f"Handling chunk: {start_in_between}-{end_in_between}")
        sorted_in_between = os.path.join(work_dir, f"{filename}_{start_in_between:0{maximum_length}d}-{end_in_between:0{maximum_length}d}.gfa")
        with open(sorted_in_between, "w") as outfile:
            subprocess.run(["vg", "chunk", "-x", sorted_gfa,
                            "-r", f"{start_in_between}:{end_in_between}",
                            "-c", "0",
                            "-O", "gfa"],
                           stdout=outfile)
        coords = get_reference_coords(sorted_in_between, reference)
        nodes = (start_in_between, end_in_between)
        chunk_files.append((sorted_in_between, coords, nodes))

    for (chunk_file, coords, nodes) in chunk_files:
        start = coords[0]
        end = coords[1]
        # Check if vg messed up path coordinates
        if start < reference_start:
            start += reference_start
            end += reference_start
        draw_plot(chunk_file, annotations, (reference_start, reference_end), nodes[0], nodes[1], reference=reference,
                  reference_start=start, reference_end=end)


if __name__ == "__main__":
    main()
