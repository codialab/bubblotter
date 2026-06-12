#!/usr/bin/env python3
import argparse
import re
import subprocess
import os
import json

def canonicalize(edge: tuple[str, bool, str, bool]) -> tuple[str, bool, str, bool]:
    if edge[0] < edge[2]:
        return edge
    if edge[0] == edge[2] and edge[1]:
        return edge
    return (edge[2], not edge[3], edge[0], not edge[1])

def edge_to_str(edge: tuple[str, bool, str, bool]) -> str:
    o1 = "+" if edge[1] else "-"
    o2 = "+" if edge[3] else "-"
    return f"L\t{edge[0]}\t{o1}\t{edge[2]}\t{o2}\t0M"

def simplify_file(filename: str) -> str:
    lines: list[str] = []
    new_s_lines: list[str] = []
    old_s_lines: list[tuple[str, str]] = []
    replacements: list[tuple[str, list[str], tuple[str, bool]]] = []
    induced_l_lines: set[tuple[str, bool, str, bool]] = set()
    used_nodes: set[str] = set()

    bubble_chains = f"{filename}_chains.json"
    subprocess.run(["BubbleGun", "-g", filename,
                    "bchains", "--only_simple", "--out_haplos", "--bubble_json", bubble_chains ], stdout=subprocess.DEVNULL)
    if os.path.exists(bubble_chains):
        chains = json.load(open(bubble_chains))
        replacements = [(key, chain['ends'],
                        (f"C{key}", False)) for (key, chain) in chains.items()]
        # Add a length for all replacements
        with open("haplotype1.fasta") as f:
            keys = []
            seqs = []
            current_seq = ""
            for line in f:
                if line[0] == ">":
                    match = re.match(r">chain_(\d+)_hap1", line)
                    if match is None:
                        print("COULD NOT PARSE BUBBLEGUN FASTA OUTPUT")
                        return ""
                    chain_id = match.group(1)
                    if current_seq != "":
                        seqs.append(current_seq)
                    keys.append(chain_id)
                else:
                    current_seq += line.strip()
            if current_seq != "":
                seqs.append(current_seq)

            for key, seq in zip(keys, seqs):
                node_id = f"C{key}"
                line = f"S\t{node_id}\t{seq}"
                new_s_lines.append(line)

    with open(filename) as f:
        for line in f:
            marker: str = line[0]
            if marker == 'P':
                fields = line.strip().split("\t")
                haplotype = list(map(lambda x: (x[:-1], True) if x[-1] == '+' else (x[:-1], False),
                                     fields[2].split(",")))
                haplotype = apply_complex_replacements(haplotype,
                                                       replacements)
                haplotype_text = ','.join(map(lambda x: f"{x[0]}+" if x[1] else f"{x[0]}-", haplotype))
                line = "\t".join([fields[0], fields[1], haplotype_text, *fields[3:]])
                lines.append(line)

                used_nodes = used_nodes.union(set(map(lambda x: x[0], haplotype)))
                for i in range(0, len(haplotype) - 1):
                    edge = (haplotype[i][0], haplotype[i][1], haplotype[i + 1][0], haplotype[i + 1][1])
                    c_edge = canonicalize(edge)
                    induced_l_lines.add(c_edge)
            elif marker == 'W':
                fields = line.strip().split("\t")
                nodes = re.split('([><][^><]*)', fields[6])[1::2]
                haplotype = [(node[1:],
                              True) if node[0] == '>' else (node[1:],
                                                            False) for node in nodes]
                haplotype = apply_complex_replacements(haplotype,
                                                       replacements)
                haplotype_text = ''.join(map(lambda x: f">{x[0]}" if x[1] else f"<{x[0]}", haplotype))
                line = "\t".join([*fields[:6], haplotype_text])
                lines.append(line)

                used_nodes = used_nodes.union(set(map(lambda x: x[0], haplotype)))
                for i in range(0, len(haplotype) - 1):
                    edge = (haplotype[i][0], haplotype[i][1], haplotype[i + 1][0], haplotype[i + 1][1])
                    c_edge = canonicalize(edge)
                    induced_l_lines.add(c_edge)
            # Skip L lines, we will use induced ones instead
            elif marker == 'L':
                continue
            elif marker == 'S':
                fields = line.strip().split()
                old_s_lines.append((fields[1], fields[2]))
            else:
                lines.append(line.strip())

    filtered_s_lines = [f"S\t{s[0]}\t{s[1]}" for s in old_s_lines if s[0] in used_nodes]
    new_l_lines = [edge_to_str(edge) for edge in induced_l_lines]
    all_lines = lines + filtered_s_lines + new_s_lines + new_l_lines

    return "\n".join(all_lines)


def apply_complex_replacements(haplotype: list[tuple[str, bool]], replacement_rules: list[tuple[str, list[str], tuple[str, bool]]]) -> list[tuple[str, bool]]:
    operations = []

    for idx, (_chain_id, borders, (new_value, already_ordered)) in enumerate(replacement_rules):
        try:
            idx1 = next(i for i, val in enumerate(haplotype) if val[0] == borders[0])
            idx2 = next(i for i, val in enumerate(haplotype) if val[0] == borders[1])
            start_idx = min(idx1, idx2)
            end_idx = max(idx1, idx2)
            if start_idx == idx1:
                if not already_ordered:
                    replacement_rules[idx] = (_chain_id, borders, (new_value, True))
                operations.append((start_idx, end_idx, (new_value, True)))
            else:
                if not already_ordered:
                    replacement_rules[idx] = (_chain_id, list(reversed(borders)), (new_value, True))
                    operations.append((start_idx, end_idx, (new_value, True)))
                else:
                    operations.append((start_idx, end_idx, (new_value, False)))
        except StopIteration:
            continue

    operations.sort(key=lambda x: x[0], reverse=True)
    new_haplotype = haplotype.copy()
    for start_idx, end_idx, new_value in operations:
        new_haplotype[start_idx:end_idx + 1] = [new_value]
    return new_haplotype


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("filename")

    args = parser.parse_args()
    filename = args.filename
    result = simplify_file(filename)
    print(result)

if __name__ == "__main__":
    main()
