# bubblotter
Simple tool for plotting bubbles within a GFA file.

## Usage
Requirements:
- python3
  - numpy
  - scipy
  - matplotlib
- [vg](https://github.com/vgteam/vg)
- [Bandage](https://github.com/rrwick/Bandage)

Subset the graph to your region of interest, e.g. a single gene. This is important as `bubblotter` will try to draw every non-simple* top-level bubble within the graph. So inputting a large file will result in many images. This can be done for example like:
```bash
vg chunk -x "my_graph.gfa" -p 'GRCh38#0#chr8:143418256-143421518' -c 20 -O gfa > roi.gfa
```
Please ensure that the GFA file contains P- or W-lines.

Now you can run `bubblotter` like:
```bash
./bubblotter.py roi.gfa work -r GRCh38 -a gencode.v50.basic.annotation.gff3 -s "chr8"
```
In this case `roi.gfa` is the region of interest graph, `work` the name of the working directory (`bubblotter` will store all of its results there).
`-r GRCh38` specifies the name of the reference, `-a gencode.v50.basic.annotation.gff3` provides a annotation of genes/exons and `-s "chr8"` tells `bubblotter` to look only for annotations in chromosome 8.

This results in multiple plots like:

![Plot showing a region of the MAFA gene in three different forms: once as a gene-arrow like haplotype visualization, once as a graph and once in a minimap of the whole inputted region with exon/gene annotations](docs/MAFA_c20.gfa_412-451.gfa.png)

Each plot shows on the left side a gene-arrow-like haplotype plot, each line consisting of all the haplotypes (paths) that have the same sequence in this bubble. Each arrow is a single node, the direction indicating the direction of the node traversal. The top-most
line is always the group of haplotypes containing the reference (if a reference is specified). The start and end nodes of the bubble are marked in bright green and red respectively. The right side contains a plot of the graph generated using Bandage, the colors matched
to the left plot. Below these plots is a minimap showing the position of this bubble in the graph given to `bubblotter` in terms of reference coordinates. If an annotation is given, then this minimap also shows the location of genes in light grey and exons in dark grey.

\* non-simple bubbles are by-default all bubbles that are not:
- bubbles containing 4 or less nodes, only of length 1 (likely some form of SNP or 1bp INS/DEL)
- bubbles containing a single insertion or deletion (can be turned back on using the parameter `--include_ins`)
