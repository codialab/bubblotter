# bubblotter
Simple tool for plotting bubbles within a GFA file.

## Usage
Requirements:
- python3
  - numpy
  - scipy
  - matplotlib
- [vg](https://github.com/vgteam/vg)

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


* non-simple bubbles are by-default all bubbles that are not:
- bubbles containing 4 or less nodes, only of length 1 (likely some form of SNP or 1bp INS/DEL)
- bubbles containing a single insertion or deletion (can be turned back on using the parameter `--include_ins`)
