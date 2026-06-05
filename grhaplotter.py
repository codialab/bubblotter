#!/usr/bin/env python3
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon


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


# --- Plot Setup ---
fig, ax = plt.subplots(figsize=(10, 3))

arrows_data = [
    {"start": 1, "end": 5, "y": 1, "color": "#4C72B0", "label": "Gene A"},
    {"start": 8, "end": 6, "y": 1, "color": "#C44E52", "label": "Gene B"},
    {"start": 6.5, "end": 9, "y": 2, "color": "#55A868", "label": "Gene C"}
]

for arr in arrows_data:
    add_custom_arrow(ax, **arr)

ax.set_xlim(0, 10)
ax.set_ylim(0, 3)
ax.set_yticks([1, 2])
ax.set_yticklabels(["Track 1", "Track 2"])
ax.set_xlabel("Genomic Position")

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.tight_layout()
plt.show()
