#!/usr/bin/env python3
import os, sys, argparse, warnings, csv
warnings.filterwarnings('ignore')
import subprocess
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib import colors
import matplotlib_venn
import math
from itertools import chain
from collections import Counter, Iterable
from tqdm import tqdm

def get_args():
    parser = argparse.ArgumentParser(description='dkato. November, 2021')
    parser.add_argument('-rps' , dest ='rps', nargs='*',
                        help = 'path to your results of rpsblast')
    parser.add_argument('-AA' , dest ='AA', nargs='*',
                        help = 'paths　to your amino acids files of genes(Venn diagram is not output if there are 6 or more files)')
    parser.add_argument('-e' , dest ='evalue', nargs='*',
                        default= ['1e-28'],  help = 'evalue in rpsblast(default:1e-28)')
    parser.add_argument('-bar' , dest ='bar_size',
                        default= 5, type = int, help = 'specify a integer value: graph size of bar plot(default:5)')
    parser.add_argument('-B', dest='n_black',
                        default=1,type = int, help = 'Number of bars dyed in black in a bar graph(default:1)')
    parser.add_argument('-PCA' , dest ='PCA_size',
                        default= 5, type = int, help = 'specify a integer value: graph size of PCA plot(default:5)')
    parser.add_argument('-G', dest='n_green',
                        default=0,type = int, help = 'Number of points dyed in green in a PCA plot(default:0)')
    parser.add_argument('-venn' , dest ='venn_size',
                        default= 7, type = int, help = 'specify a integer value: graph size of venn diagrams(default:7)')
    #parser.add_argument('-t', dest='num_threads',
     #                   default=10,type = int, help = 'num_threads(default:10)')        
                        
    parser.add_argument('-cogdb' , dest ='cogdb',
                        default= '/home/tmp/db/COG/Cog', 
                       help = 'path to your cogdb to run rpsblast(default:/home/tmp/db/COG/Cog)')    
    parser.add_argument('-cddid' , dest ='cddid',
                        default= '/home/tmp/db/COG/cdd2cog/cddid_COG.tbl',
                        help = 'path to your cddid_COG.tbl(default:/home/tmp/db/COG/cdd2cog/cddid_COG.tbl)')
    parser.add_argument('-cog', dest='cog',
                        default='/home/tmp/db/COG/cdd2cog/cog-20.def.tsv',
                        help = 'path to your cog-20.def.tsv(default:/home/tmp/db/COG/cdd2cog/cog-20.def.tsv)')

    return parser.parse_args()
#'/Users/daiki/Python/M2/rpsblast/data/cddid_COG.tbl',
#'/home/tmp/db/COG/cdd2cog/cddid_COG.tbl'
#'/Users/daiki/Python/M2/rpsblast/data/cog-20.def.tsv',
#'/home/tmp/db/COG/cdd2cog/cog-20.def.tsv'

default_colors = [
    # r, g, b, a
    [92, 192, 98, 0.5],
    [90, 155, 212, 0.5],
    [246, 236, 86, 0.6],
    [241, 90, 96, 0.4],
    [255, 117, 0, 0.3],
    [82, 82, 190, 0.2],
]
default_colors = [
    [i[0] / 255.0, i[1] / 255.0, i[2] / 255.0, i[3]]
    for i in default_colors
]

def draw_ellipse(fig, ax, x, y, w, h, a, fillcolor):
    e = patches.Ellipse(
        xy=(x, y),
        width=w,
        height=h,
        angle=a,
        color=fillcolor)
    ax.add_patch(e)

def draw_triangle(fig, ax, x1, y1, x2, y2, x3, y3, fillcolor):
    xy = [
        (x1, y1),
        (x2, y2),
        (x3, y3),
    ]
    polygon = patches.Polygon(
        xy=xy,
        closed=True,
        color=fillcolor)
    ax.add_patch(polygon)

def draw_text(fig, ax, x, y, text, color=[0, 0, 0, 1], fontsize=14, ha="center", va="center"):
    ax.text(
        x, y, text,
        horizontalalignment=ha,
        verticalalignment=va,
        fontsize=fontsize,
        color="black")

def draw_annotate(fig, ax, x, y, textx, texty, text, color=[0, 0, 0, 1], arrowcolor=[0, 0, 0, 0.3]):
    plt.annotate(
        text,
        xy=(x, y),
        xytext=(textx, texty),
        arrowprops=dict(color=arrowcolor, shrink=0, width=0.5, headwidth=8),
        fontsize=14,
        color=color,
        xycoords="data",
        textcoords="data",
        horizontalalignment='center',
        verticalalignment='center'
    )

def get_labels(data, fill=["number"]):
    """
    get a dict of labels for groups in data

    @type data: list[Iterable]
    @rtype: dict[str, str]

    input
      data: data to get label for
      fill: ["number"|"logic"|"percent"]

    return
      labels: a dict of labels for different sets

    example:
    In [12]: get_labels([range(10), range(5,15), range(3,8)], fill=["number"])
    Out[12]:
    {'001': '0',
     '010': '5',
     '011': '0',
     '100': '3',
     '101': '2',
     '110': '2',
     '111': '3'}
    """

    N = len(data)

    sets_data = [set(data[i]) for i in range(N)]  # sets for separate groups
    s_all = set(chain(*data))                     # union of all sets

    # bin(3) --> '0b11', so bin(3).split('0b')[-1] will remove "0b"
    set_collections = {}
    for n in range(1, 2**N):
        key = bin(n).split('0b')[-1].zfill(N)
        value = s_all
        sets_for_intersection = [sets_data[i] for i in range(N) if  key[i] == '1']
        sets_for_difference = [sets_data[i] for i in range(N) if  key[i] == '0']
        for s in sets_for_intersection:
            value = value & s
        for s in sets_for_difference:
            value = value - s
        set_collections[key] = value

    labels = {k: "" for k in set_collections}
    if "logic" in fill:
        for k in set_collections:
            labels[k] = k + ": "
    if "number" in fill:
        for k in set_collections:
            labels[k] += str(len(set_collections[k]))
    if "percent" in fill:
        data_size = len(s_all)
        for k in set_collections:
            labels[k] += "(%.1f%%)" % (100.0 * len(set_collections[k]) / data_size)

    return labels

def venn4(labels, ax, names=['A', 'B', 'C', 'D'], **options):
    """
    plots a 4-set Venn diagram

    @type labels: dict[str, str]
    @type names: list[str]
    @rtype: (Figure, AxesSubplot)

    input
      labels: a label dict where keys are identified via binary codes ('0001', '0010', '0100', ...),
              hence a valid set could look like: {'0001': 'text 1', '0010': 'text 2', '0100': 'text 3', ...}.
              unmentioned codes are considered as ''.
      names:  group names
      more:   colors, figsize, dpi, fontsize

    return
      pyplot Figure and AxesSubplot object
    """
    colors = options.get('colors', [default_colors[i] for i in range(4)])
    figsize = options.get('figsize', (12, 12))
    dpi = options.get('dpi', 96)
    fontsize = options.get('fontsize', 14)

    fig = plt.figure(0, figsize=figsize, dpi=dpi)
    #ax = fig.add_subplot(111, aspect='equal')
    ax.set_axis_off()
    ax.set_ylim(bottom=0.0, top=1.0)
    ax.set_xlim(left=0.0, right=1.0)

    # body
    draw_ellipse(fig, ax, 0.350, 0.400, 0.72, 0.45, 140.0, colors[0])
    draw_ellipse(fig, ax, 0.450, 0.500, 0.72, 0.45, 140.0, colors[1])
    draw_ellipse(fig, ax, 0.544, 0.500, 0.72, 0.45, 40.0, colors[2])
    draw_ellipse(fig, ax, 0.644, 0.400, 0.72, 0.45, 40.0, colors[3])
    draw_text(fig, ax, 0.85, 0.42, labels.get('0001', ''), fontsize=fontsize)
    draw_text(fig, ax, 0.68, 0.72, labels.get('0010', ''), fontsize=fontsize)
    draw_text(fig, ax, 0.77, 0.59, labels.get('0011', ''), fontsize=fontsize)
    draw_text(fig, ax, 0.32, 0.72, labels.get('0100', ''), fontsize=fontsize)
    draw_text(fig, ax, 0.71, 0.30, labels.get('0101', ''), fontsize=fontsize)
    draw_text(fig, ax, 0.50, 0.66, labels.get('0110', ''), fontsize=fontsize)
    draw_text(fig, ax, 0.65, 0.50, labels.get('0111', ''), fontsize=fontsize)
    draw_text(fig, ax, 0.14, 0.42, labels.get('1000', ''), fontsize=fontsize)
    draw_text(fig, ax, 0.50, 0.17, labels.get('1001', ''), fontsize=fontsize)
    draw_text(fig, ax, 0.29, 0.30, labels.get('1010', ''), fontsize=fontsize)
    draw_text(fig, ax, 0.39, 0.24, labels.get('1011', ''), fontsize=fontsize)
    draw_text(fig, ax, 0.23, 0.59, labels.get('1100', ''), fontsize=fontsize)
    draw_text(fig, ax, 0.61, 0.24, labels.get('1101', ''), fontsize=fontsize)
    draw_text(fig, ax, 0.35, 0.50, labels.get('1110', ''), fontsize=fontsize)
    draw_text(fig, ax, 0.50, 0.38, labels.get('1111', ''), fontsize=fontsize)

    # legend
    draw_text(fig, ax, 0.13, 0.18, names[0], colors[0], fontsize=fontsize, ha="right")
    draw_text(fig, ax, 0.18, 0.83, names[1], colors[1], fontsize=fontsize, ha="right", va="bottom")
    draw_text(fig, ax, 0.82, 0.83, names[2], colors[2], fontsize=fontsize, ha="left", va="bottom")
    draw_text(fig, ax, 0.87, 0.18, names[3], colors[3], fontsize=fontsize, ha="left", va="top")
    #leg = ax.legend(names, loc='center left', bbox_to_anchor=(1.0, 0.5), fancybox=True)
    #leg.get_frame().set_alpha(0.5)

    return fig#, ax

def venn5(labels, ax, names=['A', 'B', 'C', 'D', 'E'], **options):
    """
    plots a 5-set Venn diagram

    @type labels: dict[str, str]
    @type names: list[str]
    @rtype: (Figure, AxesSubplot)

    input
      labels: a label dict where keys are identified via binary codes ('00001', '00010', '00100', ...),
              hence a valid set could look like: {'00001': 'text 1', '00010': 'text 2', '00100': 'text 3', ...}.
              unmentioned codes are considered as ''.
      names:  group names
      more:   colors, figsize, dpi, fontsize

    return
      pyplot Figure and AxesSubplot object
    """
    colors = options.get('colors', [default_colors[i] for i in range(5)])
    figsize = options.get('figsize', (13, 13))
    dpi = options.get('dpi', 96)
    fontsize = options.get('fontsize', 14)

    fig = plt.figure(0, figsize=figsize, dpi=dpi)
    #ax = fig.add_subplot(111, aspect='equal')
    ax.set_axis_off()
    ax.set_ylim(bottom=0.0, top=1.0)
    ax.set_xlim(left=0.0, right=1.0)

    # body
    draw_ellipse(fig, ax, 0.428, 0.449, 0.87, 0.50, 155.0, colors[0])
    draw_ellipse(fig, ax, 0.469, 0.543, 0.87, 0.50, 82.0, colors[1])
    draw_ellipse(fig, ax, 0.558, 0.523, 0.87, 0.50, 10.0, colors[2])
    draw_ellipse(fig, ax, 0.578, 0.432, 0.87, 0.50, 118.0, colors[3])
    draw_ellipse(fig, ax, 0.489, 0.383, 0.87, 0.50, 46.0, colors[4])
    draw_text(fig, ax, 0.27, 0.11, labels.get('00001', ''), fontsize=fontsize)
    draw_text(fig, ax, 0.72, 0.11, labels.get('00010', ''), fontsize=fontsize)
    draw_text(fig, ax, 0.55, 0.13, labels.get('00011', ''), fontsize=fontsize)
    draw_text(fig, ax, 0.91, 0.58, labels.get('00100', ''), fontsize=fontsize)
    draw_text(fig, ax, 0.78, 0.64, labels.get('00101', ''), fontsize=fontsize)
    draw_text(fig, ax, 0.84, 0.41, labels.get('00110', ''), fontsize=fontsize)
    draw_text(fig, ax, 0.76, 0.55, labels.get('00111', ''), fontsize=fontsize)
    draw_text(fig, ax, 0.51, 0.90, labels.get('01000', ''), fontsize=fontsize)
    draw_text(fig, ax, 0.39, 0.15, labels.get('01001', ''), fontsize=fontsize)
    draw_text(fig, ax, 0.42, 0.78, labels.get('01010', ''), fontsize=fontsize)
    draw_text(fig, ax, 0.50, 0.15, labels.get('01011', ''), fontsize=fontsize)
    draw_text(fig, ax, 0.67, 0.76, labels.get('01100', ''), fontsize=fontsize)
    draw_text(fig, ax, 0.70, 0.71, labels.get('01101', ''), fontsize=fontsize)
    draw_text(fig, ax, 0.51, 0.74, labels.get('01110', ''), fontsize=fontsize)
    draw_text(fig, ax, 0.64, 0.67, labels.get('01111', ''), fontsize=fontsize)
    draw_text(fig, ax, 0.10, 0.61, labels.get('10000', ''), fontsize=fontsize)
    draw_text(fig, ax, 0.20, 0.31, labels.get('10001', ''), fontsize=fontsize)
    draw_text(fig, ax, 0.76, 0.25, labels.get('10010', ''), fontsize=fontsize)
    draw_text(fig, ax, 0.65, 0.23, labels.get('10011', ''), fontsize=fontsize)
    draw_text(fig, ax, 0.18, 0.50, labels.get('10100', ''), fontsize=fontsize)
    draw_text(fig, ax, 0.21, 0.37, labels.get('10101', ''), fontsize=fontsize)
    draw_text(fig, ax, 0.81, 0.37, labels.get('10110', ''), fontsize=fontsize)
    draw_text(fig, ax, 0.74, 0.40, labels.get('10111', ''), fontsize=fontsize)
    draw_text(fig, ax, 0.27, 0.70, labels.get('11000', ''), fontsize=fontsize)
    draw_text(fig, ax, 0.34, 0.25, labels.get('11001', ''), fontsize=fontsize)
    draw_text(fig, ax, 0.33, 0.72, labels.get('11010', ''), fontsize=fontsize)
    draw_text(fig, ax, 0.51, 0.22, labels.get('11011', ''), fontsize=fontsize)
    draw_text(fig, ax, 0.25, 0.58, labels.get('11100', ''), fontsize=fontsize)
    draw_text(fig, ax, 0.28, 0.39, labels.get('11101', ''), fontsize=fontsize)
    draw_text(fig, ax, 0.36, 0.66, labels.get('11110', ''), fontsize=fontsize)
    draw_text(fig, ax, 0.51, 0.47, labels.get('11111', ''), fontsize=fontsize)

    # legend
    draw_text(fig, ax, 0.02, 0.72, names[0], colors[0], fontsize=fontsize, ha="right")
    draw_text(fig, ax, 0.72, 0.94, names[1], colors[1], fontsize=fontsize, va="bottom")
    draw_text(fig, ax, 0.97, 0.74, names[2], colors[2], fontsize=fontsize, ha="left")
    draw_text(fig, ax, 0.88, 0.05, names[3], colors[3], fontsize=fontsize, ha="left")
    draw_text(fig, ax, 0.12, 0.05, names[4], colors[4], fontsize=fontsize, ha="right")
    #leg = ax.legend(names, loc='center left', bbox_to_anchor=(1.0, 0.5), fancybox=True)
    #leg.get_frame().set_alpha(0.5)

    return fig#, ax

def venn6(labels, ax, names=['A', 'B', 'C', 'D', 'E'], **options):
    """
    plots a 6-set Venn diagram

    @type labels: dict[str, str]
    @type names: list[str]
    @rtype: (Figure, AxesSubplot)

    input
      labels: a label dict where keys are identified via binary codes ('000001', '000010', '000100', ...),
              hence a valid set could look like: {'000001': 'text 1', '000010': 'text 2', '000100': 'text 3', ...}.
              unmentioned codes are considered as ''.
      names:  group names
      more:   colors, figsize, dpi, fontsize

    return
      pyplot Figure and AxesSubplot object
    """
    colors = options.get('colors', [default_colors[i] for i in range(6)])
    figsize = options.get('figsize', (20, 20))
    dpi = options.get('dpi', 96)
    fontsize = options.get('fontsize', 14)

    fig = plt.figure(0, figsize=figsize, dpi=dpi)
    #ax = fig.add_subplot(111, aspect='equal')
    ax.set_axis_off()
    ax.set_ylim(bottom=0.230, top=0.845)
    ax.set_xlim(left=0.173, right=0.788)

    # body
    # See https://web.archive.org/web/20040819232503/http://www.hpl.hp.com/techreports/2000/HPL-2000-73.pdf
    draw_triangle(fig, ax, 0.637, 0.921, 0.649, 0.274, 0.188, 0.667, colors[0])
    draw_triangle(fig, ax, 0.981, 0.769, 0.335, 0.191, 0.393, 0.671, colors[1])
    draw_triangle(fig, ax, 0.941, 0.397, 0.292, 0.475, 0.456, 0.747, colors[2])
    draw_triangle(fig, ax, 0.662, 0.119, 0.316, 0.548, 0.662, 0.700, colors[3])
    draw_triangle(fig, ax, 0.309, 0.081, 0.374, 0.718, 0.681, 0.488, colors[4])
    draw_triangle(fig, ax, 0.016, 0.626, 0.726, 0.687, 0.522, 0.327, colors[5])
    draw_text(fig, ax, 0.212, 0.562, labels.get('000001', ''), fontsize=fontsize)
    draw_text(fig, ax, 0.430, 0.249, labels.get('000010', ''), fontsize=fontsize)
    draw_text(fig, ax, 0.356, 0.444, labels.get('000011', ''), fontsize=fontsize)
    draw_text(fig, ax, 0.609, 0.255, labels.get('000100', ''), fontsize=fontsize)
    draw_text(fig, ax, 0.323, 0.546, labels.get('000101', ''), fontsize=fontsize)
    draw_text(fig, ax, 0.513, 0.316, labels.get('000110', ''), fontsize=fontsize)
    draw_text(fig, ax, 0.523, 0.348, labels.get('000111', ''), fontsize=fontsize)
    draw_text(fig, ax, 0.747, 0.458, labels.get('001000', ''), fontsize=fontsize)
    draw_text(fig, ax, 0.325, 0.492, labels.get('001001', ''), fontsize=fontsize)
    draw_text(fig, ax, 0.670, 0.481, labels.get('001010', ''), fontsize=fontsize)
    draw_text(fig, ax, 0.359, 0.478, labels.get('001011', ''), fontsize=fontsize)
    draw_text(fig, ax, 0.653, 0.444, labels.get('001100', ''), fontsize=fontsize)
    draw_text(fig, ax, 0.344, 0.526, labels.get('001101', ''), fontsize=fontsize)
    draw_text(fig, ax, 0.653, 0.466, labels.get('001110', ''), fontsize=fontsize)
    draw_text(fig, ax, 0.363, 0.503, labels.get('001111', ''), fontsize=fontsize)
    draw_text(fig, ax, 0.750, 0.616, labels.get('010000', ''), fontsize=fontsize)
    draw_text(fig, ax, 0.682, 0.654, labels.get('010001', ''), fontsize=fontsize)
    draw_text(fig, ax, 0.402, 0.310, labels.get('010010', ''), fontsize=fontsize)
    draw_text(fig, ax, 0.392, 0.421, labels.get('010011', ''), fontsize=fontsize)
    draw_text(fig, ax, 0.653, 0.691, labels.get('010100', ''), fontsize=fontsize)
    draw_text(fig, ax, 0.651, 0.644, labels.get('010101', ''), fontsize=fontsize)
    draw_text(fig, ax, 0.490, 0.340, labels.get('010110', ''), fontsize=fontsize)
    draw_text(fig, ax, 0.468, 0.399, labels.get('010111', ''), fontsize=fontsize)
    draw_text(fig, ax, 0.692, 0.545, labels.get('011000', ''), fontsize=fontsize)
    draw_text(fig, ax, 0.666, 0.592, labels.get('011001', ''), fontsize=fontsize)
    draw_text(fig, ax, 0.665, 0.496, labels.get('011010', ''), fontsize=fontsize)
    draw_text(fig, ax, 0.374, 0.470, labels.get('011011', ''), fontsize=fontsize)
    draw_text(fig, ax, 0.653, 0.537, labels.get('011100', ''), fontsize=fontsize)
    draw_text(fig, ax, 0.652, 0.579, labels.get('011101', ''), fontsize=fontsize)
    draw_text(fig, ax, 0.653, 0.488, labels.get('011110', ''), fontsize=fontsize)
    draw_text(fig, ax, 0.389, 0.486, labels.get('011111', ''), fontsize=fontsize)
    draw_text(fig, ax, 0.553, 0.806, labels.get('100000', ''), fontsize=fontsize)
    draw_text(fig, ax, 0.313, 0.604, labels.get('100001', ''), fontsize=fontsize)
    draw_text(fig, ax, 0.388, 0.694, labels.get('100010', ''), fontsize=fontsize)
    draw_text(fig, ax, 0.375, 0.633, labels.get('100011', ''), fontsize=fontsize)
    draw_text(fig, ax, 0.605, 0.359, labels.get('100100', ''), fontsize=fontsize)
    draw_text(fig, ax, 0.334, 0.555, labels.get('100101', ''), fontsize=fontsize)
    draw_text(fig, ax, 0.582, 0.397, labels.get('100110', ''), fontsize=fontsize)
    draw_text(fig, ax, 0.542, 0.372, labels.get('100111', ''), fontsize=fontsize)
    draw_text(fig, ax, 0.468, 0.708, labels.get('101000', ''), fontsize=fontsize)
    draw_text(fig, ax, 0.355, 0.572, labels.get('101001', ''), fontsize=fontsize)
    draw_text(fig, ax, 0.420, 0.679, labels.get('101010', ''), fontsize=fontsize)
    draw_text(fig, ax, 0.375, 0.597, labels.get('101011', ''), fontsize=fontsize)
    draw_text(fig, ax, 0.641, 0.436, labels.get('101100', ''), fontsize=fontsize)
    draw_text(fig, ax, 0.348, 0.538, labels.get('101101', ''), fontsize=fontsize)
    draw_text(fig, ax, 0.635, 0.453, labels.get('101110', ''), fontsize=fontsize)
    draw_text(fig, ax, 0.370, 0.548, labels.get('101111', ''), fontsize=fontsize)
    draw_text(fig, ax, 0.594, 0.689, labels.get('110000', ''), fontsize=fontsize)
    draw_text(fig, ax, 0.579, 0.670, labels.get('110001', ''), fontsize=fontsize)
    draw_text(fig, ax, 0.398, 0.670, labels.get('110010', ''), fontsize=fontsize)
    draw_text(fig, ax, 0.395, 0.653, labels.get('110011', ''), fontsize=fontsize)
    draw_text(fig, ax, 0.633, 0.682, labels.get('110100', ''), fontsize=fontsize)
    draw_text(fig, ax, 0.616, 0.656, labels.get('110101', ''), fontsize=fontsize)
    draw_text(fig, ax, 0.587, 0.427, labels.get('110110', ''), fontsize=fontsize)
    draw_text(fig, ax, 0.526, 0.415, labels.get('110111', ''), fontsize=fontsize)
    draw_text(fig, ax, 0.495, 0.677, labels.get('111000', ''), fontsize=fontsize)
    draw_text(fig, ax, 0.505, 0.648, labels.get('111001', ''), fontsize=fontsize)
    draw_text(fig, ax, 0.428, 0.663, labels.get('111010', ''), fontsize=fontsize)
    draw_text(fig, ax, 0.430, 0.631, labels.get('111011', ''), fontsize=fontsize)
    draw_text(fig, ax, 0.639, 0.524, labels.get('111100', ''), fontsize=fontsize)
    draw_text(fig, ax, 0.591, 0.604, labels.get('111101', ''), fontsize=fontsize)
    draw_text(fig, ax, 0.622, 0.477, labels.get('111110', ''), fontsize=fontsize)
    draw_text(fig, ax, 0.501, 0.523, labels.get('111111', ''), fontsize=fontsize)

    # legend
    draw_text(fig, ax, 0.674, 0.824, names[0], colors[0], fontsize=fontsize)
    draw_text(fig, ax, 0.747, 0.751, names[1], colors[1], fontsize=fontsize)
    draw_text(fig, ax, 0.739, 0.396, names[2], colors[2], fontsize=fontsize)
    draw_text(fig, ax, 0.700, 0.247, names[3], colors[3], fontsize=fontsize)
    draw_text(fig, ax, 0.291, 0.255, names[4], colors[4], fontsize=fontsize)
    draw_text(fig, ax, 0.203, 0.484, names[5], colors[5], fontsize=fontsize)
    #leg = ax.legend(names, loc='center left', bbox_to_anchor=(1.0, 0.5), fancybox=True)
    #leg.get_frame().set_alpha(0.5)

    return fig#, ax




def run_rpsblast(paths_to_proteins = None, 
                 path_to_cogdb = None, 
                 evalue = None):#, num_threads = None
    from subprocess import Popen
    error1 = "specify the path to your Cog database with cogdb option. (default:/home/tmp/db/COG/Cog)"
    #assert os.path.exists('/home/tmp/db/COG/Cog/'), error1
    
    if f'rps_{evalue}' not in os.listdir(path='./'):
        os.system(f'mkdir rps_{evalue}')
    
    path_to_rpsRes = []
    procs = []

    for path in paths_to_proteins:
            name = os.path.splitext(os.path.basename(path))[0]
            procs += [Popen(f"rpsblast -query {path} -db {path_to_cogdb} -out ./rps_{evalue}/{name}.txt -evalue {evalue} -outfmt 6"
                       , shell=True)]

            path_to_rpsRes.append(f"./rps_{evalue}/{name}.txt")
    
    [p.wait() for p in procs]
    return path_to_rpsRes

def preprocess(rps = None,
               cddid = None,
               cog = None):

    cddid["CDD"] = "CDD:" + cddid["CDD"].astype(str)
    _ = pd.merge(rps, cddid, on = ["CDD"]).iloc[:, [0, 1, 12]]
    _df = pd.merge(_, cog, on = ["COG"]).iloc[:, [0, 1, 2, 3, 4, 5]]
    return _df, dict(Counter("".join(_df['Group'])))

def sorter(df_i = None,
           A2Z = None):
    tmp = []
    for i in range(len(A2Z)):
        if A2Z[i] in list(df_i.keys()):
            tmp += [df_i[A2Z[i]]]
        else:
            tmp += [0]
    return tmp

def get_main_dataset(path_to_rpsRes = None,
                     path_to_cddid = None,
                     path_to_cog = None, evalue = None):
    if f'out_{evalue}' not in os.listdir(path='./'):
        os.system(f'mkdir ./out_{evalue}/')
    if 'COGdata' not in os.listdir(path=f"./out_{evalue}/"):
        os.system(f'mkdir ./out_{evalue}/COGdata/') 
    
    df_i = {}
    out_COG_i = {}
    out = pd.DataFrame()
    out_COG = pd.DataFrame()
    A2Z = [chr(i) for i in range(65, 65+26)]
    for i, path in enumerate(path_to_rpsRes):
        #load data
        rps_i = pd.read_table(path_to_rpsRes[i], names=["cdd_id", "CDD", "a", "b","c","d",
                                        "e","f","g","h","eValue","j"]).drop_duplicates(['cdd_id'])

        cddid = pd.read_table(path_to_cddid, names=["CDD", "COG", "a", "b", "c"])

        cog = pd.read_table(path_to_cog, names=["COG", "Group", "gene_name",
                                         "gene", "E3", "F3", "G3"], encoding='cp1252')
        
        col_name = os.path.splitext(os.path.basename(path))[0]#[:-len('_rpsblastout')]
        #processing
        out_COG_i[col_name], df_i[col_name] = preprocess(rps = rps_i, 
                                  cddid = cddid,
                                  cog = cog)
        
        
        COG_i = out_COG_i[col_name].rename(columns={'cdd_id': f"{col_name}"}).iloc[:, [0,1,2,3,4,5]]
        out_i = pd.DataFrame(sorter(df_i = df_i[col_name], A2Z = A2Z), columns=[f"{col_name}"])

        out = pd.concat([out, out_i], axis = 1)
        out_COG = pd.concat([out_COG, COG_i], axis = 1)
    count_data = pd.concat([pd.DataFrame(A2Z, columns=['COG']), out], axis = 1)
    

    _count_data = count_data.iloc[:,1:len(count_data.columns)]
    _ = _count_data/_count_data.sum()
 
    ratio_data = pd.concat([pd.DataFrame(A2Z, columns=['COG']), _], axis = 1)
    count_data.to_csv(f"./out_{evalue}/COGdata/COG_count.csv") ;ratio_data.to_csv(f"./out_{evalue}/COGdata/COG_ratio.csv")
    out_COG.to_csv(f"./out_{evalue}/COGdata/COG_annotation.csv")
    
    return count_data, ratio_data, out_COG_i

def plot_bar(df = None, name = None, n_black = None, size = None, evalue = None):
    if 'bar' not in os.listdir(path=f"./out_{evalue}/"):
        os.system(f'mkdir ./out_{evalue}/bar/') 
    # 棒の配置位置、ラベルを用意
    labels = list(df['COG'])
    x = np.array(range(len(labels)))

    # 各系列のデータを用意
    data, legend = [], []
    for col_name in df.columns[1:len(df.columns)]:
        data.append(df[col_name])
        legend.append(col_name) 

    # マージンを設定
    margin = 0.2  #0 <margin< 1
    totoal_width = 1 - margin
    fig = plt.figure(figsize=(3*size,2*size))
    # 棒グラフをプロット
    c1 = ['royalblue','sandybrown','yellowgreen','hotpink','0.7']*500
    for i, h in enumerate(data):
        pos = x - totoal_width *( 1- (2*i+1)/len(data) )/2
        plt.bar(pos, h, width = totoal_width/len(data), color =c1[i])
 
    plt.xticks(x, labels)
    fig.savefig(f"./out_{evalue}/bar/COG_{name}_NoLegend.pdf")
    plt.legend(legend)
    fig.savefig(f"./out_{evalue}/bar/COG_{name}.pdf")
    
    #コードが冗長
    # 棒の配置位置、ラベルを用意
    labels = list(df['COG'])
    x = np.array(range(len(labels)))

    # 各系列のデータを用意
    data, legend = [], []
    for col_name in df.columns[1:len(df.columns)]:
        data.append(df[col_name])
        legend.append(col_name) 

    # マージンを設定
    margin = 0.2  #0 <margin< 1
    totoal_width = 1 - margin
    fig = plt.figure(figsize=(3*size,2*size))
    c2 = ["0"]*n_black+["0.7"]*2500
    for i, h in enumerate(data):
        pos = x - totoal_width *( 1- (2*i+1)/len(data) )/2
        plt.bar(pos, h, width = totoal_width/len(data), color =c2[i])
    
    plt.xticks(x, labels)
    fig.savefig(f"./out_{evalue}/bar/COG_{name}_NoColor__NoLegend.pdf")
    plt.legend(legend)
    fig.savefig(f"./out_{evalue}/bar/COG_{name}_NoColor.pdf")
    
    
def CLR_PCA(df = None, size = None, delta = None, tag = None, n_green = None, CLR = None, evalue = None):#各行にCOG。
    if f'PCA_{tag}' not in os.listdir(path=f"./out_{evalue}/"):
        os.system(f'mkdir ./out_{evalue}/PCA_{tag}/') 
    if f'PCA_{tag}_{delta}' not in os.listdir(path=f"./out_{evalue}/PCA_{tag}/"):
        os.system(f'mkdir ./out_{evalue}/PCA_{tag}/PCA_{tag}_{delta}') 
        
    def Myclr(df):
        def geo_mean(iterable):
            a = np.array(iterable).astype(float)
            return a.prod()**(1.0/len(a))
        df_clr = pd.DataFrame()
        for i in range(len(df.index)):
            tmp = df.iloc[i,:]/geo_mean(df.iloc[i,:])
            df_clr = pd.concat([df_clr, tmp.map(math.log)], axis=1)
        return df_clr.T
    #ゼロ値の補完
    clr_in = df.iloc[:, 1:] + delta
    
    #CLR
    if CLR:
        df_clr = Myclr(clr_in.T).T#clr関数は行方向に和が１ものしか受け付けない
    else:
        df_clr = clr_in
    
    #PCA
    pca = PCA(n_components=2)
    _ = pca.fit_transform(df_clr.T)
    df_pca = pd.DataFrame(_, columns = ["PCA1", "PCA2"])
    
    def plot_PCA(df_pca, pca, df, evalue):
        fig = plt.figure(figsize=(size *2, size * 2))
        ax1 = fig.subplots()
        ax1.scatter(df_pca.PCA1[:n_green], df_pca.PCA2[:n_green], alpha=0.8, c='g')
        ax1.scatter(df_pca.PCA1[n_green:], df_pca.PCA2[n_green:], alpha=0.8)
        for x, y, name in zip(df_pca.PCA1, df_pca.PCA2, df.columns[1:]):
            ax1.text(x, y, name)
        
        ax1.grid()
        ax1.set_xlabel(f"PC1({(pca.explained_variance_ratio_[0]*100).round(2)}%)")
        ax1.set_ylabel(f"PC2({(pca.explained_variance_ratio_[1]*100).round(2)}%)")
        fig.savefig(f"./out_{evalue}/PCA_{tag}/PCA_{tag}_{delta}/PCA_COG_{tag}_{delta}.pdf")

        ax2 = ax1.twiny().twinx()
        for x, y, name in zip(pca.components_[0], pca.components_[1], df.COG):
            ax2.text(x, y, name)
            ax2.arrow(x=0,y=0, dx=x, dy=y,
                     width=.0001, length_includes_head=True,color='m')
        ax2.scatter(pca.components_[0],  pca.components_[1], alpha=0, color='m')
        ax1.set_title(f"PC1 Loading", fontsize=20/size*2)
        ax2.set_ylabel(f"PC2 Loading", fontsize=20/size*2)
        fig.savefig(f"./out_{evalue}/PCA_{tag}/PCA_{tag}_{delta}/PCA_COG_{tag}_{delta}_withLoadingFactor.pdf")

    plot_PCA(df_pca, pca, df, evalue)
    
    #コードが冗長
    def plot_PCA_NoName(df_pca, pca, df, evalue):
        fig = plt.figure(figsize=(size *2, size * 2))
        ax1 = fig.subplots()
        ax1.scatter(df_pca.PCA1[:n_green], df_pca.PCA2[:n_green], alpha=0.8, c='g')
        ax1.scatter(df_pca.PCA1[n_green:], df_pca.PCA2[n_green:], alpha=0.8)
        ax1.grid()
        ax1.set_xlabel(f"PC1({(pca.explained_variance_ratio_[0]*100).round(2)}%)")
        ax1.set_ylabel(f"PC2({(pca.explained_variance_ratio_[1]*100).round(2)}%)")

        ax2 = ax1.twiny().twinx()
        for x, y, name in zip(pca.components_[0], pca.components_[1], df.COG):
            ax2.text(x, y, name)
            ax2.arrow(x=0,y=0, dx=x, dy=y,
                     width=.0001, length_includes_head=True,color='m')
        ax2.scatter(pca.components_[0],  pca.components_[1], alpha=0, color='m')
        ax1.set_title(f"PC1 Loading", fontsize=20/size*2)
        ax2.set_ylabel(f"PC2 Loading", fontsize=20/size*2)
        fig.savefig(f"./out_{evalue}/PCA_{tag}/PCA_{tag}_{delta}/PCA_COG_{tag}_{delta}_NoName.pdf")

    plot_PCA_NoName(df_pca, pca, df, evalue)
    
def plot_or_not(unique_COGs):
    return sum([unique_COG==set() for unique_COG in unique_COGs]) !=len(unique_COGs)

def venn_func(dataset, unique_COG, labels, ax):
    subsets = get_labels(unique_COG, fill=['number', 'logic'])
    if len(list(dataset.keys()))==2:
        return matplotlib_venn.venn2(subsets=unique_COG, set_labels = labels)
    elif len(list(dataset.keys()))==3:
        return matplotlib_venn.venn3(subsets=unique_COG, set_labels = labels)
    elif len(list(dataset.keys()))==4:
        return venn4(subsets, ax, names = labels)
    elif len(list(dataset.keys()))==5:
        return venn5(subsets, ax, names = labels)
    elif len(list(dataset.keys()))==6:
        return venn6(subsets, ax, names = labels)
    else:
        pass
            
def plot_venn(dataset = None, size = None, evalue = None):
    if 2 <=len(dataset.keys()) <=6:
        A2Z = [chr(i) for i in range(65, 65+26)]

        #1
        fig = plt.figure(figsize=(10,10))
        ax = fig.add_subplot(111, aspect='equal')

        unique_COG = []
        for j in range(len(dataset.keys())):
            x = dataset[list(dataset.keys())[j]]
            unique_COG.append(set(x['COG'].unique()))

        venn_func(dataset, unique_COG, list(dataset.keys()), ax)
        ax.set_title('All genes')  
        fig.savefig(f"./out_{evalue}/venn{len(list(dataset.keys()))}Diagrams.pdf")

        #2
        fig = plt.figure(figsize=(size * 3, size * 4))
        for i, alphabet in enumerate(A2Z):
            unique_COG = []
            for j in range(len(dataset.keys())):
                x = dataset[list(dataset.keys())[j]]
                unique_COG.append(set(x[x['Group']==f'{alphabet}']['COG'].unique()))

            if plot_or_not(unique_COG):
                ax = fig.add_subplot(6, 5, i+1)
                venn_func(dataset, unique_COG, list(dataset.keys()), ax)
                ax.set_title(f'{alphabet}')
                plt.tight_layout()
                fig.savefig(f"./out_{evalue}/COGvenn{len(list(dataset.keys()))}Diagrams.pdf")

    #コードが冗長
    unique_COG = []
    for j in range(len(dataset.keys())):
        x = dataset[list(dataset.keys())[j]]
        unique_COG.append(set(x['COG'].unique()))
    eigengene = {}
    _ = unique_COG
    tmp = list(dataset.keys())
    eigengene[f"{tmp[0]}_eigengene"] = list(_[0])
    for i in range(1, len(tmp)):
        eigengene[f"{tmp[0]}_eigengene"] = list( set(eigengene[f"{tmp[0]}_eigengene"]) - set(list(_[i])))

    Group = []
    name, gene, gene_name =[], [], [] 
    for i in range(len(eigengene[f"{tmp[0]}_eigengene"])):
        x = dataset[f"{tmp[0]}"].COG==eigengene[f"{tmp[0]}_eigengene"][i]
        Group+=set(list(dataset[f"{tmp[0]}"]['Group'][x]))
        #assert len(set(list(dataset[f"{tmp[0]}"]['cdd_id'][x])))==1, i
        name+=set([list(dataset[f"{tmp[0]}"]['cdd_id'][x])[0]])
        gene+=set(list(dataset[f"{tmp[0]}"]['gene'][x]))
        gene_name+=set(list(dataset[f"{tmp[0]}"]['gene_name'][x]))
    pd.DataFrame([eigengene[f"{tmp[0]}_eigengene"], gene, gene_name, Group, name],
             index=[f"{tmp[0]}_eigengene", "gene", "gene name", 'Group', 'one of the names']).T.to_csv(f"./out_{evalue}/COGdata/{tmp[0]}_uniquegene.csv")
def main():
    print(f'Output directory = ', end='')
    [print(f'out_{_} ', end='') for _ in get_args().evalue]
    for e in get_args().evalue:
        if get_args().AA is not None:
            print(f'\n- rpsblast now (e-value = {e})..')
            num_files = len(get_args().AA)
            path_to_rpsRes = run_rpsblast(paths_to_proteins = get_args().AA, 
                                          path_to_cogdb = get_args().cogdb, 
                                          evalue = e)#, num_threads = get_args().num_threads
    
            count_data, ratio_data, dataset = get_main_dataset(path_to_rpsRes = path_to_rpsRes,
                                                              path_to_cddid = get_args().cddid,
                                                              path_to_cog = get_args().cog, evalue = e)
    
        elif get_args().rps is not None:
            print(f'\n- loading data..')
            num_files = len(get_args().rps)
            count_data, ratio_data, dataset = get_main_dataset(path_to_rpsRes = get_args().rps,
                                                              path_to_cddid = get_args().cddid,
                                                              path_to_cog = get_args().cog, evalue = e)
    
        if 2 <=num_files :
            print('- barplot..', end ='')
            plot_bar(df = count_data, name ='count', n_black = get_args().n_black, size = get_args().bar_size, evalue = e)
            plot_bar(df = ratio_data, name ='ratio', n_black = get_args().n_black, size = get_args().bar_size, evalue = e)
            print(f'done!')
        
        if 2 <= num_files:
            print('- PCA..', end ='')
            for i in [0]:
                CLR_PCA(df = count_data, size = get_args().PCA_size,
                        delta =i, tag = "count", n_green = get_args().n_green, CLR = False, evalue = e)
            for i in [1]:
                CLR_PCA(df = ratio_data, size = get_args().PCA_size,
                        delta =i, tag = "ratio", n_green = get_args().n_green, CLR = True, evalue = e)
            print(f'done!')
            
        if 2 <= num_files:
            if num_files <=6:
                print('- venn diagram..')
            if get_args().AA is not None:
                print(f'- finding unique genes of {os.path.splitext(os.path.basename(get_args().AA[0]))[0]}..', end ='')
            elif get_args().rps is not None:
                print(f'- finding unique genes of {os.path.splitext(os.path.basename(get_args().rps[0]))[0]}..', end ='')
            
            plot_venn(dataset = dataset, size = get_args().venn_size, evalue = e)
            print(f'done!')

if __name__ == "__main__":
    main()








