#!/usr/bin/env python3
import sys
import argparse

import numpy as np
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression
from scipy.stats import gaussian_kde

from datahandler import seqDataHandler, reader


def arg_parse(args):
    parser = argparse.ArgumentParser(description='Find linear correaltion between two BigWig data sets.')
    parser.add_argument('input_data', nargs='+', metavar='data', type=str,
                        help='List with input sequencing signals as a relative path to your current directory.'
                             'It is expected that the number of input_data values is either 2 or 4. If it is 2, the '
                             'first one is compared with the second one (e.g. wt and mutant). If it is 4, the '
                             'first two are added together and the second two are added together (e.g.'
                             'wt +strand + wt -strand and mutant +strand + mutant -strand).')
    parser.add_argument('--bed', type=str, required=True,
                        help='Annotation bed file.')
    parser.add_argument('--name', '-n', action='append', type=str, help='Names of the data sets.')
    parser.add_argument('--title', '-t', type=str, help='Plot title.')
    parsed_args = parser.parse_args(args)
    return parsed_args


def main(args):
    input_data = args.input_data
    if len(input_data) not in [2, 4]:
        raise ValueError('Number of passed input paths must be 2 or 4.')
    bed_name = args.bed
    names = args.name if args.name is not None else ['WT', 'Mutant']
    title = args.title if args.title is not None else 'Median Correlation Between WT and Mutant'
    title = title.replace('\\n', '\n')

    bed = reader.load_bam_bed_file(bed_name, rel_path='', is_abs_path=True)
    bw_files = []
    for file_path in input_data:
        bw_files.append(
            reader.load_big_file(
                name=file_path,
                rel_path='',
                is_abs_path=True
            )
        )

    all_values, chrom_start = seqDataHandler.get_values(bw_list=bw_files)
    annot, _ = seqDataHandler.annotate_all(all_values, bed, chrom_start)
    if len(input_data) == 2:
        annot_mean = [(np.median(wt), np.median(mut)) for wt, mut in zip(annot[0], annot[1])]
    else:
        annot_mean = [(np.median(wt1 + wt2), np.median(mut1 + mut2))
                      for wt1, wt2, mut1, mut2 in zip(annot[0], annot[1], annot[2], annot[3])]
    annot_wt, annot_mut = zip(*annot_mean)
    reg = LinearRegression(
        fit_intercept=True
    ).fit(np.asarray(list(annot_wt)).reshape(-1, 1), np.asarray(list(annot_mut)).reshape(-1, 1))
    r2 = reg.score(np.asarray(list(annot_wt)).reshape(-1, 1), np.asarray(list(annot_mut)).reshape(-1, 1))

    wt_mut = np.vstack([annot_wt, annot_mut])
    c_vec = gaussian_kde(wt_mut)(wt_mut)

    sc = plt.scatter(annot_wt, annot_mut, c=c_vec)
    axes = plt.gca()
    _, x_max = axes.get_xlim()
    plt.plot(np.arange(x_max), reg.coef_[0][0] * np.arange(x_max) + reg.intercept_)
    plt.plot(np.arange(x_max), np.arange(x_max), '--', color='grey')
    plt.text(
        0.1,
        0.9,
        'y = %.3fx + %.3f\nR2=%.3f' % (reg.coef_[0][0], reg.intercept_[0], r2),
        ha='left',
        va='center',
        bbox=dict(facecolor='none', edgecolor='black', boxstyle='round,pad=1'),
        transform=axes.transAxes
    )
    cbar = plt.colorbar(sc)
    cbar.ax.get_yaxis().labelpad = 15
    cbar.ax.set_ylabel('Gaussian Kernel Density', rotation=270)
    plt.xlabel(names[0])
    plt.ylabel(names[1])

    plt.title(title)
    plt.show()


if __name__ == '__main__':
    main(arg_parse(sys.argv[1:]))

