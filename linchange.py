#!/usr/bin/env python3
import sys
import argparse
import warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression
from scipy.stats import gaussian_kde

import pyBigWig
import gtfparse


def annotation_loading(path: str) -> pd.DataFrame:
    if path.split('.')[-1] == 'bed':
         annot_file = pd.read_csv(path, sep='\t', comment='t', header=None)
         header = ['chr', 'start', 'end', 'ORF', 'score', 'strand', 'thickStart', 'thickEnd',
                   'itemRgb', 'blockCount', 'blockSizes', 'blockStarts']
         annot_file.columns = header[:len(annot_file.columns)]
    elif 'gtf' in path.split('.')[-1]:
        with warnings.catch_warnings():
            warnings.filterwarnings('ignore')
            gtf = gtfparse.read_gtf(path)
            annot_file = gtf[gtf['feature'] == 'transcript'][['seqname', 'start', 'end', 'gene_name', 'strand']]
            annot_file.columns = ['chr', 'start', 'end', 'ORF', 'strand']
    else:
        raise ValueError('Annotation data file format not understood. Pass either bed or gtf.')
    return annot_file


def arg_parse(args):
    parser = argparse.ArgumentParser(description='Find linear correaltion between two BigWig data sets.')
    parser.add_argument('input_data', nargs='+', metavar='data', type=str,
                        help='List with input sequencing signals as a relative path to your current directory.'
                             'It is expected that the number of input_data values is either 2 or 4. If it is 2, the '
                             'first one is compared with the second one (e.g. wt and mutant). If it is 4, the '
                             'first two are added together and the second two are added together (e.g.'
                             'wt +strand + wt -strand and mutant +strand + mutant -strand).')
    parser.add_argument('--annotation', type=str, required=True,
                        help='Annotation file. Can be bed or gft.')
    parser.add_argument('--name', '-n', action='append', type=str, help='Names of the data sets.')
    parser.add_argument('--title', '-t', type=str, help='Plot title.')
    parser.add_argument('--distinguish_strand', action='store_true', dest='do_distinguish_strand',
                        help='If set, two consecutive strands are considered as + and - (in that order)')
    parser.add_argument('--percentile', type=float, default=100.,
                        help='Only keep the lower percentile.')
    parser.add_argument('--save_fig', action='store_true', dest='save_fig',
                        help='If set, figure is set to file at save_path')
    parser.add_argument('--save_path', type=str, default='linchange.png',
                        help='Figures is saved under this file name if save_fig is set.')
    parsed_args = parser.parse_args(args)
    return parsed_args


def main(args):
    input_data = args.input_data
    if len(input_data) not in [2, 4]:
        raise ValueError('Number of passed input paths must be 2 or 4.')
    annot_path = args.annotation
    names = args.name if args.name is not None else ['WT', 'Mutant']
    title = args.title if args.title is not None else 'Median Correlation Between WT and Mutant'
    title = title.replace('\\n', '\n')
    do_distinguish_strand = args.do_distinguish_strand
    percentile = args.percentile
    save_fig = args.save_fig
    save_path = args.save_path

    annot_file = annotation_loading(annot_path)
    if do_distinguish_strand:
        values = np.zeros((len(input_data) // 2, len(annot_file.index)))
    else:
        values = np.zeros((len(input_data), len(annot_file.index)))
    for i_sample, file_path in enumerate(input_data):
        bw_file = pyBigWig.open(file_path)
        is_minus = i_sample % 2 == 1
        for i_gene, (_, gene) in enumerate(annot_file.iterrows()):
            level = np.nanmean(bw_file.values(gene['chr'], int(gene['start']), int(gene['end'])))
            if do_distinguish_strand and is_minus and gene['strand'] == '-':
                values[i_sample // 2, i_gene] = level
            elif do_distinguish_strand and not is_minus and gene['strand'] == '+':
                values[i_sample // 2, i_gene] = level
            elif not do_distinguish_strand:
                values[i_sample, i_gene] = level
            else:
                pass

    max_value = np.percentile(values, percentile, axis=1)
    values_mask = np.all(values <= max_value.reshape(-1, 1), axis=0)
    values = values[:, values_mask]

    reg = LinearRegression(fit_intercept=True).fit(values[0].reshape(-1, 1), values[1].reshape(-1, 1))
    r2 = reg.score(values[0].reshape(-1, 1), values[1].reshape(-1, 1))
    c_vec = gaussian_kde(values)(values)
    sc = plt.scatter(*values, c=c_vec, marker='.')
    axes = plt.gca()
    x_val = np.arange(0, values[0].max(), values[0].max() * .05)
    plt.plot(x_val, reg.coef_[0][0] * x_val + reg.intercept_)
    plt.plot(x_val, x_val, '--', color='grey')
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
    if save_fig:
        plt.savefig(save_path)
        plt.close('all')
    else:
        plt.show()


if __name__ == '__main__':
    main(arg_parse(sys.argv[1:]))

