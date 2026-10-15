#!/usr/bin/env python3
import argparse
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
from ncrna_stats import save_figure


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    table = pd.read_csv(args.input, sep='\t')
    fig, ax = plt.subplots(figsize=(9, 5))
    if not table.empty:
        level = 'class' if 'class' in table else 'type'
        matrix = table.pivot_table(index='div', columns=level, values='percentage', aggfunc='sum', fill_value=0).sort_index()
        matrix.plot.bar(stacked=True, ax=ax, width=1)
    else:
        ax.text(0.5, 0.5, 'No repeat landscape annotations', ha='center', transform=ax.transAxes)
    ax.set(xlabel='Repeat divergence bin', ylabel='Genome coverage (%)', title='Repeat landscape')
    save_figure(fig, args.output)


if __name__ == '__main__':
    main()
