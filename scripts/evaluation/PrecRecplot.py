#!/usr/bin/env python
# coding: utf-8

# As usual imports first
import matplotlib.pyplot as pyplot
import pandas as pandas
import argparse

parser=argparse.ArgumentParser()

parser.add_argument('--table', help='TSV with QC_Score, Precision, Recall per line')
parser.add_argument('--title', help='Title of the plot to be drawn')

args = parser.parse_args()

def prec_rec_plot(ids, prec, rec, best, title="Precision/Recall Plot", savefig=True):
    """
    Draw a connected scatterplot (Precision-Recall plot).
    Parameters:
    ids: List of svim quality values
    prec: List of precision values
    rec: List of recall values
    title: Title of the plot(default: 'Precision/Recall Plot')
    savefig: If true, plot is saved with name %title
    """
    fig, ax = pyplot.subplots(figsize=(10,8))
    ax.set_title(title)
    ax.plot(rec, prec, alpha=0.2)
    ax.set_xlabel('Recall')
    ax.set_ylabel('Precision')
    #ax.set_xlim(0,100)
    #ax.set_ylim(0,100)
    ax.grid()
    for x,y,l in zip(rec, prec, ids):
        if(l == best["QC_Score"]):
            ax.scatter(x,y,marker='X', label=l, c='green')
        else:
            ax.scatter(x,y,marker='^', label=l, alpha=0.2, color='blue')
        ax.annotate('%s' % l, xy=(x,y))
    if(savefig):
        fig.savefig(title)


# import precision-recall files to plot
openedFile = pandas.read_csv(args.table, delimiter='\t', names=["QC_Score", "Precision", "Recall", "F1"])

# normalize values
openedFile.Precision=openedFile.Precision.multiply(100)
openedFile.Recall=openedFile.Recall.multiply(100)
openedFile.F1=openedFile.F1.multiply(100)

# compute best cutoff
bestResult = openedFile.loc[openedFile["F1"].idxmax()]

print("Best qc value: ", bestResult["QC_Score"])

# mark best result in plot

print("Best solution for ", args.title, ":\n", openedFile.loc[openedFile["F1"].idxmax()])

# plot
prec_rec_plot(openedFile["QC_Score"], openedFile["Precision"], openedFile["Recall"], bestResult, title=args.title)
