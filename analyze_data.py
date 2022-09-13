#!/usr/bin/env python3

'''Performs computing of the cond/ref counts and DE tables, produces the graphs'''

import re
import sys
import json
from glob import glob
import warnings

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt


def ranges(s=''):
    '''transforms 'x-y' into slice(x, y, None)'''
    return [
        slice(y[0], y[1] + 1) if len(y := list(map(int, x.split('-')))) > 1 else y[0]
        for x in s.split(',')
    ]


def compute_counts(seqs_series, pos='7-9'):
    slices = ranges(pos)
    from functools import reduce

    return reduce(
        lambda x, y: x + y, map(lambda x: seqs_series.str[x], slices)
    ).value_counts()


def read_file_as_series(filename, remove_stops=False):
    with open(filename) as f:
        if remove_stops:
            return pd.Series(
                [
                    line.rstrip('\n')
                    for line in f.readlines()
                    if '*' not in line.rstrip()[:-5]
                ]
            )
        return pd.Series([line.rstrip('\n') for line in f.readlines()])


def compute_all(files, names=None, pos='7-9', how='aa'):

    if how in ('aa', 'aax'):
        prefix = '_aa'
    else:
        prefix = ''

    remove_stops = how == 'aax'

    if not names:
        names = [
            re.search(f'([^/_]+){prefix}.processed.txt$', f).group(1) for f in files
        ]
    df = pd.concat(
        dict(
            zip(
                names,
                (
                    compute_counts(
                        read_file_as_series(f, remove_stops=remove_stops), pos=pos
                    )
                    for f in files
                ),
            )
        ),
        axis=1,
    )
    # df = pd.concat((compute_counts(read_file_as_series(f), pos=pos) for f in files),
    #               keys=names,
    #               axis=1)
    return df  # /df.sum()*1000


def calc_ratio(d, cond='tcx', ref='noa', norm=True, factor=1e6):
    d = d.copy()
    cols_cond = d.filter(regex=cond + r'\d+').columns
    cols_ref = d.filter(regex=ref + r'\d+').columns
    if norm:
        d[cond] = d[cols_cond].div(d[cols_cond].sum()).mean(axis=1) * factor
        d[ref] = d[cols_ref].div(d[cols_ref].sum()).mean(axis=1) * factor
    else:
        d[cond] = d[cols_cond].mean(axis=1)
        d[ref] = d[cols_ref].mean(axis=1)
    d['ratio'] = d[cond] / d[ref]
    return d.sort_values(by='ratio', ascending=False)


def DE(df, cond='tcx', ref='noa', join=False, quiet=True):

    df = df.rename_axis('id')

    if 'diffexpr.py_deseq' not in sys.modules:
        print('loading DESeq2')
    from diffexpr.py_deseq import py_DESeq2

    # get xxx000 columns
    conds = df.filter(regex=r'^[a-z]{2,5}\d+$').columns

    # split into sample (xxx) + replicate (000)
    samples, replicates = zip(*(re.search(r'(\D+)(\d+)', c).groups() for c in conds))

    # make sample_df for DESeq2
    r = re.compile(r'^(\D+)(\d+)$')
    samples, replicates = zip(*map(lambda x: r.search(x).groups(), conds))
    sample_df = pd.DataFrame(
        {'samplename': conds, 'sample': samples, 'replicate': replicates}
    ).set_index('samplename', drop=False)

    with warnings.catch_warnings():
        if quiet:
            warnings.simplefilter('ignore')
        dds = py_DESeq2(
            count_matrix=df[conds].fillna(0).reset_index(),
            design_matrix=sample_df,
            design_formula='~ replicate + sample',
            gene_column='id',
        )  # <- telling DESeq2 this should be the gene ID column

        dds.run_deseq(fitType='mean')
        dds.get_deseq_result(contrast=['sample', cond, ref])
        res = dds.deseq_result

    res['log10pvalue'] = -np.log10(res['pvalue'])
    res['log10padj'] = -np.log10(res['padj'])

    if join:
        return df.join(res)

    return res


def mk_volcano(
    df,
    query=None,
    motif=None,
    ax=None,
    x='log2FoldChange',
    y='log10pvalue',
    query_color='#BC0909',
    motif_color='#EDEA20',
    color='k',
    params={'alpha': 0.2},
    annotate=True,
    text_stroke=False,
    outfile=None,
):
    if not ax:
        f, ax = plt.subplots()

    if not annotate:
        annotate = ''

    df.plot.scatter(x=x, y=y, color=color, ax=ax, **params)

    dfs = []
    colors = []

    if query == 'auto':
        thresh = {2: (0.05, 0.2), 3: (0.05, 1), 4: (0.2, 2)}
        padj, l2fc = thresh.get(len(df.index[0]), (0.05, 1))
        query = f'(padj < {padj}) & (log2FoldChange > {l2fc})'

    if query:
        df_query = df.query(
            query, engine='python'
        )  # engine='python' to prevent bug with <NA> type
        df_query.plot.scatter(x=x, y=y, color=query_color, ax=ax, **params)
        if (annotate == True) or ('query' in annotate):
            dfs.append(df_query)
            colors.append(query_color)

    if motif:
        df_motif = df[df.index.str.match(motif)]
        df_motif.plot.scatter(x=x, y=y, color=motif_color, ax=ax, **params)
        if (annotate == True) or ('motif' in annotate):
            dfs.append(df_motif)
            colors.append(motif_color)

    if annotate:
        import matplotlib.patheffects as PathEffects

        for d, c in zip(dfs, colors):
            for name, s in d.iterrows():
                txt = ax.annotate(
                    name,
                    (s[x], s[y]),
                    xytext=(-5, 0),
                    textcoords='offset pixels',
                    ha='right',
                    va='center',
                )
                if text_stroke:
                    stroke_color = c if isinstance(text_stroke, bool) else text_stroke
                    txt.set_path_effects(
                        [PathEffects.withStroke(linewidth=2, foreground=stroke_color)]
                    )

    if outfile:
        ax.figure.savefig(outfile)

    return ax


def get_ribosome_site(pos, P_site=8):
    positions = {-1: 'E-site', 0: 'P-site', 1: 'A-site'}
    pos = int(pos) - P_site
    return positions.get(pos, f'{"+" if pos>0 else "−"}{abs(pos)}')


def from_ribosome_sites(sites, P_site=8):
    positions = {
        'E': -1,
        'P': 0,
        'A': 1,
    }
    sites = re.split(r'-(?!\d)|[,;|]', sites)
    return '-'.join(map(str, (positions.get(s, s) + P_site for s in sites)))


def mk_heatmap(table, pos='7-8', vmax=None, ax=None, how='aa'):
    if how in ('aa', 'aax'):
        # amino acids color codes
        h = '#74C170'  # hydrophobic
        s = '#E4DF51'  # special
        n = '#7094C1'  # negatively charged
        p = '#DB4755'  # positively charged
        o = '#AF70C1'  # polar

        aacolors = {
            'A': h,
            'C': s,
            'D': n,
            'E': n,
            'F': h,
            'G': s,
            'H': p,
            'I': h,
            'K': p,
            'L': h,
            'M': h,
            'N': o,
            'P': s,
            'Q': o,
            'R': p,
            'S': o,
            'T': o,
            'V': h,
            'W': h,
            'Y': h,
            '*': 'grey',
        }

        aa_order = list('HRKDESTNQCGPAVILMFYW*')
    else:
        aacolors = {'A': '#008000', 'C': '#0000FF', 'G': '#FFA806', 'T': '#FF0000'}
        aa_order = list('ACGT')

    table = table.copy()  # avoid modifying the input object
    # keep only motifs without gap
    table = table.loc[~table.index.str.contains(' ')]

    names = list(map(get_ribosome_site, pos.split('-')))

    table.index = pd.MultiIndex.from_tuples(
        [tuple(x) for x in table.index], names=names
    )

    if not ax:
        f, ax = plt.subplots()

    t = table['ratio'].unstack()
    if len(t) == 0:
        return ax

    # d = np.log2(t.loc[aa_order, aa_order]) # can't use loc in case aa are missing, must reindex
    d = np.log2(t.reindex(index=aa_order, columns=aa_order))
    # d_ref = d.copy()

    if not vmax:
        vmax = max(-d.values.min(), d.values.max())
    sns.heatmap(
        d,
        yticklabels=True,
        xticklabels=True,
        center=0,
        cmap='vlag',
        square=True,
        linecolor='#CECECE',
        linewidths=0.05,
        cbar_kws={'shrink': 0.5, 'ticks': [-0.5, -0.25, 0, 0.25, 0.5]},
        vmax=vmax,
        vmin=-vmax,
        ax=ax,
    )
    ax.figure.axes[1].set_ylabel('log2(enrichment)')
    for s in ax.spines.values():
        s.set_visible(True)

    # horizontal centered yticklabels
    ax.set_xticklabels(
        ax.get_xticklabels(),
        rotation=0,
        ha='center',
        fontdict={'fontfamily': 'monospace'},
    )
    ax.set_yticklabels(
        ax.get_yticklabels(),
        rotation=90,
        ha='center',
        fontdict={'fontfamily': 'monospace'},
    )

    ## AAs colors
    for a in [ax.xaxis, ax.yaxis]:
        for t in a.get_ticklabels():
            aa = t.get_text()
            bbox = dict(
                boxstyle='square', pad=0.15, ec="none", fc=aacolors[aa], alpha=0.5
            )
            t.set_bbox(bbox)
            t.get_bbox_patch().set_width(5)

    ax.figure.canvas.draw()

    for t in ax.get_xticklabels():
        t.get_bbox_patch().set_width(5)
    for t in ax.get_yticklabels():
        t.get_bbox_patch().set_width(5)

    ax.invert_yaxis()

    for t in ax.get_xticklabels():
        t.get_bbox_patch().set_width(5)
    for t in ax.get_yticklabels():
        t.get_bbox_patch().set_width(5)

    return ax


def read_log_json(json_files):
    if isinstance(json_files, str):
        json_files = glob(json_files)
    stats = {
        f.rsplit('/', 1)[-1].rsplit('.')[0]: json.load(open(f)) for f in json_files
    }
    return (
        pd.DataFrame(stats)
        .T.sort_index()
        .convert_dtypes()
        .select_dtypes(exclude='object')
    )
