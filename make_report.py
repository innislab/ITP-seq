#!/usr/bin/env python3

import sys, os
from datetime import datetime


def read_log_json(json_files):
    import json
    import pandas as pd

    if isinstance(json_files, str):
        from glob import glob

        json_files = glob(json_files)
    stats = {
        f.rsplit('/', 1)[-1].rsplit('.')[0]: json.load(open(f)) for f in json_files
    }
    return (
        pd.DataFrame(stats)
        .T.sort_index()
        .convert_dtypes()
        # .select_dtypes(exclude='object')
    )


def plot_to_html(fig, format='svg', size=None, outfile=None, last=' '):
    import io

    if size:
        fig.set_size_inches(*size)
    fig.tight_layout()
    if outfile:
        print(f'{last} │ ├╴saving volcano: {outfile}.svg')
        fig.savefig(outfile)
        print(f'{last} │ ╰╴saving volcano: {outfile}.png')
        fig.savefig(outfile)
    if format == 'svg':
        out = io.StringIO()
        fig.savefig(out, format='svg')
        return out.getvalue()
    if format == 'png':
        import base64

        out = io.BytesIO()
        fig.savefig(out, format='png')
        b64 = base64.b64encode(out.getvalue()).decode('utf-8')
        return f'<img src=\'data:image/png;base64,{b64}\'>'


# def compute_stats():
#     pass

# def compute_motifs(path, cond=None, ref='noa', motifs=[], how='aa', force=False):
#     dfs_motifs = {}
#     if how == 'aa':
#         prefix = '_aa'
#     else:
#         prefix = ''
#     for motif in motifs:
#         table_filename = f'{path}/table_{cond}{prefix}_{motif}.csv'
#         if not force and os.path.isfile(table_filename):
#             print(f'  ╰╴reading cached "{table_filename}"')
#             dfs_motifs[motif] = pd.read_csv(table_filename, index_col=0, na_values=[''],
#                                              keep_default_na=False).convert_dtypes()
#         else:
#             print(f'  ╰╴computing "{table_filename}"')
#             dfs_motifs[motif] = ad.calc_ratio(ad.compute_all(files, names=names,
#                                                              pos=motif, how=how),
#                                               cond=cond, ref=ref)
#             dfs_motifs[motif].to_csv(table_filename)


def make_report(
    path, cond, ref, how='aa', motifs=None, esite=7, kwargs={}, force=False, last=' '
):
    '''Generates report using all files from <cond> and <ref> for each motif in <motifs>

      path:   path of the files
      cond:   name of the condition to analyze (e.g. ab1)
       ref:   name of the reference (e.g. "noa")
    motifs:   list of motifs to analyze'''
    import pandas as pd
    import matplotlib.pyplot as plt
    import seaborn as sns
    import jinja2
    from glob import glob
    import re, os
    from itertools import combinations
    import analyze_data as ad
    from collections import defaultdict

    esite = int(esite)

    if how in ('aa', 'aax'):
        prefix = '_' + how
        if not motifs:
            motifs = [
                f'{esite-1}',
                f'{esite}',
                f'{esite+1}',
                f'{esite+2}',
                f'{esite}-{esite+1}',
                f'{esite-1}-{esite+1}',
                f'{esite-1}-{esite+2}',
            ]
            # motifs = ['7-8', '6-8', '6-9']
    else:
        prefix = ''
        if not motifs:
            # motifs = ['21-23', '24-26', '27-29', '21-29']
            motifs = [
                '18-20',
                '21-23',
                '24-26',
                '27-29',
                '21-29',
                '18-26',
                '21-26',
                '18-29',
            ]
            x = esite * 3
            motifs = [f'{x}-{x+2}', f'{x+3}-{x+5}', f'{x+6}-{x+8}', f'{x}-{x+8}']
            motifs = [
                f'{x-3}-{x-1}',
                f'{x}-{x+2}',
                f'{x+3}-{x+5}',
                f'{x+6}-{x+8}',
                f'{x}-{x+8}',
                f'{x-3}-{x+5}',
                f'{x}-{x+5}',
                f'{x-3}-{x+8}',
            ]

    template = '''
<!DOCTYPE html>
<html>
<head lang="en">
    <meta charset="UTF-8">
    <link rel="stylesheet" href="typography.css">
    <style>
            @page {
                   /*size: A4 landscape;*/
                   size: A4;
                   margin: 0in 0.44in 0.2in 0.44in;
            }
            .table-style {
                   width: 100%;
            }
            .vert {
                  transform: rotate(270deg);
            }
            .largetxt {font-size: 200%;}
            .keep-together { page-break-inside: avoid;}
            .break-before { page-break-before: always;}
            .break-after { page-break-after: always;}
    </style>
    <title>{{ title }}</title>
</head>
<body>
    <h1>Inverse Toe Printing Seq (ITP-Seq) Report</h1>

       Report generated: {{date}}
       <br/>
       Working directory: <code>{{cwd}}</code>
       <br/>
       Command: <code>{{cmd}}</code>

    {% if stats_table or stats_plot %}
    <h2>Statistics</h2>
     Sequence length {{max_len}}
     {{ stats_table }}
     {{ stats_plot }}
    {% endif %}

    {% if lengths_plot %}
     <h2>distribution of reads</h2>
     {{ lengths_plot }}
    {% endif %}

    {% if dfs_motifs %}
    <div class='keep-together'>
     <h2>Motifs</h2>
      {% for motif, df in dfs_motifs.items() %}
       <div class='keep-together'>
       <h3>positions {{ motif }}</h3>
       {{ df }}
       </div>
       {% if motif in volcano_plots %}
        {{ volcano_plots[motif] }}
       {% endif %}
      {% endfor %}
    </div>
    {% endif %}

    {% if heatmaps %}
     <div class='keep-together'>
     <h2>heatmaps</h2>
        <table>
        <thead>
        <tr>
          <th></th>
          {% for col in hm_cols %}
          <th class='largetxt'>{{ col }}</th>
          {% endfor %}
        <tr>
        </thead>
        <tbody>
        {% for row in hm_rows %}
        <tr>
            <th class='vert largetxt'>{{ row }}</th>
            {% for col in hm_cols %}
             <td>
             {% if row in heatmaps and col in heatmaps[row] %}
               {{ heatmaps[row][col] }}
             {% endif %}
             </td>
            {% endfor %}
        </tr>
        {% endfor %}
        </tbody>
        </table>
    {% endif %}
    </div>

</body>
</html>
'''

    rtemplate = jinja2.Environment().from_string(template)

    df_stats = read_log_json(path + '/*.json')

    df_html = (
        df_stats.select_dtypes(exclude='object')
        .drop('MAX_LEN', axis=1)
        .style.set_table_attributes('class="table-style"')
        .to_html()
    )

    ax1 = df_stats.select_dtypes(exclude='object').T.drop('MAX_LEN').plot.bar(rot=30)
    ax1.set_ylabel('number of reads')

    df_lengths = df_stats['lengths'].apply(pd.Series).T
    df_lengths.index = df_lengths.index.astype(int)
    ax2 = df_lengths.sort_index().loc[:90].plot()
    ax2.set_ylabel('number of reads')
    ax2.figure.set_size_inches(8, 3)

    prefix2 = '_aa' if prefix == '_aax' else prefix

    files = sorted(glob(path + f'/*_*{prefix2}.processed.txt'))
    # print(path+'/nnn15_*{prefix}.processed.txt', files)
    names = [re.search(f'([^_]+){prefix2}.processed.txt$', f).group(1) for f in files]

    ############
    #  motifs  #
    ############

    dfs_motifs = {}
    DE_dfs_motifs = {}
    for motif in motifs:
        table_filename = f'{path}/table_{cond}{prefix}_{motif}.csv'
        if not force and os.path.isfile(table_filename):
            print(f'{last} ├╴reading cached "{table_filename}"')
            dfs_motifs[motif] = pd.read_csv(
                table_filename, index_col=0, na_values=[''], keep_default_na=False
            ).convert_dtypes()
        else:
            print(f'{last} ├╴computing table "{table_filename}"')
            dfs_motifs[motif] = ad.calc_ratio(
                ad.compute_all(files, names=names, pos=motif, how=how),
                cond=cond,
                ref=ref,
            ).convert_dtypes()
            dfs_motifs[motif].to_csv(table_filename)

        DE_table_filename = f'{path}/DE_table_{cond}{prefix}_{motif}.csv'
        if not force and os.path.isfile(DE_table_filename):
            print(f'{last} ├╴reading cached "{DE_table_filename}"')
            DE_dfs_motifs[motif] = pd.read_csv(
                DE_table_filename, index_col=0, na_values=[''], keep_default_na=False
            ).convert_dtypes()
        else:
            print(f'{last} ├╴computing DE table "{DE_table_filename}"')
            DE_dfs_motifs[motif] = ad.DE(
                dfs_motifs[motif].rename_axis('id'), cond=cond, ref=ref
            )
            DE_dfs_motifs[motif].to_csv(DE_table_filename)

    ############
    # heatmaps #
    ############

    pos = list(combinations([6, 7, 8, 9], r=2))

    heatmaps = defaultdict(lambda: defaultdict(dict))

    for a, b in pos:
        sep = '-' if b == a + 1 else ','
        table_filename = f'{path}/table_{cond}{prefix}_{a}{sep}{b}.csv'
        if not force and os.path.isfile(table_filename):
            print(f'{last} ├╴reading cached "{table_filename}"')
            d = pd.read_csv(
                table_filename, index_col=0, na_values=[''], keep_default_na=False
            ).convert_dtypes()
        else:
            print(f'{last} ├╴calculating heatmap for {a},{b}')
            d = ad.calc_ratio(
                ad.compute_all(files, pos=f'{a},{b}', how=how), cond=cond, ref=ref
            )
            d = d.convert_dtypes()
            d.to_csv(table_filename)

        ax = ad.mk_heatmap(d.astype(float), pos=f'{a}-{b}', vmax=0.8, how=how)
        ax.figure.set_size_inches(4, 3)
        print(f'{last} │ ├╴saving heatmap: {path}/heatmap_{cond}{prefix}_{a}-{b}.svg')
        ax.figure.savefig(path + f'/heatmap_{cond}{prefix}_{a}-{b}.svg')
        print(f'{last} │ ╰╴saving heatmap: {path}/heatmap_{cond}{prefix}_{a}-{b}.png')
        ax.figure.savefig(path + f'/heatmap_{cond}{prefix}_{a}-{b}.png')

        a = ad.get_ribosome_site(a, P_site=esite + 1)
        b = ad.get_ribosome_site(b, P_site=esite + 1)

        heatmaps[a][b] = plot_to_html(ax.figure, 'svg')

        hm_rows = list(heatmaps.keys())
        hm_cols = list(heatmaps[hm_rows[0]].keys())

    def html_style(d):
        cmap = sns.color_palette("vlag", as_cmap=True)
        return (
            d.head(20)
            .style.format(precision=2, na_rep='−')
            .set_properties(**{'text-align': 'right'})
            .background_gradient(cmap=cmap, axis='index', vmin=0)
            .to_html()
        )

    html = rtemplate.render(
        {
            'date': datetime.today().strftime('%c'),
            'cwd': os.getcwd(),
            'cmd': ' '.join(sys.argv),
            'args': str(kwargs),
            'max_len': df_stats['MAX_LEN'].max(),
            'stats_table': df_html,
            'stats_plot': plot_to_html(ax1.figure, 'svg'),
            'lengths_plot': plot_to_html(ax2.figure, 'svg'),
            'heatmaps': heatmaps,
            'hm_rows': hm_rows,
            'hm_cols': hm_cols,
            #'dfs_motifs': {k: html_style(d) for k,d in dfs_motifs.items()},
            'dfs_motifs': {
                k: html_style(
                    d.fillna(0)
                    .convert_dtypes()
                    .join(DE_dfs_motifs[k].select_dtypes('number').astype(float))
                )
                for k, d in dfs_motifs.items()
            },
            'volcano_plots': {
                k: plot_to_html(
                    ad.mk_volcano(
                        d,
                        query='auto',
                        motif=kwargs['volcano_motif'],
                    ).figure,
                    'png',
                    outfile=f'{path}/volcano_{cond}{prefix}_{k}',
                    last=last,
                )
                for k, d in DE_dfs_motifs.items()
            },  # if how=='aa' else {}
        }
    )

    with open(f'report{prefix}_{cond}.html', 'w') as f:
        print(f'{last} ├╴writing report: "report{prefix}_{cond}.html"')
        f.write(html)
    with open(f'report{prefix}_{cond}.pdf', 'wb') as f:
        print(f'{last} ╰╴writing report: "report{prefix}_{cond}.pdf"')
        import pdfkit

        f.write(pdfkit.from_string(html))

    # from IPython.core.display import HTML
    # HTML(html)
    return heatmaps, html


def main():
    import argparse
    import logging

    parser = argparse.ArgumentParser(
        description='generate report for inverse toe-printing',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument('path', help='path of the json files', default='.', nargs='?')

    parser.add_argument(
        '-t', '--type', help='type of analysis (nuc, aa, aax)', default='aa'
    )

    parser.add_argument('-e', '--esite', help='position of E site (codons)', default=7)

    parser.add_argument('-c', '--cond', help='condition', default='')
    parser.add_argument('-r', '--ref', help='reference', default='noa')

    parser.add_argument(
        '-o',
        '--outdir',
        help='directory for the output files',
        dest='outdir',
        default='.',
    )
    parser.add_argument('-l', '--log', help='log file', dest='log', default=None)
    parser.add_argument('-v', '--verbose', action='count', default=0)

    parser.add_argument(
        '-f', '--force', help='force recomputing the tables', action='store_true'
    )

    parser.add_argument('--volcano_motif', help='highlight regex for volcano plot')

    args = parser.parse_args()
    print(args)

    # logging.basicConfig(filename='example.log', encoding='utf-8', level=logging.DEBUG)

    if not args.cond:
        from glob import glob
        import re

        conds = set(
            re.findall(r'[a-zA-Z]{2,5}(?=\d+)', f.rsplit('/', 1)[-1])[-1]
            for f in glob(args.path + '/*.json')
        ) - set([args.ref])
        print('Detected conditions:', ', '.join(sorted(conds)))
    else:
        conds = args.cond.split(',')

    for i, cond in enumerate(conds):
        last = '│'
        if i == len(conds) - 1:
            last = ' '
            print(f'╰╴computing report for "{cond}"')
        else:
            print(f'├╴computing report for "{cond}"')
        make_report(
            args.path,
            cond,
            args.ref,
            how=args.type,
            esite=args.esite,
            kwargs=vars(args),
            force=args.force,
            last=last,
        )


if __name__ == '__main__':
    main()
