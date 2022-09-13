#!/usr/bin/env python3

DEFAULTS = dict(
    a1='GTATAAGGAGGAAAAAAT',
    a2='GGTATCTCGGTGTGACTG',
)

# fmt: off
codon_table = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
    'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W',
}
# fmt: on

CONTAMINANTS = 'TCCAACATGCTGAGC|^GATCCTTTTTA'

#####

import numpy as np

to_aa = codon_table.copy()
to_aa['   '] = ' '
to_aa['\n'] = '\n'


def FastQIterator(filename):
    """
    Reads 'filename' as a fastq sequence and yield only the sequence and quality lines


    This does not perform any check on the validity od the file!
    """
    from itertools import islice

    with open(filename) as f:
        while True:
            chunk = list(islice(f, 4))
            if not chunk:
                break
            _, seq, _, score = chunk
            yield (seq[:-1], score[:-1])


def parse_trim_filter_fastq(
    filename=None,
    a1=DEFAULTS['a1'],  # adaptators
    a2=DEFAULTS['a2'],
    mm1=2,  # allowed number of mismatches
    mm2=2,
    limit=None,  # max sequences to read
    trim_start=True,
    min_seq_len=3 * 3,
    max_seq_len=3 * 15,
    round_codons=True,
    quality=30,
    contaminants=CONTAMINANTS,
    # below options for parallel processing  ## FIXME
    parallel=False,
    records=None,
    chunks=1000,  # chunks for parallel processing
    **kwargs,  # catch all other parameters
):
    '''
    Takes a 'filename' as input (fastq format), loops over the fastq using
    FastQIterator and performd several checks to extract valid ITP sequences.
    Also computes various statistics on the dataset.

    Parameters:
        filename:  fastq file to process
              a1:  left adaptator
              a2:  right adaptator
             mm1:  tolerated mismatches in a1
             mm2:  tolerated mismatches in a2
           limit:  max number of sequences to process (useful for quick tests)
      trim_start:  removes ATG
     min_seq_len:  minimum length of the matched sequence to keep
     max_seq_len:  maximum length of the matched sequence to keep
    round_codons:  trim the last nucleotides of the match to have full codons
    '''

    from collections import Counter
    import regex

    stats = Counter()
    minquals = Counter()
    lengths = Counter()

    if not min_seq_len:
        min_seq_len = 0

    if not max_seq_len:
        max_seq_len = float('inf')

    seqs = []
    MAX_LEN = 0
    extras = []

    def parse(record):
        keep = True
        seq, qual = record
        m = REG.search(seq)

        stats['total_sequences'] += 1

        if m:
            subseq = m.group()
            subqual = qual[slice(*m.span())]
            minqual = ord(min(subqual)) - 33

            # span, subsequence, quality substring, min quality
            # out.append((*m.span(), subseq, subqual, minqual))
            L = len(subseq) + (3 if trim_start else 0)
            lengths[L] += 1

            if round_codons:
                MAX = len(subseq) // 3 * 3
                extra = len(subseq) - MAX
                subseq = subseq[:MAX]
                stats[f'extra{extra}'] += 1
            else:
                MAX = len(subseq)

            # As in first 18 nucleotides
            if subseq[:18].count('A') >= 18:
                stats['A_stretch'] += 1
                keep = False

            if REG_conta.search(subseq):
                stats['contaminant'] += 1
                keep = False

            minquals[minqual] += 1

            if minqual < quality:
                stats['lowqual'] += 1
                keep = False
            if 'N' in subseq:
                stats['subseq_N'] += 1
                keep = False

            if L < min_seq_len:
                stats['tooshort'] += 1
                keep = False
            elif L > max_seq_len:
                stats['toolong'] += 1
                keep = False

            if keep:
                nonlocal MAX_LEN
                MAX_LEN = max(MAX, MAX_LEN)
                nonlocal extras
                extras.append(extra)
                return subseq

        else:
            stats['noadaptor'] += 1
        if 'N' in seq:
            stats['seq_N'] += 1

    if not filename:
        pass
    else:
        records = FastQIterator(filename)

    if limit:
        from itertools import islice

        records = islice(records, int(limit))
    if parallel:
        records = list(records)

    start = 'ATG' if trim_start else ''

    REG = regex.compile(
        f'(?<=(?:{a1}){{s<={mm1}}}{start}).+(?=(?:{a2}){{s<={mm2}}})', regex.BESTMATCH
    )
    REG_conta = regex.compile(f'(?:{contaminants}){{e<={2}}}')

    if parallel is True:  # apply by chunks             ## FIXME needs improvement
        from itertools import islice, chain
        from multiprocessing import Pool, cpu_count

        global parse2

        def parse2(x):
            print('started pool')
            return parse_trim_filter_fastq(
                records=x,
                a1=a1,
                a2=a2,
                mm1=mm1,
                mm2=mm2,
                limit=None,  # max sequences to read
                trim_start=trim_start,
                min_seq_len=min_seq_len,
                max_seq_len=max_seq_len,
                round_codons=round_codons,
                # below options for parallel processing
                parallel='raw',
            )

        NB_CPU = cpu_count() - 4

        import math

        chunks = math.ceil(len(records) / NB_CPU)
        L = len(records)

        records = [records[i : i + chunks] for i in range(0, len(records), chunks)]

        print(
            f'multiprocessing {L} records by {len(records)} chunks of {chunks} records'
        )
        with Pool(NB_CPU) as pool:
            result = pool.map(parse2, records)
        seqs, stats, MAX_LEN, dics = zip(*result)

        labels = ('min_quality_counts', 'lengths')
        seqs = list(chain.from_iterable(seqs))

        stats = sum(stats, Counter())
        stats['MAX_LEN'] = max(MAX_LEN)

        for k, d in zip(labels, dics):
            stats[k] = sum(d, Counter())
        return seqs, stats

    for record in records:
        seq = parse(record)
        if seq:
            seqs.append(seq)

    if parallel == 'raw':
        return seqs, stats, MAX_LEN, (minquals, lengths)

    stats['MAX_LEN'] = MAX_LEN
    stats['min_quality_counts'] = minquals
    stats['lengths'] = lengths

    return seqs, stats, extras


def parse_all(files=None, pattern=None, save=False, outdir=None, **kwargs):
    '''
    Apply parse_trim_filter_fastq on multiple files
    '''

    from multiprocessing import Pool, cpu_count

    if pattern:
        from glob import glob

        files = glob(pattern)

    global f

    def f(filename, *args, **kwargs):
        seqs, stats, extras = parse_trim_filter_fastq(filename, *args, **kwargs)
        if save:
            export_data(filename, seqs=seqs, stats=stats, extras=extras, outdir=outdir)
        return seqs, stats, extras

    with Pool(cpu_count() - 2) as pool:
        result = pool.map(f, files)

    short = list(map(lambda x: x.rsplit('/', 1)[-1], files))

    return dict(zip(short, result))


def simple_graph(
    counter, vmin=None, vmax=None, start='auto', trim_below=0.02, ticks=10, wrap=80
):
    import math

    bars = ' _▁▂▃▄▅▆▇█'[1:]
    L = max(counter)
    a = np.array([counter[i] for i in range(L)])
    MAX = a.max()

    # floating vmin/vmax = relative to MAX
    if isinstance(vmin, float):
        vmin = int(vmin * MAX)
    if isinstance(vmax, float):
        vmax = int(vmax * MAX)

    # trim below percent MAX
    if isinstance(trim_below, float):
        trim_below = int(trim_below * MAX)

    if trim_below:
        l = len(a)
        a[a < trim_below] = 0
        a = np.trim_zeros(a, trim='b')
        L -= l - len(a)

    if start == 'auto':
        start = np.min(np.nonzero(np.hstack((a, 1)))) // ticks * ticks
    L -= start
    a = a[start:]

    if not vmin:
        vmin = 0
    if not vmax:
        vmax = MAX
    # handle custom vmax
    if vmax < MAX:
        bars += '▒'
        vmax += 1

    b = (
        ((a - vmin) / (vmax - vmin) * (len(bars) - 1))
        .astype(int)
        .clip(0, len(bars) - 1)
        .tolist()
    )

    graph = ''.join(bars[i] for i in b)
    xticks = (f'╹{" "*(ticks-1)}' * math.ceil((L + start) / ticks))[:L]

    # left-aligned tick labels
    # labels = ''.join(f'{i: <{ticks}}' for i in range(start, L+start, ticks))

    # right-aligned tick labels
    labels = ''.join(f'{i: >{ticks}}' for i in range(start + ticks, L + start, ticks))
    S = str(start)
    labels = S + labels[len(S) - 1 :]

    out = [graph, xticks, labels]

    if wrap:
        import textwrap
        from itertools import chain

        wrapper = textwrap.TextWrapper(width=wrap)
        out = list(chain.from_iterable(zip(*(wrapper.wrap(s) for s in out))))

    return '\n'.join(out)


def seq2aa(seq):
    return ''.join([to_aa.get(seq[i : i + 3], '???') for i in range(0, len(seq), 3)])


def export_data(filename, seqs=None, stats=None, extras=None, outdir=None, MAX=None):

    if '/' in filename:
        path, fname = filename.rsplit('/', 1)
    else:
        fname = filename
    if outdir:
        path = outdir.rstrip('/') + '/'
    if not path:
        path = '.'

    import re

    base_fname = re.sub(r'(\.assembled)?\.fastq$', '', fname)
    f_log = base_fname + '.processed.log'
    f_json = base_fname + '.processed.json'
    f_extra = base_fname + '.processed.extra'
    f_seq = base_fname + '.processed.txt'
    f_seq_aa = base_fname + '_aa.processed.txt'

    if stats:
        with open(f'{path}/{f_json}', 'w') as f_json:
            import json

            f_json.write(json.dumps(stats))

        with open(f'{path}/{f_log}', 'w') as f_log:
            f_log.write('')
            N = len(str(stats['total_sequences'])) + 2
            f_log.write(
                f"Total sequences read: {stats['total_sequences']: >{N}}\n\n"
                f"Rejected sequences\n"
                f"    no adaptor found: {stats['noadaptor']: >{N}}\n"
                f"         contaminant: {stats['contaminant']: >{N}}\n"
                f"      stretches of A: {stats['A_stretch']: >{N}}\n"
                f"      seq contains N: {stats['seq_N']: >{N}}\n"
                f"   subseq contains N: {stats['subseq_N']: >{N}}\n"
                f"         low quality: {stats['lowqual']: >{N}}\n"
                f"    subseq too short: {stats['tooshort']: >{N}}\n"
                f"     subseq too long: {stats['toolong']: >{N}}\n\n"
            )
            if stats['lengths'].values():
                f_log.write(
                    'Distribution of subsequences length:\n'
                    f"scale: {min(stats['lengths'].values())}–{max(stats['lengths'].values())}\n"
                )
            else:
                f_log.write('Distribution of subsequences length:\nno sequences\n')
            f_log.write(simple_graph(stats['lengths'], start=1) + '\n')
    if extras:
        np.array(extras, dtype=np.int8).tofile(f'{path}/{f_extra}')
    if seqs:
        with open(f'{path}/{f_seq}', 'w') as f, open(f'{path}/{f_seq_aa}', 'w') as f_aa:
            if not MAX:
                MAX = stats['MAX_LEN'] if stats else max(map(len, seqs))
            for seq in seqs:
                s = f'{seq: >{MAX}}\n'
                f.write(s)
                f_aa.write(seq2aa(s))


def export_all(results_all, outdir=None):
    for filename, (seqs, stats) in results_all.items():
        print(k)
        print(simple_graph(stats['lengths'], start=0, wrap=81))
        export_data(filename, seqs=seqs, stats=stats, outdir=outdir, MAX=None)


def main(out=False):
    import argparse
    import datetime

    def ArgRange(string):
        from argparse import ArgumentTypeError
        import re

        m = re.match(r'^(?:(\d+)?-)?(\d+)?$', string)
        if not m:
            raise ArgumentTypeError(
                f'invalid format "{string}", use 1-10 or 1- or -10 or 10'
            )
        return tuple(int(x) if x else 0 for x in m.groups())

    parser = argparse.ArgumentParser(
        description='parse fastq files and extract peptides sequence for inverse toe-printing',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument('files', help='input files', nargs='+')

    parser.add_argument(
        '-o',
        '--outdir',
        help='directory for the output files',
        dest='outdir',
        default='.',
    )

    parser.add_argument(
        '-a1',
        '--left-adaptator',
        help='Sequence (5ʹ→3ʹ) of the left adaptator',
        dest='a1',
        default=DEFAULTS['a1'],
    )
    parser.add_argument(
        '-a2',
        '--right-adaptator',
        help='Sequence (5ʹ→3ʹ) of the right adaptator',
        dest='a2',
        default=DEFAULTS['a2'],
    )

    parser.add_argument(
        '-s',
        '--peptide-size',
        help='Allowed length range of peptide considered (min-max)',
        dest='range',
        default='3-14',
        type=ArgRange,
    )
    parser.add_argument(
        '-q',
        '--quality',
        help='Threshold for the PHRED quality cutoff',
        dest='quality',
        default=30,
    )

    parser.add_argument(
        '--limit',
        help='Max number of sequences to process (useful for quick tests)',
        dest='limit',
        type=int,
        default=None,
    )

    parser.add_argument('--trim_start', action='store_true')
    parser.add_argument('--no-trim_start', dest='trim_start', action='store_false')
    parser.set_defaults(trim_start=True)

    args = parser.parse_args()

    args_dic = vars(args)
    print(args_dic)
    kwargs = args_dic.copy()

    kwargs['min_seq_len'] = args_dic['range'][0] * 3
    kwargs['max_seq_len'] = (args_dic['range'][1] + 1) * 3

    del kwargs['files']
    del kwargs['outdir']
    print(datetime.datetime.now().strftime('Start time: %Y-%m-%d %H:%M'))
    results_all = parse_all(
        args_dic['files'], outdir=args_dic['outdir'], save=True, **kwargs
    )
    print(datetime.datetime.now().strftime('  End time: %Y-%m-%d %H:%M'))
    if out:
        return results_all


if __name__ == '__main__':
    main()
