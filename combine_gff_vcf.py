import os
import gzip
import argparse


parsers = {
    'MDO': {
        'gff_chr': lambda ll: ll[0].split(':')[-1],
        'gene_id': lambda ll: list(filter(lambda d: d[0:7] == 'Parent=',
                                          ll[8].split(';')))[0][7:].split(',')[0].split(':')[-1]
        },
    'ATH': {
        'gff_chr': lambda ll: ll[0],
        'gene_id': lambda ll: list(filter(lambda d: d[0:7] == 'Parent=',
                                          ll[8].split(';')))[0][7:].split(',')[0]
        }
    }

def get_seq(my_seq_file):
    seq = ''
    chrm = ''
    tmp = []
    with open(my_seq_file) as f:
        # fasta format. Header has to be the first line
        chrm = next(f).strip().split()[0]
        # return sequence for all the other lines.
        # quicker if the whole the seq is in a single line.
        for l in f:
            seq += l.rstrip()
    f.close()
    # remove >
    chrm = chrm[1:]
    return (chrm, seq)


# Return a dict using cds id as key and the gff3 features of cds as list of dict.
# Field names for gff3 features are took from sequenceontology.org
# def parse_gff(my_gff_file, my_chr):
#     my_gff = open(my_gff_file)
#     my_dict = {}
#     for l in my_gff:
#         # tmp={}
#         l = l.rstrip()
#         l = l.split("\t")
#         if l[2] != 'CDS': continue
#         if l[0] == my_chr:
#             desc = l[8].split(';')
#             parent = ''
#             for d in desc:
#                 if d[0:7] == 'Parent=':
#                     parent = d[7:]
#             tmp = {'start': int(l[3]),
#                    'end': int(l[4]),
#                    'strand': l[6],
#                    'phase': int(l[7]),
#                    'cds_id': parent}
#             if parent in my_dict:
#                 my_dict[parent].append(tmp)
#             else:
#                 my_dict[parent] = [tmp]
#     my_gff.close()
#     return (my_dict)

def parse_gff_new(my_gff_file, my_chr):
    position_map = {}
    cds = {}

    if args.organism:
        gffchrparser = parsers[args.organism]['gff_chr']
        geneidparser = parsers[args.organism]['gene_id']
    else:
        gffchrparser = lambda ll: ll[0]
        geneidparser = lambda ll: list(filter(lambda d: d[0:7] == 'Parent=',
                                              ll[8].split(';')))[0][7:]


    with open(my_gff_file) as my_gff:
        for l in my_gff:

            l = l.rstrip()
            l = l.split("\t")
            if l[2] != 'CDS': continue
            # chr = l[0].split(':')[-1]
            chr = gffchrparser(l)
            if chr == my_chr:
                desc = l[8].split(';')
                parent = ''
                parent = geneidparser(l)
                # for d in desc:
                #     if d[0:7] == 'Parent=':
                #         parent = d[7:].split(',')[0]

                tmp = {
                    'start': int(l[3]),
                    'end': int(l[4]),
                    'strand': l[6],
                    'phase': int(l[7]),
                    'cds_id': parent}

                for i in xrange(int(l[3]), int(l[4]) + 1):
                    position_map[i] = tmp

                cds.setdefault(parent, list()).append(tmp)

    return position_map, cds


def parse_vcf_header(header_file):
    with open(header_file) as f:
        return next(f).strip('#').strip().split()


# While reading the vcf variants it keeps only those ones that overlaps with a cds
def parse_vcf(my_vcf_file, my_cds, vcf_header=None):
    my_keys = parse_vcf_header(vcf_header) if vcf_header else list()
    my_sample_name = []

    for i in my_keys[::-1]:
        if i == 'FORMAT':
            break
        my_sample_name.append(i)

    openfunc = gzip.open if 'gz' in os.path.splitext(my_vcf_file)[1] else open

    with openfunc(my_vcf_file) as my_vcf:
        for l in my_vcf:
            l = l.rstrip()
            if l[0] == '#':
                if l[1: 6] == 'CHROM':
                    my_keys = l.split("\t")
                    my_keys[0] = my_keys[0][1:]
                continue
            l = l.split("\t")
            tmp = dict(zip(my_keys, l))
            tmp['POS'] = int(tmp['POS'])

            mut = combine_cds_vcf_new(position_map, tmp)
            if mut:

                # whole_mut.append(mut)
                yield mut, my_sample_name

    # get sample names
    # In VCF, samples are the last columns of header line, after the "FORMAT" column


# Check the overlap btw the variant and the cds
def combine_cds_vcf(my_cds, my_vcf):
    my_dict = {}
    #
    rv = False
    for _, m in my_cds.iteritems():
        for i in m:
            # print(type(my_vcf['POS']),type(i['start']),type(my_vcf['POS']),type(i['end']))
            if my_vcf['POS'] >= i['start'] and my_vcf['POS'] <= i['end']:
                my_dict.update(i)
                my_dict.update(my_vcf)
                rv = True
                break
            if rv:
                break
    return (my_dict)


def combine_cds_vcf_new(position_map, my_vcf):
    my_dict = None
    parent = position_map.get(my_vcf['POS'])
    if parent:
        my_dict = my_vcf
        my_dict.update(parent)

    return my_dict


# print results in a big table. positions are protein based
def print_table(my_seq):
    first_iteration = True
    header = ['CDS_ID', 'ALT']

    for m, my_sample_name in parse_vcf(args.vcf, cds, vcf_header=args.vcfHeader):
        if first_iteration is True:
            header += my_sample_name
            print("\t".join(header))
            first_iteration = False
        # temporary
        # Currently only SNP variant are take into account
        if len(m['REF']) == 1 and len(m['ALT']) == 1:
            (codon, alt_codon) = get_snp_mutation(m, my_seq)

            # Mutation format is encoded as original aa + protein pos + mutated aa
            # Eg. F426L
            # means Phe at postion 426 is mutated in Leu
            variant = get_aa(codon, genetic_table)
            variant += str(get_alt_pos_on_cds(m, cds[m['cds_id']]))
            variant += get_aa(alt_codon, genetic_table)
            to_print = [m['cds_id'], variant]
            to_print += get_cultivar(m, my_sample_name)
            print("\t".join(to_print))


# Get codon and alternative codon from postion and sequence
def get_snp_mutation(my_mut, my_seq):
    # python string coords are zero based
    my_pos = my_mut['POS'] - 1
    # use % func to get the nucleotide position of the codon take into account exon phase and strand
    # eg. strand + ATG->210, strand - ATG->012
    #
    if my_mut['strand'] == '+':
        tmp = my_mut['POS'] - (my_mut['start'] + my_mut['phase']) + 1
    else:
        if my_mut['phase'] > 0:
            tmp = (my_mut['end'] + my_mut['phase']) - my_mut['POS']
        else:
            tmp = (my_mut['end'] + my_mut['phase']) - my_mut['POS'] + 1
    rest = tmp % 3
    # once that the codon pos is known get the codon
    my_codon = get_codon(my_pos, rest, my_seq, my_mut['strand'])
    # create the seq with the snp and get the mutated codon
    my_alt_seq = my_seq[:my_pos] + my_mut['ALT'] + my_seq[my_pos + 1:]
    my_alt_codon = get_codon(my_pos, rest, my_alt_seq, my_mut['strand'])
    return (my_codon, my_alt_codon)


# Adjust the phase oc cds
def get_codon(my_pos, my_rest, my_seq, my_strand):
    my_codon = 'NNN'
    #
    # for sure there is a shorter and maybe easier way to do that ... but it works
    if my_strand == '-':
        if my_rest == 0:
            st = my_pos
            ed = my_pos + 3
        elif my_rest == 2:
            st = my_pos - 1
            ed = my_pos + 2
        else:
            st = my_pos - 2
            ed = my_pos + 1
    else:
        if my_rest == 0:
            st = my_pos - 2
            ed = my_pos + 1
        elif my_rest == 2:
            st = my_pos - 1
            ed = my_pos + 2
        else:
            st = my_pos
            ed = my_pos + 3
    #
    my_codon = my_seq[st:ed]
    if my_strand == '-':
        my_codon = reverse_complement(my_codon)
    return (my_codon)


# Get the reverse complement
def reverse_complement(my_str):
    rc = ''
    my_str = my_str.upper()
    for i in my_str[::-1]:
        if i == 'A':
            rc += 'T'
        elif i == 'T':
            rc += 'A'
        elif i == 'G':
            rc += 'C'
        elif i == 'C':
            rc += 'G'
        else:  # It loose eventual iupac codes
            rc += 'N'
    return (rc)


# Built a useful dictionary
def get_genetic_table():
    bases = ['T', 'C', 'A', 'G']
    amino_acid = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
    codons = [a + b + c for a in bases for b in bases for c in bases]
    codon_table = dict(zip(codons, amino_acid))
    return (codon_table)


# Return the aminoacid encoded as sinlge letter or "Stop" for \*
def get_aa(my_codon, my_codon_table):
    my_aa = 'X'
    if my_codon in my_codon_table:
        my_aa = my_codon_table[my_codon]
    if my_aa == '*':
        my_aa = 'Stop'
    return (my_aa)


# Since cds is divided into exons; it gets the right exons
def get_alt_pos_on_cds(my_mut, my_cds):
    start = []
    my_cds_pos = 0
    # get the exons
    # for i in my_cds :
    #    exons[i['start']]=i['end']
    exons = {i['start']: i['end'] for i in my_cds}
    # order the exons by start and strand
    if my_mut['strand'] == '-':
        start = reversed(sorted(exons.keys()))
    else:
        start = sorted(exons.keys())
    # get the length from stat to the snp
    cds_len = 0
    for s in start:
        # print my_mut['cds_id'],s,exons[s],cds_len
        if my_mut['POS'] > s and my_mut['POS'] < exons[s]:
            if my_mut['strand'] == '-':
                cds_len += exons[s] - my_mut['POS']
            else:
                cds_len += my_mut['POS'] - s + 1
            break
        else:
            cds_len += exons[s] - s + 1
    my_cds_pos = cds_len // 3
    if cds_len % 3 > 0:
        my_cds_pos += 1
    return (my_cds_pos)


def get_cultivar(my_alt, my_sample_name):
    alleles = []
    for s in my_sample_name:
        # split the sample line. and get different code for the predicted ploidity:
        # _N_o variation for 0/0
        # _H_eterozygous var for 0/1
        # homoz_Y_gotic var 1/1
        # if my_alt['POS'] == 3767:
        #     print(s)
        #     break
        tmp = my_alt[s].split(':')
        if tmp[0] in {'1|1', '1/1'}:
            alleles.append('Y')
        elif tmp[0] in {'1|0', '0|1', '1/0', '0/1'}:
            alleles.append('H')
        else:
            alleles.append('N')
    return (alleles)


def parse_args():
    parser = argparse.ArgumentParser(prog='combine_vcf_gff')
    parser.add_argument('vcf', help='vcf file')
    parser.add_argument('gff', help='gff file')
    parser.add_argument('seq', help='chromosome sequence file')
    parser.add_argument('-vh', '--vcfHeader', help='header of the vcf file in case the latter '
                                                   'doesn\'t have one')
    parser.add_argument('-org', '--organism', default=None)
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()

    # Store name and sequence of reference sequence.
    (chrm, seq) = get_seq(args.seq)
    # Store cds positions in the reference sequence
    position_map, cds = parse_gff_new(args.gff, chrm)
    # Predicted variants in a vcf format overlapping with cds
    # (mut_cds, vcf_sample_name) = parse_vcf(args.vcf, position_map, vcf_header=args.vcfHeader)
    # To produce a suitable genetic table
    genetic_table = get_genetic_table()
    # print result
    # print_table(vcf_sample_name, mut_cds, seq)
    print_table(seq)
