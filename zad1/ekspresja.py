import numpy as np
import csv
import sys

def read_data(filename):
    data = np.genfromtxt(filename, skip_header=1, dtype='string')
    expr = {}
    GENE_COL = 6
    FPKM_COL = 9
    for row in data:
        chrom, rangee = row[GENE_COL].split(':')
        if chrom not in expr:
            expr[chrom] = {}
        (start, end) = (int(rangee.split('-')[0]), int(rangee.split('-')[1]))
        expr[chrom][(start, end)] = float(row[FPKM_COL])
    return expr

def bucket_genes(expr):
    result = {}
    for chrom, e_data in expr.iteritems():
        result[chrom] = {
            25: set(),
            50: set(),
            75: set(),
            100: set(),
            1000000000: set()
        }

        for gene, fpkm in e_data.iteritems():
            for b_range in sorted(result[chrom]):
                if (fpkm < b_range):
                    result[chrom][b_range].add(gene)
                    break
    return result

def export_buckets_to_csv(data, filename):
    with open(filename, "w") as output:
        writer = csv.writer(output, lineterminator='\n')
        writer.writerow(['chromosome', '<25', '25-49', '50-74', '75-100', '>100'])
        for chrom, buckets in data.iteritems():
            writer.writerow([chrom] + [len(buckets[k]) for k in sorted(buckets)])

def export_correlation_to_csv(data, filename):
    with open(filename, "w") as output:
        writer = csv.writer(output, lineterminator='\n')
        writer.writerow(['chromosome', '<25', '25-49', '50-74', '75-100', '>100'])
        for chrom, buckets in data.iteritems():
            writer.writerow([chrom] + [buckets[k] for k in sorted(buckets)])

def read_MACS(filename):
    MACS_peaks = {}
    file = open(filename, "r")
    iterfile = iter(file)
    next(iterfile)
    for line in iterfile:
        spl = line.split('\t')
        shortLine = "\t".join([spl[0], spl[1], spl[2]])
        if (spl[0] not in MACS_peaks):
            MACS_peaks[spl[0]] = [(spl[1], spl[2])]
        else:
            MACS_peaks[spl[0]].append((spl[1], spl[2]))

    return MACS_peaks

def intersect(l1, l2):
    l1_i = 0
    l2_i = 0
    res = []
    while (l1_i < len(l1) and l2_i < len(l2)):
        l1_el = (int(l1[l1_i][0]), int(l1[l1_i][1]))
        l2_el = (int(l2[l2_i][0]), int(l2[l2_i][1]))
        if (l1_el[0] > l2_el[1]):
            l2_i += 1
        elif (l2_el[0] > l1_el[1]):
            l1_i += 1
        elif (l1_el[0] >= l2_el[0] and l1_el[1] <= l2_el[1]):
            res.append(l1_el)
            l1_i += 1
        elif (l1_el[0] <= l2_el[0] and l1_el[1] >= l2_el[1]):
            res.append(l2_el)
            l2_i += 1
        elif (l1_el[0] >= l2_el[0] and l1_el[1] >= l2_el[1]):
            res.append((l1_el[0], l2_el[1]))
            l2_i += 1
        elif (l2_el[0] >= l1_el[0] and l2_el[1] >= l1_el[1]):
            res.append((l2_el[0], l1_el[1]))
            l1_i += 1
    return res

def find_correlation(genes, peaks):
    result = {}
    for chrom, buckets in genes.iteritems():
        if chrom in peaks:
            result[chrom] = {}
            for b_range, genes in buckets.iteritems():
                result[chrom][b_range] = len(intersect(sorted(list(genes)), peaks[chrom])) / float(len(list(genes)))
    return result

if __name__ == '__main__':
    filename_expr = sys.argv[1]
    filename_macs = sys.argv[2]
    expr = read_data(filename_expr)
    genes_in_buckets = bucket_genes(expr)
    export_buckets_to_csv(genes_in_buckets, 'C_genes_expr_buckets.csv')

    macs_peaks = read_MACS(filename_macs)
    corr = find_correlation(genes_in_buckets, macs_peaks)
    export_correlation_to_csv(corr, "C_correlation.csv")