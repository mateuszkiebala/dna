import sys

#l1 = [(1,2), (4,7), (11,11), (12,13)]
#l2 = [(3,4), (6,8), (9,11)]

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

# od l1 odejmij l2
def subtract(l1, l2):
    l1_i = 0
    l2_i = 0
    res = []
    while (l1_i < len(l1) and l2_i < len(l2)):
        #print l1
        #print l1_i
        #print l2
        #print l2_i
        l1_el = (int(l1[l1_i][0]), int(l1[l1_i][1]))
        l2_el = (int(l2[l2_i][0]), int(l2[l2_i][1]))
        if (l1_el[0] > l2_el[1]):
            l2_i += 1
        elif (l2_el[0] > l1_el[1]):
            res.append(l1_el)
            l1_i += 1
        elif (l1_el[0] >= l2_el[0] and l1_el[1] <= l2_el[1]):
            l1_i += 1
        elif (l1_el[0] < l2_el[0]):
            res.append((l1_el[0], l2_el[0]-1))
            if (l1_el[1] > l2_el[1]):
                l1[l1_i] = (l2_el[1]+1, l1_el[1])
                l2_i += 1
            else:
                l1_i += 1
        elif (l2_el[0] <= l1_el[0]):
            if (l1_el[1] > l2_el[1]):
                l1[l1_i] = (l2_el[1]+1, l1_el[1])
                l2_i += 1
            else:
                l1_i += 1
        #print res
        #print '\n'
    while (l1_i < len(l1)):
        l1_el = (int(l1[l1_i][0]), int(l1[l1_i][1]))
        res.append(l1_el)
        l1_i += 1
    return res

#print subtract(l2, l1)

def intervals_len(l):
    return len(l)
    # odkomentować w przypadku analizy obszarowej
    #res = 0
    #for p in l:
    #    res += (int(p[1]) - int(p[0]))
    #return res

transcripts = {}
exons = {}
genes = {}
file = open("genes.gtf", "r")
for line in file:
    spl = line.split('\t')
    shortLine = "\t".join([spl[0], spl[2], spl[3], spl[4]])
    if (spl[2] == 'transcript'):
        if (spl[0] not in transcripts):
            transcripts[spl[0]] = [(spl[3], spl[4])]
        else:
            transcripts[spl[0]].append((spl[3], spl[4]))
    elif (spl[2] == 'exon'):
        if (spl[0] not in exons):
            exons[spl[0]] = [(spl[3], spl[4])]
        else:
            exons[spl[0]].append((spl[3], spl[4]))
    if (spl[0] not in genes):
        genes[spl[0]] = [(spl[3], spl[4])]
    else:
        genes[spl[0]].append((spl[3], spl[4]))
        
introns = {}
for chromosome, transcript_list in transcripts.items():
    introns[chromosome] = subtract(transcript_list, exons[chromosome])
    #print val
    #print exons[key]
#print transcripts['AB325691']
#print exons['AB325691']
#print introns['AB325691']
        
filename = sys.argv[1]
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

# wyliczenie procentowego udzialu obszaru zajmowanego przez peaki w obszarze zajmowanym przez exony, introny i obszary miedzygenowe
# dla kazdego z chromosomow (I, II, III) oraz laczne (total)
exons_percentage  = {}
introns_percentage = {}
intergene_percentage = {}
total_peaks_len = 0
total_exons_len = 0
total_introns_len = 0
for chromosome, peaks in MACS_peaks.items():
    peaks_len = intervals_len(peaks)
    exons_len = intervals_len(intersect(exons[chromosome], peaks))
    exons_percentage[chromosome] = float((exons_len) * 100)/float(peaks_len)
    introns_len = intervals_len(intersect(introns[chromosome], peaks))
    introns_percentage[chromosome] = float((introns_len) * 100)/float(peaks_len)
    intergene_percentage[chromosome] = 100 - exons_percentage[chromosome] - introns_percentage[chromosome]
    total_peaks_len += peaks_len
    total_exons_len += exons_len
    total_introns_len += introns_len
exons_percentage['total'] = float((total_exons_len) * 100)/float(total_peaks_len)
introns_percentage['total'] = float((total_introns_len) * 100)/float(total_peaks_len)
intergene_percentage['total'] = 100 - exons_percentage['total'] - introns_percentage['total']

print "exons coverage (in %):"
print exons_percentage
print "introns coverage (in %):"
print introns_percentage
print "intergene coverage (in %):"
print intergene_percentage


# koordynaty centromerow i telomerow ze strony www.pombase.org
centromeres = {'I': [(3753687, 3789421)],'II': [(1602264, 1644747)],'III': [(1070904, 1137003)]}
telomeres = {'I': [(1, 29663), (5554844, 5579133)],'II': [(1, 39186), (4500619, 4539800)], 'III': []}

# wyliczenie procentowego udzialu obszaru zajmowanego przez peaki w obszarze zajmowanym przez centromery i telomery
# dla kazdego z chromosomow (I, II, III) oraz laczne (total)
centromeres_percentage  = {}
telomeres_percentage = {}
total_peaks_len = 0
total_centromeres_len = 0
total_telomeres_len = 0
for chromosome, peaks in MACS_peaks.items():
    peaks_len = intervals_len(peaks)
    centromeres_len = intervals_len(intersect(centromeres[chromosome], peaks))
    centromeres_percentage[chromosome] = float((centromeres_len) * 100)/float(peaks_len)
    telomeres_len = intervals_len(intersect(telomeres[chromosome], peaks))
    telomeres_percentage[chromosome] = float((telomeres_len) * 100)/float(peaks_len)
    total_peaks_len += peaks_len
    total_centromeres_len += centromeres_len
    total_telomeres_len += telomeres_len
centromeres_percentage['total'] = float((total_centromeres_len) * 100)/float(total_peaks_len)
telomeres_percentage['total'] = float((total_telomeres_len) * 100)/float(total_peaks_len)

print "centromeres coverage (in %):"
print centromeres_percentage
print "telomeres coverage (in %):"
print telomeres_percentage
