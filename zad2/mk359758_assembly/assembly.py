import sys

################################ Wziete z cwiczen ##############################

class DeBruijnGraph:
    """ A de Bruijn multigraph built from a collection of strings.
        User supplies strings and k-mer length k.  Nodes of the de
        Bruijn graph are k-1-mers and edges correspond to the k-mer
        that joins a left k-1-mer to a right k-1-mer. """

    @staticmethod
    def chop(st, k):
        """ Chop a string up into k mers of given length """
        for i in xrange(0, len(st)-(k-1)):
            yield (st[i:i+k], st[i:i+k-1], st[i+1:i+k])
    
    class Node:
        """ Node in a de Bruijn graph, representing a k-1 mer.  We keep
            track of # of incoming/outgoing edges so it's easy to check
            for balanced, semi-balanced. """
        
        def __init__(self, km1mer):
            self.km1mer = km1mer
        
        def __hash__(self):
            return hash(self.km1mer)
        
        def __str__(self):
            return self.km1mer
    
    def __init__(self, strIter, k):
        """ Build de Bruijn multigraph given string iterator and k-mer
            length k """
        self.G = {}     # multimap from nodes to neighbors
        self.nodes = {} # maps k-1-mers to Node objects
        for st in strIter:
            for kmer, km1L, km1R in self.chop(st, k):
                nodeL, nodeR = None, None
                if km1L in self.nodes:
                    nodeL = self.nodes[km1L]
                else:
                    nodeL = self.nodes[km1L] = self.Node(km1L)
                if km1R in self.nodes:
                    nodeR = self.nodes[km1R]
                else:
                    nodeR = self.nodes[km1R] = self.Node(km1R)
                self.G.setdefault(nodeL, []).append(nodeR)

def neighbors1mm(kmer, alpha):
    ''' Generate all neighbors at Hamming distance 1 from kmer '''
    neighbors = []
    for j in xrange(len(kmer)-1, -1, -1):
        oldc = kmer[j]
        for c in alpha:
            if c == oldc: continue
            neighbors.append(kmer[:j] + c + kmer[j+1:])
    return neighbors

def kmerHist(reads, k):
    ''' Return k-mer histogram and average # k-mer occurrences '''
    kmerhist = {}
    for read in reads:
        for kmer in [ read[i:i+k] for i in xrange(0, len(read)-(k-1)) ]:
            kmerhist[kmer] = kmerhist.get(kmer, 0) + 1
    return kmerhist

def correct1mm(read, k, kmerhist, alpha, thresh):
    ''' Return an error-corrected version of read.  k = k-mer length.
        kmerhist is kmer count map.  alpha is alphabet.  thresh is
        count threshold above which k-mer is considered correct. '''
    # Iterate over k-mers in read
    for i in xrange(0, len(read)-(k-1)):
        kmer = read[i:i+k]
        # If k-mer is infrequent...
        if kmerhist.get(kmer, 0) <= thresh:
            # Look for a frequent neighbor
            for newkmer in neighbors1mm(kmer, alpha):
                if kmerhist.get(newkmer, 0) > thresh:
                    # Found a frequent neighbor; replace old kmer
                    # with neighbor
                    read = read[:i] + newkmer + read[i+k:]
                    break
    # Return possibly-corrected read
    return read

################################################################################

# Tworzymy zachlannie graf zbudowany tylko z najlepszych przejsc
def buildGreedyGraph(deBruijn):
    graph = {}
    for src, dsts in deBruijn.G.iteritems():
        srclab = src.km1mer
        weightmap = {}
        for dst in dsts:
            weightmap[dst] = weightmap.get(dst, 0) + 1
        best_dst = ""
        best_weight = 0
        for dst, v in weightmap.iteritems():
            dstlab = dst.km1mer
            if best_weight < v:
                best_dst = dstlab
                best_weight = v
        graph[srclab] = best_dst
    return graph

# Przchodzimy graf uproszczonym BFS-em, poniewaz z jednego wierzcholka
# wychodzi maksymalnie jedna krawedz. Pojedyncza sciezka tworzy contig.
def traverseForContigs(graph):
    contigs = []
    visited = set()
    for src in graph:
        if src not in visited:
            v = src
            contig = v[:-1]
            while (v and v not in visited):
                visited.add(v)
                contig = contig + v[-1:]
                v = graph.get(v, None)
            contigs.append(contig)
    return contigs

def correctReads(reads, k, threshold):
    hist = kmerHist(reads, K)
    corrected_reads = []
    for i in range(len(reads)):
        corrected_reads.append(correct1mm(reads[i], k, hist, "ACTG", threshold))
    return corrected_reads

def importReads(filename):
    reads = []
    with open(filename, "r") as f:
        is_header = True
        for read in f:
            if is_header:
                is_header = False
            else:
                reads = reads + [read.rstrip()]
                is_header = True
    return reads

def exportContigs(filename, contigs):
  with open(filename, "w") as f:
    for i, contig in enumerate(contigs):
      f.write('>' + str(i) + '\n')
      f.write(contig + '\n')

if __name__ == '__main__':
    filename_in = sys.argv[1]
    filename_out = sys.argv[2]
    reads = importReads(filename_in)
    K = 17
    ERR_THRESHOLD = 1
    corr_reads = correctReads(reads, K, ERR_THRESHOLD)
    graph = DeBruijnGraph(corr_reads, K)
    greedy_graph = buildGreedyGraph(graph)
    exportContigs(filename_out, traverseForContigs(greedy_graph))
