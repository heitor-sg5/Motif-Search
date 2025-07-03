
def build_profile(motifs):
    k = len(motifs[0])
    profile = {nuc: [1] * k for nuc in "ACGT"}

    for motif in motifs:
        for i, nucleotide in enumerate(motif):
            profile[nucleotide][i] += 1

    for i in range(k):
        total = sum(profile[nuc][i] for nuc in "ACGT")
        for nuc in "ACGT":
            profile[nuc][i] /= total

    return profile

def score(motifs):
    k = len(motifs[0])
    t = len(motifs)
    score = 0

    for i in range(k):
        counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
        for motif in motifs:
            counts[motif[i]] += 1
        max_count = max(counts.values())
        score += (t - max_count) 

    return score

def most_probable_kmer(text, k, profile):
    max_prob = -1
    best_kmer = text[0:k]

    for i in range(len(text) - k + 1):
        kmer = text[i:i + k]
        prob = 1
        for j, nucleotide in enumerate(kmer):
            prob *= profile[nucleotide][j]
        if prob > max_prob:
            max_prob = prob
            best_kmer = kmer

    return best_kmer

def consensus(motifs):
    k = len(motifs[0])
    consensus = ""
    for i in range(k):
        counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
        for motif in motifs:
            counts[motif[i]] += 1
        consensus += max(counts, key=counts.get)
    return consensus

def greedy_motif_search(Dna, k, t):
    best_motifs = [dna[:k] for dna in Dna] 
    first_string = Dna[0]

    for i in range(len(first_string) - k + 1):
        motifs = [first_string[i:i + k]]
        for j in range(1, t):
            profile = build_profile(motifs)
            next_motif = most_probable_kmer(Dna[j], k, profile)
            motifs.append(next_motif)

        if score(motifs) < score(best_motifs):
            best_motifs = motifs
    
    print("Consensus:", consensus(best_motifs))
    print("Score:", score(best_motifs))

    return best_motifs

def read_sequences_from_file():
    sequences = []
    with open('dna_list_example.txt', 'r') as file:
        lines = file.readlines()
        for line in lines:
            line = line.strip()
            if line and not line.lower().startswith("seq"):
                sequences.append(line)
    return sequences

best_motifs = greedy_motif_search(read_sequences_from_file(), k=15, t=6)
print("Best Motifs:", best_motifs)
