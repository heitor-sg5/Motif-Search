import random

def random_motif_search(Dna, k, t):
    motifs = [random_kmer(dna, k) for dna in Dna]
    best_motifs = motifs[:]

    while True:
        profile = build_profile(motifs)
        motifs = [most_probable_kmer(dna, k, profile) for dna in Dna]
        
        if score(motifs) < score(best_motifs):
            best_motifs = motifs[:]
        else:
            return best_motifs

def random_kmer(dna, k):
    start = random.randint(0, len(dna) - k)
    return dna[start:start + k]

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

def consensus(motifs):
    k = len(motifs[0])
    consensus = ""
    for i in range(k):
        counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
        for motif in motifs:
            counts[motif[i]] += 1
        consensus += max(counts, key=counts.get)
    return consensus

def repeat_search(Dna, k, t, n):
    best_motifs = None
    best_score = float('inf')

    for _ in range(n):
        motifs = random_motif_search(Dna, k, t)
        current_score = score(motifs)
        if current_score < best_score:
            best_score = current_score
            best_motifs = motifs[:]

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

best_motifs = repeat_search(read_sequences_from_file(), k=15, t=6, n=100)
print("Best Motifs:", best_motifs)
