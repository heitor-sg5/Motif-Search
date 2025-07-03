import random

def gibbs_sampler(Dna, k, t, N):
    motifs = [random_kmer(dna, k) for dna in Dna]
    best_motifs = motifs[:]

    for _ in range(N):
        i = random.randint(0, t - 1)
        motifs_except_i = motifs[:i] + motifs[i+1:]
        profile = build_profile(motifs_except_i)
        
        motifs[i] = profile_random_kmer(Dna[i], k, profile)
        
        if score(motifs) < score(best_motifs):
            best_motifs = motifs[:]

    print("Consensus:", consensus(best_motifs))
    print("Score:", score(best_motifs))

    return best_motifs

def random_kmer(dna, k):
    start = random.randint(0, len(dna) - k)
    return dna[start:start + k]

def build_profile(motifs):
    k = len(motifs[0])
    t = len(motifs)
    profile = {nuc: [1] * k for nuc in "ACGT"}

    for motif in motifs:
        for i, nucleotide in enumerate(motif):
            profile[nucleotide][i] += 1

    for i in range(k):
        total = sum(profile[nuc][i] for nuc in "ACGT")
        for nuc in "ACGT":
            profile[nuc][i] /= total

    return profile

def profile_random_kmer(text, k, profile):
    n = len(text)
    probabilities = []

    for i in range(n - k + 1):
        kmer = text[i:i+k]
        prob = 1
        for j, nucleotide in enumerate(kmer):
            prob *= profile[nucleotide][j]
        probabilities.append(prob)

    total_prob = sum(probabilities)
    if total_prob == 0:
        return text[random.randint(0, n - k): random.randint(0, n - k) + k]
    probabilities = [p / total_prob for p in probabilities]
    chosen_index = random.choices(range(n - k + 1), weights=probabilities)[0]

    return text[chosen_index:chosen_index + k]

def consensus(motifs):
    k = len(motifs[0])
    consensus = ""
    for i in range(k):
        counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
        for motif in motifs:
            counts[motif[i]] += 1
        consensus += max(counts, key=counts.get)
    return consensus

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

def read_sequences_from_file():
    sequences = []
    with open('dna_list_example.txt', 'r') as file:
        lines = file.readlines()
        for line in lines:
            line = line.strip()
            if line and not line.lower().startswith("seq"):
                sequences.append(line)
    return sequences

best_motifs = gibbs_sampler(read_sequences_from_file(), k=15, t=6, N=1000)
print("Best Motifs:", best_motifs)
