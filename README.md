# Motif Finding Algorithms

This repository contains implementations of three classic motif-finding algorithms in bioinformatics:

- **Greedy Motif Search**
- **Randomized Motif Search**
- **Randomized Motif Search with Gibbs Sampling**

These algorithms aim to identify shared motifs (subsequences) in a given set of DNA sequences.

---

## 🧬 What is a Motif?

In bioinformatics, a **motif** is a short, recurring sequence pattern that is presumed to have a biological function. Motif-finding is crucial in identifying regulatory elements such as promoters or enhancers in DNA.

---

## 📁 Files in This Repository

- `greedy_motif_search_with_pseudocounts.py`: Implements the **Greedy Motif Search** algorithm.
- `randomized_motif_search.py`: Implements the **Randomized Motif Search** algorithm.
- `randomized_motif_search_with_gibbs_sampler.py`: Implements the **Gibbs Sampling Motif Search** algorithm.
- `dna_list_example.txt`: An example input file with DNA sequences.

---

## ⚙️ How to Use

### 1. Prepare Input

Ensure your DNA sequences are listed line-by-line in a plain text file (e.g., `dna_list_example.txt`). Each line should contain a reference ID (e.g. seq1) followed by a single DNA string. For example:

  seq1
  ACATCGATCATGCTGACTGA

  seq2
  ACAGCTTTTACGGAGCGTTA

### 2. Run the Algorithms

Each script reads from the input file and prints:

- The **consensus motif**
- The **score** (lower score = better consensus)
- The list of **best motifs** found

---

#### Greedy Motif Search

  bash
```python greedy_motif_search_with_pseudocounts.py```

#### Randomized Motif Search

  bash
```python randomized_motif_search.py```

#### Gibbs Sampler Motif Search 

  bash
```python randomized_motif_search_with_gibbs_sampler.py```

You can change the default parameters (k, t, and N) directly in each script's function call at the bottom.

---

## 🧠 Algorithm Overviews

### Greedy Motif Search

- Iteratively builds a motif set by choosing the most probable k-mer using a profile.
- Fast and deterministic, but may get stuck in local optima.

### Randomized Motif Search

- Starts with random k-mers, updates motifs using profiles until convergence.
- Repeated multiple times to improve chances of a better result.

### Gibbs Sampling

- Iteratively refines motifs by randomly sampling one at a time based on profiles.
- Often yields better motifs with enough iterations.

---

## 👤 Author

Heitor Gelain do Nascimento
Email: heitorgelain@outlook.com
GitHub: @heitor-sg5

---

## 🧪 Example Output

Consensus: GAAAAAAATTTTTTT

Score: 2

Best Motifs: 
'CAAAAAAATTTTTTT', 
'GAAAAAAATTTTTTT', 
'GAAAAAAATTTTTTT', 
'CAAAAAAATTTTTTT', 
'GAAAAAAATTTTTTT', 
'GAAAAAAATTTTTTT'
