#!/usr/bin/env python3
import argparse
import math
import sys
import matplotlib.pyplot as plt


class Sequence:
    nucleotides = ['A', 'C', 'T', 'G']
    genetics_code = {
        'AAA': 'K', 'AAC': 'N', 'AAT': 'N', 'AAG': 'K',
        'ACA': 'T', 'ACC': 'T', 'ACT': 'T', 'ACG': 'T',
        'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
        'AGA': 'R', 'AGC': 'S', 'AGT': 'S', 'AGG': 'R',
        'CAA': 'Q', 'CAC': 'H', 'CAT': 'H', 'CAG': 'Q',
        'CCA': 'P', 'CCC': 'P', 'CCT': 'P', 'CCG': 'P',
        'CTA': 'L', 'CTC': 'L', 'CTT': 'L', 'CTG': 'L',
        'CGA': 'R', 'CGC': 'R', 'CGT': 'R', 'CGG': 'R',
        'TAA': '*', 'TAC': 'Y', 'TAT': 'Y', 'TAG': '*',
        'TCA': 'S', 'TCC': 'S', 'TCT': 'S', 'TCG': 'S',
        'TTA': 'L', 'TTC': 'F', 'TTT': 'F', 'TTG': 'L',
        'TGA': '*', 'TGC': 'C', 'TGT': 'C', 'TGG': 'W',
        'GAA': 'E', 'GAC': 'D', 'GAT': 'D', 'GAG': 'E',
        'GCA': 'A', 'GCC': 'A', 'GCT': 'A', 'GCG': 'A',
        'GTA': 'V', 'GTC': 'V', 'GTT': 'V', 'GTG': 'V',
        'GGA': 'G', 'GGC': 'G', 'GGT': 'G', 'GGG': 'G'
    }
    amino_code = {

    }

    def __init__(self, nsequence: str, index: int = 0):
        self.index = index
        if nsequence.startswith('>'):  # fasta
            lines = nsequence.splitlines()
            nsequence = ''.join(lines[1:])
        self.nsequence = nsequence
        self.psequence = ''.join([codon for codon in self])
        self.proteins = [protein for protein in self.psequence.split('*') if protein]
        self.proteins_amino_num = list(map(len, self.proteins))
        self.max_protein_length = max(self.proteins_amino_num)

    def __iter__(self):
        while self.index + 3 <= len(self.nsequence) and len(self.nsequence) > 0:
            codon = self.nsequence[self.index:self.index + 3]
            yield self.genetics_code[codon]
            self.index = self.index + 3
        self.index = 0


def main():
    parser = argparse.ArgumentParser(description='Reading proteins from nucleotide sequence')
    parser.add_argument('sequence', nargs='?', default='-', type=str, help='nucleotide sequence')
    parser.add_argument('-s', '--sort', action='store_true', default=False, help='sort proteins sequences by length')
    args = parser.parse_args()
    if args.sequence == '-':
        args.sequence = ''.join([line for line in sys.stdin])
    seq = Sequence(args.sequence)
    print('\n\n'.join(seq.proteins if not args.sort else sorted(seq.proteins, key=len)))
    print('\n')
    print(f'Proteins found = {len(seq.proteins)}')
    nums = []
    for i in range(1, seq.max_protein_length + 1):
        proteins_num_by_length = len(list(filter(lambda x: x == i, seq.proteins_amino_num)))
        if proteins_num_by_length > 0:
            print(f'\t length {i} = {proteins_num_by_length}')
        nums.append(proteins_num_by_length)
    xaxis = list(range(1, seq.max_protein_length + 1))
    plt.plot(xaxis, nums)
    plt.grid(True)
    plt.show()


if __name__ == "__main__":
    main()
