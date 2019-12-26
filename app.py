#!/usr/bin/env python3
import argparse
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
        'A': {'one': 'A', 'tri': 'Ala', 'full': 'Alanine '},
        'C': {'one': 'C', 'tri': 'Cys', 'full': 'Cysteine '},
        'D': {'one': 'D', 'tri': 'Asp', 'full': 'Aspartic acid '},
        'E': {'one': 'E', 'tri': 'Glu', 'full': 'Glutamic acid '},
        'F': {'one': 'F', 'tri': 'Phe', 'full': 'Phenylalanine '},
        'G': {'one': 'G', 'tri': 'Gly', 'full': 'Glycine '},
        'H': {'one': 'H', 'tri': 'His', 'full': 'Histidine '},
        'I': {'one': 'I', 'tri': 'Ile', 'full': 'Isoleucine '},
        'K': {'one': 'K', 'tri': 'Lys', 'full': 'Lysine '},
        'L': {'one': 'L', 'tri': 'Leu', 'full': 'Leucine '},
        'M': {'one': 'M', 'tri': 'Met', 'full': 'Methionine '},
        'N': {'one': 'N', 'tri': 'Asn', 'full': 'Asparagine '},
        'P': {'one': 'P', 'tri': 'Pro', 'full': 'Proline '},
        'Q': {'one': 'Q', 'tri': 'Gln', 'full': 'Glutamine '},
        'R': {'one': 'R', 'tri': 'Arg', 'full': 'Arginine '},
        'S': {'one': 'S', 'tri': 'Ser', 'full': 'Serine '},
        'T': {'one': 'T', 'tri': 'Thr', 'full': 'Threonine '},
        'V': {'one': 'V', 'tri': 'Val', 'full': 'Valine '},
        'W': {'one': 'W', 'tri': 'Trp', 'full': 'Tryptophan '},
        'Y': {'one': 'Y', 'tri': 'Tyr', 'full': 'Tyrosine '}
    }
    amino_modes = ['one', 'tri', 'full']

    def convert_codon(self, protein, mode='one'):
        return ''.join([self.amino_code[codon][mode] for codon in protein if mode in self.amino_code[codon]])

    def convert_protein(self, proteins, mode='one'):
        return [self.convert_codon(protein, mode) for protein in proteins]

    def __init__(self, nsequence: str, index: int = 0):
        self.index = index
        if nsequence.startswith('>'):  # fasta
            lines = nsequence.splitlines()
            nsequence = ''.join(lines[1:])
        self.nsequence = nsequence
        self.psequence = ''.join([codon for codon in self])
        proteins = [protein for protein in self.psequence.split('*') if protein]
        self.proteins = {key: self.convert_protein(proteins, key) for key in self.amino_modes}
        self.proteins_amino_num = list(map(len, self.proteins['one']))
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
    parser.add_argument('-l', '--length', action='store_true', default=False, help='show sequences length analysis')
    parser.add_argument('-p', '--plot', action='store_true', default=False, help='show proteins sequences length plot')
    parser.add_argument('-a', '--amino', type=str, default='one', choices=['one', 'tri', 'full'],
                        help='amino acids representation')
    args = parser.parse_args()
    if args.sequence == '-':
        args.sequence = ''.join([line for line in sys.stdin])
    seq = Sequence(args.sequence)
    print('\n\n'.join(seq.proteins[args.amino] if not args.sort else sorted(seq.proteins[args.amino], key=len)))
    print('\n')
    print(f'Proteins found = {len(seq.proteins["one"])}')
    if args.length:
        nums = []
        for i in range(1, seq.max_protein_length + 1):
            proteins_num_by_length = len(list(filter(lambda x: x == i, seq.proteins_amino_num)))
            if proteins_num_by_length > 0:
                print(f'\t length {i} = {proteins_num_by_length}')
            nums.append(proteins_num_by_length)
        if args.plot:
            xaxis = list(range(1, seq.max_protein_length + 1))
            plt.plot(xaxis, nums)
            plt.grid(True)
            plt.show()


if __name__ == "__main__":
    main()
