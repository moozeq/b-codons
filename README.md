# Description
Simple program which converts nucleotide sequence to amino acid sequences and does some length frequency analysis. 

# Installation

```bash
git clone https://github.com/moozeq/B_Codons.git

cd B_Codons
pip3 install -r requirements.txt
```

# Usage

## Help
```bash
./app.py -h
```

## In pipeline
```bash
cat seq.fasta | ./app.py
```

## Input sequence
```bash
./app.py <sequence>
```