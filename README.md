# Positional-Burrows-Wheeler-Transform-PBWT-Haplotype Matching

This repository provides **Python and C++ implementations** of the  
**Positional Burrows–Wheeler Transform (PBWT)** algorithm for detecting
long shared haplotype segments.

The implementation follows the original algorithm introduced by  
**Richard Durbin** and is intended for:

- Learning and understanding PBWT
- Algorithmic bioinformatics demonstrations
- Academic projects and portfolios
- Testing on small to medium haplotype datasets

---

## Clone the Repository

Clone this repository to your local machine:

```bash
git clone https://github.com/SherazAhmadd/PBWT-a-fast-Haplotypes-Matching-Algorithm.git
cd PBWT-a-fast-Haplotypes-Matching-Algorithm
```

## Repository Structure

```text
PBWT.py                  # Python implementation of PBWT
PBWT.cpp                 # C++ implementation of PBWT
test_data/
 └── example_haplotypes.csv # Example haplotype dataset
README.md                # Project documentation
```

## Input File Format

Input data must be provided as a CSV file.
Format Rules

- The first column is a haplotype identifier (ignored by the algorithm)
- The remaining columns are haplotype values (0 or 1)
- All rows must have the same number of sites
```
id,s1,s2,s3,s4,s5
h0,0,1,1,0,1
h1,0,1,1,0,1
h2,1,0,0,1,0
h3,1,0,0,1,0
```
## Python Installation
Requirements
The Python version requires:
- Python 3.8+
- NumPy
- Pandas

```
pip install numpy pandas
```
- Usage: PBWT.py <haplotypes.csv> <min_match_length>
```
python PBWT.py test_data/example_haplotypes.csv 25
```
## C++ Implementation
- Usage: ./pbwt <haplotypes.csv> <min_match_length>

Compile
```
g++ -O3 PBWT.cpp -o pbwt
./pbwt test_data/example_haplotypes.csv 25
```

## Output Format
Output file: matched_haplotypes.csv
```
hap1,hap2,start,end,length
3,7,15,42,27
```
| Column | Description                              |
| ------ | ---------------------------------------- |
| hap1   | Index of first haplotype                 |
| hap2   | Index of second haplotype                |
| start  | Start site of shared segment             |
| end    | End site of shared segment               |
| length | Length of shared segment (`end - start`) |

## PBWT Variable Glossary

| Symbol  | Description                           |
| ------- | ------------------------------------- |
| N       | Number of haplotypes (rows)           |
| M       | Number of sites (columns)             |
| X       | Haplotype matrix                      |
| a[k]    | PBWT ordering at site `k`             |
| d[k][i] | Divergence value for haplotype `i`    |
| k       | Site (column) index                   |
| i       | Index within PBWT ordering            |
| p, q    | Divergence trackers for allele groups |


## Contributor / Author / Creator

**Name:** Rana Sheraz Ahmad  
**Role:** Implementator  
**GitHub:** [SherazAhmadd](https://github.com/SherazAhmadd)


## Contact & Issues
Email: ranasheraz.202101902@gcuf.edu.pk 

Issues: https://github.com/SherazAhmadd/PBWT-a-fast-Haplotypes-Matching-Algorithm/issues

## Acknowledgment
Based on the PBWT algorithm by:  
**Richard Durbin**  
*[Efficient haplotype matching and storage using the positional Burrows–Wheeler transform (Bioinformatics, 2014)](https://doi.org/10.1093/bioinformatics/btu014)*

