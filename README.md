# Bioinformatic toolbox

Bioinformatic toolbox contains basic functions, which help to handle basic bioinformatic tasks.
User can retrieve some information using commandline or importing functions in own python program and work further with obtained objects and information.

##### Steps to get bioinformatic toolbox running:
```console
$ git clone git@gitlab.mff.cuni.cz:dorcakot/bioinformatic-toolbox.git
$ cd bioinformatic-toolbox
$ python setup.py install
``` 

##### Measuring sequence similarity using Hamming distance
Function **hd(seq1, seq2)** prints Hamming distance of two given sequences of identical length. In case of different length, error message is returned.

*Commandline access:* 
```console
$ binf hd seq1 seq2
```
Example use: 
```console
$ binf hd ATCG ATGC
2
```
```console
$ binf hd ATCG ATG 
Exception: Sequences have different length!
```
*Python access:*
```
>>> from binf.ham import * 
>>> hd(seq1, seq2)
```
Example use: 
```
>>> from binf.fasta import *
>>> file = Fasta("inputs/ls_orchid.fasta") 
>>> a = file.get_subsequence(0, 0, 5) 
>>> b = file.get_subsequence(0, 0, 5)
>>> print(hd(a, b))
0
```
```
>>> file = fasta.Fasta("inputs/ls_orchid.fasta") 
>>> a = file.get_subsequence(0, 0, 5) 
>>> b = file.get_subsequence(0, 5, 10) 
>>> print(hd(a, b))
4
```
```
>>> file = fasta.Fasta("inputs/ls_orchid.fasta") 
>>> a = file.get_subsequence(0, 0, 5) 
>>> b = file.get_subsequence(0, 0, 10) 
>>> print(hd(a, b))
Exception: Sequences have different length!
```
##### Sequence alignment using edit distance
Function **edit_distance(seq1, seq2)** returns value of edit distance of the two sequences. To get list of all alignments use function **all_alignments(seq1, seq2)**.

*Commandline access:*
```console
$ binf ed seq1 seq2 [--align]
```
Example use: 
```console
$ binf ed abeceda abecede 
1
```
```console
$ binf ed writers vintner --align 
5 
wri-t-ers 
-vintner-

wri-t-ers 
v-intner-

writ-ers 
vintner-
```
*Python access:*
```
>>> from binf.edit import *
>>> edit_distance(seq1, seq2) 
>>> all_alignments(seq1, seq2)
```
Example use: 
```
>>> print(edit_distance('stall', 'table')) 
3
```
```
>>> for alignment in all_alignments('stall', 'table'): 
>>>     print(alignment + "\n") 
stall- 
-table

sta-ll 
-table
```
##### Processing FASTA files
Parser of FASTA files is implemented as a class **Fasta**. An instance is initialized by file name of data stored in fasta format. An array of molecules contained in this file is created. Such instance contains an array of molecules from this file accesible via property molecules or function get_molecules. To get information (description, sequence, sequence length) about individual molecules, use the number of the molecule and corresponding function (**get_description**, **get_sequence**...).

*Commandline access:*
```console
$ binf fasta file_name [--description molecule_num] [--sequence molecule_num] [--length molecule_num] [--sub molecule_num start end]
```
Example use:
```console
$ binf fasta inputs/ls_orchid.fasta --description 0 
gi|2765658|emb|Z78533.1|CIZ78533 C.irapeanum 5.8S rRNA gene and ITS1 and ITS2 DNA
```
```console
$ binf fasta inputs/1tup.fasta --length 1
21
```
*Python access:*
```
>>> from binf.fasta import *
>>> molecules = Fasta(file_name)
>>> molecules.get_molecules()
>>> molecules.get_description(seq_num)
>>> molecules.get_sequence(seq_num)
>>> molecules.get_sequence_length(seq_num)
>>> molecules.get_subsequence(seq_num, start, end)
```
Example use:
```
>>> molecules = Fasta('inputs/1tup.fasta') 
>>> print(molecules.get_subsequence(1, 0, 15))
ATAATTGGGCAAGTC
```
```
>>> molecules = Fasta('inputs/ls_orchid.fasta')
>>> print(molecules.get_sequence(10)
CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGCAGAATATATGATCGAGTGAATCCGGTGGACTTGTGGTTACTCAGCTCGACATAGGCTTTGCTTTTGCGGTGACCCTAATTTGTCATTGGGCCTCCCCCCAAGCTTTCCTTGTGGGTTTGAACCTCTAGCACGGTGCAGTATGCGCCAAGTCATATGAAGCATCACTGATGAATGACATTATTGTCCAAAAAGTTGGAGTGGAAGCGTGCTATTGCATACATGCAAATGAATTTTTTATGACTCTCGACATATCGTGGTGTGATCGCAGAATCCCGTGAACCATCGAGTCTTTGAACGCAAGTTGCGCCCGATGCCATCAGGCTAAGGGAACGCCTGCCTGGGCGTCGTGTGCTGCGTCTCTCCTGTCAATGGTTATACGTCATATAGACAGGTTTGCATTGCGTGGATGTGAAAGATTGGCCCCTTGTGCCTAGGTGCGGTGGGTCTAAGGACTAGTGTTTTGATGGTTCGAAACCTGGCAGGAGGTGGAGGATGTTGGCAGCTATAAGGCTATCATTTGAATCCCCCATATTGTCGTGTTTGTCGGACCTAGAGAAGAACATGTTTGAATCCCAATGGGGGCAAACAACCCTCGGGCGGTTGATTGCCATTCATATGCGACCCCAGGTCAGGCGGCCACCCGCTGAG
```
##### Processing multiple sequence alignment
Parser of multiple sequence alignment data stored in clustal format implemented as a class **MSA**. To create an instance use name of file. Functionalities of the parser include obtaining sequence by position/id, column by its position number, score of such column or score of the whole msa with respect to given scoring matrix passed as argument (default matrix is Blosum62, default gap penalty is 1 (passed as positive number and subtracted from the score)). To use in other ways, user is able to access alignment via class property align.

*Commandline access:*
```console
$ binf msa file_name [--position num] [--id id] [--column column_num] [--column_score column_num] [--score]
```
Example use: 
```console
$ binf msa inputs/p53.clustal --position 17 
ID: UniRef90_UPI000 Name: Description: UniRef90_UPI000 Number of features: 0 Seq('MLGNSPARVECEPGGSGQGSRGIRVLPEARIGETVSFSCSCGNGVLPEATMTDP...DSD')
```
```console
$ binf msa inputs/p53.clustal --score
2050433.0
```
```console
$ binf msa inputs/p53.clustal --column_score 100 
1333.0
```
*Python access:* 
```
>>> from binf.msa import * 
>>> data = MSA(file_name)
>>> data.get_seq_by_position(position_num)
>>> data.get_seq_by_id(position_id)
>>> data.get_column(col_num)
>>> data.get_column_score(col_num, matrix=Align.substitution_matrices.load('PAM250'), gap=1) 
>>> data.get_score(matrix=Align.substitution_matrices.load('PAM250'), gap=1)
```
Example use:
```
>>> data = MSA('inputs/p53.clustal')
>>> print(data.get_column(300))
TNNNATNNNNVINSNNNNHNNNNNNNNNTNNNNTTNQQHQNHNNQQQQNIVNQKT
```
```
>>> data = MSA('inputs/p53.clustal')
>>> print(data.get_seq_by_id('UniRef90_H3B1Z4'))
ID: UniRef90_H3B1Z4 Name: Description: UniRef90_H3B1Z4 Number of features: 0 Seq('--------------------------------------------------MTDP...GGG')
```

##### Conservation determination from multiple aligned sequences
Extension of msa library. Specify sequence in msa and retrieve its conservation score or display N top scoring positions of this alignment.

*Commandline access:*
```console
$ binf msa file [--top N] [--conservation seq_num]
```
Example use: 
```console
$ binf msa inputs/p53.clustal --top 9 
Position 74 has score 25192.0 
Position 370 has score 18480.0 
Position 368 has score 18480.0 
Position 332 has score 18480.0 
Position 328 has score 18480.0 
Position 259 has score 18480.0
Position 218 has score 18480.0 
Position 224 has score 17830.0 
Position 310 has score 15400.0
```
*Python access:*
```
>>> from binf.msa import * 
>>> data = MSA(file_name) 
>>> data.get_conservation_scores(seq_num, matrix=Align.substitution_matrices.load('BLOSUM62'), gap=1)
>>> data.get_top_scoring_position(N, matrix=Align.substitution_matrices.load('BLOSUM62'), gap=1)
```
Example use:
```
>>> data = MSA('inputs/p53.clustal')
>>> print(data.get_conservation_scores(17)
```

##### Processing PDB files
PDB files parser class. Initialize with pdb file name and easily retrieve structural information. User is able to obtain an object respresenting a model, chain, residue, atom; display information about the macromolecule or compute its width. Another useful functionality is to obtain list of all atoms(residues) at given distance from a heteroatom.

*Commandline access:*
```console
$ binf pdb file_name [--model] [--chain id] [--residue num] [--atom num][--width] [--distance distance] [-r] [-a]
```
Example use: 
```console
$ binf pdb inputs/1TUP.pdb
number of models: 1
number of chains: 5
number of residues: 1014
number of atoms: 5828
```
```console
$ binf pdb inputs/1TUP.pdb --width
87.83
```
```
$ binf pdb inputs/1TUP.pdb --residue 951 --atom ZN --distance 4 -a
('1tup', 0, 'A', (' ', 176, ' '), ('N', ' ')) ('1tup', 0, 'A', (' ', 239, ' '), ('N', ' '))
```
*Python access:*
```
>>> from binf.pdb import *
>>> data = PDBMolecule(file_name)
>>> data.get_chain(chain_id)
>>> data.get_residue(res_id)
>>> data.get_atom(atom_id, residue)
>>> data.get_information()
>>> data.get_width()
>>> data.get_atoms_at_distance(hetatm, distance, round_digits=1)
>>> data.get_residues_at_distance(hetatm, distance, round_digits=1)
```
Example use:
```
>>> data = PDBMolecule('inputs/1TUP.pdb')
>>> chain_A = data.get_chain('A')
>>> residue_951 = data.get_residue(951)
>>> zn_hetatm = data.get_atom('ZN', residue_951)
>>> for residue in data.get_residues_at_distance(zn_hetatm, 3, round_digits=0): >>>     print(residue)
<Residue CYS het= resseq=176 icode= >
<Residue HIS het= resseq=179 icode= > 
<Residue CYS het= resseq=242 icode= >
```

##### Computing structure-related properties
Extension of PDB parser. These functions compute the ratio of surface and buried amino acids; quantify portion of polar amino acids in the core and on the surface of the protein. Furthermore, user is able to display/save a histogram showing amino acid composition of buried and exposed amino acids.

*Commandline access:*
```console
$ binf pdb file_name [--show_histogram] [--store_histogram HIS_FILE] [-e] [--polar] [--exposed]
```
Example use: 
```console
$ binf pdb inputs/1b0b.pdb --polar
28% of buried residues is polar
72% of buried residues is not polar
49% of exposed residues is polar
51% of exposed residues is not polar
```
*Python access:*
```
>>> from binf.pdb import *
>>> data = PDBMolecule(file_name)
>>> data.show_histogram(exposed=False)
>>> data.save_histogram(file='histogram', exposed='False')
>>> data.ratio_exposed()
>>> data.polarity_ratio()
>>> data.ratio_exposed()
```
Example use:
```
>>> data = PDBMolecule('inputs/a2a.pdb')
>>> print(data.ratio_exposed())
56% of residues is buried
44% of residues is exposed
```

###### A2a vs 1b0b
A2a

56% of residues is buried
44% of residues is exposed

29% of buried residues is polar
71% of buried residues is not polar
49% of exposed residues is polar
51% of exposed residues is not polar

1b0b

45% of residues is buried
55% of residues is exposed


28% of buried residues is polar
72% of buried residues is not polar
49% of exposed residues is polar
51% of exposed residues is not polar

A2a molecule has more than half of the residues buried, while in case of 1b0b molecule, more than half of the residues are exposed. Considering their structure, this ration difference makes sense. Molecule 1b0b is smaller and more flattened - molecule has bigger proportion of surface. A2a receptor is "thicker", it contains more buried molecules in the core.
When comparing the portion of polar and nonpolar amino acids of these molecules, there is a significant trend of having more nonpolar molecules "hidden" in the core from the polar environment.
