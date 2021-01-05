import argparse
import sys
from .msa import *
from .pdb import *
from .fasta import *
from .edit import *
from .ham import *

def main():
    """
    Build commandline argument parser and act upon given arguments.
    """
    arguments = argparse.ArgumentParser()
    arguments_sub = arguments.add_subparsers()

    arguments_ed = arguments_sub.add_parser('ed', help='Count edit distance and align of two sequences.')
    arguments_ed.set_defaults(action='ed')
    arguments_ed.add_argument('seq1')
    arguments_ed.add_argument('seq2')
    arguments_ed.add_argument('--align', action='store_true')

    arguments_fasta = arguments_sub.add_parser('fasta', help='Parse fasta file.')
    arguments_fasta.set_defaults(action='fasta')
    arguments_fasta.add_argument('file')
    arguments_fasta.add_argument('--description', dest='description', type=int)
    arguments_fasta.add_argument('--sequence', dest='sequence', type=int)
    arguments_fasta.add_argument('--length', dest='length', type=int)
    arguments_fasta.add_argument('--sub', dest='sub', type=int, nargs=3)

    arguments_hd = arguments_sub.add_parser('hd', help='Count hamming distance of two sequences of identical length.')
    arguments_hd.set_defaults(action='ham')
    arguments_hd.add_argument('seq1')
    arguments_hd.add_argument('seq2')

    arguments_msa = arguments_sub.add_parser('msa', help='Parse msa file in clustal format.')
    arguments_msa.set_defaults(action='msa')
    arguments_msa.add_argument('file')
    arguments_msa.add_argument('--position', dest='position', type=int)
    arguments_msa.add_argument('--id', dest='id')
    arguments_msa.add_argument('--column', dest='column', type=int)
    arguments_msa.add_argument('--column_score', dest='cs', type=int)
    arguments_msa.add_argument('--score', action='store_true')
    arguments_msa.add_argument('--top', dest='top', type=int)
    arguments_msa.add_argument('--conservation', dest='cons', type=int)

    arguments_pdb = arguments_sub.add_parser('pdb', help='Parse pdb files.')
    arguments_pdb.set_defaults(action='pdb')
    arguments_pdb.add_argument('file')
    arguments_pdb.add_argument('--model', action='store_true')
    arguments_pdb.add_argument('--chain', dest='chain')
    arguments_pdb.add_argument('--residue', dest='residue')
    arguments_pdb.add_argument('--atom', dest='atom')
    arguments_pdb.add_argument('--width', action='store_true')
    arguments_pdb.add_argument('--distance', dest='dist', type=int)
    arguments_pdb.add_argument('-r', action='store_true')
    arguments_pdb.add_argument('-a', action='store_true')
    arguments_pdb.add_argument('--show_histogram', action='store_true')
    arguments_pdb.add_argument('--store_histogram', dest='his_file')
    arguments_pdb.add_argument('-e', action='store_true')
    arguments_pdb.add_argument('--polar', action='store_true')
    arguments_pdb.add_argument('--exposed', action='store_true')


    if len(sys.argv)<2:
        class HelpConfig:
            def __init__(self):
                self.action = 'help'
        config = HelpConfig()
    else:
        config = arguments.parse_args()

    if config.action == 'help':
        arguments.print_usage()
        arguments.print_help()
    elif config.action == 'ed':
        print(edit_distance(config.seq1, config.seq2))
        if config.align:
            for align in all_alignments(config.seq1, config.seq2):
                print(align+"\n")
    elif config.action == 'fasta':
        file = Fasta(config.file)
        if config.description is not None:
            print(file.get_description(config.description))
        elif config.sequence is not None:
            print(file.get_sequence(config.sequence))
        elif config.length is not None:
            print(file.get_sequence_length(config.length))
        elif config.sub is not None:
            print(file.get_subsequence(config.sub[0], config.sub[1], config.sub[2]))
        else:
            print(file.get_molecules())
    elif config.action == 'ham':
        print(hd(config.seq1, config.seq2))
    elif config.action == 'msa':
        file = MSA(config.file)
        if config.position is not None:
            print(file.get_seq_by_position(config.position))
        elif config.id is not None:
            print(file.get_seq_by_id(config.id))
        elif config.column is not None:
            print(file.get_column(config.column))
        elif config.cs is not None:
            print(file.get_column_score(config.cs))
        elif config.score:
            print(file.get_score())
        elif config.top is not None:
            for t in (file.get_top_scoring_positions(config.top)):
                print('Position '+str(t[1])+' has score '+str(t[0]))
        elif config.cons is not None:
            print(file.get_conservation_scores(file.get_seq_by_position(config.cons)))
    elif config.action == 'pdb':
        file = PDBMolecule(config.file)
        if config.model:
            for chain in file.structure.get_chains():
                print(chain.id)
        elif config.chain is not None:
            ch = file.get_chain(config.chain)
            for residue in ch:
                print(residue.get_resname(), residue.full_id)
        elif config.dist:
            if config.a:
                for atom in file.get_atoms_at_distance(file.get_atom(config.atom, file.get_residue(config.residue)),config.dist):
                    print(atom.full_id)
            elif config.r:
                for residue in file.get_residues_at_distance(file.get_atom(config.atom, file.get_residue(config.residue)),config.dist):
                    print(residue.full_id)
        elif config.atom is not None:
            res = file.get_residue(config.residue)
            atom = file.get_atom(config.atom, res)
            print(atom.full_id)
        elif config.residue is not None:
            res = file.get_residue(config.residue)
            for atom in res:
                print(atom.full_id)
        elif config.width:
            print(file.get_width())
        elif config.show_histogram:
            if config.e:
                file.show_histogram(exposed=True)
            else:
                file.show_histogram()
        elif config.his_file is not None:
            if config.e:
                file.save_histogram(config.his_file, exposed=True)
            else:
                file.save_histogram(config.his_file)
        elif config.polar:
            print(file.polarity_ratio())
        elif config.exposed:
            print(file.ratio_exposed())
        else:
            print(file.get_information())





if __name__ == '__main__':
    main()
