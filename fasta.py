from Bio import SeqIO

class Fasta:
        """
        Class representing a fasta file parser.
        An instance of this class is initialized by the name of a fasta file.
        An array of molecules contained in this file is created and store as 'molecules' property.
        """
        def __init__(self, file_name):
                """
                Initialize Fasta parser using a name of a fasta file. During initialization, fasta file is read
                and molecules are stored in 'molecules' property.
                :param file_name: name of a file containing data in fasta format
                """
                self.molecules = []
                for record in SeqIO.parse(file_name, "fasta"):
                        self.molecules.append(record)

        def get_molecules(self):
                """
                Get list of all the molecules in this file.
                :return: list of molecules
                """
                molecules = []
                for i in range(len(self.molecules)):
                        molecules.append(('[', i, '] ', self.molecules[i].id))
                return molecules

        def get_description(self, num):
                """
                Function serves to access information about individual molecules.
                :param num: number of a molecule
                :return: description of given molecule
                """
                return(self.molecules[num].description)

        def get_sequence(self, num):
                """
                Function retrieves sequence of a molecule.
                :param num: number of a molecule
                :return: sequence of given molecule
                """
                return(self.molecules[num].seq)

        def get_sequence_length(self, num):
                """
                Function provides information about length of a sequence of a given molecule.
                :param num: number of a molecule
                :return: length of a sequence of given molecule
                """
                return(len(self.molecules[num].seq))

        def get_subsequence(self, num, start, end):
                """
                Function serves to obtain a subsequence of a given molecule.
                :param num: number of a molecule
                :param start: start index of subsequence
                :param end: end index of subsequence
                :return: subsequence of a given molecule
                """
                return(self.molecules[num].seq[start:end])
