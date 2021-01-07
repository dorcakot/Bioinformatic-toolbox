from Bio import AlignIO, Align, SeqRecord


class MSA:
    """
    Class MSA is a representation of multiple sequence alignment data stored in clustal format.
    """
    def __init__(self, file: str) -> None:
        """
        Initialize class instance with name of a file containing data in clustal format.
        Function reads this file and store information accessible as 'align' property.
        :param file: name of a file containing data
        """
        self.align = AlignIO.read(file, "clustal")

    def get_seq_by_position(self, position: int) -> Align.SeqRecord:
        """
        To obtain a sequence by its position, use this function.
        :param position: position of a sequence in file
        :return: sequence
        """
        return self.align[position]

    def get_seq_by_id(self, id: str) -> SeqRecord.SeqRecord:
        """
        Use this function to retrieve sequence by its id.
        :param id: sequence id
        :return: sequence
        """
        for seq in self.align:
            if seq.id == id:
                print(type(seq))
                return seq

    def get_column(self, num: int) -> str:
        """
        Use to get column of multiple sequence alignment.
        :param num: number of a column
        :return: column from multiple sequence alignment
        """
        return self.align[:, num]

    def get_column_score(self, column: int, matrix: Align.substitution_matrices
                        = Align.substitution_matrices.load('PAM250'), gap: int = 1) -> int:
        """
        Function calculates score of a column with respect to given scoring matrix and gap penalty value.
        :param column: number of a column
        :param matrix: scoring matrix
        :param gap: score of a gap (passed as positive number, when evaluating, subtracted from score)
        :return: score of a column
        """
        score = 0
        rows = len(self.align[:, 0])
        for i in range(rows):
            for j in range(i, rows):
                if self.align[i, column] == '-' and self.align[j, column] == '-':
                    score += gap
                elif self.align[i, column] == '-' or self.align[j, column] == '-':
                    score -= gap
                else:
                    score += matrix[self.align[i, column], self.align[j, column]]
        return score

    def get_score(self, matrix: Align.substitution_matrices = Align.substitution_matrices.load('BLOSUM62'), gap: int = 1
                  ) -> int:
        """
        Function calculates score of the whole alignment with respect to given scoring matrix and gap penalty value.
        :param matrix: scoring matrix
        :param gap:  score of a gap (passed as positive number, when evaluating, subtracted from score)
        :return: score of the whole multiple sequence alignment
        """
        score = 0
        columns = len(self.align[0])
        for i in range(columns):
            score += self.get_column_score(i, matrix, gap)
        return score

    def get_conservation_scores(self, sequence: int, matrix: Align.substitution_matrices = Align.substitution_matrices
                                .load('BLOSUM62')) -> int:
        """
        Function calculates score of each position within a sequence with respect to given scoring matrix.
        :param sequence: number of a sequence within the alignment
        :param matrix: scoring matrix
        :return: score of each position in a sequence
        """
        columns = len(self.align[0])
        scores = [(float, int)]
        for i in range(columns):
            score = 0
            for j in range(len(self.align[:, i])):
                if self.align[j, i] != '-' and self.align[sequence, i] != '-':
                    score += matrix[self.align[j, i], self.align[sequence, i]]
            scores.append((score , i))
        return scores

    def get_top_scoring_positions(self, N: int, matrix: Align.substitution_matrices = Align.substitution_matrices.load
                                ('BLOSUM62'), gap: int = 1) -> list:
        """
        Function first calculates scores of all the positions with respect to a given scoring matrix. Then it returns
        the first top N scoring ones.
        :param N: number of positions to be returned
        :param matrix: scoring matrix
        :param gap: score of a gap (passed as positive number, when evaluating, subtracted from score)
        :return: top N scoring positions with values of their score
        """
        columns = len(self.align[0])
        scores = []
        for i in range(columns):
            scores.append((self.get_column_score(i, matrix, gap), i))
        scores.sort(reverse=True)
        return scores[:N]
