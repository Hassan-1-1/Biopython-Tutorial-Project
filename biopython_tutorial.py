import os
from Bio import AlignIO
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Align import PairwiseAligner
from Bio.Align import substitution_matrices
from Bio.Align.Applications import ClustalOmegaCommandline
from Bio.Align.Applications import ClustalwCommandline
from Bio.Align.Applications import MuscleCommandline

# -------------------functions for complement and reverse complement------------------

# --------------complement---------------

# when the sequence is uploaded directly
def seq_complement(sequence):
    sequence = sequence.upper()
    string = 'BDEFHIJKLMNOPQRSVWXYZ'
    for letter in string:
        if letter in sequence:
            print('Please input a valid DNA/RNA sequence')
            break
    else:
        sequence = Seq(sequence)
        print('Sequence:', sequence)
        print('Complement:', sequence.complement())
    return()

# When a sequnce file is uploaded
def file_complement(file_name, file_formate):
    for seq_record in SeqIO.parse(file_name, file_formate):
        print(seq_record.id)
        print("Sequence:\n", seq_record.seq)
        print("Complement:\n", seq_record.seq.complement())
        print('Seq length: ', len(seq_record))
        print(' ')
    return()

# -------------reverse complement----------------

# when the sequence is uploaded directly
def seq_reverse_complement(sequence):
    sequence = sequence.upper()
    string = 'BDEFHIJKLMNOPQRSVWXYZ'
    for letter in string:
        if letter in sequence:
            print('Please input a valid DNA/RNA sequence')
            break
    else:
        sequence = Seq(sequence)
        print('Sequence: ', sequence)
        print('Reverse Complement: ', sequence.reverse_complement())
    return()

# When a sequnce file is uploaded
def file_reverse_complement(file_name, file_formate):
    for seq_record in SeqIO.parse(file_name, file_formate):
        print(seq_record.id)
        print("Sequence:\n", seq_record.seq)
        print("Reverse Complement:\n", seq_record.seq.reverse_complement())
        print('Seq length: ', len(seq_record))
        print(' ')
    return()


# ----------------------functions for pairwise sequence alignment--------------------

# When the sequence is uploaded directly
def seq_pairwise_alignment(seq1, seq2, mode, match = 1, missmatch = 0, gap_open = 0, gap_extend = 0):
    aligner = PairwiseAligner()
    aligner.mode = mode # specifying mode (global, local)

    # defining alignment score values
    aligner.match_score = match
    aligner.mismatch_score = missmatch
    aligner.query_internal_open_gap_score = gap_open
    aligner.query_internal_extend_gap_score = gap_extend

    # converting seq in uppercase
    seq1 = seq1.upper()
    seq2 = seq2.upper()

    #Performing Alignment
    string = 'BDEFHIJKLMNOPQRSVWXYZ'
    for letter in string:
        if letter not in seq1 and seq2: # for DNA/RNA sequences
            alignments = aligner.align(seq1, seq2)
        else: # for protein alignment
            aligner.match_score = 1
            aligner.mismatch_score = 0
            aligner.query_internal_open_gap_score = 0
            aligner.query_internal_extend_gap_score = 0
            aligner.substitution_matrix =  substitution_matrices.load("BLOSUM62")
            alignments = aligner.align(seq1, seq2)

    # Printing alignment
    for alignment in alignments:
        print(alignment)
        print('Score: ', alignment.score)
        break
    return()

# When sequence file is uploaded
def file_pairwise_alignment(seq1_file, seq2_file, fileformat, mode = 'global', match = 1, missmatch = 0, gap_open = 0, gap_extend = 0):
    aligner = PairwiseAligner()
    aligner.mode = mode

    # defining alignment score values
    aligner.match_score = match
    aligner.mismatch_score = missmatch
    aligner.query_internal_open_gap_score = gap_open
    aligner.query_internal_extend_gap_score = gap_extend

    # Extracting sequence from file
    for record in SeqIO.parse(seq1_file, fileformat):
        seq1 = record.seq
        break
    for record in SeqIO.parse(seq2_file, fileformat):
        seq2 = record.seq
        break

    #Performing Alignment
    string = 'BDEFHIJKLMNOPQRSVWXYZ'
    for letter in string:
        if letter not in seq1 and seq2: # for DNA/RNA sequences
            alignments = aligner.align(seq1, seq2)
        else: # for protein alignment
            aligner.match_score = 1
            aligner.mismatch_score = 0
            aligner.query_internal_open_gap_score = 0
            aligner.query_internal_extend_gap_score = 0
            aligner.substitution_matrix =  substitution_matrices.load("BLOSUM62")
            alignments = aligner.align(seq1, seq2)

    for alignment in alignments:
        print(alignment)
        print('Score: ', alignment.score)
        break
    return()

# ----------------------functions for multiple sequence alignment--------------------

#---------Display alignment function()
def display_alignment(align_file, fileType):
    aligned = AlignIO.read(align_file, fileType)
    for record in aligned:
        print(f">{record.id}")
        print(record.seq)
        print()

    return()

 # -----Perform multiple sequence alignment using Clustal Omega---------
def multiple_seq_alignment_omega(input_file, output_file):

    clustalomega_exe = r"C:\Program Files (x86)\clustal-omega-1.2.2-win64\clustal-omega-1.2.2-win64\clustalo.exe"
    clustalomega_cline = ClustalOmegaCommandline(clustalomega_exe, infile=input_file, outfile=output_file, verbose=True, auto=True, force = True)
    assert os.path.isfile(clustalomega_exe), "Clustal omega executable missing"
    stdout, stderr = clustalomega_cline()

    print(f"Alignement file '{output_file}' saved in the current directory\n")
    alignment = AlignIO.read(output_file, 'fasta')
    print(alignment)

    choose = input('\nDisplay complete Alignment? [Y / N]: ').upper()
    if choose == 'Y':
        display_alignment(output_file, 'fasta')
    else:
        return()

# -------Perform multiple sequence alignment using Clustalw-------
def multiple_seq_alignment_clustalW(input_file, output_file):

    clustalW_exe = r"C:\Program Files (x86)\ClustalW2\clustalw2.exe"
    clustalW_cline = ClustalwCommandline(clustalW_exe, infile=input_file, outfile=output_file)
    assert os.path.isfile(clustalW_exe), "Clustal W executable missing"
    stdout, stderr = clustalW_cline()

    print(f"Alignement file '{output_file}' saved in the current directory\n")
    alignment = AlignIO.read(output_file, 'clustal')
    print(alignment)

    choose = input('\nDisplay complete Alignment? [Y / N]: ').upper()
    if choose == 'Y':
        display_alignment(output_file, 'clustal')
    else:
        return()

# -------Perform multiple sequence alignment using Muscles--------
def multiple_seq_alignment_muscles(input_file, output_file):

    muscle_exe = r"C:\Program Files (x86)\muscle-commandline\muscle5.1.win64.exe"
    muscle_cline = MuscleCommandline(muscle_exe, infile=input_file, outfile=output_file)
    assert os.path.isfile(muscle_exe), "muscle executable missing"
    stdout, stderr = muscle_cline()

    # Display output alignment file
    alignment = AlignIO.read(output_file, 'fasta')
    print(alignment)
    return()






def main():
    print('\n---------This is your biopython program, what would you like to perform?------------')

    choice = input('''
    Select any of the following options;

        1. Complement
        2. Reverse Complement
        3. Pairwise Sequence Alignment
        4. Multiple Sequence Alignment\n
        Enter your choice [1, 2, 3, 4]: ''')

    if choice not in '1234':
        print('Choice is invalid')
    elif choice == '1':
        print('\n----------------------------------')
        choice1 = input('''Complement:
        Do you have;
            a. raw sequence
            b. sequence file
            \n Enter choice [a or b]: ''').upper()

        if choice1 == 'A':
            print('\n----------------------------------')
            sequence = input('Input DNA/RNA sequence: ')
            seq_complement(sequence)
            print('\n----------------------------------')

        elif choice1 == 'B':
            print('\n----------------------------------')
            filename = input('Enter path to your file: ')
            fileformat = input('Enter file format: ')

            file_complement(filename, fileformat)
            print('\n----------------------------------')
        else:
            print('\n----------------------------------')
            print('Invalid input')
            print('\n----------------------------------')

    elif choice == '2':
        print('\n----------------------------------')
        choice1 = input('''Reverse Complement:
        Do you have;
            a. raw sequence
            b. sequence file
            \n Enter choice [a or b]: ''')
        if choice1 in 'Aa':
            print('\n----------------------------------')
            sequence = input('Input DNA/RNA sequence: ')
            seq_reverse_complement(sequence)
            print('\n----------------------------------')
        elif choice1 in 'Bb':
            print('\n----------------------------------')
            filename = input('Enter path to your file: ')
            fileformat = input('Enter file format: ')
            file_reverse_complement(filename, fileformat)
            print('\n----------------------------------')
        else:
            print('\n----------------------------------')
            print('Invalid input')
            print('\n----------------------------------')

    elif choice == '3':
        print('\n----------------------------------')
        print('------Pairwise Alignment------')
        choice1 = input('''Do you want to;
                1. Manualy enter match/missmatch and gap penalties
                2. Use default values.
                (default: match=1, missmatch=0, gap-penalties=0 | for protein seq: BLOSUM62 matrix)
            \nEnter choice [1 or 2]: ''')
        if choice1 == '1':
            print('\n----------------------------------')
            match_score = float(input('Match score: '))
            missmatch_score = float(input('Miss-match score: '))
            gap_opening = float(input('Gap_opening_penalty: '))
            gap_extend = float(input('Gap_extending_penalty: '))
            print('\n----------------------------------')
            choice2 = input('''Do you have;
            a. raw sequence
            b. sequence file
            \n Enter choice [a or b]: ''')
            if choice2 in 'Aa':
                print('\n----------------------------------')
                sequenceA = input('Input 1st sequence: ')
                sequenceB = input('Input 2nd sequence: ')
                mode = input('''Enter mode [global, local] *must check spelling*: ''')
                seq_pairwise_alignment(seq1=sequenceA, seq2=sequenceB, mode=mode,match=match_score, missmatch=missmatch_score, gap_open=gap_opening, gap_extend=gap_extend)
                print('\n----------------------------------')
            elif choice2 in 'Bb':
                print('\n----------------------------------')
                seq1 = input('Enter path to 1st sequence file: ')
                seq2 = input('Enter path to 2nd sequence file: ')
                mode = input('''Enter mode [global, local] *must check spelling*: ''')
                file_pairwise_alignment(seq1_file=seq1, seq2_file=seq2, mode=mode,match=match_score, missmatch=missmatch_score, gap_open=gap_opening, gap_extend=gap_extend)
                print('\n----------------------------------')
            else:
                print('\n----------------------------------')
                print('Invalid input')
                print('\n----------------------------------')
        elif choice1 == '2':
            print('\n----------------------------------')
            choice2 = input('''Do you have;
            a. raw sequence
            b. sequence file
            \n Enter choice [a or b]: ''')
            if choice2 in 'Aa':
                print('\n----------------------------------')
                sequenceA = input('Input 1st sequence: ')
                sequenceB = input('Input 2nd sequence: ')
                mode = input('''Enter mode [global, local] *must check spelling*: ''')
                seq_pairwise_alignment(seq1=sequenceA, seq2=sequenceB, mode=mode)
                print('\n----------------------------------')
            elif choice2 in 'Bb':
                print('\n----------------------------------')
                seq1 = input('Enter path to 1st sequence file: ')
                seq2 = input('Enter path to 2nd sequence file: ')
                mode = input('''Enter mode [global, local] *must check spelling*: ''')
                file_pairwise_alignment(seq1_file=seq1, seq2_file=seq2, mode=mode)
                print('\n----------------------------------')
            else:
                print('\n----------------------------------')
                print('Invalid input')
                print('\n----------------------------------')
        else:
            print('\n----------------------------------')
            print('Invalid input')
            print('\n----------------------------------')


    elif choice == '4':
        print('\n----------------------------------')
        inputFile = input("Enter path to the file containing all the sequences to be aligned in FASTA format: ")
        outFile = input("Enter name of your output file: ")
        print('\n----------------------------------')
        choice1 = input('''\nChoose the method for multiple sequence alignment:
            a. Clustal Omega
            b. Clustal W
            c. Muscles
        \nEnter your choice [a, b, c]: ''').upper()
        print('\n----------------------------------')
        if choice1 == 'A':
            multiple_seq_alignment_omega(inputFile, outFile)

        elif choice1 == 'B':
            multiple_seq_alignment_clustalW(inputFile, outFile)

        elif choice1 == 'C':
            multiple_seq_alignment_muscles(inputFile, outFile)

        else:
            print('Invalid Input......')
            print('\n----------------------------------')

    else:
        print('Invalid Input......')
        print('\n----------------------------------')


    run = input('''Try again?
        \nEnter Yes or No [Y / N]: ''').upper()
    if run == 'Y':
        os.system('cls')
        main()
    else:
        print('\n----------------------------------')
        print('Thank you for using my program😇')
        print('\n----------------------------------')
        exit()


# run main()
main()


