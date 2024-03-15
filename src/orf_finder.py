# Author: Brandon Lee De Bruyn
# Date: 12/04/2024
# Description: A class for finding open reading frames (ORFs) in a DNA or RNA sequence.

from helpers import Strand, StartCodon

class ORFFinder:
    """
        Class for finding open reading frames (ORFs) in a DNA or RNA sequence. 

        Attributes:
        - sequence: str, the DNA or RNA sequence to search for ORFs.
        - orf_data: dict, a dictionary containing the ORF data.

        Methods:
        - __init__(sequence: str) -> None: Initializes the ORFFinder object with a given sequence input stream.
        - _extract_descriptor_and_sequence(fasta_content: str) -> tuple[str, str]: Extracts the descriptor and sequence from a FASTA/TXT formatted content.
        - _get_dna_complement_sequence(sequence: str) -> str: Returns the complement sequence of a given DNA sequence.
        - _is_valid_sequence(sequence: str) -> tuple[bool, str]: Checks if a given complete sequence is valid for ORF search.
        - _translate_to_amino_acid(sequence: str) -> str: Translates a given DNA or RNA sequence to an amino acid sequence.
        - __get_correct_nucleotide_position(position: int, sequence_length: int) -> int: Maps the position in the reverse complement back to the original sequence's position.
        - _find_orfs_in_frame_sequence(sequence: str, strand: str): Finds ORFs in a given frame sequence.
        - find_orfs() -> dict: Finds open reading frames (ORFs) in a DNA or RNA sequence.
    """

    # algorithm const definitions
    # ---------------------------
    MIN_SEQUENCE_LENGTH = 30 # minimum acceptable sequence length 
    MIN_ORF_LENGTH = 30 # minimum acceptable ORF length

    VALID_NUCLEOTIDES_DNA = {'A', 'T', 'C', 'G', 'N'}  # for DNA
    VALID_NUCLEOTIDES_RNA = {'A', 'U', 'C', 'G', 'N'}  # for RNA

    START_CODON_DNA = "ATG" # start codon for DNA

    STOP_CODONS_DNA = {"TAA", "TAG", "TGA"} # stop codons DNA
    STOP_CODONS_RNA = {"UAA", "UAG", "UGA"} # stop codons RNA

    IS_DNA = True # default sequence type is DNA

    # map codons to amino acids
    DNA_CODON_TO_AA_MAP = {
        "ATA":"I", "ATC":"I", "ATT":"I", "ATG":"M",
        "ACA":"T", "ACC":"T", "ACG":"T", "ACT":"T",
        "AAC":"N", "AAT":"N", "AAA":"K", "AAG":"K",
        "AGC":"S", "AGT":"S", "AGA":"R", "AGG":"R",
        "CTA":"L", "CTC":"L", "CTG":"L", "CTT":"L",
        "CCA":"P", "CCC":"P", "CCG":"P", "CCT":"P",
        "CAC":"H", "CAT":"H", "CAA":"Q", "CAG":"Q",
        "CGA":"R", "CGC":"R", "CGG":"R", "CGT":"R",
        "GTA":"V", "GTC":"V", "GTG":"V", "GTT":"V",
        "GCA":"A", "GCC":"A", "GCG":"A", "GCT":"A",
        "GAC":"D", "GAT":"D", "GAA":"E", "GAG":"E",
        "GGA":"G", "GGC":"G", "GGG":"G", "GGT":"G",
        "TCA":"S", "TCC":"S", "TCG":"S", "TCT":"S",
        "TTC":"F", "TTT":"F", "TTA":"L", "TTG":"L",
        "TAC":"Y", "TAT":"Y", "TAA":"*", "TAG":"*",
        "TGC":"C", "TGT":"C", "TGA":"*", "TGG":"W",
    }

    # complement nucleotide map for DNA and ambiguous nucleotides complement map 
    COMPLEMENT_MAP = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

    # all based on the IUPAC (International Union of Pure and Applied Chemistry) nucleotide code
    # Johnson AD. An extended IUPAC nomenclature code for polymorphic nucleic acids. Bioinformatics. 2010 May 15;26(10):1386-9. 
    # doi: 10.1093/bioinformatics/btq098. Epub 2010 Mar 3. PMID: 20202974; PMCID: PMC2865858.
    AMBIGUOUS_NUCLEOTIDES_DNA = {'W', 'Y', 'M', 'B', 'K', 'S', 'R', 'N', 'V', 'H', 'D', 'U'}
    AMBIGUOUS_COMPLEMENT_MAP = { 
        'R': 'Y', 'Y': 'R',  # purine and pyrimidine
        'K': 'M', 'M': 'K',  # keto and amino
        'S': 'S', 'W': 'W',  # strong and weak
        'B': 'V', 'V': 'B',  # not A and not T
        'D': 'H', 'H': 'D',  # not C and not G
        'N': 'N'             # any nucleotide
    }

    def __init__(self, sequence, min_orf_length=MIN_ORF_LENGTH, start_codon_option=StartCodon.AUTOMATIC) -> None:
        self.descriptor, self.sequence = self._extract_descriptor_and_sequence(sequence)
        self.orf_data = {
            'sequence_type': None,
            'complete_sequence_length': None,
            'complete_sequence': self.sequence,
            'descriptor': None,
            'no_orfs': 0,
            'orfs': []
        }
        self.min_orf_length = min_orf_length
        self.start_codon_option = start_codon_option

        # set the sequence type based on the start codon option
        if (self.start_codon_option == StartCodon.AUG):
            self.IS_DNA = False

    def _extract_descriptor_and_sequence(self, fasta_content) -> tuple[str, str]:
        """
            Extracts the descriptor and sequence from a FASTA/TXT formatted content.

            Parameters:
            - fasta_content: str, the entire content of a FASTA file or entry as a string.

            Returns:
            - descriptor: str, the descriptor line from the FASTA/TXT content, without the '>' prefix.
            - sequence: str, the sequence from the FASTA/TXT content.
        """

        # split the content into lines and extract the descriptor line, if present otherwise return not found
        lines = fasta_content.split('\n')
        descriptor_line = "Descriptor not found."
        sequence = []

        # iterate over the lines and extract the descriptor line and sequence
        for line in lines:
            if line.startswith('>'):
                descriptor_line = line[1:].strip() # remove the '>' prefix and strip whitespace
            else:
                sequence.append(line.strip()) # append the sequence line to the sequence list

        # return the descriptor line and sequence
        return descriptor_line, ''.join(sequence) if sequence else None

    def _get_dna_complement_sequence(self, sequence) -> str:
        """
            Returns the complement sequence of a given DNA sequence.

            Parameters:
                - sequence: str, a DNA sequence.
            
            Returns:
                - complement_sequence: str, the complement sequence of the input sequence.
        """
        # get the full complement map by combining the complement map and ambiguous complement dicts/maps (curious unpacking notation **)
        full_complement_map = {
            **self.COMPLEMENT_MAP,
            **self.AMBIGUOUS_COMPLEMENT_MAP
        }

        # get the complement sequence by mapping each nucleotide to its complement
        complement_sequence = ''.join([full_complement_map.get(nucleotide, 'N') for nucleotide in sequence.upper()])
        reverse_complement_sequence = complement_sequence[::-1] # reverse the complement sequence

        # return the reverse complement sequence
        return reverse_complement_sequence
            
    def __get_rna_codon_to_aa_map(self) -> dict:
        """
            Returns a map of RNA codons to amino acids by converting the existing DNA codon to amino acid map.

            Returns:
                - rna_codon_map: dict, a map of RNA codons to amino acids based on the DNA existing map.
        """
        # convert the DNA codon to amino acid map to RNA
        rna_codon_map = {}
        for codon, aa in self.DNA_CODON_TO_AA_MAP.items():
            rna_codon = codon.replace('T', 'U')
            rna_codon_map[rna_codon] = aa
        return rna_codon_map
    
    def _is_valid_sequence(self, sequence) -> tuple[bool, str]:
        """
            Checks if a given complete sequence is valid for ORF search.

            Parameters:
                - sequence: str, a DNA or RNA sequence.
            
            Returns:
                - is_valid: bool, True if the sequence is valid, False otherwise.
                - message: str, a message indicating the reason for invalidity, if applicable.
        """
        message = ""

        # check if the start codon option is set to automatic
        if self.start_codon_option == StartCodon.AUTOMATIC:
            # determine if sequence is DNA or RNA
            if 'U' in sequence and 'T' not in sequence:
                self.IS_DNA = False
            elif 'T' in sequence and 'U' not in sequence:
                self.IS_DNA = True
            else:
                # handle ambiguous sequences by assuming DNA
                message += "Warning: Ambiguous sequence type. Assuming DNA sequence, ambiguious nucleotides will be assumed as X in AA sequence\n\n\n"
                self.IS_DNA = True
                print(message)
        
        # set valid nucleotides, start, and stop codons dependent on sequence type
        valid_nucleotides = self.VALID_NUCLEOTIDES_DNA if self.IS_DNA else self.VALID_NUCLEOTIDES_RNA
        start_codon = self.START_CODON_DNA if self.IS_DNA else self.START_CODON_DNA.replace('T', 'U')
        stop_codons = self.STOP_CODONS_DNA if self.IS_DNA else self.STOP_CODONS_RNA

         # check minimum length
        if len(sequence) < self.MIN_SEQUENCE_LENGTH:
            message += "Error: Sequence is too short to contain a meaningful ORF.\n\n"
            return False, message

        # check for invalid characters
        invalid_chars = set(sequence) - valid_nucleotides
        if invalid_chars:
            message += "Warning: Ambiguous nucleotides detected, ambiguious codons with ambiguous nucleotides will be assumed as X in AA sequence translation.\n\n"
            message += f"Warning: Invalid nucleotides found in the sequence: {', '.join(invalid_chars)}\n\n"

        # check for start codon presence
        if start_codon not in sequence:
            message += "Error: No start codon found in the sequence.\n\n"
            return False, message

        # optionally check for stop codon presence
        if not any(stop_codon in sequence for stop_codon in stop_codons):
            message += "Error: No stop codon found in the sequence.\n\n"
            return False, message

        # check divisibility by three
        if len(sequence) % 3 != 0:
            message += "Warning: Sequence length is not a multiple of three. Trailing nucleotides will be ignored.\n\n"

        message += "Sequence passed all checks.\n\n"
        return True, message
    
    def _translate_to_amino_acid(self, sequence) -> str:
        """
            Translates a given DNA or RNA sequence to an amino acid sequence using the codon to amino acid map, setting ambiguous nucleotides dependent on sequence type 
            denoted by X.

            Parameters:
                - sequence: str, a DNA or RNA sequence.
            
            Returns:
                - amino_acid_sequence: str, the amino acid sequence.
        """
        # get the codon to amino acid map based on DNA or RNA and set ambiguous nucleotides dependent on sequence type
        codon_to_aa_map = self.DNA_CODON_TO_AA_MAP if self.IS_DNA else self.__get_rna_codon_to_aa_map()
        ambiguous_nucleotides = self.AMBIGUOUS_NUCLEOTIDES_DNA if self.IS_DNA else {n if n != 'U' else 'T' for n in self.AMBIGUOUS_NUCLEOTIDES_DNA}

        # translate the sequence to an amino acid sequence
        amino_acid_sequence = ""
        for i in range(0, len(sequence) - 2, 3):
            codon = sequence[i : i + 3]

            # check for any nucleotide in the current codon that is ambiguous and mark as X otherwise lookup the amino acid (easier as iterable)
            if any(nucleotide in ambiguous_nucleotides for nucleotide in codon):
                amino_acid_sequence += 'X'
            else:
                # translate the codon to an amino acid, defaulting to X for unknown codons
                amino_acid = codon_to_aa_map.get(codon, 'X')
                if amino_acid == "*":  # Stop codon
                    break

                amino_acid_sequence += amino_acid

        return amino_acid_sequence
    
    def __get_correct_nucleotide_position(self, position, sequence_length) -> int:
        """
            Maps the position in the reverse complement back to the original sequence's position

            Parameters:
                - position: int, the position in the reverse complement sequence.
                - sequence_length: int, the length of the complete sequence.
            
            Returns:    
                - int, the position in the original sequence.
        """
        return sequence_length - position + 1

    def _find_orfs_in_frame_sequence(self, sequence, strand):
        """
            Finds ORFs in a given frame sequence.

            Parameters:
                - sequence: str, a DNA or RNA sequence.
                - strand: str, the strand of the sequence (positive or negative).
            
            Returns:
                - orf_data: dict, a dictionary containing the ORF data.
        """
        # set start and stop codons dependent on sequence type
        start_codon = self.START_CODON_DNA if self.IS_DNA else self.START_CODON_DNA.replace('T', 'U')
        stop_codons = self.STOP_CODONS_DNA if self.IS_DNA else self.STOP_CODONS_RNA

        # three reading frames per sequence
        for frame in range(3):
            i = frame
            while i + 2 <= len(sequence): # for each full codon in the nucleotide framed sequence 
                # get the current codon and check if it is a start codon
                codon = sequence[i : i + 3]
                if codon == start_codon:
                    for j in range(i + 3, len(sequence), 3): # iterate over all full codons untill the next stop codon
                        if sequence[j : j + 3] in stop_codons:  # if a stop codon is found, store the ORF and move to the next codon

                            # calculate the ORF length and if the ORF falls within the minimum length store the ORF 
                            orf_length = j - i + 3
                            if orf_length >= self.min_orf_length: 
                                start_nucleotide = i + 1
                                stop_nucleotide = j + 3
                                if strand == Strand.NEGATIVE: # if the strand is negative, get the correct nucleotide position
                                    start_nucleotide = self.__get_correct_nucleotide_position(start_nucleotide, self.orf_data['complete_sequence_length'])
                                    stop_nucleotide = self.__get_correct_nucleotide_position(stop_nucleotide, self.orf_data['complete_sequence_length'])

                                self.orf_data['no_orfs'] += 1                            
                                self.orf_data['orfs'].append({
                                    'strand': strand.value, # store the strand +/-
                                    'frame': frame + 1,
                                    'start_nucleotide': start_nucleotide,
                                    'stop_nucleotide': stop_nucleotide,
                                    'length_nt': orf_length,
                                    'length_aa': (orf_length - 3) // 3,
                                    'nucleotide_sequence': sequence[i : j + 3], # store the ORF nucleotide sequence
                                    'aa_sequence': self._translate_to_amino_acid(sequence[i : j + 3]) # translate the ORF to an amino acid sequence
                                })
                                i = j + 3
                                break
                                
                    else:
                        i += 3
                        continue
                i += 3

    def find_orfs(self) -> dict:
        """
            Finds open reading frames (ORFs) in a DNA or RNA sequence.

            Returns:
                - orf_data: dict, a dictionary containing the ORF data.
        """
        # check if sequence is valid
        is_valid, message = self._is_valid_sequence(self.sequence)
        if (is_valid):
            print(f"Sequence is valid: {message}")

            # find ORFs in the positive strand both RNA and DNA
            strands_to_search = [Strand.POSITIVE]

            # if sequence is DNA, also search the negative strand
            if (self.IS_DNA):
                strands_to_search.append(Strand.NEGATIVE)

            # build ORF data for UI
            self.orf_data['sequence_type'] = 'DNA' if self.IS_DNA else 'RNA'
            self.orf_data['complete_sequence_length'] = len(self.sequence)
            self.orf_data['complete_sequence'] = self.sequence
            self.orf_data['descriptor'] = self.descriptor

            # find ORFs in the specified strands
            for strand in strands_to_search:
                    sequence = self.sequence 

                    # if current strand is negative then get complement sequence
                    if (strand == Strand.NEGATIVE):
                        sequence = self._get_dna_complement_sequence(self.sequence) 

                    # find ORFs in the given sequence
                    self._find_orfs_in_frame_sequence(sequence, strand)

            # display ORF data to the console and return the data to a caller parent
            display_data = message
            if self.orf_data['no_orfs'] > 0:

                # overview of findings
                display_data += f"{self.orf_data['no_orfs']} ORFs were found for < {self.orf_data['descriptor']} > in a sequence of length {self.orf_data['complete_sequence_length']} nucleotides.\n\n"

                # displaying the initial part of the complete nucleotide sequence
                display_data += "Complete Nucleotide Sequence (first 50 nt only): {}\n\n".format(self.orf_data['complete_sequence'][:50])
                                
                # adding the header for the summary table
                display_data += f"{'Label':<6} {'Strand':<8} {'Frame':<6} {'Start':<6} {'Stop':<6} {'Length (nt | aa)':<18} {'Summary'}\n"

                # adding ORF summary information to the table
                for idx, orf in enumerate(self.orf_data['orfs'], start=1):
                    summary_info = f"ORF{idx} ({orf['length_nt']} nt | {orf['length_aa']} aa)"
                    display_data += f"ORF{idx:<5} {orf['strand']:<8} {orf['frame']:<6} {orf['start_nucleotide']:<6} {orf['stop_nucleotide']:<6} {orf['length_nt']} | {orf['length_aa']:<15} {summary_info}\n"
            
                # display detailed ORF data
                for idx, orf in enumerate(self.orf_data['orfs'], start=1):
                    # Add nucleotide and amino acid sequences for each ORF
                    display_data += f"ORF{idx} details:\n"
                    display_data += f"  Nucleotide Sequence: {orf['nucleotide_sequence']}\n"
                    display_data += f"  Amino Acid Sequence: {orf['aa_sequence']}\n\n"

                print(display_data)
                
            else:
                message += "No ORFs were found in the sequence.\n\n"
                display_data = message
                print(display_data)

            # return the ORF data and display data
            return self.orf_data, display_data
            
        # if sequence is invalid return None     
        else:
            print(f"Sequence is invalid:\n\n {message}")
            return None, message