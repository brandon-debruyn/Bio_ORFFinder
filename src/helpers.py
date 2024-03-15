# Author: Brandon Lee De Bruyn
# Date: 12/04/2024
# Description: A helper module for the ORF Finder application.

from enum import Enum

class Strand(Enum):
    # enum for sequence types
    POSITIVE = '+'
    NEGATIVE = '-'

class StartCodon(Enum):
    # enum for start codon types
    AUTOMATIC = 'Detect Automatically'
    ATG = 'ATG'
    AUG = 'AUG'