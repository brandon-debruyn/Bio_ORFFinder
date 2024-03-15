# Author: Brandon Lee De Bruyn
# Date: 12/04/2024
# Description: A test suite for the ORFFinder class

import unittest
from src.orf_finder import ORFFinder

class TestORFFinder(unittest.TestCase):
    """
        Test the ORFFinder class using the unittest module (extremely familiar and popular module for unit testing in C#, Python).
        
        The ORFFinder class is responsible for finding open reading frames in a DNA sequence.
        
        The class has the following methods:
        - _extract_descriptor_and_sequence: Extracts the descriptor and sequence from the fasta content input stream
        - _is_valid_sequence: Validates the sequence
        - find_orf: Finds the open reading frames in the sequence

        Simple tests are implemented for now...
    """
    
    def test_descriptor_and_sequence_extraction(self):
        """
            Test that the descriptor and sequence are extracted correctly from fasta content input stream
        """

        fasta_content = """>Example Descriptor
        ATGAAATGA
        """
        orf_finder = ORFFinder(fasta_content)
        descriptor, sequence = orf_finder._extract_descriptor_and_sequence(fasta_content)
        self.assertEqual(descriptor, "Example Descriptor")
        self.assertEqual(sequence, "ATGAAATGA")

    def test_sequence_validation(self):
        """
            Test that the sequence is validated correctly
        """
        valid_sequence = "ATGAAATGA"  
        invalid_sequence = "ATGXXXTGA" 
        orf_finder = ORFFinder(valid_sequence)
        is_valid, message = orf_finder._is_valid_sequence(valid_sequence)
        self.assertTrue(is_valid)
        
        orf_finder = ORFFinder(invalid_sequence)
        is_valid, message = orf_finder._is_valid_sequence(invalid_sequence)
        self.assertFalse(is_valid)

    def test_orf_finding(self):
        """
            Test that ORFs are found correctly
        """

        sequence = "ATGAAATGA"
        orf_finder = ORFFinder(sequence)
        orf_finder.find_orf()  
        self.assertEqual(len(orf_finder.orf_data['orfs']), 1)  

# run the tests
if __name__ == '__main__':
    unittest.main()