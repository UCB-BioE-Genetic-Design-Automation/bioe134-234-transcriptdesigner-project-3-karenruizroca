import unittest
from genedesign.transcript_designer import TranscriptDesigner
from genedesign.rbs_chooser import RBSChooser
from unittest.mock import MagicMock

class TestTranscriptDesigner(unittest.TestCase):
    def setUp(self):
        self.designer = TranscriptDesigner()
        # Mock RBSChooser to avoid dependency on its functionality
        self.designer.rbsChooser = MagicMock()
        self.designer.rbsChooser.initiate = MagicMock()
        self.designer.rbsChooser.run = MagicMock(return_value="MockedRBS")
        self.designer.initiate()
    
    def test_generate_codon(self):
        """Test that generate_codon produces valid codons for a given amino acid."""
        valid_amino_acids = "ACDEFGHIKLMNPQRSTVWY"
        for aa in valid_amino_acids:
            codon = self.designer.generate_codon(aa)
            self.assertEqual(len(codon), 3)  # Codon must be 3 bases
            self.assertIn(codon, self.designer.aminoAcidToCodon.values())

    def test_rbs_chooser_initiation(self):
        """Test that RBSChooser is properly initiated."""
        self.designer.rbsChooser.initiate.assert_called_once()  # Check initiate was called
        self.assertIsInstance(self.designer.rbsChooser, MagicMock)  # Mocked instance

    def test_sliding_window(self):
        """Test the sliding window functionality."""
        peptide = "MYPF"
        optimized_transcript = self.designer.sliding_window_optimization(peptide)
        # Ensure the final transcript contains the stop codon
        self.assertTrue(optimized_transcript.endswith("TAA"))
        # Ensure the length matches the peptide (3 bases per AA + stop codon)
        expected_length = len(peptide) * 3 + 3  # 3 bases per AA + 3 for stop codon
        self.assertEqual(len(optimized_transcript), expected_length)

    def test_run(self):
        """Test the full run method."""
        peptide = "MYPF"
        ignores = set()
        transcript = self.designer.run(peptide, ignores)
        self.assertEqual(transcript.peptide, peptide)  # Peptide should match input
        self.assertTrue(len(transcript.codons) > 0)  # Ensure codons are generated
        self.assertEqual(transcript.codons[-1], "TAA")  # Ensure stop codon is added

    def test_error_handling_invalid_amino_acid(self):
        """Test that invalid amino acids raise an error."""
        with self.assertRaises(KeyError):  # Expect KeyError for unknown amino acids
            self.designer.generate_codon("Z")  # 'Z' is not a valid amino acid

if __name__ == "__main__":
    unittest.main()
