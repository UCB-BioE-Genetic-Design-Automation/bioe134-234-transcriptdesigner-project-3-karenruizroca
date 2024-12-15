from genedesign.rbs_chooser import RBSChooser
from genedesign.models.transcript import Transcript
from genedesign.checkers.codon_checker import CodonChecker
from genedesign.checkers.forbidden_sequence_checker import ForbiddenSequenceChecker
from genedesign.checkers.hairpin_checker import hairpin_checker
from genedesign.checkers.internal_promoter_checker import PromoterChecker
import random
from concurrent.futures import ThreadPoolExecutor

class TranscriptDesigner:
    def __init__(self):
        self.aminoAcidToCodonWeighted = {}
        self.rbsChooser = None
        self.codonChecker = CodonChecker()
        self.forbiddenChecker = ForbiddenSequenceChecker()
        self.promoterChecker = PromoterChecker()

    def initiate(self) -> None:
        """
        Initializes the codon table, RBS chooser, and all checkers.
        """
         # Use an empty set of RBS options initially
        empty_rbs_options = set()
        self.rbsChooser = RBSChooser()
        self.rbsChooser.initiate(empty_rbs_options)

        # Initialize checkers
        self.codonChecker.initiate()
        self.forbiddenChecker.initiate()
        self.promoterChecker.initiate()

        # Weighted codon table (simplified example; replace with full weights)
        self.aminoAcidToCodonWeighted = {
            'A': ['GCG'] * 40 + ['GCC'] * 35 + ['GCA'] * 20 + ['GCT'] * 5,
            'C': ['TGC'] * 55 + ['TGT'] * 45,
            'D': ['GAT'] * 60 + ['GAC'] * 40,
            'E': ['GAA'] * 70 + ['GAG'] * 30,
            'F': ['TTC'] * 80 + ['TTT'] * 20,
            'G': ['GGT'] * 35 + ['GGC'] * 30 + ['GGA'] * 25 + ['GGG'] * 10,
            'H': ['CAC'] * 65 + ['CAT'] * 35,
            'I': ['ATC'] * 50 + ['ATT'] * 40 + ['ATA'] * 10,
            'K': ['AAA'] * 75 + ['AAG'] * 25,
            'L': ['CTG'] * 49 + ['TTA'] * 13 + ['TTG'] * 13 + ['CTT'] * 11 + ['CTC'] * 10 + ['CTA'] * 4,
            'M': ['ATG'] * 100,
            'N': ['AAC'] * 60 + ['AAT'] * 40,
            'P': ['CCG'] * 50 + ['CCC'] * 20 + ['CCA'] * 20 + ['CCT'] * 10,
            'Q': ['CAG'] * 70 + ['CAA'] * 30,
            'R': ['CGT'] * 40 + ['CGC'] * 30 + ['CGA'] * 15 + ['CGG'] * 10 + ['AGA'] * 5,
            'S': ['TCT'] * 30 + ['TCC'] * 25 + ['TCA'] * 20 + ['TCG'] * 15 + ['AGT'] * 5 + ['AGC'] * 5,
            'T': ['ACC'] * 50 + ['ACT'] * 25 + ['ACA'] * 15 + ['ACG'] * 10,
            'V': ['GTT'] * 45 + ['GTC'] * 25 + ['GTA'] * 20 + ['GTG'] * 10,
            'W': ['TGG'] * 100,
            'Y': ['TAC'] * 60 + ['TAT'] * 40,
        }

    def generate_codon(self, amino_acid: str) -> str:
        """
        Randomly selects a codon for the given amino acid based on usage frequency.
        """
        codon_list = self.aminoAcidToCodonWeighted[amino_acid]
        return random.choice(codon_list)

    def evaluate_candidate(self, candidate: str) -> float:
        """
        Evaluates a candidate DNA sequence using all four checkers.
        Returns a cumulative penalty score (lower is better).
        """
        penalty = 0.0

        # Check codon properties
        codons = [candidate[i:i + 3] for i in range(0, len(candidate), 3)]
        codon_pass, diversity, rare_count, cai = self.codonChecker.run(codons)
        if not codon_pass:
            penalty += 50  # Penalty for failing thresholds
        penalty += 1.0 - cai  # Penalize low CAI
        penalty += rare_count  # Penalize rare codons

        # Check forbidden sequences
        forbidden_pass, forbidden_seq = self.forbiddenChecker.run(candidate)
        if not forbidden_pass:
            penalty += 100  # Heavy penalty for forbidden sequences

        # Check hairpins
        hairpin_pass, hairpin_seq = hairpin_checker(candidate)
        if not hairpin_pass:
            penalty += 50  # Penalty for problematic hairpins

        # Check internal promoters
        promoter_pass, promoter_seq = self.promoterChecker.run(candidate)
        if not promoter_pass:
            penalty += 100  # Heavy penalty for promoters

        return penalty

    def monte_carlo_optimize(self, window: str, num_candidates: int = 10) -> str:
        """
        Optimize a DNA window using Monte Carlo sampling with weighted codon selection.
        """
        candidates = [
            ''.join(self.generate_codon(aa) for aa in window)
            for _ in range(num_candidates)
        ]

        with ThreadPoolExecutor() as executor:
            scores = list(executor.map(self.evaluate_candidate, candidates))

        # Return the candidate with the lowest score
        best_index = scores.index(min(scores))
        return candidates[best_index]

    def sliding_window_optimization(self, peptide: str, num_candidates: int = 10) -> str:
        """
        Optimize the peptide sequence using sliding windows and Monte Carlo sampling.
        """
        peptide_codons = [self.generate_codon(aa) for aa in peptide]
        preamble = ""  # Initialize with an empty preamble
        optimized_transcript = preamble
        window_size = 9  # 3 codons

        for i in range(0, len(peptide_codons) - 3 + 1, 3):  # Slide by one codon
            window = optimized_transcript[-3:] + ''.join(peptide_codons[i:i + 3])
            best_window = self.monte_carlo_optimize(window, num_candidates)
            optimized_transcript += best_window[3:6]  # Append only the middle codon

        # Handle boundary case (end of sequence)
        final_window = optimized_transcript[-3:] + ''.join(peptide_codons[-3:])
        optimized_transcript += self.monte_carlo_optimize(final_window, num_candidates)[3:6]

        # Append stop codon
        optimized_transcript += "TAA"

        return optimized_transcript

    def run(self, peptide: str, ignores: set) -> Transcript:
        """
        Optimized transcript design using sliding windows and Monte Carlo sampling.

        Parameters:
            peptide (str): The protein sequence to translate.
            ignores (set): RBS options to ignore.

        Returns:
            Transcript: The transcript object with the selected RBS and translated codons.
        """
        optimized_cds = self.sliding_window_optimization(peptide)

        # Choose an RBS
        selectedRBS = self.rbsChooser.run(optimized_cds, ignores)

        # Convert to codons
        codons = [optimized_cds[i:i + 3] for i in range(0, len(optimized_cds), 3)]

        return Transcript(selectedRBS, peptide, codons)

if __name__ == "__main__":
    # Example usage of TranscriptDesigner
    peptide = "MYPFIRTARMTV"  # Example peptide sequence
    designer = TranscriptDesigner()
    designer.initiate()

    ignores = set()  # Example ignores set
    transcript = designer.run(peptide, ignores)

    # Print out the transcript information
    print(transcript)