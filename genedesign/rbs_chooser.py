from genedesign.models.rbs_option import RBSOption
from genedesign.seq_utils.hairpin_counter import hairpin_counter
from genedesign.seq_utils.Translate import Translate
from genedesign.seq_utils.calc_edit_distance import calculate_edit_distance

from typing import Set

class RBSChooser:
    """
    A class to choose the best RBS for a given CDS sequence.
    """
    # Class variable to store RBS options
    rbs_options: Set[RBSOption] = set()  # Initialize with an empty set

    def initiate(self, rbs_option_list: Set[RBSOption]) -> None:
        """
        Initialization method for RBSChooser.
        It populates the class variable with RBS options.

        Parameters:
            rbs_option_list (Set[RBSOption]): A set of RBSOption instances.
        """
        self.rbs_options = rbs_option_list

    def run(self, cds: str, ignores: Set[RBSOption]) -> RBSOption:
        """
        A method to choose the best RBS for a given CDS sequence.

        Parameters:
            cds (str): The coding sequence to analyze.
            ignores (Set[RBSOption]): A set of RBSOption instances to exclude from consideration.

        Returns:
            RBSOption: The best RBSOption based on the specified criteria.
        """
        # Step 1: Check for valid length
        if len(cds) % 3 != 0:
            print(f"Skipping invalid CDS sequence (length not multiple of 3): {cds}")
            return None  # Skip this sequence or raise a ValueError if needed

        #Step 2: Exclude RBSOptions based on ignores
        candidates = [option for option in self.rbs_options if option not in ignores]

        print(f"Number of candidates after exclusion: {len(candidates)}")
        # Prepare a list to hold viable candidates (those that do not form hairpins)
        viable_candidates = []
        for option in candidates:
            # Combine RBS UTR and CDS for analysis
            combined_sequence = option.utr + cds

            # Check for hairpins using hairpin_counter
            hairpin_count = hairpin_counter(combined_sequence)

            # Debugging: Print hairpin counts for each option
            print(f"RBS: {option.utr}, Hairpin Count: {hairpin_count}")

            # Only keep candidates with no hairpins (or adjust logic as needed)
            if hairpin_count == 0:
              viable_candidates.append(option)

        # Debugging: Print the number of viable candidates
        print(f"Number of viable candidates after hairpin check: {len(viable_candidates)}")
        # Step 3: Translate CDS to peptide
        translator = Translate()  # Create an instance of the Translate class
        translator.initiate()  # Initialize the translator
        protein_sequence = translator.run(cds)  # Get the protein sequence from the CDS
        cds_peptide = protein_sequence[:6]  # Get the first six amino acids

        best_rbs = None
        best_edit_distance = float('inf')  # Initialize with infinity for comparison

        # Step 4: Compare the input peptide to the RBS source gene's peptide
        for rbs in viable_candidates:
            edit_distance = calculate_edit_distance(cds_peptide, rbs.first_six_aas)
            print(f"RBS: {rbs.utr}, Edit Distance: {edit_distance}")

            # Select the RBS with the smallest edit distance
            if edit_distance < best_edit_distance:
                best_edit_distance = edit_distance
                best_rbs = rbs

        print(f"Chosen RBS Option: {best_rbs}")
        return best_rbs

if __name__ == "__main__":
    # Example usage of RBSChooser
    cds = "ATGGTAAGAAAACAGTTGCAGAGAGTTGAATT..."

    # Initialize the chooser
    chooser = RBSChooser()
    chooser.initiate()

    # Choose RBS with no ignores
    ignores = set()
    selected1 = chooser.run(cds, ignores)
    
    # Add the first selection to the ignore list
    ignores.add(selected1)
    
    # Choose another RBS option after ignoring the first
    selected2 = chooser.run(cds, ignores)

    # Print the selected RBS options
    print("Selected1:", selected1)
    print("Selected2:", selected2)
