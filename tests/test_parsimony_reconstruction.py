import unittest
import sys
import os
from typing import Dict, List

# Adjust path to import script from parent directory
# This assumes tests are run from the project root or 'tests' directory
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from scripts.parsimony_reconstruction import (
    PhyloNode,
    parse_newick,
    calculate_initial_sequences,
    reconstruct_ancestral_sequences_parsimony,
    apply_mutations
)

def get_all_node_sequences(node: PhyloNode) -> Dict[str, str]:
    """Helper to collect all inferred sequences from a tree into a dictionary."""
    sequences: Dict[str, str] = {}

    q: List[PhyloNode] = [node]
    while q:
        curr = q.pop(0)
        sequences[curr.name] = "".join(curr.inferred_sequence)
        for child in curr.children:
            q.append(child)
    return sequences

class TestParsimonyReconstruction(unittest.TestCase):

    def test_simple_tree(self):
        newick_string = "(Leaf1:[0A>T],Leaf2:[0A>C])Anc1;"
        root_sequence_str = "AGTC"

        expected_sequences = {
            "Anc1": "YGTC", # Root after parsimony
            "Leaf1": "TGTC",
            "Leaf2": "CGTC"
        }

        tree_root = parse_newick(newick_string)
        calculate_initial_sequences(tree_root, list(root_sequence_str))
        reconstruct_ancestral_sequences_parsimony(tree_root)

        actual_sequences = get_all_node_sequences(tree_root)
        self.assertEqual(actual_sequences, expected_sequences)

    def test_complex_tree(self):
        newick_string = "((LeafA:[1G>T],LeafB:[2T>A])Anc1:[0C>G],(LeafC:[3A>G],LeafD:[0C>T,4G>C])Anc2)Root;"
        root_sequence_str = "CGTAG"

        # Revised expected sequences after careful re-tracing
        expected_sequences = {
            "Root": "BGTAG",
            "Anc1": "GKWAG",
            "LeafA": "GTTAG",
            "LeafB": "GGAAG",
            "Anc2": "YGTRS",
            "LeafC": "CGTGG",
            "LeafD": "TGTAC"
        }

        tree_root = parse_newick(newick_string)
        calculate_initial_sequences(tree_root, list(root_sequence_str))
        reconstruct_ancestral_sequences_parsimony(tree_root)

        actual_sequences = get_all_node_sequences(tree_root)
        self.assertEqual(actual_sequences, expected_sequences)

    def test_parse_newick_simple_leaf(self):
        node = parse_newick("LeafA;")
        self.assertEqual(node.name, "LeafA")
        self.assertIsNone(node.branch_info)
        self.assertEqual(len(node.children), 0)

    def test_parse_newick_leaf_with_branch_info(self):
        node = parse_newick("LeafB:[1A>T];")
        self.assertEqual(node.name, "LeafB")
        self.assertEqual(node.branch_info, "[1A>T]")
        self.assertEqual(len(node.children), 0)

    def test_parse_newick_internal_node_no_branch_info(self):
        node = parse_newick("(L1,L2)Anc;")
        self.assertEqual(node.name, "Anc")
        self.assertIsNone(node.branch_info)
        self.assertEqual(len(node.children), 2)
        self.assertEqual(node.children[0].name, "L1")
        self.assertEqual(node.children[1].name, "L2")

    def test_parse_newick_internal_node_with_branch_info(self):
        node = parse_newick("(L1,L2)Anc:[0C>G];")
        self.assertEqual(node.name, "Anc")
        self.assertEqual(node.branch_info, "[0C>G]")
        self.assertEqual(len(node.children), 2)

    def test_parse_newick_nested(self):
        node = parse_newick("((L1A,L1B)Anc1,(L2A,L2B)Anc2)Root;");
        self.assertEqual(node.name, "Root")
        self.assertEqual(len(node.children), 2)
        self.assertEqual(node.children[0].name, "Anc1")
        self.assertEqual(node.children[1].name, "Anc2")
        self.assertEqual(len(node.children[0].children), 2)
        self.assertEqual(node.children[0].children[0].name, "L1A")

    def test_apply_mutations_no_mutations(self):
        seq = list("ACGT")
        self.assertEqual(apply_mutations(seq, None), list("ACGT"))
        self.assertEqual(apply_mutations(seq, "[]"), list("ACGT"))
        self.assertEqual(apply_mutations(seq, ""), list("ACGT"))

    def test_apply_mutations_single(self):
        seq = list("ACGT")
        self.assertEqual(apply_mutations(seq, "[0A>T]"), list("TCGT"))

    def test_apply_mutations_multiple(self):
        seq = list("ACGT")
        self.assertEqual(apply_mutations(seq, "[0A>T,2G>C]"), list("TCCT"))

    def test_apply_mutations_out_of_bounds(self):
        # Should ignore out-of-bounds mutations and print warning (not tested here)
        seq = list("ACGT")
        self.assertEqual(apply_mutations(seq, "[5A>T,0A>G]"), list("GCGT"))

    def test_apply_mutations_no_ref_in_branch_info(self):
        seq = list("ACGT")
        self.assertEqual(apply_mutations(seq, "[0>T,2>C]"), list("TCCT"))


if __name__ == '__main__':
    unittest.main()
