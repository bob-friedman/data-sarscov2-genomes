"""
Parses a Newick formatted tree file with annotations on the nodes for
mutational changes, along with a reference sequence at the root node.
Outputs the biological sequences at the internal nodes and leaves by a
parsimony method.

Input:
  A text file where:
  - The first line is a Newick tree string with mutations annotated
    on branches. E.g., ((Leaf1:[1C>A,3T>A],Leaf2:[1C>T])Anc1)Root;
  - The second line is the reference sequence string. E.g., CGTGA

Output:
  Prints the inferred sequences for all nodes in the tree to standard output.
"""
import re
import sys
import argparse
from typing import List, Dict, Set, Optional

# --- IUPAC Data Structures ---
IUPAC_EXPANSION: Dict[str, Set[str]] = {
    'A': {'A'}, 'C': {'C'}, 'G': {'G'}, 'T': {'T'}, 'R': {'A', 'G'},
    'Y': {'C', 'T'}, 'M': {'A', 'C'}, 'K': {'G', 'T'}, 'W': {'A', 'T'},
    'S': {'C', 'G'}, 'H': {'A', 'C', 'T'}, 'D': {'A', 'G', 'T'},
    'V': {'A', 'C', 'G'}, 'B': {'C', 'G', 'T'}, 'N': {'A', 'C', 'G', 'T'}
}
IUPAC_MAP: Dict[str, str] = {
    'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T', 'AG': 'R', 'CT': 'Y',
    'AC': 'M', 'GT': 'K', 'AT': 'W', 'CG': 'S', 'ACT': 'H',
    'AGT': 'D', 'ACG': 'V', 'CGT': 'B', 'ACGT': 'N'
}

class PhyloNode:
    """Represents a node in a phylogenetic tree."""
    def __init__(self, name: str, branch_info: Optional[str] = None):
        self.name: str = name
        self.branch_info: Optional[str] = branch_info
        self.children: List[PhyloNode] = []
        self.inferred_sequence: List[str] = []

    def __repr__(self) -> str:
        return f"PhyloNode(name='{self.name}')"

def _split_top_level(s: str, delimiter: str) -> List[str]:
    """
    Splits a string by a delimiter, but only at the top level
    (not within parentheses or brackets).
    """
    parts: List[str] = []
    buf: List[str] = []
    p_depth: int = 0  # Parentheses depth
    b_depth: int = 0  # Brackets depth
    for ch in s:
        if ch == '(': p_depth += 1
        elif ch == ')': p_depth -= 1
        elif ch == '[': b_depth += 1
        elif ch == ']': b_depth -= 1

        if ch == delimiter and p_depth == 0 and b_depth == 0:
            parts.append(''.join(buf))
            buf = []
        else:
            buf.append(ch)
    parts.append(''.join(buf))
    return [p.strip() for p in parts if p.strip()]

def parse_newick(s: str) -> PhyloNode:
    """
    Parses a Newick string (potentially with branch annotations) into a PhyloNode tree.
    Example: ((Leaf1:[1C>A],Leaf2)Anc1)Root;
    """
    s = s.strip().rstrip(';')

    def _recursive_parse(substring: str) -> PhyloNode:
        substring = substring.strip()
        # Base case: Leaf node or node without explicit children string
        if not substring.startswith('('):
            parts = substring.split(':', 1)
            node_name = parts[0]
            branch_info = parts[1] if len(parts) > 1 else None
            return PhyloNode(name=node_name, branch_info=branch_info)

        # Recursive case: Internal node
        # Find the matching parenthesis for the children block
        depth = 0
        match_idx = -1
        for i, char in enumerate(substring):
            if char == '(':
                depth += 1
            elif char == ')':
                depth -= 1
                if depth == 0:
                    match_idx = i
                    break

        if match_idx == -1:
            raise ValueError(f"Invalid Newick substring: Unmatched parentheses in '{substring}'")

        children_str = substring[1:match_idx]
        parent_info_str = substring[match_idx+1:]

        parent_parts = parent_info_str.split(':', 1)
        parent_name = parent_parts[0]
        parent_branch_info = parent_parts[1] if len(parent_parts) > 1 else None

        node = PhyloNode(name=parent_name, branch_info=parent_branch_info)

        child_substrings = _split_top_level(children_str, ',')
        for child_substring in child_substrings:
            node.children.append(_recursive_parse(child_substring))

        return node

    return _recursive_parse(s)


def apply_mutations(base_seq: List[str], branch_info: Optional[str]) -> List[str]:
    """
    Applies mutations from branch_info to a base sequence.
    Mutations are expected in a format like "[1C>A,3T>G]", where '1' is the
    0-indexed position, 'C' is the original nucleotide (optional in parsing but good for clarity),
    and 'A' is the new nucleotide.
    """
    new_seq = list(base_seq) # Make a copy
    if not branch_info:
        return new_seq

    # Regex finds all occurrences of "POS_REF_NUC > NEW_NUC" or "POS > NEW_NUC"
    # e.g., "1C>A" or "0>T". Captures position and new nucleotide.
    # The reference nucleotide in the middle is descriptive but not used by this function.
    mutation_pattern = re.compile(r'(\d+)[A-Z]?>([A-Z])')
    mutations_found = mutation_pattern.findall(branch_info)

    for pos_str, new_nucleotide in mutations_found:
        pos = int(pos_str)
        if 0 <= pos < len(new_seq):
            new_seq[pos] = new_nucleotide
        else:
            print(f"Warning: Mutation '{pos_str}>{new_nucleotide}' refers to position {pos}, "
                  f"which is out of bounds for sequence of length {len(new_seq)}. Mutation skipped.", file=sys.stderr)
    return new_seq

def calculate_initial_sequences(node: PhyloNode, parent_seq: List[str]):
    """
    Recursively calculates the initial sequence for each node by applying mutations
    down the tree from the parent. This is the first pass, typically from root to tips.
    """
    node.inferred_sequence = apply_mutations(parent_seq, node.branch_info)
    for child in node.children:
        calculate_initial_sequences(child, node.inferred_sequence)

def reconstruct_ancestral_sequences_parsimony(node: PhyloNode):
    """
    Recursively reconstructs ancestral sequences using a parsimony approach (Fitch's algorithm).
    This is the second pass, typically from tips up to the root.
    Assumes initial sequences (possibly from tips or a previous pass) are present.
    """
    if not node.children:  # Leaf node
        # For leaf nodes, the sequence is already determined (or given).
        # If they were ambiguous, this is where they'd be resolved from IUPAC_EXPANSION.
        # For this script, apply_mutations already set them to specific bases.
        return

    # Post-order traversal: process children first
    for child in node.children:
        reconstruct_ancestral_sequences_parsimony(child)

    # Determine sequence for the current internal node based on its children
    if not node.children: # Should not happen due to the check above, but defensive
        return

    # All children must have sequences of the same length
    first_child_seq_len = len(node.children[0].inferred_sequence)
    if not all(len(child.inferred_sequence) == first_child_seq_len for child in node.children):
        raise ValueError(f"Child sequences of node {node.name} have inconsistent lengths.")

    if first_child_seq_len == 0: # Can happen if root sequence was empty
        node.inferred_sequence = []
        return

    new_parent_seq: List[str] = []
    for i in range(first_child_seq_len):  # Iterate over each site in the alignment

        def get_expanded_bases(base_char: str, node_name: str, site_idx: int) -> Set[str]:
            """Safely get expanded bases, defaulting to 'N' for unknowns."""
            if base_char not in IUPAC_EXPANSION:
                print(f"Warning: Unknown base '{base_char}' in sequence for node {node_name} "
                      f"at site {site_idx}. Treating as 'N'.", file=sys.stderr)
                return set(IUPAC_EXPANSION['N'])
            return set(IUPAC_EXPANSION[base_char])

        # Initialize intersection with bases from the first child
        intersected_bases = get_expanded_bases(
            node.children[0].inferred_sequence[i],
            node.children[0].name,
            i
        )

        # Intersect with bases from other children
        for child_idx in range(1, len(node.children)):
            child_node = node.children[child_idx]
            current_child_bases = get_expanded_bases(
                child_node.inferred_sequence[i],
                child_node.name,
                i
            )
            intersected_bases.intersection_update(current_child_bases)

        # If intersection is empty, take the union
        if not intersected_bases:
            # This union logic is part of Fitch's algorithm for parsimony
            final_bases: Set[str] = set()
            for child_node in node.children:
                # Reuse get_expanded_bases for consistent handling of unknowns
                child_bases_for_union = get_expanded_bases(
                    child_node.inferred_sequence[i],
                    child_node.name,
                    i
                )
                final_bases.update(child_bases_for_union)
        else:
            final_bases = intersected_bases

        # Convert set of bases back to IUPAC code
        # Sort for consistent key generation for IUPAC_MAP
        sorted_bases_key = "".join(sorted(list(final_bases)))
        # Default to 'N' if the specific combination of bases isn't in IUPAC_MAP
        # (e.g., if final_bases was empty, though logic should prevent that for valid inputs)
        new_parent_seq.append(IUPAC_MAP.get(sorted_bases_key, 'N'))

    node.inferred_sequence = new_parent_seq


def print_final_tree_sequences(node: PhyloNode, depth: int = 0):
    """Recursively prints the name and inferred sequence of each node in the tree."""
    indent = "  " * depth
    node_type = "Internal" if node.children else "Leaf"
    sequence_str = ''.join(node.inferred_sequence) if node.inferred_sequence else "N/A"
    print(f"{indent}- {node_type} Node: {node.name}")
    print(f"{indent}  - Final Sequence: {sequence_str}")
    for child in node.children:
        print_final_tree_sequences(child, depth + 1)

def main():
    """
    Main function to parse arguments, read input file, run parsimony reconstruction,
    and print results.
    """
    parser = argparse.ArgumentParser(
        description="Reconstructs ancestral sequences in a phylogenetic tree using parsimony."
    )
    parser.add_argument(
        "input_file_path",
        type=str,
        help="Path to the input file (first line: Newick string, second line: root sequence)."
    )
    args = parser.parse_args()

    try:
        with open(args.input_file_path, 'r') as f:
            newick_string = f.readline().strip()
            root_sequence_str = f.readline().strip().upper() # Ensure uppercase
            if not newick_string :
                raise ValueError("Input file is missing the Newick string on the first line.")
            if not root_sequence_str:
                 raise ValueError("Input file is missing the root sequence on the second line.")
            if not re.match(r"^[ACGTNRYMKSWBDHV]+$", root_sequence_str):
                raise ValueError(f"Invalid characters in root sequence: {root_sequence_str}. Only ACGTNRYMKSWBDHV allowed.")

    except FileNotFoundError:
        print(f"Error: The file '{args.input_file_path}' was not found.", file=sys.stderr)
        sys.exit(1)
    except ValueError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

    print(f"Reading from: {args.input_file_path}")
    print(f"Using Newick string: {newick_string}")
    print(f"Using known root sequence: {root_sequence_str}\n")

    try:
        tree_root = parse_newick(newick_string)
    except ValueError as e:
        print(f"Error parsing Newick string: {e}", file=sys.stderr)
        sys.exit(1)

    # First pass: Calculate initial sequences from root down to tips by applying mutations
    calculate_initial_sequences(tree_root, list(root_sequence_str))

    # Second pass: Reconstruct ancestral sequences using parsimony from tips up to root
    # This will refine internal node sequences based on their children.
    # The root's sequence will also be re-evaluated based on its children after this step,
    # effectively using the provided root sequence as a strong prior for the first pass.
    reconstruct_ancestral_sequences_parsimony(tree_root)

    print("--- Final Inferred Sequences for All Nodes (with Parsimony) ---")
    print_final_tree_sequences(tree_root)

if __name__ == '__main__':
    main()
