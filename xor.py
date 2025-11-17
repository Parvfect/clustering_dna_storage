
from crc_encoding import \
    convert_dna_to_binary_string, convert_binary_string_to_dna
from itertools import cycle
import random
from utils import read_synthesized_strands_from_file
from tqdm import tqdm

strands = read_synthesized_strands_from_file(
    'data/final_run/bird_strands.fasta')[0]

def check_obeys_dna_constraints(
        strand, max_hp_length=4, max_gc_content=0.6):
    """
    Checks if homopolymers are less than 3 and 
    GC content is less than 45 %
    """

    run = 1
    gc = 0
    prev_nt = ""
    violations = 0
    for nt in strand:
        if nt == prev_nt:
            run += 1
        else:
            prev_nt = nt
            run = 1

        if nt == 'G' or nt == 'C':
            gc += 1
        
        if run >= max_hp_length + 1:
            return False
        
        
        if run >= max_hp_length:
            violations += 1
    
    gc_content = gc / len(strand)

    if gc_content > max_gc_content:
        return False
    
    print(f"GC content is {gc_content}")
    print(f"{violations} Violations")
    
    return True

def add_xor_mask(bitstream, seed='01001'):
    return ''.join(
        '1' if b != mask else '0' for b, mask in zip(
            bitstream, cycle(seed)))

def generate_random_xor_seed(seed_length):
    return "".join([
        str(random.randint(0, 1)) for i in range(seed_length)])

def get_valid_xor_seed_and_strand(strand, xor_seed_length=20):
    bin_strand = convert_dna_to_binary_string(strand)
    xor_seed_dna = ""
    final_strand = ""
    while True:
        xor_seed = generate_random_xor_seed(xor_seed_length)
        xored_bin_string = add_xor_mask(bin_strand, xor_seed)
        new_strand = convert_binary_string_to_dna(xored_bin_string)
        if check_obeys_dna_constraints(new_strand):
            final_strand = new_strand
            xor_seed_dna = convert_binary_string_to_dna(xor_seed)
            break

    return final_strand, xor_seed_dna



for i in strands:
    print(get_valid_xor_seed_and_strand(i))











