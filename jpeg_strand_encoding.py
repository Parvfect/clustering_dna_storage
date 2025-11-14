
import numpy as np
from textwrap import wrap
import math
import random
import os
from itertools import cycle
from datetime import datetime
import json
from crc_encoding import get_crc_strand


bits_to_base = {'00': 'A', '01': 'C', '10': 'G', '11': 'T'}
strand_ids = [
        "CGTCTCGCGCCGGACCGAGC", "GGTAGGCCTGGCTAGTTGCT", "TTTGCGGCAGTGTCTGCGAT",
        "GGGCACAAATGGCTAGCCAC", "CGTCTTTGCCAAGGAGGTGT", "ACGACGCTGAGACAGAGATG"]

def add_xor_mask(bitstream, seed='01001'):
    return ''.join(
        '1' if b != mask else '0' for b, mask in zip(
            bitstream, cycle(seed)))

def index_to_bases(idx: int, length: int = 8) -> str:
    bits = format(idx, f'0{2*length}b')  # 2 bits per base
    return ''.join(bits_to_base[bits[i:i+2]] for i in range(0, len(bits), 2))


def encode_strands(
        dir, id_length, strand_length, fix_strand_ids=False,
        crc_encoding=False):
    """
    Encodes strands by adding XOR mask of 0101 and basic bit to dna conversion
    Splits to length such that its equally balanced
    """

    input_file = os.path.join(dir, "bird.jpg")
    fasta_file = os.path.join(dir, "bird_strands.fasta")
    cfg_file = os.path.join(dir, "cfg.json")

    #input_file = "bird.jpg"
    fasta_file = "bird_strands.fasta"
    cfg_file = "cfg.json"

    with open(input_file, "rb") as f:
        data = np.frombuffer(f.read(), dtype=np.uint8)

    bitstream = ''.join(format(b, '08b') for b in data)
    bitstream = add_xor_mask(bitstream)

    dna = ''.join(bits_to_base[bitstream[i:i+2]] for i in range(
        0, len(bitstream), 2))

    payload_length = strand_length
    total_bases = len(dna)
    n_strands = math.ceil(total_bases / strand_length)

    # Adjust if too long for the payload section
    bases_per_strand = min(n_strands, payload_length)
    chunks = wrap(dna, strand_length)
    
    if fix_strand_ids:
        strand_ids = [
            "CGTCTCGCGCCGGACCGAGC", "GGTAGGCCTGGCTAGTTGCT", "TTTGCGGCAGTGTCTGCGAT",
            "GGGCACAAATGGCTAGCCAC", "CGTCTTTGCCAAGGAGGTGT", "ACGACGCTGAGACAGAGATG"]
    else:
        strand_ids = [
        ''.join(random.choice(['A', 'C', 'T', 'G']) for _ in range(id_length))
        for _ in range(n_strands)
        ]


    # Add primers (truncate or pad as needed)
    strands = [
        f"{strand_ids[i]}{chunk}" for i, chunk in enumerate(
            chunks[:n_strands])]
    
    padding = np.zeros(n_strands)
    padded_strands = []
    strand_length = strand_length + id_length

    for ind, strand in enumerate(strands):
        padding_length = 0
        if len(strand) < strand_length: # Padding
            padding_length = strand_length - len(strand)
            strand += "A" * (padding_length)

        padded_strands.append(strand)
        padding[ind] = padding_length
    
    strands = padded_strands

    if crc_encoding:
        strands = [
            get_crc_strand(strand) for strand in strands]

    # ==== WRITE TO FASTA ====
    with open(fasta_file, "w") as f:
        for i, strand in enumerate(strands):
            f.write(f">strand_{strand_ids[i]}\n{strand}\n")

    cfg_dict = {
        'date': str(datetime.today()),
        'image_encoded': input_file,
        'fasta_encoded': fasta_file,
        'XOR seed': xor_seed,
        'strand_ids': strand_ids,
        'padding': list(padding)
    }

     # ==== WRITE config file ====
    with open(cfg_file, "w") as f:
        json.dump(cfg_dict, f)

    print(f"Encoded {len(strands)} strands written to {fasta_file}")
    return strands, strand_ids, total_bases, n_strands, padding


def save_partially_decoded_jpeg(
        recovered_strands, total_bases, n_strands,
        id_length, strand_ids, filename="decoded_partial.jpg",
        maintain_order=True, crc_encoding=True,
        padding=[0, 0, 0, 0, 0, 2]):

    # Remove primers and join payloads safely
    recovered_payloads = []
    strand_order = []
    for ind, s in enumerate(recovered_strands):
        if len(s) > id_length:
            recovered_strand = s[id_length:]
            recovered_id = s[:id_length]

            if crc_encoding:
                recovered_strand = recovered_strand[:-16]
            
            recovered_payloads.append(recovered_strand)

            if recovered_id in strand_ids:
                strand_order.append(strand_ids.index(recovered_id))
        else:
            print("Warning: strand too short, skipping:", s[:20], "…")

    # Sort recovered_strands by strand_order - reject if not strictly correct
    if sorted(strand_order) != list(np.arange(n_strands))[:len(recovered_payloads)]:
        print('Strands are in incorrect order! Decoding incomplete')
        return False
    
    if not recovered_payloads:
        raise ValueError("No valid strand payloads found")
    
    recovered_payloads = np.array(recovered_payloads)[
        np.argsort(strand_order)]
    
    if padding is not None:
        recovered_payloads_ = []

        for ind, i in enumerate(recovered_payloads):
            strand = i
            if padding[ind] > 0:
                strand = i[:-int(padding[ind])]
            recovered_payloads_.append(strand)

        recovered_payloads = recovered_payloads_

    recovered_payload = ''.join(recovered_payloads)
    #recovered_payload = recovered_payload[:total_bases]  # trim to max length if needed

    base_to_bits = {'A': '00', 'C': '01', 'G': '10', 'T': '11'}

    filtered_payload = ''.join(
        b for b in recovered_payload if b in base_to_bits)

    recovered_bits = ''.join(
        base_to_bits[b] for b in filtered_payload)
    
    recovered_bits = add_xor_mask(recovered_bits)

    bit_len = len(recovered_bits) - (
        len(recovered_bits) % 8)
    if bit_len == 0:
        raise ValueError("Recovered bits too short for even one byte")

    recovered_bits = recovered_bits[:bit_len]

    byte_data = int(recovered_bits, 2).to_bytes(bit_len // 8, byteorder='big')

    decoded_file = filename
    print(decoded_file)
    with open(decoded_file, "wb") as f:
        f.write(byte_data)

    return True  # Would like a more robust test here


# ==== PARAMETERS ====

id_length = 20
dir = os.path.join("data", "bird")

decoded_file = "data/bird/bird_decoded.jpg"
cfg_file = "data/bird/cfg.json"
xor_seed = '01001'

strand_length = 1093          # total bases per strand (including primer)
primer_prefix = "ACGTACGTACGT" # 20nt primer (easy to read/recognize)
crc_encoding=True

strands, strand_ids, total_bases, n_strands, padding = encode_strands(
    dir=dir, id_length=id_length,
    strand_length=strand_length, crc_encoding=crc_encoding)

for i in range(n_strands):
    save_partially_decoded_jpeg(
        recovered_strands=strands[:i+1], total_bases=total_bases,
        n_strands=n_strands, id_length=id_length,
        strand_ids=strand_ids, filename=f"decoded_partial_{i}.jpg",
        crc_encoding=crc_encoding, padding=padding)