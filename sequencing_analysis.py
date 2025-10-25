
import argparse

parser = argparse.ArgumentParser(
                    prog='Sequencing data analysis',
                    description='Automated sequencing data analysis')

parser.add_argument('--fasta_filepath', type=str, help="Filepath with all the fasta files")
parser.add_argument('--forward_adapter', type=str)
parser.add_argument('--reverse_adapter', type=str)
parser.add_argument('--reference_filepath', type=str)

parser.set_defaults()

args = parser.parse_args()


if __name__ == '__main__':
    epochs = args.epochs



# Reference strands

# Remove primers

# Get reference markers using edit distance

# Align to reference

# Get IDS rates

# Junk reads

# Position of error
