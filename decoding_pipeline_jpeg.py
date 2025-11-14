
from utils import read_synthesized_strands_from_file
import Levenshtein
from clustering import Clustering
from utils import get_fastq_records
from strand_reconstruction import make_prediction
from tqdm import tqdm
from utils import reverse_complement
from crc_encoding import get_crc_strand
from jpeg_strand_encoding import save_partially_decoded_jpeg

def validate_crc(strand, info_length=1113):
    return get_crc_strand(strand[:info_length]) == strand

original_strands = read_synthesized_strands_from_file('bird_strands.fasta')[0]

records = get_fastq_records(r'birding.fastq')
strands = [str(i.seq) for i in records]
ids = [i.id for i in records]
strand_length = 1129

strands_length_filtered = [i for i in strands if len(i) > strand_length - 5 and len(i) < strand_length + 10]


clustering_obj = Clustering(strand_pool=strands_length_filtered, reference_length=1128, n_reference_strands=6, distance_threshold=100)

clustering_obj.run_pipeline(fix_orientation=True)


def generate_candidates_crc_validated(
        clustered_seqs, n_clusters=6, n_attempts=5,
        strand_length=200, ma_sample_size=10):
    
    validated_strands = []
    for ind, i in tqdm(enumerate(clustered_seqs[:n_clusters]), total=n_clusters):  # Iterate through clusters
        # Can be RC remember!
        for k in range(n_attempts): # Repeat n_attempts time
            candidate = make_prediction(i, sample_size=ma_sample_size)  # Make candidate prediction
            rev = reverse_complement(candidate)  # Obtain RC vector
            if validate_crc(candidate):  # Validate CRC code for forward and rc prediction
                validated_strands.append(candidate)
                break
            elif validate_crc(rev):
                validated_strands.append(rev)
                break
            else:
                continue
    
    validated_strands = list(set(validated_strands))
    print(f"{len(validated_strands)} valid strands found")

    return validated_strands

validated_strands = generate_candidates_crc_validated(clustering_obj.clustered_seqs)

no_crc_original_strands = [i[:-16] for i in original_strands]
print(len(set(validated_strands).intersection(
    original_strands)))