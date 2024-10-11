from Bio import SeqIO
import gzip
from Bio.Seq import Seq
import os
import argparse

current_dir= os.getcwd()
parser = argparse.ArgumentParser()

parser.add_argument('-p','--path', type=str, default=current_dir, help="path to working dir")
parser.add_argument('-i', '--infile', type=str, help="Input fastq, can be gzipped")
parser.add_argument('-o', '--outfile', type=str, help="output file name, can output gz file with additional argument --gzd ")
parser.add_argument('-t', '--trimlen', type=int, default=60, help="the value supplied with with the length of nt cut at both sides of the reads")
parser.add_argument('-l', '--minlen', type=int, default=0, help="remove reads that are below the supplied length")
parser.add_argument('-g', '--gzd', action='store_true')
args = parser.parse_args()


os.chdir(args.path)

def trim_fastq(input_file, output_file, trim_length, min_length):
    with open(output_file, "w") as output_handle:
        for record in SeqIO.parse(input_file, "fastq"):
            if len(record.seq) > 2*trim_length:
                trimmed_seq = record.seq[trim_length:-trim_length]  # Trim 60 from each side
                ano= record.letter_annotations["phred_quality"][trim_length:-trim_length]
                record.annotations.clear()
                record.letter_annotations.clear()
                record.seq =trimmed_seq
                record.letter_annotations={"phred_quality":ano}
            if len(record.seq) > min_length:
                SeqIO.write(record, output_handle, "fastq")

input_fastq = args.infile
output_fastq = args.outfile
trim_length = args.trimlen
min_length= args.minlen


if input_fastq.endswith(".gz"):
    with gzip.open(input_fastq, "rt") as file:
        trim_fastq(file, output_fastq, trim_length, min_length)
else:
    trim_fastq(input_fastq, output_fastq, trim_length, min_length)

if args.gzd:
    with open(output_fastq, 'rb') as fq_file:
        with gzip.open(f"{output_fastq}.gz", 'wb') as comped_file:
            comped_file.writelines(fq_file)
    os.remove(output_fastq)

