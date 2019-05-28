import pysam
import argparse

parser = argparse.ArgumentParser(description='make fastq from possorted_genome_bam.bam from cellranger')

parser.add_argument('-f', '--bam', required=True, help="cellranger bam")
parser.add_argument('-b', '--barcodes', required=True, help="cellranger barcodes.tsv")
parser.add_argument('-o', '--out', required=True, help="output fastq name")
args = parser.parse_args()

fn = args.bam#"possorted_genome_bam.bam"#files[0]
bam = pysam.AlignmentFile(fn, "rb")

cell_barcodes = set([])
with open(args.barcodes) as barcodes:
    for line in barcodes:
        tokens=line.strip().split()
        cell_barcodes.add(tokens[0])

recent_umis = {}
with open(args.out,'w') as fastq:
    for (index,read) in enumerate(bam):
        if not read.has_tag("CB"):
            continue
        cell_barcode = read.get_tag("CB")
        if read.is_secondary or read.is_supplementary:
            continue
        if not read.has_tag("UB"):
            continue
        UMI = read.get_tag("UB")
        pos = read.pos
        full_umi = cell_barcode + UMI + str(pos)
        if full_umi in recent_umis:
            continue
        recent_umis[full_umi] = pos
        if index % 10000 == 0:
            keys_to_remove = []
            for (key, val) in recent_umis.items():
                if val - pos > 2000:
                    keys_to_remove.append(key)
            for key in keys_to_remove:
                del recent_umis[key]
    
        readname = read.qname
        if read.has_tag("CB") and read.get_tag("CB") in cell_barcodes:
            fastq.write("@"+read.qname+";"+cell_barcode+";"+UMI+"\n")
            fastq.write(read.seq+"\n")
            fastq.write("+\n")
            fastq.write(read.qual+"\n")

