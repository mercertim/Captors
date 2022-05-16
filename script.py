#
# python3 script.py <FASTA FILE> <FASTQ FILE> <BED FILE> <BAM FILE>
#
#    bedtools intersect -a bc10long.BRCAcon.sort.bam -b BA_BA.bed > BA_BA_bc10long.BRCAcon.sort.bam
#    samtools index BA_BA_bc10long.BRCAcon.sort.bam
#    /home/linuxbrew/.linuxbrew/bin/python3.7 script.py BRCAconjoined.fa barcode10.long.fastq BA_BA.bed BA_BA_bc10long.BRCAcon.sort.bam > perBase.tsv
#    /home/linuxbrew/.linuxbrew/bin/python3.7 merge.py perBase.tsv > perRead.tsv
#

import os
import sys
import pysam
from Bio import SeqIO

faA  = sys.argv[1]
fqA  = sys.argv[2]
bedA = sys.argv[3]
bamA = sys.argv[4]

bed = {}
with open(bedA, "r") as r:
    for line in r:
        if line == "\n":
            continue
        toks = line.strip().split("\t")
        assert(len(toks) >= 3)
        bed[toks[0]] = { "start":int(toks[1]), "end":int(toks[2]) }
assert(len(bed) == 1)
#assert(len(bed) == 1 and bed["BA_BA"] is not None) # Only BA_BA...

# Eg: read=8415 ch=77 start_time=2020-06-01T04:55:27Z flow_cell_id=ACR542 protocol_group_id=FLFL073983 sample_id=FLFL073983 barcode=barcode10
fastq = {}
for r in SeqIO.parse(fqA, "fastq"):
    toks = r.description.split(" ")
    for tok in toks:
        if not "=" in tok:
            continue
        name = str(r.name)
        if not name in fastq:
            fastq[name] = {}
        fastq[name][tok.split("=")[0]] = tok.split("=")[1]
assert(len(fastq) > 0)

bam = pysam.AlignmentFile(bamA, "rb")

# Create a BAM with only the read
def createBAM(r):
    name = "/tmp/A.bam"
    tmp = pysam.AlignmentFile(name, "wb", template=bam)
    tmp.write(r)
    tmp.close()
    os.system("samtools index " + name)
    return name

def pysamstats(bam):
    name = "/tmp/A.bed"
    os.system("pysamstats --fasta " + faA + " --type variation " + bam + " > " + name)
    return name

def bedtools(name):
    out = "/tmp/B.bed"
    w = open(out, "w")
    n = 0 # Number of matching lines
    with open(name) as r:
        for l in r:
            l = l.strip()
            toks = l.split("\t")
            if toks[0] == "chrom":
                w.write(str(l) + "\n")
                continue
            elif not toks[0] in bed:
                continue

            n1  = bed[toks[0]]["start"]
            n2  = bed[toks[0]]["end"]            
            pos = int(toks[1])

            # Only if Tim explicily asked for it
            if pos >= n1 and pos <= n2:
                n += 1
                w.write(str(l) + "\n")
    w.close()
    return (n, out)

os.system("rm -rf data")
os.system("mkdir -p data")

# For each read...
for r in bam.fetch():
    # Create a BAM file with only the read so anything calculated will also be per-read based
    name = createBAM(r)
    
    # Running statistical computation
    name = pysamstats(name)
    
    # Apply intersection
    (n, name) = bedtools(name)
    
    if n == 0:
        continue
    assert(n > 0)
    
    # Per-base performance for this read only
    os.system("python3 analyzePile.py " + name + " > /tmp/analyzePile.txt")
    
    with open("/tmp/analyzePile.txt") as x:
        assert(r.qname in fastq)
        i = fastq[r.qname]
        keys = ["runid", "read", "ch", "start_time", "flow_cell_id", "protocol_group_id", "sample_id", "barcode"]
        
        print("ReadName\trunid\tread\tch\tstart_time\tflow_cell_id\tprotocol_group_id\tsample_id\tbarcode\tChrom\tNumber\tREF_NT\tA\tT\tC\tG\tIns\tDel\tCoverage\tError")
        for l in x:
            toks = l.strip().split("\t")
            if toks[0] == "Chrom":
                continue
            assert(len(toks) == 11)

            # Single line on all information
            print(r.qname + "\t" + str(i["runid"]) + "\t" + str(i["read"]) + "\t" + i["ch"] + "\t" + i["start_time"] + "\t" +
                  i["flow_cell_id"] + "\t" + i["protocol_group_id"] + "\t" + i["sample_id"] + "\t" + i["barcode"] + "\t", end='')
            
            # Randomly pick insertion to check
            Ins = toks[7]

            # We are always reporting counts now
            assert(float(Ins).is_integer())
            
            # Single line on all information
            print(toks[0]  + "\t" + \
                  toks[1]  + "\t" + \
                  toks[2]  + "\t" + \
                  toks[3]  + "\t" + \
                  toks[4]  + "\t" + \
                  toks[5]  + "\t" + \
                  toks[6]  + "\t" + \
                  toks[7]  + "\t" + \
                  toks[8]  + "\t" + \
                  toks[9]  + "\t" + \
                  toks[10])
    
    # Generate a data directory for each indivudal read (one for each read)
    os.system("mv /tmp/analyzePile.txt data/" + r.qname + ".txt")
