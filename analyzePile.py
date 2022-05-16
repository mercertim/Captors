#
# python3 analyzePile.py <FROM pysamstats> <Optional BED File>
#

import sys

x = {} # Indexed by position
with open(sys.argv[1]) as r:
    for line in r:
        if "reads_all" in line:
            continue
        toks = line.strip().split("\t")
        assert(len(toks) == 23)
        
        chr = toks[0]      # Sequence name
        pos = int(toks[1]) # 1-based coordinate
        ref = toks[2]      # Reference

        A = int(toks[13]) + int(toks[14])
        T = int(toks[17]) + int(toks[18])
        C = int(toks[15]) + int(toks[16])
        G = int(toks[19]) + int(toks[20])
        D = int(toks[8]) + int(toks[9])
        I = int(toks[10]) + int(toks[11])
        dp = A + T + C + G
        
        if not chr in x:
            x[chr] = {}
        if ref == "N":
            continue
        assert(ref == 'A' or ref == 'C' or ref == 'T' or ref == 'G')
        x[chr][pos] = { "R":ref, "A":A, "C":C, "T":T, "G":G, "N":0, "I":I, "D":D, "DP":dp }

bed = {}
if len(sys.argv) > 2:
    with open(sys.argv[2]) as r:
        for line in r:
            toks = line.strip().split("\t")
            bed[toks[0]] = { "start":int(toks[1]), "end":int(toks[2]) }
    assert(len(bed) > 0)

print("Chrom\tNumber\tREF_NT\tA\tT\tC\tG\tIns\tDel\tCoverage\tError")
for chr in x:
    for pos in x[chr]:
        if chr in bed:
            if pos < bed[chr]["start"] or pos > bed[chr]["end"]:
                continue
        
        A = x[chr][pos]["A"]
        T = x[chr][pos]["T"]
        C = x[chr][pos]["C"]
        G = x[chr][pos]["G"]
        I = x[chr][pos]["I"]
        D = x[chr][pos]["D"]
        
        assert(float(A).is_integer())
        assert(float(T).is_integer())
        assert(float(C).is_integer())
        assert(float(G).is_integer())
        assert(float(I).is_integer())
        assert(float(D).is_integer())

        dp = A + T + C + G        
        total = A + T + C + G + I + D
    
        data = {}
        data["A"] = A #(1.0 * A) / total # Relative fraction
        data["T"] = T #(1.0 * T) / total
        data["C"] = C #(1.0 * C) / total
        data["G"] = G #(1.0 * G) / total
        data["I"] = I #(1.0 * I) / total
        data["D"] = D #(1.0 * D) / total
        
        error = data["A"] + data["T"] + data["C"] + data["G"] + data["I"] + data["D"] # Sum of everything
        error -= data[x[chr][pos]["R"]]        
        assert(float(error).is_integer())
        
        print(chr + "\t" + str(pos) + "\t" +
              x[chr][pos]["R"] + "\t" +
              str(data["A"]) + "\t" +
              str(data["T"]) + "\t" +
              str(data["C"]) + "\t" +
              str(data["G"]) + "\t" +
              str(data["I"]) + "\t" +
              str(data["D"]) + "\t" +
              str(dp) + "\t" + str(error))
