#
# Merge performance per base to per read
#
#    python3 merge.py <Per_Base.tsv>
#

import sys

keys = ["ReadName", "runid", "read", "ch", "start_time", "flow_cell_id", \
        "protocol_group_id", "sample_id", "barcode", "Chrom", "Number",  \
        "REF_NT", "A", "T", "C", "G", "Ins", "Del", "Coverage", "Error"]

data = {}

# Read everything from script.py
with open(sys.argv[1]) as r:
    for line in r:
        toks = line.strip().split("\t")
        assert(len(toks) == 20)
        if "ReadName" in line:
            continue
        name = toks[0]
        if not name in data:
            data[name] = { "runid":[], "read":[], "ch":[], "start_time":[], "flow_cell_id":[], "protocol_group_id":[],
                           "sample_id":[], "barcode":[], "Chrom":[], "Number":[], "REF_NT":[],
                           "A":[], "T":[], "C":[], "G":[], "Ins":[], "Del":[], "Coverage":[], "Error":[] }
        data[name]["runid"].append(toks[1])
        data[name]["read"].append(toks[2])
        data[name]["ch"].append(toks[3])
        data[name]["start_time"].append(toks[4])
        data[name]["flow_cell_id"].append(toks[5])
        data[name]["protocol_group_id"].append(toks[6])
        data[name]["sample_id"].append(toks[7])
        data[name]["barcode"].append(toks[8])
        data[name]["Chrom"].append(toks[9])
        data[name]["Number"].append(toks[10])
        data[name]["REF_NT"].append(toks[11])
        data[name]["A"].append(toks[12])
        data[name]["T"].append(toks[13])
        data[name]["C"].append(toks[14])
        data[name]["G"].append(toks[15])
        data[name]["Ins"].append(toks[16])
        data[name]["Del"].append(toks[17])
        data[name]["Coverage"].append(toks[18])
        data[name]["Error"].append(toks[19])

assert(len(data) > 10)
print("ReadName; read; ch; start_time; flow_cell_id; protocol_group_id; sample_id; barcode; Chrom; Length; Ref; Mismatch; Ins; Del; Error; A; T; C; G")

for x in data:
    #
    # ReadName; read; ch; start_time; flow_cell_id; protocol_group_id; sample_id; barcode; Chrom; Length; Ref; Mismatch; Ins; Del; Error)
    #
    
    # How many we need to divide
    n = len(data[x]["read"])
    
    REF = 0
    Mismatch = 0
    Ins = 0
    Del = 0
    
    A_n = 0
    C_n = 0
    T_n = 0
    G_n = 0

    # Loop through all the nucleotides
    for i in range(n):
        A = float(data[x]["A"][i]) # Coverage
        T = float(data[x]["T"][i]) # Coverage
        C = float(data[x]["C"][i]) # Coverage
        G = float(data[x]["G"][i]) # Coverage
        R = data[x]["REF_NT"][i]   # Reference allele
        
        A_n += A
        C_n += C
        T_n += T
        G_n += G
        
        if R == "A":
            REF += A
            Mismatch += (T + C + G)
        elif R == "T":
            REF += T
            Mismatch += (A + C + G)
        elif R == "C":
            REF += C
            Mismatch += (A + T + G)
        elif R == "G":
            REF += G
            Mismatch += (A + T + C)
        else:
            raise Exception("????")
        
        Ins += float(data[x]["Ins"][i])
        Del += float(data[x]["Del"][i])        
        assert(float(Ins).is_integer())
        assert(float(Del).is_integer())

    Error = A_n + T_n + C_n + G_n + Ins + Del # Sum of everything
    Error -= REF
    assert(Error >= 0)
    
    print(x + ";" + data[x]["read"][0] + ";" + data[x]["ch"][0] + ";" + data[x]["start_time"][0] + ";" + \
          data[x]["flow_cell_id"][0] + ";" + data[x]["protocol_group_id"][0] + ";" + data[x]["sample_id"][0] + ";" + \
          data[x]["barcode"][0] + ";" + data[x]["Chrom"][0] + ";" + str(n) + ";" +
          str(int(REF)) + ";" + 
          str(int(Mismatch)) + ";" + 
          str(int(Ins)) + ";" + 
          str(int(Del)) + ";" +
          str(int(Error)) + ";" +
          str(int(A_n)) + ";" +
          str(int(T_n)) + ";" +
          str(int(C_n)) + ";" +
          str(int(G_n)))