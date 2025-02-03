import json
import pandas as pd
import sys
import os

### USAGE ###
#python 20250124_go_term_id_mapping_hj.py [*.json] [id mapping file name]

if len(sys.argv) != 3:
    print("please, check the arguments")
else:
    inputfile = str(sys.argv[1])
    outputfile = str(sys.argv[2])

with open(inputfile, 'r') as f:
    json_data = json.load(f)
GOTERM = []
GOID = []
term = list(json_data.keys())
for key, val in json_data.items():
    GOTERM.append(key)
    GOID.append(val["exactSource"])
    
term_id_mapping = pd.DataFrame({"GOTERM" : GOTERM, "GOID" : GOID})

term_id_mapping.to_csv(outputfile, index = False, mode = 'w', header = True, sep = "\t")
