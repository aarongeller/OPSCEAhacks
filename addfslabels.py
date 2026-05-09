#!/usr/bin/python

import os, sys, glob, re
import pandas as pd

# add freesurfer labels to the TSV produced by exporting channel file from brainstorm.

# Usage: python addlabels.py SUBJ SUBJ.tsv

if len(sys.argv)<3:
    print("Usage: python addlabels.py SUBJ SUBJ.tsv")
    sys.exit()
    
subj = sys.argv[1]
channeltsv = sys.argv[2]
ofname = channeltsv.split(".")[0] + "fs.tsv"

opsceapath = "/Users/aaron/Documents/MATLAB/OPSCEA-main/OPSCEADATA"
opsceasubjpath = os.path.join(opsceapath, subj)
opsceasubjlabelspath = os.path.join(opsceasubjpath, "Imaging/Recon/labels")

labelfiles = glob.glob(os.path.join(opsceasubjlabelspath, '*.tex'))

df = pd.read_csv(channeltsv, sep="\t")

newlead = 1
fslabels = []

def channel_key(ch): # thanks Claude!
    m = re.match(r'^([A-Za-z]+)(\d+)$', ch)
    return (m.group(1), int(m.group(2))) if m else (ch, 0)

df.sort_values(by=["Channel"], key=lambda x: x.map(channel_key), inplace=True)

for contact in df["Channel"]:
    if contact[-1]=="1" and not contact[-2].isdigit():
        # it's a new lead
        thislead = contact[:-1]

        # get the corresponding tex file
        texf = open(os.path.join(opsceasubjlabelspath, thislead + "_labels.tex"))

        texlines = texf.readlines()
        
        for ll in texlines:
            llparts = ll.split("&")
            fslabels.append(llparts[1].replace("\\", "").strip())
    
df["Freesurfer Labels"] = fslabels

df.to_csv(ofname, sep="\t")
