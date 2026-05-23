#!/usr/bin/python

import os, sys, glob, re, math
import pandas as pd

# for electrodes.tsv, brainstorm exports columns: name, x, y, z, size (empty), group, type

# ofname = channeltsv.split(".")[0] + "_edited.tsv"
bsprefix = "/Users/aaron/Documents/brainstorm_db/IEEG_visualization/data/"

def driver(channeltsv):
    df = pd.read_csv(channeltsv, sep="\t")
    
    # sort rows
    def channel_key(ch): # thanks Claude!
        m = re.match(r'^([A-Za-z]+)(\d+)$', ch)
        return (m.group(1), int(m.group(2))) if m else (ch, 0)

    df.sort_values(by=["name"], key=lambda x: x.map(channel_key), inplace=True)

    max_contact_dict = {}
    stem = ''
    num = ''
    for i,c in enumerate(df["name"]):
        thisstem = re.match(r'^([A-Za-z]+)', c).group()
        if thisstem!=stem: # it's a new stem
            if i>0:
                max_contact_dict[stem] = num
            stem = thisstem
        else:
            num = int(re.findall(r'(\d+)$', c)[0])

    # drop rows starting with XXX
    df = df[~df["name"].str.startswith("XXX")]

    # add size value- for now, assuming only DIXI D08 microdeep
    contact_length = 2 # mm
    electrode_diameter = 0.8 # mm
    surface_area = 2*math.pi*electrode_diameter*0.5*contact_length
    dixi_contact_nums = [8, 10, 12, 15, 18]

    for k in max_contact_dict.keys():
        if max_contact_dict[k] not in dixi_contact_nums:
            df.loc[df["group"]==k, "size"] = "n/a"
        else:
            df.loc[df["group"]==k, "size"] = str(f"{surface_area:.3f}")

    df.to_csv(channeltsv, sep="\t", index=False)

if __name__=="__main__":
    if len(sys.argv)<2:
        subjlist = ["UCHAK240403", "UCHAM250108", "UCHDR220801",
                    "UCHDR240313", "UCHGG230823", "UCHJR250122", "UCHSM240205",
                    "UCHSN230406", "UCHTD250331", "UCHVG230719"]
        hasshort = ['UCHSN230406', 'UCHGG230823', 'UCHVG230719', 'UCHDR220801', 'UCHAK240403']

        for s in subjlist:
            if s in hasshort:
                s = s[:5]
            fpath = os.path.join(bsprefix, s, s + "_electrodes.tsv")
            print("Editing " + fpath)
            driver(fpath)

    else:
        s = sys.argv[1]
        fpath = os.path.join(bsprefix, s, s + "_electrodes.tsv")
        print("Editing " + fpath)
        driver(fpath)
