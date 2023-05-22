#!/usr/bin/python

import os, sys, shutil, scipy.io, math, string, re
from img_pipe import img_pipe
import pandas as pd
import numpy as np
from glob import glob
from create_opscea_subj import OpsceaMaker
import nibabel as nib

# script requirements:
# - freesurfer dir under /Applications/freesurfer/subjects
# - nerdco dir under the freesurfer dir with
# --- SUBJ_SEEGext.mat
# --- SUBJ_CT_ESec.nii

class NerdcoOpsceaMaker(OpsceaMaker):
    def __init__(self):
        super().__init__()
        self.fs_nerdco_dir = os.path.join(self.freesurfer_subjdir, 'nerdco')
        self.lr_borderval = 130
        
    def do_imaging_elecs(self):
        elecmatrix = self.build_elecmatrix()
        print("====> Selected %d channels." % len(self.contiguous_labels))
        
        elecmatrix = self.do_mri_coordinate_transformation(elecmatrix)
            
        if not os.path.isdir(os.path.join(self.freesurfer_subjdir, "elecs")):
            os.makedirs(os.path.join(self.freesurfer_subjdir, "elecs"))
        output_path1 = os.path.join(self.freesurfer_subjdir, "elecs", "elecs_all.mat")
        if self.do_label:
            print("====> Getting FreeSurfer labels...")
            self.fslabels = self.patient.label_elecs(elecfile_prefix = "elecs_all", all_depth = True, quietmode = True)
            # self.process_labels()
        scipy.io.savemat(output_path1, {'eleclabels': self.make_cell_array(self.contiguous_labels), 'elecmatrix': elecmatrix})
        print("====> Extreme contacts:")
        self.print_contactdict(self.contactdict)
        output_path2 = os.path.join(self.opscea_subjdir, 'Imaging/Elecs', 'Electrodefile.mat')
        os.system('ln -sf %s %s' % (output_path1, output_path2))

    def build_elecmatrix(self):
        wirelist_fname = os.path.join(self.fs_nerdco_dir, 'WireEndpoints.csv')
        if not os.path.isfile(wirelist_fname):
            matlab_command = "/Applications/MATLAB_R2021a.app/bin/matlab -nodisplay -nodesktop -r \"load(\\\"%s\\\"); writetable(WireEndpoints, \\\"%s\\\"); quit;\"" % (os.path.join(fs_nerdco_dir, self.subjname + '_SEEGext.mat'), wirelist_fname)
            print(matlab_command)
            os.system(matlab_command)
        df = pd.read_csv(wirelist_fname)
        allwires = df['WIREnum']
    
        wirelist = list(set(allwires)) # list of unique wire nums

        points_dict = {}
        for w in wirelist:
            rowsel = df['WIREnum']==w
            points = df[rowsel]['PointID']
            ptx = df[rowsel]['PTX']
            pty = df[rowsel]['PTY']
            ptz = df[rowsel]['PTZ']

            for k in points.keys():
                try:
                    # note dimension swap
                    points_dict[string.ascii_uppercase[w-1]].append([points[k], pty[k], ptx[k], ptz[k]])
                except KeyError:
                    # note dimension swap
                    points_dict[string.ascii_uppercase[w-1]] = [[points[k], pty[k], ptx[k], ptz[k]]]

        self.sorted_dict = {}
        self.contactdict = {}
        contactcount = 0
    
        for w in points_dict.keys():
            sortdim = self.get_max_gradient(points_dict[w])
            maxgradient_list = []
            for p in points_dict[w]:
                maxgradient_list.append(p[sortdim])
            sortinds = np.argsort(maxgradient_list)
        
            for i,si in enumerate(sortinds):
                try:
                    self.sorted_dict[w].append([points_dict[w][si][1], points_dict[w][si][2], points_dict[w][si][3]])
                    contactcount += 1
                except KeyError:
                    self.sorted_dict[w] = [[points_dict[w][si][1], points_dict[w][si][2], points_dict[w][si][3]]]
                    contactcount += 1
                    first = contactcount
            self.contactdict[w] = [first, contactcount]

            # want low-order contacts to be the most mesial (need to flip if on the left)
            if sortdim==1 and self.sorted_dict[w][-1][0] < self.lr_borderval:
                self.sorted_dict[w].reverse()
            # elif sortdim==2 and self.sorted_dict[w][-1][1] < SOMEVAL:
            #     self.sorted_dict[w].reverse()
            # elif sortdim==3 and self.sorted_dict[w][-1][2] < SOMEVAL:
            #     self.sorted_dict[w].reverse()
            
        self.contiguous_labels = []
        for k in self.sorted_dict.keys():
            for i,p in enumerate(self.sorted_dict[k]):
                self.contiguous_labels.append(k + str(i+1))

        elecdict = {}
        for k in self.sorted_dict.keys():
            for i,p in enumerate(self.sorted_dict[k]):
                contactlabel = k + str(i+1)
                elecdict[contactlabel] = p

        for i,cl in enumerate(self.contiguous_labels):
            if i==0:
                elecmatrix = np.array([elecdict[cl]])
            else:
                elecmatrix = np.append(elecmatrix, [elecdict[cl]], axis=0)

        self.included_channels = range(len(self.contiguous_labels)) # assuming we won't be using actual EEG for this set

        return elecmatrix*1e-3

    def get_max_gradient(self, pointlist):
        # for each wire, need to identify dimension of greatest change
        # sort the contacts for each wire according to the dimension of
        # greatest change, because as it is they are not sorted!
        
        maxdx = -math.inf
        maxdy = -math.inf
        maxdz = -math.inf

        for i,p in enumerate(pointlist):
            if i==0:
                x0 = p[1]
                y0 = p[2]
                z0 = p[3]
            else:
                thisdx = abs(p[1] - x0)
                if thisdx > maxdx:
                    maxdx = thisdx

                thisdy = abs(p[2] - y0)
                if thisdy > maxdy:
                    maxdy = thisdy
                    
                thisdz = abs(p[3] - z0)
                if thisdz > maxdz:
                    maxdz = thisdz
                    
        if maxdx > maxdy and maxdx > maxdz:
            return 1 # dx is greatest
        elif maxdy > maxdx and maxdy > maxdz:
            return 2 # dy is greatest
        else:
            return 3 # dz is greatest

    def do_mri_coordinate_transformation(self, elecmatrix):
        # transform nifti coordinates to RAS
        
        scale_factor = 1000
        
        fs_mri_dir = os.path.join(self.freesurfer_subjdir, 'mri')
        fs_nerdco_dir = os.path.join(self.freesurfer_subjdir, 'nerdco')

        orig_file = os.path.join(fs_mri_dir, 'orig.mgz')
        voxsize_file = os.path.join(fs_nerdco_dir, 'voxsize_subj.txt')
        os.system('mri_info --cres %s > %s' % (orig_file, voxsize_file))
        # because we only use volumetric scans, only need 1 dimension
        voxsize = np.loadtxt(voxsize_file)
        
        # os.system('mri_info --vox2ras %s > %s' % (orig_file, affine_file))
        # vox2ras = np.loadtxt(affine_file)
        
        nifile = glob(os.path.join(fs_nerdco_dir, self.subjname + '_CT_*.nii*'))[0]
        myni = nib.load(nifile)
        vox2ras = myni.affine
        
        # see /Applications/brainstorm3/toolbox/anatomy/cs_convert.m
        # cs_convert l.83
        # Convert from MRI(mm) to voxels
        vox2ras = np.dot(vox2ras, np.diag(1./np.append([voxsize, voxsize, voxsize], 1)))
        # Convert to meters: transformation MRI=>WORLD (in meters)
        scalematrix = np.vstack((np.hstack((np.ones([3,3]), 1e-3*np.ones([3,1]))), np.ones([1,4])))
        mri2world = np.multiply(vox2ras, scalematrix)
        
        eleconesmatrix = np.hstack((elecmatrix, np.ones([elecmatrix.shape[0],1])))
        eleconesmatrix = np.transpose(np.dot(mri2world, np.transpose(eleconesmatrix)))
        return scale_factor*eleconesmatrix[:, 0:3]

    def process_labels(self):
        label_spreadsheet = os.path.join(self.fs_nerdco_dir, 'macroChanInfo.xlsx')
        labeldf = pd.read_excel(label_spreadsheet)
        # print(labeldf)

        self.letter_location_dict = {}
        # lateralize every letter
        for k in self.sorted_dict.keys():
            self.letter_location_dict[k] = [self.get_wire_hemisphere(self.sorted_dict[k])]

        for k in self.sorted_dict.keys():
            self.letter_location_dict[k].extend(self.get_wire_locations(k))

        self.print_contactdict(self.letter_location_dict)
        standard_targets = ["hippocampus", "amygdala", "insula", "cingulate"]

        self.sidedict = {}
        self.sidedict["L"] = "left"
        self.sidedict["R"] = "right"
        
        for hem in ["L", "R"]:
            rowsel = labeldf['hemi']==hem
            hem_longN = labeldf[rowsel]['longN']
            
            for k in hem_longN.keys():
                nameparts = hem_longN[k].split()
                struct = ""
                # for every standard structure in labeldf,
                if len(nameparts)==1 and nameparts[0] in standard_targets:
                    struct = nameparts[0]
                elif len(nameparts)==2 and nameparts[1] in standard_targets:
                    struct = nameparts[1]
                if len(struct)>0:
                    # make a list of lettered wires that could match
                    print(hem, struct, self.get_possible_letters(hem, struct))
        
        # for every ant/mod/post,

        # use coordinates to assign a single letter to each

    def get_wire_hemisphere(self, lis):
        coords_mat = np.array(lis)
        if np.mean(coords_mat[:,0]) < self.lr_borderval:
            return "L"
        else:
            return "R"

    def get_wire_locations(self, elecname):
        numcontacts = self.contactdict[elecname][1] - self.contactdict[elecname][0] + 1
        locations = []
        for i in range(numcontacts):
            thiscontactname = elecname + str(i+1)
            for contactlabel in self.fslabels:
                if contactlabel[0]==thiscontactname:
                    locations.append(contactlabel[3])
        return locations

    def get_possible_letters(self, hem, location):
        possible_letters = []
        for k in self.letter_location_dict.keys():
            thishem = self.letter_location_dict[k][0]
            if hem==thishem:
                thisloc = self.letter_location_dict[k][1:]
                for loc in thisloc:
                    locparts = list(map(lambda x: x.lower(), re.split('_|-', loc)))
                    if location in locparts:
                        possible_letters.append(k)
        return possible_letters
                            
    def do_seizure_dirs(self):
        badch = np.zeros((len(self.included_channels), 1))

        sfx = 256
        tim_sec = 30
        d = np.zeros((len(self.included_channels), sfx*tim_sec))

        szstring = self.subjname + "_01"
        szdir = os.path.join(self.opscea_subjdir, szstring)
        if not os.path.isdir(szdir):
            os.makedirs(szdir)

        szmatname = szstring + "_sz.mat"
        output_path = os.path.join(szdir, szmatname)
        scipy.io.savemat(output_path, {'d': d, 'sfx': sfx})
    
        badch_matname = szstring + "_badch.mat"
        badch_output_path = os.path.join(szdir, badch_matname)
        scipy.io.savemat(badch_output_path, {'badch': badch})
     
if __name__=="__main__":
    om = NerdcoOpsceaMaker()
    om.driver()

