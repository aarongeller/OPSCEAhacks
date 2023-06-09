#!/usr/bin/python

# The script assumes recon-all has been done for SUBJ and electrodes
# have been localized, with EEG saved as edf files under SUBJ/eeg.

import os, sys, shutil, scipy.io, mne
import numpy as np
import numpy.matlib
import pandas as pd
from glob import glob
from abc import ABC, abstractmethod

class OpsceaMaker(ABC):
    def __init__(self):
        self.handle_args()
        os.environ['SUBJECTS_DIR'] = os.path.join(os.environ['FREESURFER_HOME'], "subjects")
        self.freesurfer_subjdir = os.path.join(os.environ['SUBJECTS_DIR'], self.subjname)
        self.fs_eeg_dir = os.path.join(self.freesurfer_subjdir, 'eeg')
        self.all_eeg_files = glob(os.path.join(self.fs_eeg_dir, '*.edf'))
        self.all_eeg_files.sort()

        path_to_opsceadata = os.path.join(os.environ['HOME'], "Documents/MATLAB/OPSCEA-main/OPSCEADATA")
        self.opscea_subjdir = os.path.join(path_to_opsceadata, self.subjname)

        # path_to_imgpipe = "/opt/local/Library/Frameworks/Python.framework/Versions/3.8/lib/python3.8/site-packages/img_pipe"

    def handle_args(self):
        if len(sys.argv) < 2:
            print("Usage: python3.8 create_opscea_subj.py SUBJ <rh/lh/stereo> <0/1> <0/1> <lowpass> <number-elecs>")
            print("\twhere SUBJ is a directory in /Applications/freesurfer/subjects")
            print("\toptional args: hemisphere, do_subcort, do_label, numberlabels, low pass freq, and list of electrodes ending with a digit.")
            quit()
            
        self.subjname = sys.argv[1]

        if len(sys.argv)>2:
            self.hem = sys.argv[2]
        else:
            self.hem = "stereo"

        if len(sys.argv)>3:
            self.do_subcort = int(sys.argv[3])
        else:
            self.do_subcort = 1

        if len(sys.argv)>4:
            self.do_label = int(sys.argv[4])
        else:
            self.do_label = 1

        if len(sys.argv)>5:
            self.lowpass = int(sys.argv[5])
        else:
            self.lowpass = -1

        if len(sys.argv)>6:
            self.numberlabels = sys.argv[6:]
        else:
            self.numberlabels = []

    def driver(self):
        self.do_imaging_dirs()
        self.do_seizure_dirs()
        print("Done.")

    def do_imaging_dirs(self):
        chosen_structures = ['Hipp', 'Amgd']
        
        targets = ['Imaging/Elecs',
                   'Imaging/Meshes',
                   'Imaging/Meshes/Subcortical',
                   'Imaging/MRI',
                   'Imaging/Recon/figs',
                   'Imaging/Recon/labels']
        
        mritargets = ['aparc+aseg', 'brain']

        fs_mesh_dir = os.path.join(self.freesurfer_subjdir, 'Meshes')
        fs_surf_dir = os.path.join(self.freesurfer_subjdir, 'surf')
    
        from AG_img_pipe import freeCoG # wait to import until after environment variables are set
        self.patient = freeCoG(subj = self.subjname, hem = self.hem)
        self.patient.convert_fsmesh2mlab()
        if self.do_subcort:
            self.patient.get_subcort()

        for t in targets:
            if not os.path.isdir(os.path.join(self.opscea_subjdir, t)):
                os.makedirs(os.path.join(self.opscea_subjdir, t))

            if t=='Imaging/Elecs':
                self.do_imaging_elecs() # needs to be overridden
                  
            elif t=='Imaging/Meshes':
                # copy required meshes to opscea
                mesh_name = 'pial'
                for h in ['lh', 'rh']:
                    out_file_struct = os.path.join(fs_mesh_dir, '%s_%s_%s.mat'%(self.subjname, h, mesh_name))
                    shutil.copy2(out_file_struct, os.path.join(self.opscea_subjdir, t))
        
            elif t=='Imaging/Meshes/Subcortical':
                subcort_dir = os.path.join(fs_mesh_dir,'subcortical')     
                for side in ['l','r']:
                    for s in chosen_structures:
                        structname = os.path.join(subcort_dir, side + s + '_subcort.mat')
                        shutil.copy2(structname, os.path.join(self.opscea_subjdir, t))
                
            elif t=='Imaging/MRI':
                for mt in mritargets:
                    shutil.copy2(os.path.join(self.freesurfer_subjdir, 'mri', mt + '.mgz'), 
                                 os.path.join(self.opscea_subjdir, t))

            elif t=='Imaging/Recon/figs':
                self.do_imaging_recon_figs()
            
            else: # t=='Imaging/Recon/labels'
                self.do_imaging_recon_labels()

    @abstractmethod
    def do_imaging_elecs(self):
        pass
    
    def do_imaging_recon_figs(self):
        pass
    
    def do_imaging_recon_labels(self):
        if self.do_label:
            self.do_label_output(os.path.join(self.opscea_subjdir, 'Imaging/Recon/labels'))
        
    def do_seizure_dirs(self):
        # process bad channel list if it's there
        fs_badch_path = os.path.join(self.fs_eeg_dir, 'badch.xls')
        try:
            df = pd.read_excel(fs_badch_path, header=None)
            bad_chans_list = list(df[0])
            print("Processed bad channels from " + fs_badch_path + ".") 
        except FileNotFoundError:
            bad_chans_list = []

        badch = []
        numbad = 0
        for cl in self.contiguous_labels:
            if cl in bad_chans_list:
                badch.append(1)
                numbad += 1
            else:
                badch.append(0)

        print(str(numbad) + " channels marked bad.")    
        badch = np.transpose([np.array(badch)])

        # create seizure directories
        for i,f in enumerate(self.all_eeg_files):
            print("====> Importing seizure %d from %s..." % (int(i+1), f))
            # load edf
            raw = mne.io.read_raw_edf(f, preload=True, verbose=False) 
            # get SR
            sfx = raw.info.get('sfreq')

            d_unfiltered = raw.get_data(picks=self.included_channels) # this reorders channels according to contiguous_labels
            if self.lowpass > 0:
                d = mne.filter.filter_data(d_unfiltered, sfx, l_freq=None, h_freq=self.lowpass, verbose=False)
                if self.lowpass > 60:
                    d = mne.filter.notch_filter(d, sfx, 60, verbose=False)
            else:
                d = mne.filter.notch_filter(d_unfiltered, sfx, 60, verbose=False) # mne 60 hz notch with default settings

            sznumstr = "_%02d" % int(i+1)
            szstring = self.subjname + sznumstr
            szdir = os.path.join(self.opscea_subjdir, szstring)
            if not os.path.isdir(szdir):
                os.makedirs(szdir)
            
            fstem = f.split('/')[-1][:-4]
            fstemparts = fstem.split("_")

            # handle patient labels with underscore(s)
            if len(fstemparts)==8:
                ptlabel_limit = 1
            elif len(fstemparts)==9:
                ptlabel_limit = 2
            else:
                raise EEGFilenameException("Unexpected format for patient label.")

            szmatname = "_".join(fstemparts[:ptlabel_limit]) + sznumstr + "__" + "_".join(fstemparts[ptlabel_limit:]) + "_sz.mat"
            output_path = os.path.join(szdir, szmatname)
            scipy.io.savemat(output_path, {'d': d, 'sfx': sfx})

            badch_matname = szstring + "_badch.mat"
            badch_output_path = os.path.join(szdir, badch_matname)
            scipy.io.savemat(badch_output_path, {'badch': badch})

    def do_label_output(self, outputdir):
        stem = ''
        for ll in self.fslabels:
            thisstem = self.get_stem(ll[1], self.numberlabels)
            if thisstem!=stem: # starting new lead
                if len(stem)>0:
                    outf.close()
                stem = thisstem
                ofname = thisstem + "_labels.tex"
                outf = open(os.path.join(outputdir, ofname), 'w')
            outf.write(ll[1] + ' & ' + ll[3].replace('_', '\_') + ' \\\\\n')
        outf.close()

    def make_cell_array(self, lis):
        ca = np.array([lis], dtype=object)
        return np.transpose(ca)

    def print_contactdict(self):
        for k,v in self.contactdict.items():
            print(k + ": " + str(v))

    def save_contactdict(self):
        output_path = os.path.join(self.opscea_subjdir, 'Imaging/Elecs', 'contactdict.mat')
        wirelabels = list(self.contactdict.keys())
        extremecontacts = self.matrix_from_dict_vals(self.contactdict)
        allcolors = np.array([[1,0,0],
                          [0,1,0],
                          [0,0,1],
                          [1,0,1],
                          [1,0.43921568627451,0],
                          [0.6,0,1],
                          [0.325490196078431,0.568627450980392,0.992156862745098],
                          [0.0313725490196078,0.815686274509804,0.976470588235294],
                          [0.403921568627451,0.466666666666667,0.0980392156862745],
                          [0.588235294117647,0.588235294117647,0.588235294117647],
                          [0.180392156862745,0.545098039215686,0.341176470588235],
                          [0.545098039215686,0.266666666666667,0.0745098039215686],
                          [1,0.549019607843137,0.411764705882353],
                          [0.419607843137255,0.352941176470588,0.803921568627451],
                          [1,0,0],
                           [0,1,0]])
        colors = allcolors[:len(wirelabels), :]
        scipy.io.savemat(output_path, {'wirelabels': self.make_cell_array(wirelabels), 'extremecontacts': extremecontacts, 'colors':colors})

    def matrix_from_dict_vals(self, d, l=[]):
        if len(l)==0:
            l = list(d.keys())
            
        for i,cl in enumerate(l):
            if i==0:
                m = np.array([d[cl]])
            else:
                m = np.append(m, [d[cl]], axis=0)

        return m

    def remove_spacers(self, s):
        # assumes you can have spacing with either a space or underscore
        # but not both
        sparts = s.split()
        if len(sparts)==2:
            return sparts[0] + sparts[1]
        sparts = s.split("_")
        if len(sparts)==2:
            return sparts[0] + sparts[1]
        return s

    def get_stem(self, st, numberlabels):
        for nl in numberlabels:
            if st[:len(nl)]==nl:
                return nl
        return self.stripdigits(st)

    def stripdigits(self, s):
        s1 = ''
        for c in s:
            if not c.isalpha() and not c.isspace():
                break
            else:
                s1 += c
        return s1

class BrainstormOpsceaMaker(OpsceaMaker):
    def do_imaging_elecs(self):
        path_to_brainstormdb =  os.path.join(os.environ['HOME'], "Documents/brainstorm_db/IEEG_visualization/")
        brainstorm_channel_data_path = os.path.join(path_to_brainstormdb, "data", self.subjname, "*", "channel.mat") 
        # necessary because importing any edf to brainstorm creates a
        # channel.mat with null coordinates for each electrode, so we
        # need to find the one with non-null coordinates
        all_channel_mats = glob(brainstorm_channel_data_path)

        print("====> Choosing channel.mat... ")
        found_good_elecmatrix = False
        for f in all_channel_mats:
            try:
                elecmatrix = self.build_elecmatrix(f)
                if self.is_good_elecmatrix(elecmatrix):
                    print("====> Selected channel.mat: " + f)
                    found_good_elecmatrix = True
                    break
            except DollarException:
                print("WARNING: Skipping channel mat " + f + " because of $ in channel name.")
            except ChanfileException:
                print("WARNING: Skipping channel mat " + f + ".")

        if not found_good_elecmatrix:
            raise Exception("No valid channel.mat was found. Please find a valid channel.mat.")
        print("====> Selected %d EDF channels." % len(self.contiguous_labels))

        brainstorm_anat_path = os.path.join(path_to_brainstormdb, "anat", self.subjname, "subjectimage_MRI.mat")
        elecmatrix = self.do_mri_coordinate_transformation(elecmatrix, brainstorm_anat_path)

        if not os.path.isdir(os.path.join(self.freesurfer_subjdir, "elecs")):
            os.makedirs(os.path.join(self.freesurfer_subjdir, "elecs"))
        output_path1 = os.path.join(self.freesurfer_subjdir, "elecs", "elecs_all.mat")
        scipy.io.savemat(output_path1, {'eleclabels': self.make_cell_array(self.contiguous_labels), 'elecmatrix': elecmatrix})

        if self.do_label:
            print("====> Getting FreeSurfer labels...")
            self.fslabels = self.patient.label_elecs(elecfile_prefix = "elecs_all", all_depth = True, quietmode = True)

        print("====> Extreme contacts:")
        self.print_contactdict()
        self.save_contactdict()
        
        output_path2 = os.path.join(self.opscea_subjdir, 'Imaging/Elecs', 'Electrodefile.mat')
        os.system('ln -sf %s %s' % (output_path1, output_path2))

    def build_elecmatrix(self, channelmat):
        # inputs:
        # channelmat- matrix listing 3D CTF coords for every contact as listed in brainstorm,
        #             can be discontiguous, may include contacts not in edffile
        # edffile- edf file with EEG
        # numberlabels- optional list of electrodes that end with a digit
        # 
        # outputs:
        # clabels: alphabetized, contiguous edf-based contact listing with any whitespage removed    
        # elecmatrix- matrix listing 3D coords (not transformed) for every contact as listed in clabels
        # included_channels- list of indices of included edf channels
        # contact_dict- dict mapping each electrode to extreme contacts
        channelinfo = scipy.io.loadmat(channelmat)
        df = pd.DataFrame(channelinfo['Channel'][0]) # convert to pandas
        rowsel = (df['Type']=='SEEG') | (df['Type']=='ECOG')
        if not rowsel.any(): 
            raise ChanfileException("No valid channels in " + channelmat + ".")
        bslabels = df[rowsel]['Name']
        self.bslabels = np.transpose(np.array(bslabels))
        coords_list = df[rowsel]['Loc']
        
        raw = mne.io.read_raw_edf(self.all_eeg_files[0], verbose=False)
        self.edflabels = raw.info.get('ch_names')
        self.get_contiguous_labels()
        
        elecdict = {} # map contacts to 3D coordinates
        for i,c in enumerate(coords_list):
            elecdict[self.remove_spacers(self.bslabels[i][0])] = [c[0][0], c[1][0], c[2][0]]
        
        elecmatrix = self.matrix_from_dict_vals(elecdict, self.contiguous_labels)
        
        return elecmatrix

    def get_contiguous_labels(self):
        # first, choose the channels; use bslabels to exclude unwanted
        # ones like $LOF11
        
        non_seeg_channels = ["SpO2", "EtCO2", "Pulse", "CO2Wave", "EKG", "EKG1", "EKG2", "BP1", "Annotations"]
        included_labels = []
        edfdict = {} # map label to edf channel #
        for e in self.bslabels:
            label = self.remove_spacers(e[0]) # handle underscore/space
            if label not in non_seeg_channels:
                for k,edl in enumerate(self.edflabels):
                    if edl[0:4]=='POL ' and self.remove_spacers(edl[4:])==label:
                        included_labels.append(label)
                        edfdict[label] = k
                        break

        # make a dict to record the extreme contacts for each lead
        included_labels.sort()
        self.check_dollar_channel(included_labels)
        self.contactdict = {} # map label to extreme contacts
        for i,l in enumerate(included_labels):
            if i==0:
                this_stem = self.get_stem(l, self.numberlabels)
                first = 1
                current = 1
            else:
                new_stem = self.get_stem(l, self.numberlabels)
                current += 1
                if i==len(included_labels)-1:
                    self.contactdict[this_stem] = [first, current]
                elif new_stem!=this_stem:
                    self.contactdict[this_stem] = [first, current - 1]
                    this_stem = new_stem
                    first = current

        # use contactdict to make a list of all contacts for all leads
        self.contiguous_labels = []
        for k,v in self.contactdict.items():
            numcontacts = v[1] - v[0] + 1
            for i in range(numcontacts):
                self.contiguous_labels.append(k + str(int(i+1)))

        # make a list of indices of channels we will select from the edf
        self.included_channels = []
        for cl in self.contiguous_labels:
            self.included_channels.append(edfdict[cl])

    def do_mri_coordinate_transformation(self, bsmatrix, bs_anat_path):
        # transform brainstorm (SCS) coordinates to World coordinates (via
        # MRI coodinates)
        #
        # SCS: subject based, with origin at midpoint of line from LPA to RPA
        # 
        # MRI coordinates: mm, indexed from corner of volume; can convert
        # to voxel index by dividing by voxel size
    
        scale_factor = 1000
        mriinfo = scipy.io.loadmat(bs_anat_path)
    
        # see /Applications/brainstorm3/toolbox/anatomy/cs_convert.m
        # cs_convert l.79
        vox2ras = mriinfo['InitTransf'][0][1]
        # Convert from MRI(mm) to voxels
        vox2ras = np.dot(vox2ras, np.diag(1./np.append(mriinfo['Voxsize'], 1)))
        # Convert to meters: transformation MRI=>WORLD (in meters)
        scalematrix = np.vstack((np.hstack((np.ones([3,3]), 1e-3*np.ones([3,1]))), np.ones([1,4])))
        mri2world = np.multiply(vox2ras, scalematrix)
        
        # cs_convert l.112
        # scsmap maps MRI (mm) coordinates to SCS (brainstorm) coordinates
        rot = mriinfo['SCS']['R'][0][0] # 3 x 3
        tra = mriinfo['SCS']['T'][0][0] # 3 x 1
        scsmap = np.vstack((np.hstack((rot, 1e-3*tra)), np.array([0, 0, 0, 1])))
        # Transf maps SCS ==> MRI ==> World
        Transf = np.dot(mri2world, np.linalg.inv(scsmap))
        eleconesmatrix = np.hstack((bsmatrix, np.ones([bsmatrix.shape[0],1])))
        eleconesmatrix = np.transpose(np.dot(Transf, np.transpose(eleconesmatrix)))
        return scale_factor*eleconesmatrix[:, 0:3]

    def check_dollar_channel(self, lis):
        for c in lis:
            if c[0]=="$":
                raise DollarException("Channel listing contains invalid name: " + c + ".  Please correct channel file in Brainstorm.")
    
    def is_good_elecmatrix(self, elecmatrix):
        min_size = 1
        return np.linalg.norm(np.sum(abs(elecmatrix), axis=0)) >= min_size

    def do_imaging_recon_figs(self):
        bsimages = glob(os.path.join(self.freesurfer_subjdir, 'ct', '*.png'))
        for bsi in bsimages:
            shutil.copy2(bsi, os.path.join(self.opscea_subjdir, 'Imaging/Recon/figs'))

class ChanfileException(Exception):
    pass

class DollarException(ChanfileException):
    pass

class EEGfilenameException(Exception):
    pass

if __name__=="__main__":
    om = BrainstormOpsceaMaker()
    om.driver()
