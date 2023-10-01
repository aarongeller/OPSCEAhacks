# Selected functions from the UCSF img_pipe module
# See https://github.com/ChangLabUcsf/img_pipe

import os, glob
import nibabel as nib
import numpy as np
import scipy.io

class freeCoG:
    ''' This defines the class freeCoG, which creates a patient object      
    for use in creating brain surface reconstructions, electrode placement,     
    and warping.        
    
    To initialize a patient, you must provide the subject ID, hemisphere,       
    freesurfer subjects directory, and (optionally) the freesurfer      
    executable directory. If these aren't specified they will default to the environment
    variables $SUBJECTS_DIR and $FREESURFER_HOME.
    
    For example:        
    
    >>> subj = 'EC1'     
    >>> subj_dir = '/usr/local/freesurfer/subjects'      
    >>> hem = 'rh'       
    >>> fs_dir = '/usr/local/freesurfer'     
    >>> patient = img_pipe.freeCoG(subj = subj, subj_dir = subj_dir, hem = hem, fs_dir = fs_dir)
    
    Parameters
    ----------       
    subj : str 
           The subject ID      
    hem : {'lh', 'rh', 'stereo'}
          The hemisphere of implanation. Can be 'lh', 'rh', or 'stereo'.        
    zero_indexed_electrodes: bool, optional
                             Whether or not to use zero-indexing for the electrode numbers (default: True)
    fs_dir : str, optional
             Path to the freesurfer install (default is $FREESURFER_HOME environmental variable) 
    subj_dir : str, optional
               the freesurfer subjects directory (e.g. /usr/local/freesurfer/subjects). Default
               is the $SUBJECTS_DIR environmental variable set by Freesurfer.
    
    Attributes
    ----------
    subj : str
        The subject ID.
    subj_dir : str
        The subjects directory where the data live.
    patient_dir : str
        The directory containing data for this particular patient.
        Usually [subj_dir]/[subj]
    hem : str
        The hemisphere of implantation
    img_pipe_dir : str
        The path to img_pipe code.
    zero_indexed_electrodes : bool
        Whether zero-indexed electrode numbers are used (True)
        or not (False)
    fs_dir : str
        Path to the freesurfer installation
    CT_dir : str
        Path to the CT scan
    elecs_dir : str
        Path to the marked electrode file locations
    mesh_dir : str
        Path to the generated triangle-mesh surfaces (pial, subcortical, etc.)
    pial_surf_file : dict
        Dictionary containing full file with path for the left and right pial surfaces.
        Left pial surface is pial_surf_file['lh'] and right pial surface is
        pial_surf_file['rh']
    surf_dir : str
        The freesurfer surf directory for this patient.
    mri_dir : str
        The freesurfer MRI directory for this patient.

    Returns
    ----------
    patient : <img_pipe.freeCoG instance>
        patient object, including information about subject ID, hemisphere, where data live, 
        and related functions for creating surface reconstructions and/or plotting.    

    '''

    def __init__(self, subj, hem, zero_indexed_electrodes=True, fs_dir=os.environ['FREESURFER_HOME'], subj_dir=os.environ['SUBJECTS_DIR']):
        '''
        Initializes the patient object.

        Parameters
        ----------       
        subj : str 
            The subject ID (e.g. 'SUBJ_25')     
        hem : {'lh', 'rh', 'stereo'}
            The hemisphere of implanation. Can be 'lh', 'rh', or 'stereo'.        
        zero_indexed_electrodes: bool, optional
            Whether or not to use zero-indexing for the electrode numbers (default: True)
        fs_dir : str, optional
            Path to the freesurfer install (default is $FREESURFER_HOME environmental variable) 
        subj_dir : str, optional
            Path to the freesurfer subjects (default is $SUBJECTS_DIR environmental variable)
    
        Returns
        ----------   
        patient : <img_pipe.freeCoG instance>
            patient object, including information about subject ID, hemisphere, where data live, 
            and related functions for creating surface reconstructions and/or plotting.    
        '''
        # Check if hem is valid
        if not hem in ['rh', 'lh', 'stereo']:
            raise NameError('Invalid hem for freeCoG')
        
        self.subj = subj
        self.subj_dir = subj_dir
        self.patient_dir = os.path.join(self.subj_dir, self.subj)
        self.hem = hem
        self.img_pipe_dir = os.path.dirname(os.path.realpath(__file__))
        self.zero_indexed_electrodes = zero_indexed_electrodes

        # Freesurfer home directory
        self.fs_dir = fs_dir

        # acpc_dir: dir for acpc MRIs
        self.acpc_dir = os.path.join(self.subj_dir, self.subj, 'acpc')

        # CT_dir: dir for CT img data
        self.CT_dir = os.path.join(self.subj_dir, self.subj, 'CT')

        # dicom_dir: dir for storing dicoms
        self.dicom_dir = os.path.join(self.subj_dir, self.subj, 'dicom')

        # elecs_dir: dir for elecs coordinates
        self.elecs_dir = os.path.join(self.subj_dir, self.subj, 'elecs')

        # Meshes directory for matlab/python meshes
        self.mesh_dir = os.path.join(self.subj_dir, self.subj, 'Meshes')
        if self.hem == 'stereo':
            surf_file = os.path.join(self.subj_dir, self.subj, 'Meshes', 'lh_pial_trivert.mat')
        else:
            surf_file = os.path.join(self.subj_dir, self.subj, 'Meshes', self.hem + '_pial_trivert.mat')
        if os.path.isfile(surf_file):
            self.pial_surf_file = dict()
            self.pial_surf_file['lh'] = os.path.join(self.subj_dir, self.subj, 'Meshes', 'lh_pial_trivert.mat')
            self.pial_surf_file['rh'] = os.path.join(self.subj_dir, self.subj, 'Meshes', 'rh_pial_trivert.mat')

        # surf directory
        self.surf_dir = os.path.join(self.subj_dir, self.subj, 'surf')

        # mri directory
        self.mri_dir = os.path.join(self.subj_dir, self.subj, 'mri')

        # Check if subj is valid
        if not os.path.isdir(os.path.join(self.subj_dir, self.subj)):
            print('No such subject directory as %s' % (os.path.join(self.subj_dir, self.subj)))
            ans = raw_input('Would you like to create a new subject directory? (y/N): ')
            if ans.upper() == 'Y' or ans.upper() == 'YES':
                os.chdir(self.subj_dir)
                print("Making subject directory")
                os.mkdir(self.subj)
                os.mkdir(self.acpc_dir)
                os.mkdir(self.CT_dir)
                os.mkdir(self.dicom_dir)
                os.mkdir(self.elecs_dir)
                os.mkdir(self.mesh_dir)
                os.mkdir(self.mri_dir)
            else:
                raise NameError('No such subject directory as %s' % (os.path.join(subj_dir, subj)))

        #if paths are not the default paths in the shell environment:
        os.environ['FREESURFER_HOME'] = fs_dir
        os.environ['SUBJECTS_DIR'] = subj_dir

    def convert_fsmesh2mlab(self, mesh_name='pial'):
        '''Creates surface mesh triangle and vertex .mat files
        If no argument for mesh_name is given, lh.pial and rh.pial
        are converted into lh_pial_trivert.mat and rh_pial_trivert.mat
        in the Meshes directory (for use in python) and *_lh_pial.mat
        and *_rh_pial.mat for use in MATLAB.

        Parameters
        ----------
        mesh_name : {'pial', 'white', 'inflated'}
        
        '''

        hems = ['lh', 'rh']

        if not os.path.isdir(self.mesh_dir):
            print('Making Meshes Directory')
            # Make the Meshes directory in subj_dir if it does not yet exist
            os.mkdir(self.mesh_dir)

        # Loop through hemispheres for this mesh, create one .mat file for each
        for h in hems:
            print("Making %s mesh"%(h))
            mesh_surf = os.path.join(self.surf_dir, h+'.'+mesh_name)
            vert, tri = nib.freesurfer.read_geometry(mesh_surf)
            out_file = os.path.join(self.mesh_dir, '%s_%s_trivert.mat'%(h, mesh_name))
            out_file_struct = os.path.join(self.mesh_dir, '%s_%s_%s.mat'%(self.subj, h, mesh_name))
            scipy.io.savemat(out_file, {'tri': tri, 'vert': vert})

            cortex = {'tri': tri+1, 'vert': vert}
            scipy.io.savemat(out_file_struct, {'cortex': cortex})

        if mesh_name=='pial':
            self.pial_surf_file = dict()
            self.pial_surf_file['lh'] = os.path.join(self.subj_dir, self.subj, 'Meshes', 'lh_pial_trivert.mat')
            self.pial_surf_file['rh'] = os.path.join(self.subj_dir, self.subj, 'Meshes', 'rh_pial_trivert.mat')
        else:
            setattr(self, mesh_name+'_surf_file', out_file)
  


    def get_subcort(self):
        '''Obtains .mat files for vertex and triangle
           coords of all subcortical freesurfer segmented meshes'''

        # set ascii dir name
        subjAscii_dir = os.path.join(self.subj_dir, self.subj, 'ascii')
        if not os.path.isdir(subjAscii_dir):
            os.mkdir(subjAscii_dir)

        # tessellate all subjects freesurfer subcortical segmentations
        print('::: Tesselating freesurfer subcortical segmentations from aseg using aseg2srf... :::')
        print(os.path.join(self.img_pipe_dir, 'aseg2srf.sh'))
        os.system(os.path.join(self.img_pipe_dir, 'aseg2srf.sh') + ' -s "%s" -l "4 5 10 11 12 13 17 18 26 \
                 28 43 44  49 50 51 52 53 54 58 60 14 15 16" -d' % (self.subj))

        # get list of all .srf files and change fname to .asc
        srf_list = list(set([fname for fname in os.listdir(subjAscii_dir)]))
        asc_list = list(set([fname.replace('.srf', '.asc') for fname in srf_list]))
        asc_list.sort()
        for fname in srf_list:
            new_fname = fname.replace('.srf', '.asc')
            os.system('mv %s %s'%(os.path.join(subjAscii_dir,fname), os.path.join(subjAscii_dir,new_fname)))

        # convert all ascii subcortical meshes to matlab vert, tri coords
        subcort_list = ['aseg_058.asc', 'aseg_054.asc', 'aseg_050.asc',
                        'aseg_052.asc', 'aseg_053.asc', 'aseg_051.asc', 'aseg_049.asc',
                        'aseg_043.asc', 'aseg_044.asc', 'aseg_060.asc', 'aseg_004.asc',
                        'aseg_005.asc', 'aseg_010.asc', 'aseg_011.asc', 'aseg_012.asc',
                        'aseg_013.asc', 'aseg_017.asc', 'aseg_018.asc', 'aseg_026.asc',
                        'aseg_028.asc', 'aseg_014.asc', 'aseg_015.asc', 'aseg_016.asc']

        nuc_list = ['rAcumb', 'rAmgd', 'rCaud', 'rGP', 'rHipp', 'rPut', 'rThal',
                    'rLatVent', 'rInfLatVent', 'rVentDienceph', 'lLatVent', 'lInfLatVent',
                    'lThal', 'lCaud', 'lPut',  'lGP', 'lHipp', 'lAmgd', 'lAcumb', 'lVentDienceph',
                    'lThirdVent', 'lFourthVent', 'lBrainStem']

        subcort_dir = os.path.join(self.mesh_dir,'subcortical')     
        if not os.path.isdir(subcort_dir):      
            print('Creating directory %s'%(subcort_dir))        
            os.mkdir(subcort_dir)

        print('::: Converting all ascii segmentations to matlab tri-vert :::')
        for i in range(len(subcort_list)):
            subcort = os.path.join(subjAscii_dir, subcort_list[i])
            nuc = os.path.join(subcort_dir, nuc_list[i])
            self.subcortFs2mlab(subcort, nuc)

    def subcortFs2mlab(self, subcort, nuc):
        '''Function to convert freesurfer ascii subcort segmentations
           to triangular mesh array .mat style. This helper function is
           called by get_subcort().
        
        Parameters
        ----------
        subcort : str
            Name of the subcortical mesh ascii file (e.g. aseg_058.asc).
            See get_subcort().
        nuc : str
            Name of the subcortical nucleus (e.g. 'rAcumb')
        
        '''

        # use freesurfer mris_convert to get ascii subcortical surface
        subcort_ascii = subcort

        # clean up ascii file and extract matrix dimensions from header
        subcort = open(subcort_ascii, 'r')
        subcort_mat = subcort.readlines()
        subcort.close()
        subcort_mat.pop(0)  # get rid of comments in header
        # get rid of new line char
        subcort_mat = [item.strip('\n') for item in subcort_mat]

        # extract inds for vert and tri
        subcort_inds = subcort_mat.pop(0)
        # seperate inds into two strings
        subcort_inds = subcort_inds.split(' ')
        subcort_inds = [int(i) for i in subcort_inds]  # convert string to ints

        # get rows for vertices only, strip 0 column, and split into seperate
        # strings
        subcort_vert = [item.strip(' 0')
                        for item in subcort_mat[:subcort_inds[0]]]
        subcort_vert = [item.split('  ')
                        for item in subcort_vert]  # seperate strings

        # Convert into an array of floats
        subcort_vert = np.array(np.vstack((subcort_vert)), dtype=float)

        # get rows for triangles only, strip 0 column, and split into seperate
        # strings
        subcort_tri = [item[:-2] for item in subcort_mat[subcort_inds[0] + 1:]]
        subcort_tri = [item.split(' ')
                       for item in subcort_tri]  # seperate strings
        subcort_tri = np.array(np.vstack((subcort_tri)), dtype=int)

        outfile = '%s_subcort_trivert.mat' % (nuc)
        scipy.io.savemat(outfile, {'tri': subcort_tri, 'vert': subcort_vert})  # save tri/vert matrix

        # convert inds to scipy mat
        subcort_inds = scipy.mat(subcort_inds)
        scipy.io.savemat('%s_subcort_inds.mat' %
                         (nuc), {'inds': subcort_inds})  # save inds
        
        out_file_struct = '%s_subcort.mat' % (nuc)
        
        cortex = {'tri': subcort_tri+1, 'vert': subcort_vert}
        scipy.io.savemat(out_file_struct, {'cortex': cortex})

    def nearest_electrode_vert(self, cortex_verts, elecmatrix):
        ''' Find the vertex on a mesh that is closest to the given electrode
        coordinates.
        
        Parameters
        ----------
        cortex_verts : array-like
            [nvertices x 3] matrix of vertices on the cortical surface mesh
        elecmatrix : array-like
            [nchans x 3] matrix of 3D electrode coordinates 

        Returns
        -------
        vert_inds : array-like
            Array of vertex indices that are closest to each of the 
            electrode 
        nearest_verts : array-like
            Coordinates for the nearest cortical vertices
        '''

        nchans = elecmatrix.shape[0]
        d = np.zeros((nchans, cortex_verts.shape[0]))

        # Find the distance between each electrode and all possible vertices
        # on the surface mesh
        for chan in np.arange(nchans):
            d[chan,:] = np.sqrt((elecmatrix[chan,0] - cortex_verts[:,0])**2 + \
                                (elecmatrix[chan,1] - cortex_verts[:,1])**2 + \
                                (elecmatrix[chan,2] - cortex_verts[:,2])**2)

        # Find the index of the vertex nearest to each electrode
        vert_inds = np.argmin(d, axis = 1)
        nearest_verts = cortex_verts[vert_inds,:]

        return vert_inds, nearest_verts

    def label_elecs(self, elecfile_prefix='TDT_elecs_all', atlas_surf='desikan-killiany', atlas_depth='destrieux', elecs_all=True,all_depth=False,quietmode=False):
        ''' Automatically labels electrodes based on the freesurfer annotation file.
        Assumes TDT_elecs_all.mat or clinical_elecs_all.mat files
        Uses both the Desikan-Killiany Atlas and the Destrieux Atlas, as described 
        here: https://surfer.nmr.mgh.harvard.edu/fswiki/CorticalParcellation
        
        Parameters
        ----------
        elecfile_prefix : str, optional
            prefix of the .mat file with the electrode coordinates matrix
        atlas_surf : {'desikan-killiany', 'destrieux'}
            The atlas to use for labeling of surface electrodes.
        atlas_depth : {'destrieux', 'desikan-killiany'}
            The atlas to use for labeling of depth electrodes.
        elecs_all : bool
            Label all electrodes
        
        Returns
        -------
        elec_labels : array-like
            [nchans x 4] matrix of electrode labels. Columns include short name,
            long name, 'grid'/'depth'/'strip' label, and assigned anatomical label.
        '''

        if atlas_surf == 'desikan-killiany':
            surf_atlas_flag = ''
        elif atlas_surf == 'destrieux':
            surf_atlas_flag = '--a2009s'
        else:
            surf_atlas_flag = ''

        print(self.subj_dir)
        print('Creating labels from the freesurfer annotation file for use in automated electrode labeling')
        gyri_labels_dir = os.path.join(self.subj_dir, self.subj, 'label', 'gyri')
        if not os.path.isdir(gyri_labels_dir):
            os.mkdir(gyri_labels_dir)
         
        # This version of mri_annotation2label uses the coarse labels from the Desikan-Killiany Atlas, unless
        # atlas_surf is 'destrieux', in which case the more detailed labels are used
        if self.hem in ['lh','rh']:
            os.system('mri_annotation2label --subject %s --hemi %s --surface pial %s --outdir %s'%(self.subj, self.hem, surf_atlas_flag, gyri_labels_dir))
        else: # it's stereo
            for h in ['lh', 'rh']:
                os.system('mri_annotation2label --subject %s --hemi %s --surface pial %s --outdir %s'%(self.subj, h, surf_atlas_flag, gyri_labels_dir))
        print('Loading electrode matrix')
        elecfile = os.path.join(self.elecs_dir, elecfile_prefix+'.mat')
        elecmatrix = scipy.io.loadmat(elecfile)['elecmatrix']
        
        # Initialize empty variable for indices of grid and strip electrodes
        isnotdepth = []
        
        # Choose only the surface or grid electrodes (if not using hd_grid.mat)
        if elecfile_prefix == 'TDT_elecs_all' or elecfile_prefix == 'clinical_elecs_all' or elecs_all:
            elecmontage = scipy.io.loadmat(elecfile)['eleclabels']
            # Make the cell array into something more usable by python
            short_label = []
            long_label = []
            grid_or_depth = []

            for r in elecmontage:
                short_label.append(r[0][0]) # This is the shortened electrode montage label
                if all_depth: # added by AG
                    long_label.append(r[0][0]) 
                    grid_or_depth.append('depth') 
                else:                    
                    long_label.append(r[1][0]) # This is the long form electrode montage label
                    grid_or_depth.append(r[2][0]) # This is the label for grid, depth, or strip
            
            # These are the indices that won't be used for labeling
            #dont_label = ['EOG','ECG','ROC','LOC','EEG','EKG','NaN','EMG','scalpEEG']
            indices = [i for i, x in enumerate(long_label) if ('EOG' in x or 'ECG' in x or 'ROC' in x or 'LOC' in x or 'EEG' in x or 'EKG' in x or 'NaN' in x or 'EMG' in x or x==np.nan or 'scalpEEG' in x)]
            indices.extend([i for i, x in enumerate(short_label) if ('EOG' in x or 'ECG' in x or 'ROC' in x or 'LOC' in x or 'EEG' in x or 'EKG' in x or 'NaN' in x or 'EMG' in x or x==np.nan or 'scalpEEG' in x)])
            indices.extend([i for i, x in enumerate(grid_or_depth) if ('EOG' in x or 'ECG' in x or 'ROC' in x or 'LOC' in x or 'EEG' in x or 'EKG' in x or 'NaN' in x or 'EMG' in x or x==np.nan or 'scalpEEG' in x)])
            indices.extend(np.where(np.isnan(elecmatrix)==True)[0])
            indices = list(set(indices))
            indices_to_use = list(set(range(len(long_label))) - set(indices))

            # Initialize the cell array that we'll store electrode labels in later
            elec_labels_orig = np.empty((len(long_label),4),dtype=object)
            elec_labels_orig[:,0] = short_label
            elec_labels_orig[:,1] = long_label
            elec_labels_orig[:,2] = grid_or_depth 
            elec_labels = np.empty((len(indices_to_use),4), dtype = object)
            elecmatrix_orig = elecmatrix
            elecmatrix = elecmatrix[indices_to_use,:]
            
            short_label_orig,long_label_orig,grid_or_depth_orig = short_label,long_label,grid_or_depth
            short_label = [i for j, i in enumerate(short_label) if j not in indices]
            long_label = [i for j, i in enumerate(long_label) if j not in indices]
            grid_or_depth = [i for j, i in enumerate(grid_or_depth) if j not in indices]
            elec_labels[:,0] = short_label
            elec_labels[:,1] = long_label
            elec_labels[:,2] = grid_or_depth
            
            # Find the non depth electrodes
            isnotdepth = np.array([r!='depth' for r in grid_or_depth])
            
        # Use the surface label files to get which label goes with each surface vertex
        label_files = glob.glob(os.path.join(gyri_labels_dir, '%s.*.label'%(self.hem)))
        vert_label = {}
        for label in label_files:
            label_name = label.split('.')[1]
            print('Loading label %s'%label_name)
            fid = open(label,'r')
            d = np.genfromtxt(fid, delimiter=' ', \
                              skip_header=2)
            vertnum, x, y, z, junk=d[~np.isnan(d)].reshape((-1,5)).T
            for v in vertnum:
                vert_label[np.int(v)] = label_name.strip()
            fid.close()

        if self.hem in ['lh','rh']:
            trivert_file = os.path.join(self.mesh_dir, '%s_pial_trivert.mat'%(self.hem))
            cortex_verts = scipy.io.loadmat(trivert_file)['vert']
        else: # it's stereo
            for h in ['lh', 'rh']:
                trivert_file = os.path.join(self.mesh_dir, '%s_pial_trivert.mat' % h)
                cortex_verts = scipy.io.loadmat(trivert_file)['vert']

        # Only use electrodes that are grid or strips
        if len(isnotdepth)>0:
            elecmatrix_new = elecmatrix[isnotdepth,:]
        else:
            elecmatrix_new = elecmatrix

        print('Finding nearest mesh vertex for each electrode')
        vert_inds, nearest_verts = self.nearest_electrode_vert(cortex_verts, elecmatrix_new)

        ## Now make a dictionary of the label for each electrode
        elec_labels_notdepth=[]
        for v in range(len(vert_inds)):
            if vert_inds[v] in vert_label:
                elec_labels_notdepth.append(vert_label[vert_inds[v]].strip())
            else:
                elec_labels_notdepth.append('Unknown')

        if elecfile_prefix == 'TDT_elecs_all' or elecfile_prefix == 'clinical_elecs_all' or elecs_all:
            elec_labels[isnotdepth,3] = elec_labels_notdepth
            elec_labels[np.invert(isnotdepth),3] = '' # Set these to an empty string instead of None type
        else:
            elec_labels = np.array(elec_labels_notdepth, dtype = np.object)
        print('Saving electrode labels for surface electrodes to %s'%(elecfile_prefix))
        ## added by BKD so that elec_mat_grid='hd_grid' works. It does not contain elecmontage
        save_dict = {'elecmatrix': elecmatrix, 'anatomy': elec_labels}
        if 'elecmontage' in locals():
            save_dict['eleclabels'] = elecmontage
        else:
            print('electmontage does not exist')
        #scipy.io.savemat('%s/%s/elecs/%s'%(self.subj_dir, self.subj, elecfile_prefix), save_dict)

        if np.any(np.invert(isnotdepth)): # If there are depth electrodes, run this part
            print('*************************************************')
            print('Now doing the depth electrodes')

            # Get the volume corresponding to the labels from the Destrieux atlas, which is more 
            # detailed than Desikan-Killiany (https://surfer.nmr.mgh.harvard.edu/fswiki/CorticalParcellation)
            if atlas_depth == 'desikan-killiany':
                depth_atlas_nm = ''
            elif atlas_depth == 'destrieux':
                depth_atlas_nm = '.a2009s'
            else:
                depth_atlas_nm = '.a2009s'

            aseg_file = os.path.join(self.subj_dir, self.subj, 'mri', 'aparc%s+aseg.mgz'%(depth_atlas_nm))
            dat = nib.freesurfer.load(aseg_file)
            aparc_dat = dat.get_fdata() #UPDATED DEPRECATED CALL was dat.get_data() now dat.get_fdata()
             
            # Define the affine transform to go from surface coordinates to volume coordinates (as CRS, which is
            # the slice *number* as x,y,z in the 3D volume. That is, if there are 256 x 256 x 256 voxels, the
            # CRS coordinate will go from 0 to 255.)
            affine = np.array([[  -1.,    0.,    0.,  128.],
                               [   0.,    0.,    1., -128.],
                               [   0.,   -1.,    0.,  128.],
                               [   0.,    0.,    0.,    1.]])

            elecs_depths = elecmatrix[np.invert(isnotdepth),:]
            intercept = np.ones(len(elecs_depths))
            elecs_ones = np.column_stack((elecs_depths,intercept))

            # Find voxel CRS
            VoxCRS = np.dot(np.linalg.inv(affine), elecs_ones.transpose()).transpose().astype(int)
            # Make meshgrid the same size as aparc_dat (only for gaussian blob version), ignore
            #xx, yy, zz = np.mgrid[0:aparc_dat.shape[0], 0:aparc_dat.shape[1], 0:aparc_dat.shape[2]]
            #unique_labels = np.unique(aparc_dat)
            #unique_labels = unique_labels[unique_labels>0]

            # Get the names of these labels using Freesurfer's lookup table (LUT)
            print("Loading lookup table for freesurfer labels")
            fid = open(os.path.join(self.fs_dir,'FreeSurferColorLUT.txt'))
            LUT = fid.readlines()
            fid.close()

            # Make dictionary of labels
            LUT = [row.split() for row in LUT]
            lab = {}
            for row in LUT:
                if len(row)>1 and row[0][0]!='#' and row[0][0]!='\\': # Get rid of the comments
                    lname = row[1]
                    lab[int(row[0])] = lname

            # Label the electrodes according to the aseg volume
            nchans = VoxCRS.shape[0]
            anatomy = np.empty((nchans,), dtype=object)
            print("Labeling electrodes...")

            for elec in np.arange(nchans):
                anatomy[elec] = lab[aparc_dat[VoxCRS[elec,0], VoxCRS[elec,1], VoxCRS[elec,2]]]
                if not quietmode:
                    print("E%d, Vox CRS: [%d, %d, %d], Label #%d = %s"%(elec, VoxCRS[elec,0], VoxCRS[elec,1], VoxCRS[elec,2], 
                                                                        aparc_dat[VoxCRS[elec,0], VoxCRS[elec,1], VoxCRS[elec,2]], 
                                                                        anatomy[elec]))

            elec_labels[np.invert(isnotdepth),3] = anatomy
            
            #make some corrections b/c of NaNs in elecmatrix
        elec_labels_orig[:,3] = ''
        elec_labels_orig[indices_to_use,3] = elec_labels[:,3] 
        
        print('Saving electrode labels to %s'%(elecfile_prefix))
        scipy.io.savemat(os.path.join(self.elecs_dir, elecfile_prefix+'.mat'), {'elecmatrix': elecmatrix_orig, 
                                                                                'anatomy': elec_labels_orig, 
                                                                                'eleclabels': elecmontage})

        return elec_labels
