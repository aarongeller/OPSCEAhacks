#!/usr/bin/python

import os, sys
from glob import glob

subjname = sys.argv[1]

path_to_opsceadata = "/Users/aaron/Documents/MATLAB/OPSCEA-main/OPSCEADATA"
opscea_subjdir = os.path.join(path_to_opsceadata, subjname)
opscea_texdir = os.path.join(opscea_subjdir, "Imaging", "Recon")
opscea_slicedir = os.path.join(opscea_texdir, "figs")
opscea_labeldir = os.path.join(opscea_texdir, "labels")

texfname = os.path.join(opscea_texdir, subjname + "_recon.tex")
texf = open(texfname, 'w')

subjname_parts = subjname.split("_")
if len(subjname_parts) > 1:
    subjname = "\_".join(subjname_parts)

preamble = "\\documentclass[12pt]{article}\n\\usepackage{graphicx}\n\\usepackage{hyperref}\n\\hypersetup{colorlinks=true,linkcolor=blue}\n\\renewcommand{\\familydefault}{\\sfdefault}\n\\title{" + subjname + " ICEEG Implant Reconstruction}\n\\begin{document}\n\\maketitle\n\n\\tableofcontents\n\clearpage\n\n"
texf.write(preamble)

surfsize_a = 0.4
surfsize_c = 0.3
slicesize_a = 0.7
slicesize_c = 0.9
wholebrainsize = 0.6

imgfiles = glob(os.path.join(opscea_slicedir, '*.png'))
imgfiles.sort()

wholebrains = []
for f in imgfiles:
    if os.path.basename(f)[:4]=="SEEG":
        wholebrains.append(f)

if len(wholebrains)==4:
    texlines = "\\setcounter{section}{-1}\\section{Whole Brain}\n\\begin{tabular}{cc}\n\\includegraphics[width=" + str(wholebrainsize) + "\\textwidth]{" + wholebrains[1] + "} &\n\\includegraphics[width=" + str(wholebrainsize) + "\\textwidth]{" + wholebrains[0] + "}\\\\\n\\includegraphics[width=" + str(wholebrainsize) + "\\textwidth]{" + wholebrains[2] + "} &\n\\includegraphics[width=" + str(wholebrainsize) + "\\textwidth]{" + wholebrains[3] + "}\\\\\n\\end{tabular}\n\\clearpage\n\n"
    texf.write(texlines)
    
for f in imgfiles:
    fparts = os.path.basename(f).split('_')
    if os.path.basename(f)[:4]=="SEEG" or len(fparts)<4 or fparts[3]=='a.png':
        continue
    axial_surf = os.path.join(os.path.dirname(f), fparts[0] + "_" + fparts[1] + "_surf_a.png")
    coronal_surf = os.path.join(os.path.dirname(f), fparts[0] + "_" + fparts[1] + "_surf_c.png")
    axial_nonsurf = os.path.join(os.path.dirname(f), fparts[0] + "_" + fparts[1] + "_a.png")
    coronal_nonsurf = os.path.join(os.path.dirname(f), fparts[0] + "_" + fparts[1] + "_c.png")
    elecname = fparts[1]
    labelfile = os.path.join(opscea_labeldir, elecname + "_labels")
    texlines = "\\section{" + elecname + "}\n\\subsection{" + elecname + " Axial View}\n\\includegraphics[width=" + str(surfsize_a) + "\\textwidth]{" + axial_surf + "}\\\\\n\\includegraphics[width=" + str(slicesize_a) + "\\textwidth]{" + axial_nonsurf + "}\n\n\\subsection{" + elecname + " Coronal View}\n\\includegraphics[width=" + str(surfsize_c) + "\\textwidth]{" + coronal_surf + "}\\\\\n\\includegraphics[width=" + str(slicesize_c) + "\\textwidth]{" + coronal_nonsurf + "}\n\n\\subsection{" + elecname + " FreeSurfer Labels}\n\\begin{tabular}{ll}\n\\input{" + labelfile + "}\n\\end{tabular}\n\\clearpage\n\n"
    texf.write(texlines)

postamble = "\\end{document}\n"
texf.write(postamble)

texf.close()

os.system("pdflatex -output-directory " + opscea_texdir + " " + texfname)
os.system("pdflatex -output-directory " + opscea_texdir + " " + texfname)
