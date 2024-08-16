.. _`data-preprocess`:
========================================
Sampling Design and Data Organization
========================================

Stereo-seq sequencing data were preprocessed using 

**SAW** (https://github.com/STOmics/SAW) 

to generate spatial gene expression matrices in 

**GEM** format (https://stereopy.readthedocs.io/en/latest/Tutorials/IO.html#GEM).

Image files, which are usually in **TIFF** format, generated during the stereo-seq library construction process should be ready for further process together with the GEM files.


Sampling
---------------------------------
<img src='../_static/sampling design.png' width=300 align='right' hspace='50' />Planarian animals of one uninjured individual and 16 regenerative individuals from eight distinct stages were sampled.


Sample fixation
---------------------------------
<img src='../_static/oct embedding.png' width=300 align='right' hspace='50' /> Individuals were embedded into two OCT blocks, of which relevant data were organized as two samples in [STOmicsDB](https://db.cngb.org/stomics/) database. 
- [STSA0000254](https://db.cngb.org/stomics/sample/STSA0000254/)
- [STSA0000255](https://db.cngb.org/stomics/sample/STSA0000255/)


Stereo-seq
---------------------------------
<img src='../_static/section.png' width=300 align='right' hspace='50' /> Each sample block was serially sectioned at 10 µm intervals throughtout the entire animal body. Each section was individually mounted onto the Stereo-seq chip. 

Sections on the chip were stained with nucleic acid dye for single-stranded DNA (ssDNA) visualization and further processed throughout library construction and sequencing. 


Data preprocess
---------------------------------
Stereo-seq sequencing data were preprocessed using [SAW](https://github.com/STOmics/SAW) to generation spatial gene expression matrices in [GEM format](https://stereopy.readthedocs.io/en/latest/Tutorials/IO.html#GEM)

Data of each section were packed with a sinlge **Tissue Section** ID in STOmicsDB database, including staining image *(.tif)*, GEM file *(.gem)* as well as relevant annotation information *(.txt)*.
- [STTS0000461 - 515](https://db.cngb.org/stomics/project/STT0000028)


Data availability
---------------------------------

All data generated in this study were deposited at CNGB Nucleotide Sequence Archive (accession code: STT0000028).

Processed data can be interactively explored from our PRISTA4D database (https://db.cngb.org/stomics/prista4d). 

