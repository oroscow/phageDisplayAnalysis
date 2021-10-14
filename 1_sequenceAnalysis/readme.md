# Sequence Analysis
 
Analyses the University of Guelph's Advanced Analytics Centre (AAC) Genomics Facility Sanger sequencing output. Original
files are left unaltered and copies are renamed, reorganised, converted to fasta, trimmed, translated, and aligned.
Final amino acid/nucleotide alignments are in fasta, clustal, and xlsx formats.

![gui](https://user-images.githubusercontent.com/55511532/137400899-cea27a53-72dc-4130-bcb7-eda1edcd1a98.png)

## Details:

**Input:** Raw DNA Sanger sequencing data from the University of Guelph's Advanced Analytics Centre (AAC) Genomics Facility
(seq and ab1).

**Output:** Raw DNA sequence data (original seq and ab1 files) are sorted into their own folders, copied and converted
into fasta format, trimmed at 5' and 3' cut sites, translated to amino acid sequences (fasta), and both amino acid and
nucleotide sequences are aligned (fasta, clustal, xlsx). The Excel alignment (xlsx) has four worksheets containing all
amino acid sequences, unique amino acid sequences, all nucleotide sequences, and unique nucleotide sequences,
respectively. Frequency and corresponding sequence IDs are generated for each unique sequence. Sequences that cannot be
trimmed are moved to their own folder for independent manual analysis.

* **Final alignments exclude:**

    a) Sequences that are not full length.

    b) Sequences that have premature stop codons.
* **Final alignments include:**

    a) Sequences with uncalled residues/base pairs (Xs and Ns, respectively).
* Amino acid/nucleotide sequences that cannot be trimmed at the 5' or 3' ends are excluded from the alignment and moved
to a separate folder for manual analysis.
* Truncated amino acid/nucleotide sequences that have been excluded from alignment are still included in batch files.

## Compatibility:
* PyCharm is the recommended IDE to use for running terminal scripts. If using Spyder, avoid version 5 as this version
for has conflicts with the xlsxwriter package and may get stuck on importing modules.
* Confirmed to work with Python 3.9. Later/earlier versions may work but have not been verified.
* Confirmed to work in Windows and unconfirmed in Macs and Linux. Path names may need to be changed to suit Macs
and Linux' formats.

## Usage:

### GUI/executable

1. Go to _phageDisplayAnalysis/1_sequenceAnalysis/gui/executable_ and download the entire folder as a zip file.
2. Extract to _C:\Program Files_ (give administrator permission if necessary).
3. Right click '_gui_sequenceAnalysis.exe_' and click 'Create shortcut'. The shortcut will likely relocate to your
desktop.
4. Open the program with the shortcut.
5. The first time the program is run, Windows may try to prevent you from opening the file in order to protect your PC
(see below). Click on 'More info' and then click on 'Run anyway'. Once you've done this once, it won't ask again.

![warning](https://external-content.duckduckgo.com/iu/?u=https%3A%2F%2Fmonstersocial.net%2Fwp-content%2Fuploads%2F2015%2F08%2Fwindowsprotectedyourpc.jpg&f=1&nofb=1)

###### Image credit to [monstersocial.net](https://monstersocial.net/).

Note: Only do this for trusted sources. Be very careful about opening executables sent by strangers and make sure to do
your due diligence and carefully examine all sources online.

### Terminal

1. Go to _phageDisplayAnalysis/1_sequenceAnalysis/terminal/_ and download '_terminal_sequenceAnalysis.py_'.
2. Run in IDE of choice (PyCharm recommended).
3. Follow the prompts.
