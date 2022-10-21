# !python3

##################
#    MODULES
##################

import warnings
import glob
import logging
import os
import re
import shutil
import xlsxwriter
from collections import Counter, OrderedDict
from shutil import copyfile
from Bio import AlignIO, BiopythonWarning
from Bio.Seq import Seq


##################
#    FUNCTIONS
##################

def cyanprint(text):
    """Print cyan coloured text in the console."""
    print('\033[0;36m' + text + '\033[0m')


def greenprint(text):
    """Print green coloured text in the console."""
    print('\033[0;32m' + text + '\033[0m')


##################
#    CLASSES
##################

# Ordered list of counts.
class OrderedCounter(Counter, OrderedDict):
    pass


##################
#    WARNINGS
##################

# Suppresses BiopythonWarning: 'Partial codon, len(sequence) not a multiple of three. Explicitly trim the sequence or
# add trailing N before translation. This may become an error in future.' When troubleshooting, comment out this line.
warnings.filterwarnings('ignore',
                        category=BiopythonWarning
                        )

##################
#    MAIN
##################

##################
# Set up working directory and logging file.
##################

# Working directory setup. User prompt.
greenprint('\nProgram started.')
cyanprint('''\nEnter parent folder location/path where files are located:
* This will also be the location for the output files.''')
while True:
    path = input()
    if os.path.exists(path):
        path = path.replace('\\',
                            '/'
                            )
        os.chdir(path)
        if len(glob.glob('*.seq')) > 0:
            break
        # Redirect user to input the path again if no seq files are found.
        else:
            cyanprint('''\nInvalid input.
No seq files were found in this directory
Please try again.'''
                      )
            continue
    # Redirect user to input the path again if path doesn't exist.
    else:
        cyanprint('''\nInvalid input.
The entered path does not exist.
Please try again.'''
                  )
        continue

# Logging setup.
pathRegex = re.compile('[^/]*$')
pathList = pathRegex.findall(path)
folderName = str(pathList[0])
logging.basicConfig(filename=path + '/' + folderName + '.log',
                    level=logging.INFO,
                    format='%(asctime)s - %(message)s',
                    filemode='w'
                    )
logging.info('Working directory changed to %s.' % path)

##################
# Rename seq files to remove well designations used by the AAC, add a leading zero to well IDs, and move to raw
# sequences folder.
##################

# Create raw sequences folder.
rawSeqFolder = '01_rawSeq'
rawSeqPath = path + '/' + rawSeqFolder
if not os.path.exists(rawSeqPath):
    os.makedirs(rawSeqPath)
    logging.info('%s directory created.' % rawSeqFolder)
for fileName in glob.glob('*.seq'):
    # Remove AAC well designation.
    shortName = re.sub(r'(_[A-H]\d\d[.])',
                       '.',
                       fileName
                       )
    # Add leading zero to well IDs.
    if len(re.findall(r'(\d{2,}_\s*\w*)', shortName)) < 1:
        newName = re.sub(r'(\d_\s*\w*)',
                         r'0\g<1>',
                         shortName
                         )
    else:
        newName = shortName
    # Move to raw sequences folder.
    shutil.move(path + '/' + fileName,
                rawSeqPath + '/' + newName
                )
logging.info('Sequences renamed and moved to %s.' % rawSeqFolder)

##################
# Rename ab1 files to remove well designations used by the AAC, add a leading zero to well IDs, and move to raw ab1
# folder. If no ab1 files are found, a text file stating such will be created.
##################

# TODO: Analyse edited ab1 files.
# Create raw ab1 folder.
rawAb1Folder = '02_rawAb1'
rawAb1Path = path + '/' + rawAb1Folder
if not os.path.exists(rawAb1Path):
    os.makedirs(rawAb1Path)
    logging.info('%s directory created.' % rawAb1Folder)
if len(glob.glob('*.ab1')) > 0:
    for fileName in glob.glob('*.ab1'):
        # Remove AAC well designation.
        shortName = re.sub(r'(_[A-H]\d\d[.])',
                           '.',
                           fileName
                           )
        # Add leading zero to well IDs.
        if len(re.findall(r'(\d{2,}_\s*\w*)', shortName)) < 1:
            newName = re.sub(r'(\d_\s*\w*)',
                             r'0\g<1>',
                             shortName
                             )
        else:
            newName = shortName
        # Move to raw ab1 folder.
        shutil.move(path + '/' + fileName,
                    rawAb1Path + '/' + newName
                    )
    logging.info('Chromatograms renamed and moved to %s.' % rawAb1Path)
else:
    # Create a text file stating that no ab1 files were found.
    noAb1FileName = 'noAb1FilesFound.txt'
    noAb1File = open(noAb1FileName, 'x')
    noAb1File.write('No ab1 files found.')
    noAb1File.close()
    shutil.move(path + '/' + noAb1FileName,
                rawAb1Path + '/' + noAb1FileName
                )
    logging.info('No ab1 files found so %s was created.' % noAb1FileName)
greenprint('\nFile renaming finished.')
logging.info('File renaming finished.')

##################
# Convert seq files to fasta, move to nucleotide fasta folder, and combine all sequences into a single batch file.
##################

# Create nucleotide fasta folder.
ntFastaFolder = '03_ntFasta'
ntFastaPath = path + '/' + ntFastaFolder
if not os.path.exists(ntFastaFolder):
    os.makedirs(ntFastaFolder)
    logging.info('%s directory created.' % ntFastaFolder)
os.chdir(rawSeqPath)
for fileName in glob.glob('*.seq'):
    with open(fileName, 'r') as preConvFile:
        # Retrieve sequence and remove unnecessary formatting.
        rawSeq = preConvFile.read().replace('\n',
                                            ''
                                            )
        # Retrieve sequence name to write in file.
        seqName = fileName.replace('.seq',
                                   ''
                                   )
        # Assemble name and sequence in fasta format.
        fastaSeq = '>' + seqName + '\n' + rawSeq + '\n'
        # Write output fasta file.
        outputFileName = seqName + '.fasta'
        with open(outputFileName, 'w') as postConvFile:
            postConvFile.write(fastaSeq)
        preConvFile.close()
greenprint('\nSeq to fasta conversion finished.')
logging.info('Seq to fasta conversion finished.')

# Move fasta files to nucleotide fasta folder.
for fileName in glob.glob('*.fasta'):
    shutil.move(rawSeqPath + '/' + fileName,
                ntFastaPath + '/' + fileName)
logging.info('Fasta files moved to %s.' % ntFastaFolder)

# Combine all fasta sequences together into a single batch file.
os.chdir(ntFastaPath)
fastaBatchName = folderName + '_batch.fasta'
fastaFiles = glob.glob('*.fasta')
with open(fastaBatchName, 'w') as batchFile:
    # Write the content of each fasta file to the new batch file.
    for file in fastaFiles:
        with open(file) as inFile:
            batchFile.write(inFile.read() + '\n')
            inFile.close()
greenprint('\nFasta batch file created.')
logging.info('Fasta batch file created.')

##################
# Trim fasta sequences at a 5' motif.
##################

# User prompt.
cyanprint('''\nEnter the sequence to begin the 5' trim at:
* Not case-sensitive.
* Must be at least six nucleotides long.
* Use conserved nucleotides to prevent trimming at multiple sites or lack of trimming.
* Do not enter nucleotides that contribute to the codons of diversified residues, where possible. 
    E.g. For UbVs, the three nucleotides succeeding the start codon are diversified, so use the preceding three '''
          '''nucleotides to ensure specificity (i.e. AAA ATG; the last FLAG tag codon & the UbV start codon).'''
          )

# Create no trim folder.
noTrimFolder = 'noTrim'
noTrimPath = path + '/' + noTrimFolder
if not os.path.exists(noTrimPath):
    os.makedirs(noTrimPath)
    logging.info('%s directory created.' % noTrimFolder)

fivePrimeRegex = re.compile('[agtcurynwsmkbhdv]{6,}',
                            re.IGNORECASE)

while True:
    startSite = input()
    startSite = startSite.upper()
    logging.info('%s chosen as the start site sequence.' % startSite)
    # Find the first instance of the 5' motif.
    if fivePrimeRegex.match(startSite):
        fiveFileList = glob.glob('*.fasta')
        for fileName in fiveFileList:
            # Change back to the directory containing the files.
            os.chdir(ntFastaPath)
            index = 0
            with open(fileName, 'r') as tempFile:
                # Retrieve ID and sequence data.
                lines = tempFile.readlines()
                seq = str(lines[1])
                index = seq.find(startSite,
                                 index
                                 )
                # If the motif can be found, trim the sequence and write it to a new file.
                if index != -1:
                    trimmedSeq = seq[index:]
                    trimName = fileName.replace('.fasta',
                                                '_ntTrimmed.fasta'
                                                )
                    with open(trimName, 'w') as trimFile:
                        trimFasta = trimName.replace('.fasta',
                                                     ''
                                                     )
                        trimFile.write('>' + trimFasta + '\n' + trimmedSeq)
                    index += 1
                # If the motif cannot be found in the sequence, move the sequence to the untrimmed folder.
                elif index == -1:
                    os.chdir(noTrimPath)
                    untrimName = fileName.replace('.fasta',
                                                  '_5trimfail.fasta'
                                                  )
                    with open(untrimName, 'w') as unTrimFile:
                        untrimFasta = untrimName.replace('.fasta',
                                                         ''
                                                         )
                        unTrimFile.write('>' + untrimFasta + '\n' + seq)
        break
    # Redirect user to input the motif again if regex doesn't match anything.
    else:
        cyanprint('''\nInvalid input.
Make sure the sequence is typed correctly, is at least six nucleotides long, and contains only IUPAC approved letters.
Please try again.'''
                  )
greenprint('''\n5' trim finished.''')
logging.info('''5' trim finished.''')

##################
# Create trimmed nucleotide folder, trim fasta sequences at the 3' end by length, and move trimmed sequences to trimmed
# nucleotide folder.
##################

# Trim 3' end of sequence. User prompt.
cyanprint('''\nEnter the number of nucleotides downstream of the 5' trim site to trim at:
* Must be at least ten nucleotides long.
  E.g. For UbVs, use 237 (nucleotides downstream of the 5' trim site).'''
          )

# Create trimmed nucleotide folder.
ntTrimmedFolder = '04_ntTrimmed'
ntTrimmedPath = path + '/' + ntTrimmedFolder
if not os.path.exists(ntTrimmedPath):
    os.makedirs(ntTrimmedPath)
    logging.info('%s directory created.' % ntTrimmedFolder)

threePrimeRegex = re.compile('[0-9]{2,}')

while True:
    endSite = input()
    logging.info('%s chosen as the sequence length.' % endSite)
    if threePrimeRegex.match(endSite):
        endSite = int(endSite)
        os.chdir(ntFastaPath)
        fiveFileList = glob.glob('*trimmed.fasta')
        for fileName in fiveFileList:
            os.chdir(ntFastaPath)
            with open(fileName, 'r') as tempFile1:
                # Retrieve ID and sequence data.
                lines = tempFile1.readlines()
                fastaName = str(lines[0])
                seq = str(lines[1])
                trimSeq = seq[0:endSite]
                tempFile1.close()
                # If the sequence cannot be trimmed at the input length, move the sequence to the untrimmed folder.
                if len(trimSeq) < endSite:
                    if not os.path.exists(noTrimPath):
                        os.remove(fileName)
                        os.makedirs(noTrimPath)
                        logging.info('%s directory created.' % noTrimFolder)
                    os.chdir(noTrimPath)
                    untrimName = fileName.replace('_ntTrimmed.fasta',
                                                  '_3trimfail.fasta'
                                                  )
                    with open(untrimName, 'w') as untrimFile:
                        newFasta = untrimName.replace('.fasta',
                                                      ''
                                                      )
                        untrimFile.write('>' + newFasta + '\n' + seq)
                    os.chdir(ntFastaPath)
                    os.remove(fileName)
                # If the sequence can be trimmed at the input length, trim the sequence, write it to the trimmed file,
                # and move the sequence to the trimmed nucleotide folder.
                else:
                    with open(fileName, 'r') as tempFile2:
                        lines = tempFile2.readlines()
                        fastaName = str(lines[0])
                        seq = str(lines[1])
                        trimSeq = seq[0:endSite]
                        with open(fileName, 'w') as trimFile:
                            trimFile.write(fastaName + trimSeq)
                            trimFile.close()
                        tempFile2.close()
                        shutil.move(ntFastaPath + '/' + fileName,
                                    ntTrimmedPath + '/' + fileName
                                    )
        break
    # Redirect user to input the length again if regex doesn't match anything.
    else:
        cyanprint('''\nInvalid input.
Make sure the number is typed correctly, greater than or equal to ten, and doesn't contain non-numeric characters.
Please try again.'''
                  )
greenprint('''\n3' trim finished.''')
logging.info('''3' trim finished.''')

##################
# Delete faulty/incorrect files and combine all trimmed sequences into a single batch file.
##################

# Delete unnecessary/faulty files.
os.chdir(ntTrimmedPath)
faultyNtTrimBatchNames = [ntTrimmedPath + '/' + folderName + '_batch_ntTrimmed.fasta',
                          noTrimPath + '/' + folderName + '_batch_5trimfail.fasta',
                          noTrimPath + '/' + folderName + '_batch_3trimfail.fasta'
                          ]
for faultyFiles in faultyNtTrimBatchNames:
    if os.path.exists(faultyFiles):
        os.remove(faultyFiles)
        logging.info('Deleted %s (faulty batch file).' % faultyFiles)
    else:
        continue

# Combine all trimmed nucleotide sequences together into a single batch file.
correctNtTrimBatchName = folderName + '_ntTrimmed_batch.fasta'
correctNtTrimBatchPath = ntTrimmedPath + '/' + correctNtTrimBatchName
ntTrimFileNames = glob.glob('*trimmed.fasta')
with open(correctNtTrimBatchName, 'w') as ntTrimBatchFile:
    # Write the content of each trimmed file to the new batch file.
    for file in ntTrimFileNames:
        with open(file) as ntTrimFile:
            ntTrimBatchFile.write(ntTrimFile.read() + '\n\n')
greenprint('\nTrimmed nucleotide batch file created.')
logging.info('Trimmed nucleotide batch file created.')

##################
# Create trimmed amino acids folder, make copies of the trimmed nucleotides in this folder, and overwrite the file
# contents with the translated translated sequences.
##################

# Create a folder for trimmed amino acid sequences.
os.chdir(ntTrimmedPath)
aaTrimmedFolder = '05_aaTrimmed'
aaTrimmedPath = path + '/' + aaTrimmedFolder
if not os.path.exists(aaTrimmedPath):
    os.makedirs(aaTrimmedPath)
    logging.info('%s directory created.' % aaTrimmedFolder)

# Make copies of trimmed nucleotide sequences in the trimmed amino acids folder.
ntTrimFileList = glob.glob('*.fasta')
for fileName in ntTrimFileList:
    copyfile(ntTrimmedPath + '/' + fileName,
             aaTrimmedPath + '/' + re.sub('_ntTrimmed',
                                          '_aaTrimmed',
                                          fileName)
             )
logging.info('Copies of nucleotide fasta files made in  %s.' % aaTrimmedFolder)

# Translate trimmed nucleotides to amino acids.
os.chdir(aaTrimmedPath)
aaTrimFileList = glob.glob('*.fasta')
for file in aaTrimFileList:
    with open(file, 'r+') as aaTrimFile:
        # Retrieve data.
        lines = aaTrimFile.readlines()
        fastaName = lines[0]
        fastaName = re.sub('_ntTrimmed',
                           '_aaTrimmed',
                           fastaName
                           )
        ntSeq = lines[1]
        # Remove previous file contents.
        aaTrimFile.seek(0)
        aaTrimFile.truncate()
        # Translate nucleotide codons to amino acids.
        gene = Seq(ntSeq)
        aaSeq = gene.translate(table='Standard',
                               stop_symbol='*'
                               )
        aaSeq = str(aaSeq)
        # Write new contents to same file.
        aaTrimFile.write(fastaName + aaSeq)
greenprint('''\nTranslation finished.''')
logging.info('Nucleotide sequences translated into amino acid sequences.')

##################
# Delete faulty/incorrect files and combine all translated sequences into a single batch file.
##################

# Delete unnecessary/faulty files.
faultyAaTrimBatchPath = aaTrimmedPath + '/' + folderName + '_aaTrimmed_batch.fasta'
if os.path.exists(faultyAaTrimBatchPath):
    os.remove(faultyAaTrimBatchPath)
    logging.info('Deleted %s (faulty batch file).' % faultyAaTrimBatchPath)
else:
    pass

# Combine all trimmed amino acid sequences together into a single batch file.
correctAaTrimBatchPath = faultyAaTrimBatchPath
aaTrimFileNames = glob.glob('*.fasta')
aaTrimData = ''
for file in aaTrimFileNames:
    with open(file) as aaTrimFile:
        aaTrimData += str(aaTrimFile.read() + '\n\n')
with open(correctAaTrimBatchPath, 'a') as aaTrimBatchFile:
    aaTrimBatchFile.write(aaTrimData)

##################
# Create a new batch file where truncated and non-full length amino acid sequences are removed.
##################

# Remove truncated amino acid sequences.
truncAaTrimBatchPath = re.sub('.fasta',
                              '_noTrunc.fasta',
                              correctAaTrimBatchPath
                              )
# with open(truncAaTrimBatchPath, 'x') as truncBatchFile:
with open(correctAaTrimBatchPath, 'r') as fullBatchFile:
    batchData = fullBatchFile.readlines()
seqBatchList = []
nameBatchList = []
# Retrieve amino acid sequences and names in batch file.
for entry in batchData:
    # Retrieve amino acid sequence.
    if re.match(r'^\w+', entry):
        seqBatchList += [entry]
    # Retrieve sequence name.
    else:
        nameBatchList += [entry]
# Remove unnecessary formatting from sequence list.
seqBatchList = [entry[:-1] for entry in seqBatchList]
# Remove unnecessary formatting from name list.
nameBatchList = [entry[:-1] for entry in nameBatchList]
# Remove blank entries from name list.
nameBatchList = list(filter(None, nameBatchList))
# Remove sequences that are not full length and put in a temporary dictionary.
maxLen = max([len(entry) for entry in seqBatchList])
aaAlignList = []
for entry in seqBatchList:
    # Put full length sequences in the alignment pool.
    if len(entry) == maxLen:
        aaAlignList.append(entry)
    # Skip sequences that are not full length.
    elif len(entry) < maxLen:
        pass
logging.info('Non-full length amino acid sequences removed.')
aaAlignDict = dict(zip(nameBatchList,
                       aaAlignList)
                   )
# Remove sequences with at least one stop codon and write all remaining sequences to a batch file.
with open(truncAaTrimBatchPath, 'x') as truncAaTrimBatchFile:
    for name, seq in aaAlignDict.items():
        if seq.find('*') == -1:
            truncAaTrimBatchFile.write(name + '\n' + seq + '\n\n')
        else:
            pass
logging.info('Non-full length and truncated amino acid sequences removed from alignment pool.')
greenprint('Non-full length and truncated amino acid sequences removed from alignment pool.')

##################
# Create clustal and fasta alignments for both nucleotide and amino acid sequences.
##################

# Create alignment folder.
alignmentFolder = '06_alignment'
alignmentPath = path + '/' + alignmentFolder
if not os.path.exists(alignmentPath):
    os.makedirs(alignmentPath)
    logging.info('%s directory created.' % alignmentFolder)
os.chdir(alignmentPath)

#########
# Amino acids
#########

# Align amino acid sequences.
aaAlignment = AlignIO.read(open(truncAaTrimBatchPath),
                           'fasta'
                           )
aaAlignLen = aaAlignment.get_alignment_length()
# Amino acid fasta alignment.
aaAlignFasta = format(aaAlignment,
                      'fasta'
                      )
aaFastaName = folderName + '_aaTrimmed_alignment.fasta'
with open(aaFastaName, 'w') as aaFastaFile:
    aaFastaFile.write(aaAlignFasta)
    logging.info('Amino acid fasta alignment created.')
# Amino acid clustal alignment.
aaAlignClustal = format(aaAlignment,
                        'clustal'
                        )
aaClustalName = folderName + '_aaTrimmed_alignment.clustal'
with open(aaClustalName, 'w') as aaClustalFile:
    aaClustalFile.write(aaAlignClustal)
    logging.info('Amino acid clustal alignment created.')

#########
# Nucleotides
#########

# Align nucleotide sequences.
ntAlignment = AlignIO.read(open(correctNtTrimBatchPath),
                           'fasta'
                           )
ntAlignLen = ntAlignment.get_alignment_length()
# Nucleotide fasta alignment.
ntAlignFasta = format(ntAlignment,
                      'fasta'
                      )
ntFastaName = folderName + '_ntTrimmed_alignment.fasta'
with open(ntFastaName, 'w') as ntFastaFile:
    ntFastaFile.write(ntAlignFasta)
    logging.info('Nucleotide fasta alignment created.')
# Nucleotide clustal alignment.
ntAlignClustal = format(ntAlignment,
                        'clustal'
                        )
ntClustalName = folderName + '_ntTrimmed_alignment.clustal'
with open(ntClustalName, 'w') as ntClustalFile:
    ntClustalFile.write(ntAlignClustal)
    logging.info('Nucleotide clustal alignment created.')

##################
# Extract data for writing to xlsx file.
##################

#########
# Amino acids
#########

# Amino acid alignment data extraction.
fastaNameRegex = re.compile('>(.*)')
seqRegex = re.compile(r'(?<!>)([A-Z]{5,})(?!\d|[a-z]|_)')
stopCodonRegex = re.compile('([*]+[A-Z]*)')
nameTrimRegex = re.compile(r'(_M\w*)')
alignSource = folderName + '_aaTrimmed_alignment.fasta'
with open(alignSource, 'r') as alignFile:
    allDataAa = alignFile.read()
    # Remove '_aaTrimmed' from names.
    shortNamesAa = re.sub('_aaTrimmed',
                          '',
                          allDataAa)
    # Retrieve names.
    aaNameList = fastaNameRegex.findall(shortNamesAa)
    # Remove text formatting.
    formattedDataAa = allDataAa.replace('\n',
                                        ''
                                        )
    # Remove truncated sequences.
    noTruncSeq = stopCodonRegex.sub('',
                                    formattedDataAa
                                    )
    # Retrieve sequences.
    aaList = seqRegex.findall(noTruncSeq)
    # Remove '_<primer>_aaTrimmed' from names.
    shorterNamesAa = nameTrimRegex.sub('',
                                       allDataAa
                                       )
    # Retrieve well IDs.
    shortWellListAa = fastaNameRegex.findall(shorterNamesAa)

# Create list of unique amino acid sequences ordered by frequency.
aaUnique = OrderedCounter(aaList)
aaUnique = aaUnique.most_common()
aaUniqueDict = dict(aaUnique)
logging.info('Dictionary of unique amino acid sequences created.')

#########
# Nucleotides
#########

# Nucleotide alignment data extraction.
alignFile2 = folderName + '_ntTrimmed_alignment.fasta'
with open(alignFile2, 'r') as alignFile:
    allDataNt = alignFile.read()
    # Remove '_ntTrimmed' from names.
    shortNamesNt = re.sub('_ntTrimmed',
                          '',
                          allDataNt
                          )
    # Retrieve names.
    nameListNt = fastaNameRegex.findall(shortNamesNt)
    # Remove text formatting.
    formattedDataNt = allDataNt.replace('\n',
                                        ''
                                        )
    ntList = seqRegex.findall(formattedDataNt)
    # Remove '_<primer>_ntTrimmed' from names.
    shorterNamesNt = nameTrimRegex.sub('',
                                       allDataNt
                                       )
    # Retrieve well IDs.
    shortWellListNt = fastaNameRegex.findall(shorterNamesNt)

# Create list of unique nucleotide sequences ordered by frequency.
uniqueNt = OrderedCounter(ntList)
uniqueNt = uniqueNt.most_common()
ntUniqueDict = dict(uniqueNt)
logging.info('Dictionary of unique nucleotide sequences created.')

##################
# Order sequences and wells so they can be attributed to unique sequences. Necessary for subsequent statistics.
##################

#########
# Amino acids
#########

# Get ordered list of wells that correspond to unique sequences; necessary for assigning wells to unique IDs.
aaOrderedSeq = [key for key in aaUniqueDict.keys()]
aaDict = dict(zip(shortWellListAa,
                  aaList)
              )
aaOrderedIndex = []
for seq in aaOrderedSeq:
    for key, value in aaDict.items():
        if seq in value:
            aaOrderedIndex.append(key)
logging.info('Ordered index of amino acid well IDs created.')

# Associate specific well IDs with corresponding unique sequence (amino acids).
aaCountID = []
begin = 0
for uniqueSeq, count in aaUniqueDict.items():
    end = int(count) + begin
    aaCountID.append(aaOrderedIndex[begin:end])
    begin += count
logging.info('List of specific well IDs associated with amino acid sequences created.')

#########
# Nucleotides
#########

# Get ordered list of wells that correspond to unique sequences; necessary for attributing wells to unique sequences.
newNtDict = dict(zip(shortWellListNt,
                     ntList)
                 )
ntOrderedSeq = [key for key in ntUniqueDict.keys()]
ntOrderedIndex = []
for seq in ntOrderedSeq:
    for key, value in newNtDict.items():
        if seq in value:
            ntOrderedIndex.append(key)
logging.info('Ordered index of nucleotide well IDs created.')
# Associate specific well IDs with corresponding unique sequence (nucleotides).
ntCountID = []
begin = 0
for uniqueSeq, count in ntUniqueDict.items():
    end = int(count) + begin
    ntCountID.append(ntOrderedIndex[begin:end])
    begin += count

##################
# Export alignments as a single xlsx file.
##################

# Create workbook.
workbook = xlsxwriter.Workbook(alignmentPath + '/' + folderName + '_alignment.xlsx')

# TODO: Find a way to clean up this section's formatting.
#########
# Cell formatting rules.
#########

# General.
general_format = workbook.add_format()
general_format.set_align('center')
general_format.set_align('vcenter')
# Titles.
title_format = workbook.add_format({'bold': True,
                                    'font_size': 12
                                    }
                                   )
title_format.set_align('center')
title_format.set_align('vcenter')
wellTitle_format = workbook.add_format({'bold': True,
                                        'font_size': 12
                                        }
                                       )
wellTitle_format.set_align('left')
wellTitle_format.set_align('vcenter')
# Wells.
wellList_format = workbook.add_format({'font_size': 11})
wellID_format = workbook.add_format({'font_size': 12})
wellID_format = workbook.add_format({'font_size': 12})
wellID_format.set_align('center')
wellID_format.set_align('vcenter')
# Residue numbers.
residue_format = workbook.add_format({'font_size': 10})
residue_format.set_align('center')
residue_format.set_align('vcenter')
# Sequences.
sequence_format = workbook.add_format({'font_size': 10})
sequence_format.set_align('center')
sequence_format.set_align('vcenter')
sequence_format.set_font_name('Lucida Console')
logging.info('Cell formatting rules set.')

##################
# Create worksheet for all amino acid sequences.
##################

# Initialize sheet characteristics.
worksheet1Name = 'All AA Seq'
worksheet1 = workbook.add_worksheet(worksheet1Name)
worksheet1.hide_gridlines(option=2)
idColWidth = round(len(aaNameList[0]) * 1.4)
worksheet1.set_column(0, 0, idColWidth)
worksheet1.set_column(1, aaAlignLen, 3)
worksheet1.freeze_panes(0, 1)
logging.info('%s worksheet created.' % worksheet1Name)

# Write well IDs.
worksheet1.merge_range(0, 0, 1, 0, 'ID', title_format)
idRow = 2
for name in aaNameList:
    worksheet1.write(idRow, 0, name, wellID_format)
    idRow += 1
logging.info('Amino acid well IDs written to %s worksheet.' % worksheet1Name)

# Write all amino acid sequences.
worksheet1.merge_range(0, 1, 0, 9, 'Amino Acid Sequence', title_format)
seqRow = 2
seqCol = 1
for aa in aaList:
    letterList = list(aa)
    for letter in letterList:
        worksheet1.write(seqRow, seqCol, letter, sequence_format)
        seqCol += 1
    seqRow += 1
    seqCol = 1
logging.info('All Amino acid sequences written to %s worksheet.' % worksheet1Name)

# Write amino acid residue numbers above sequences.
aaResList = list(range(1,
                       aaAlignLen + 1)
                 )
residueCol = 1
for number in aaResList:
    worksheet1.write(1, residueCol, number, residue_format)
    residueCol += 1
logging.info('Residue numbers written to %s worksheet.' % worksheet1Name)

##################
# Create worksheet for unique amino acid sequences.
##################

# Initialize sheet characteristics.
worksheet2Name = 'Unique AA Seq'
worksheet2 = workbook.add_worksheet(worksheet2Name)
worksheet2.hide_gridlines(option=2)
worksheet2.set_column(0, 0, 10)
worksheet2.set_column(1, aaAlignLen, 2)
worksheet2.set_column(1, aaAlignLen, 2)
worksheet2.freeze_panes(2, 1)
logging.info('%s worksheet created.' % worksheet2Name)

# Write unique amino acid sequences.
worksheet2.merge_range(0, 1, 0, 9, 'Amino Acid Sequence', title_format)
seqRow = 2
seqCol = 1
for seq in aaUniqueDict.keys():
    letterList = list(seq)
    for letter in letterList:
        worksheet2.write(seqRow, seqCol, letter, sequence_format)
        seqCol += 1
    seqRow += 1
    seqCol = 1
logging.info('Unique amino acid sequences written to %s worksheet.' % worksheet2Name)

# Write sequence counts.
count = list(aaUniqueDict.values())
countRow = 2
countCol = aaAlignLen + 1
worksheet2.merge_range(0, aaAlignLen + 1, 1, aaAlignLen + 1, 'Count', title_format)
for number in count:
    worksheet2.write_number(countRow, countCol, number, general_format)
    countRow += 1
logging.info('Amino acid sequence counts written to %s worksheet.' % worksheet2Name)

# Assign arbitrary IDs to each unique amino acid sequence.
worksheet2.merge_range(0, 0, 1, 0, 'ID', title_format)
aaIdList = list(range(1,
                      len(aaUnique) + 1)
                )
idRow = 2
for number in aaIdList:
    worksheet2.write_number(idRow, 0, number, general_format)
    idRow += 1
logging.info('Arbitrary unique amino acid sequence IDs written to %s worksheet.' % worksheet2Name)

# Write amino acid residue numbers above sequences.
residueCol = 1
for number in aaResList:
    worksheet2.write(1, residueCol, number, residue_format)
    residueCol += 1
logging.info('Residue numbers written to %s worksheet.' % worksheet2Name)

# Change column width to fit all IDs.
wellColWidth = round((len(aaCountID[0]) * 3) * 1.4)
worksheet2.set_column(aaAlignLen + 2, aaAlignLen + 2, wellColWidth)
# Write specific IDs to worksheet.
worksheet2.merge_range(0, aaAlignLen + 2, 1, aaAlignLen + 2, 'Wells', wellTitle_format)
wellRow = 2
wellCol = aaAlignLen + 2
wellRegex = re.compile(r'([A-Z][0-1][0-9])')
sep = ', '
for wellList in aaCountID:
    wellList = wellRegex.findall(str(wellList))
    wellList = sep.join(wellList)
    worksheet2.write(wellRow, wellCol, wellList, wellList_format)
    wellRow += 1
logging.info('Amino acid sequence-specific well IDs written to %s worksheet.' % worksheet2Name)

##################
# Create worksheet for all nucleotide sequences.
##################

# Initialize sheet characteristics.
worksheet3Name = 'All NT Seq'
worksheet3 = workbook.add_worksheet(worksheet3Name)
worksheet3.hide_gridlines(option=2)
worksheet3.set_column(0, 0, idColWidth)
worksheet3.set_column(1, ntAlignLen, 3)
worksheet3.freeze_panes(2, 1)
logging.info('%s worksheet created.' % worksheet3Name)

# Write well IDs.
worksheet3.merge_range(0, 0, 1, 0, 'ID', title_format)
idRow = 2
for name in nameListNt:
    worksheet3.write(idRow, 0, name, wellID_format)
    idRow += 1
logging.info('Nucleotide well IDs written to %s worksheet.' % worksheet3Name)

# Write all nucleotide sequences.
worksheet3.merge_range(0, 1, 0, 7, 'Nucleotide Sequence', title_format)
seqRow = 2
seqCol = 1
for nt in ntList:
    letterList = list(nt)
    for letter in letterList:
        worksheet3.write(seqRow, seqCol, letter, sequence_format)
        seqCol += 1
    seqRow += 1
    seqCol = 1
logging.info('All nucleotide sequences written to %s worksheet.' % worksheet3Name)

# Write nucleotide base pair numbers above sequences.
ntBpList = list(range(1,
                      ntAlignLen + 1)
                )
bpCol = 1
for number in ntBpList:
    worksheet3.write(1, bpCol, number, residue_format)
    bpCol += 1
logging.info('Base pair numbers written to %s worksheet.' % worksheet3Name)

##################
# Create worksheet for unique nucleotide sequences.
##################

# Initialize sheet characteristics.
worksheet4Name = 'Unique NT Seq'
worksheet4 = workbook.add_worksheet(worksheet4Name)
worksheet4.hide_gridlines(option=2)
worksheet4.set_column(0, 0, 10)
worksheet4.set_column(1, ntAlignLen, 3)
worksheet4.freeze_panes(2, 1)
logging.info('%s worksheet created.' % worksheet4Name)

# Write unique nucleotide sequences.
worksheet4.merge_range(0, 1, 0, 7, 'Nucleotide Sequence', title_format)
seqRow = 2
seqCol = 1
for seq in ntUniqueDict.keys():
    letterList = list(seq)
    for letter in letterList:
        worksheet4.write(seqRow, seqCol, letter, sequence_format)
        seqCol += 1
    seqRow += 1
    seqCol = 1
logging.info('Unique nucleotide sequences written to %s worksheet.' % worksheet4Name)

# Write sequence counts.
count = list(ntUniqueDict.values())
countRow = 2
countCol = ntAlignLen + 1
worksheet4.merge_range(0, ntAlignLen + 1, 1, ntAlignLen + 1, 'Count', title_format)
for number in count:
    worksheet4.write_number(countRow, countCol, number, general_format)
    countRow += 1
logging.info('Nucleotide sequence counts written to %s worksheet.' % worksheet4Name)

# Assign arbitrary IDs to each unique nucleotide sequence.
worksheet4.merge_range(0, 0, 1, 0, 'ID', title_format)
ntIdList = list(range(1,
                      len(uniqueNt) + 1)
                )
idRow = 2
for number in ntIdList:
    worksheet4.write_number(idRow, 0, number, general_format)
    idRow += 1
logging.info('Arbitrary unique nucleotide sequence IDs written to %s worksheet.' % worksheet4Name)

# Write nucleotide base pair numbers above sequences.
bpCol = 1
for number in ntBpList:
    worksheet4.write(1, bpCol, number, residue_format)
    bpCol += 1

# Change column width to fit all IDs.
wellColWidth = round((len(ntCountID[0]) * 3) * 1.4)
worksheet4.set_column(ntAlignLen + 2, ntAlignLen + 2, wellColWidth)
# Write specific IDs to worksheet.
worksheet4.write(0, ntAlignLen + 2, 'Wells', wellTitle_format)
wellRow = 2
wellCol = ntAlignLen + 2
sep = ', '
for wellList in ntCountID:
    wellList = wellRegex.findall(str(wellList))
    wellList = sep.join(wellList)
    worksheet4.write(wellRow, wellCol, wellList, wellList_format)
    wellRow += 1
logging.info('Well IDs grouped to unique nucleotide sequences and written to %s worksheet.' % worksheet4Name)

##################
# Transform data in every sheet into proper Excel-formatted tables without any design style applied.
##################

worksheet1.add_table(2, 0, len(aaNameList) + 1, aaAlignLen,
                     {'header_row': False,
                      'style': None
                      }
                     )
worksheet2.add_table(2, 0, len(aaUnique) + 1, aaAlignLen + 2,
                     {'header_row': False,
                      'style': None
                      }
                     )
worksheet3.add_table(2, 0, len(nameListNt) + 1, ntAlignLen,
                     {'header_row': False,
                      'style': None
                      }
                     )
worksheet4.add_table(2, 0, len(uniqueNt) + 1, ntAlignLen + 2,
                     {'header_row': False,
                      'style': None
                      }
                     )

##################
# Close workbook and logging and give post-analysis advice.
##################

workbook.close()
greenprint('\nExcel alignment saved as %s_alignment.xlsx.' % folderName)
logging.info('Excel alignment exported as %s_alignment.xlsx.' % folderName)
cyanprint('''\nSequence Analysis program finished running, see log file for details.
\n\nPost-analysis help:
\nNon-trimmed files are in the 'noTrim' folder and couldn't be trimmed because of one of the following reasons:
\n      a) There's no sequence.
In this case, either the sequence was not good enough quality/quantity to be sequenced or no sequence was 
present in the sample.
\n      b) The 5' trim site couldn't be found.       
In this case, look at the corresponding ab1 file and assess whether or not the bases were called correctly. If
not, create a copy of the ab1 file, edit the base calls (type manually called bases in lowercase to distinguish
them from bases originally called by the sequencing equipment), move the original unedited ab1 file to a separate'''
          '''folder, replace its original position with the new edited ab1 file, and run this script again.'''
          )
logging.info('Sequence Analysis program finished running.')
logging.shutdown()
