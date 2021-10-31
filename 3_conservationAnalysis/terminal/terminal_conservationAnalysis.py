#! python3

##################
#    MODULES
##################

import os
import re
import logging
import xlsxwriter
import pandas
from collections import Counter, OrderedDict


##################
#    FUNCTIONS
##################

def cyanprint(text):
    """Print cyan coloured text in the console."""
    print('\033[0;36m' + text + '\033[0m')


def greenprint(text):
    """Print green coloured text in the console."""
    print('\033[0;32m' + text + '\033[0m')


def comparestrings(a, b):
    """Return a string containing only the characters in string b that differ from string a."""
    diff = ''
    for x, y in zip(a, b):
        if x == y:
            diff += '-'
        else:
            diff += y
    return diff


##################
#    CLASSES
##################

# Ordered list of counts.
class OrderedCounter(Counter, OrderedDict):
    pass


##################
#    MAIN
##################

##################
# Set up working directory.
##################

greenprint('\nProgram started.')
cyanprint('''\n\nEnter parent folder location/path where files are located:
* This will also be the location for the output files.'''
          )
while True:
    path = input()
    if os.path.exists(path):
        path = path.replace('\\',
                            '/'
                            )
        os.chdir(path)
        break
    else:
        cyanprint('''\nInvalid input.
The entered path does not exist.
Please try again.'''
                  )

##################
# Select data input format.
##################

# Choose whether to analyse ELISA and sequencing data or just sequencing data. User prompt.
cyanprint('''\nChoose input format by typing the corresponding number:
\n[1] ELISA and sequencing data
* Requires xlsx output from Binding Analysis program.
\n[2] Sequencing data only
* Requires amino acid alignment fasta output from Sequence Analysis program.'''
          )
inputOptions = {'1': 'ELISA and sequencing',
                '2': 'sequencing'
                }
while True:
    inputFormat = input()
    if inputFormat in inputOptions:
        cyanprint('\nOption %s chosen, analysing %s data.' % (inputFormat,
                                                              inputOptions[inputFormat]
                                                              )
                  )
        break
    else:
        cyanprint('\nInvalid option.'
                  'Please try again.')

##################
# Select input files.
##################
# TODO: Put while loop here to catch mistakes.
# ELISA input format.
if inputFormat == '1':
    # Select ELISA data source. User prompt.
    cyanprint('''\nEnter the analysed ELISA data file name:
* Must be in xlsx format.
* Include the file extension in the name.'''
              )
    inFileName = input()

# Alignment input format.
elif inputFormat == '2':
    # Select amino acid alignment file. User prompt.
    cyanprint('''\nEnter amino acid alignment file name:
    * Must be in fasta format.
    * Include the file extension in the name.'''
              )
    inFileName = input()

##################
# Setup logging file.
##################

inFileNameShort = re.sub(r'_a.*[.].*',
                         '_conservation',
                         inFileName
                         )

logging.basicConfig(filename=path + '/' + inFileNameShort + '.log',
                    level=logging.INFO,
                    format='%(asctime)s - %(message)s',
                    filemode='w'
                    )

logging.info('Working directory changed to %s.' % path)
logging.info('Option %s chosen, analysing %s data.' % (inputFormat,
                                                       inputOptions[inputFormat]
                                                       )
             )
logging.info('%s chosen as the data source.' % inFileName)

##################
# Retrieve and parse data.
##################

if inputFormat == '1':
    # Read ELISA file.
    allCells = pandas.read_excel(inFileName,
                                 sheet_name=1
                                 )
    # Remove rows with more than eight NaN values.
    allCells = allCells.dropna(axis=0, thresh=8)
    logging.info('%s data read.' % inFileNameShort)

    # Retrieve statistical data.
    countListFloat = list(allCells['Count'])
    countListFloat = countListFloat[1:]
    countList = []
    for count in countListFloat:
        countList.append(int(count))
    logging.info('Count values extracted from %s.' % inFileNameShort)
    maxList = list(allCells['Max.'])
    maxList = maxList[1:]
    logging.info('Maximum values extracted from %s.' % inFileNameShort)
    minList = list(allCells['Min.'])
    minList = minList[1:]
    logging.info('Minimum values extracted from %s.' % inFileNameShort)
    medianList = list(allCells['Median'])
    medianList = medianList[1:]
    logging.info('Median values extracted from %s.' % inFileNameShort)
    meanList = list(allCells['Mean'])
    meanList = meanList[1:]
    logging.info('Mean values extracted from %s.' % inFileNameShort)
    devList = list(allCells['St. Dev.'])
    devList = devList[1:]
    logging.info('Standard deviation values extracted from %s.' % inFileNameShort)

    # Retrieve well data.
    countID = list(allCells['Wells'])
    countID = countID[1:]
    logging.info('Wells extracted from %s.' % inFileNameShort)

    # Retrieve amino acid sequences from ELISA file.
    seqCells = allCells.iloc[1:, 1:]
    seqCells = seqCells.to_string(index=False)
    aaList = seqCells.replace(' ', '')
    seqRegex = re.compile(r'[ARNDCEQGHILKMFPSTWYVX]{10,}')
    aaList = seqRegex.findall(aaList)

    # Create list of unique nucleotide sequences ordered by frequency.
    uniqueDict = dict(zip(aaList, countList))

elif inputFormat == '2':
    # Read alignment file.
    seqRegex = re.compile(r'[ARNDCEQGHILKMFPSTWYVX]{10,}')
    stopRegex = re.compile(r'([*]+[A-Z]*)')
    with open(inFileName, 'r') as alignFile:
        allData = alignFile.read()
        # Retrieve amino acid sequences.
        seqClean = allData.replace('\n',
                                   ''
                                   )
        seqClean = stopRegex.sub('',
                                 seqClean
                                 )
        aaList = seqRegex.findall(seqClean)
        logging.info('Amino acid sequences retrieved from %s.' % inFileName)
        alignFile.close()

    # Retrieve well data.
    wellList = re.findall(r'([A-H][0-1][0-9])',
                          allData
                          )
    # Create list of unique nucleotide sequences ordered by frequency.
    unique = OrderedCounter(aaList)
    unique = unique.most_common()
    uniqueDict = dict(unique)
    # Create ordered list of wells.
    orderedSeq = []
    for key in uniqueDict.keys():
        orderedSeq.append(key)
    aaDict = dict(zip(wellList,
                      aaList)
                  )
    orderedIndex = []
    for seq in orderedSeq:
        for key, value in aaDict.items():
            if seq in value:
                orderedIndex.append(key)
    logging.info('Ordered index of amino acid well IDs created.')
    wellList = []
    begin = 0
    for uniqueSeq, count in uniqueDict.items():
        end = int(count) + begin
        wellList.append(orderedIndex[begin:end])
        begin += count
    sep = ', '
    countID = []
    for wells in wellList:
        wellList = sep.join(wells)
        countID.append(wellList)
    logging.info('List of specific well IDs associated with amino acid sequences created.')
    # Create new list of unique amino acids.
    aaList = []
    for uniqueSeq, count in uniqueDict.items():
        aaList.append(uniqueSeq)

##################
# For each amino acid sequence, replace non-diversified regions with dashes.
##################

# Remove amino acid prior to start codon.
aaShortList = []
for seq in aaList:
    aaShortSeq = seq[1:]
    aaShortList.append(aaShortSeq)
logging.info('Initial non-ubiquitin amino acid residue removed from all sequences.')

# Compare UbV sequences against a consensus sequence and replace conserved amino acids with dashes.
consensusSeq = 'MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGGGG'
logging.info('''Consensus sequence set as '%s'.''' % consensusSeq)
consensusLen = len(consensusSeq)
conservedList = []
for ubvSeq in aaShortList:
    conservedList.append(comparestrings(consensusSeq,
                                        ubvSeq
                                        )
                         )
logging.info('List of conserved sequences created.')

##################
# Create lists of arbitrary IDs and residue numbers.
##################

IDlist = list(range(1,
                    len(uniqueDict) + 1)
              )

residueList = list(range(1,
                         consensusLen + 1)
                   )

##################
# Export data as a single xlsx file.
##################

# Create workbook.
workbook = xlsxwriter.Workbook(path + '/' + inFileNameShort + '.xlsx')
logging.info('''Excel spreadsheet created as '%s.xlsx'.''' % inFileNameShort)

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
title_format.set_text_wrap()
wellTitle_format = workbook.add_format({'bold': True,
                                        'font_size': 12
                                        }
                                       )
wellTitle_format.set_align('left')
wellTitle_format.set_align('vcenter')
library_format = workbook.add_format({'font_size': 12})
library_format.set_align('left')
library_format.set_align('vcenter')
# Numbers.
stats_format = workbook.add_format({'num_format': '#,##0.000'})
stats_format.set_align('center')
stats_format.set_align('vcenter')
residue_format = workbook.add_format({'font_size': 10})
residue_format.set_align('center')
residue_format.set_align('vcenter')
integer_format = workbook.add_format({'num_format': '#,##0'})
integer_format.set_align('center')
integer_format.set_align('vcenter')
percent_format = workbook.add_format({'num_format': '#,##0.0%'})
percent_format.set_align('center')
percent_format.set_align('vcenter')
# Wells.
wellList_format = workbook.add_format({'font_size': 11})
wellID_format = workbook.add_format({'font_size': 12})
wellID_format.set_align('center')
wellID_format.set_align('vcenter')
# Sequences.
sequence_format = workbook.add_format({'font_size': 10})
sequence_format.set_align('center')
sequence_format.set_align('vcenter')
sequence_format.set_font_name('Lucida Console')
# Region 1.
region1_format = workbook.add_format()
region1_format.set_bg_color('#BD7191')
# Region 2.
region2_format = workbook.add_format()
region2_format.set_bg_color('#8FC1C0')
# Region 3.
region3_format = workbook.add_format()
region3_format.set_bg_color('#DCA16A')
logging.info('Cell formatting rules set.')

##################
# Create worksheet for unique amino acid sequence conservation.
##################

worksheet1Name = 'Conserved AA'
worksheet1 = workbook.add_worksheet(worksheet1Name)
worksheet1.hide_gridlines(option=2)
worksheet1.set_column(0, 0, 10)
worksheet1.set_column(1, consensusLen, 2)
worksheet1.set_column(consensusLen + 1, consensusLen + 5, 8)
worksheet1.freeze_panes(3, 1)
logging.info('%s worksheet created.' % worksheet1Name)

# Assign IDs to each unique amino acid sequence.
worksheet1.merge_range(0, 0, 2, 0, 'ID', title_format)
idRow = 3
for ID in IDlist:
    worksheet1.write(idRow, 0, ID, general_format)
    idRow += 1
logging.info('IDs written to %s worksheet.' % worksheet1Name)

# Write amino acid residue numbers above sequences.
residueCol = 1
for residue in residueList:
    worksheet1.write(2, residueCol, residue, residue_format)
    residueCol += 1

# Write unique amino acid sequences.
worksheet1.merge_range(0, 1, 0, 9, 'Amino Acid Sequence', title_format)
seqRow = 3
seqCol = 1
for seq in conservedList:
    letterList = list(seq)
    for letter in letterList:
        worksheet1.write(seqRow, seqCol, letter, sequence_format)
        seqCol += 1
    seqRow += 1
    seqCol = 1
logging.info('Unique conserved sequences written to %s worksheet.' % worksheet1Name)

# Write counts for each unique amino acid sequence.
worksheet1.merge_range(0, consensusLen + 1, 2, consensusLen + 1, 'Count', title_format)
count = list(uniqueDict.values())
countRow = 3
countCol = consensusLen + 1
for number in count:
    worksheet1.write(countRow, countCol, number, general_format)
    countRow += 1
logging.info('Counts written to %s worksheet.' % worksheet1Name)

# Write statistics to worksheet.
if inputFormat == '1':
    # Max.
    maxRow = 3
    maxCol = consensusLen + 2
    worksheet1.merge_range(0, consensusLen + 2, 2, consensusLen + 2, 'Max.', title_format)
    for number in maxList:
        worksheet1.write(maxRow, maxCol, number, stats_format)
        maxRow += 1
    logging.info('Maximum values written to %s worksheet.' % worksheet1Name)

    # Min.
    minRow = 3
    minCol = consensusLen + 3
    worksheet1.merge_range(0, consensusLen + 3, 2, consensusLen + 3, 'Min.', title_format)
    for number in minList:
        worksheet1.write(minRow, minCol, number, stats_format)
        minRow += 1
    logging.info('Minimum values written to %s worksheet.' % worksheet1Name)

    # Median.
    medianRow = 3
    medianCol = consensusLen + 4
    worksheet1.merge_range(0, consensusLen + 4, 2, consensusLen + 4, 'Median', title_format)
    for number in medianList:
        worksheet1.write(medianRow, medianCol, number, stats_format)
        medianRow += 1
    logging.info('Median values written to %s worksheet.' % worksheet1Name)

    # Mean.
    meanRow = 3
    meanCol = consensusLen + 5
    worksheet1.merge_range(0, consensusLen + 5, 2, consensusLen + 5, 'Mean', title_format)
    for number in meanList:
        worksheet1.write(meanRow, meanCol, number, stats_format)
        meanRow += 1
    logging.info('Mean values written to %s worksheet.' % worksheet1Name)

    # St. dev.
    stdevRow = 3
    stdevCol = consensusLen + 6
    worksheet1.merge_range(0, consensusLen + 6, 2, consensusLen + 6, 'St. Dev.', title_format)
    for number in devList:
        worksheet1.write(stdevRow, stdevCol, number, stats_format)
        stdevRow += 1
    logging.info('Standard deviation values written to %s worksheet.' % worksheet1Name)

    # Wells.
    worksheet1.merge_range(0, consensusLen + 7, 2, consensusLen + 7, 'Wells', wellTitle_format)
    wellRow = 3
    wellCol = consensusLen + 7
    # Change column width to fit all IDs.
    wellColWidth = round((len(countID[0]) / 1.16))
    worksheet1.set_column(consensusLen + 7, consensusLen + 7, wellColWidth)
    # Write specific IDs to worksheet.
    for wellList in countID:
        worksheet1.write(wellRow, wellCol, wellList, wellList_format)
        wellRow += 1
    logging.info('Wells written to %s worksheet.' % worksheet1Name)

    # Conditional formatting for statistics.
    worksheet1.conditional_format(2,
                                  consensusLen + 2,
                                  len(conservedList) + 2,
                                  consensusLen + 6,
                                  {'type': '2_color_scale', 'min_color': '#FAFAFA', 'max_color': '#008000'})
    logging.info('Conditional formatting applied to statistics.')

elif inputFormat == '2':
    # Wells.
    worksheet1.write(1, consensusLen + 2, 'Wells', wellTitle_format)
    wellRow = 3
    wellCol = consensusLen + 2
    # Change column width to fit all IDs.
    wellColWidth = round((len(countID[0]) / 1.16))
    worksheet1.set_column(consensusLen + 2, consensusLen + 2, wellColWidth)
    # Write specific IDs to worksheet.
    for wellList in countID:
        worksheet1.write(wellRow, wellCol, wellList, wellList_format)
        wellRow += 1
    logging.info('Amino acid sequence-specific well IDs written to %s worksheet.' % worksheet1Name)

##################
# Choose phage display library to overlay.
##################

# User prompt.
cyanprint('''\nChoose which library design to use (type corresponding number in square brackets):

[1] Library 1 (Ernst et al., 2013).
Diversified residues:
(Region 1) 2, 4, 6, 8-12, 14
(Region 2) 35, 37, 39-40, 42, 44, 46-49
(Region 3) 62-64, 66, 68, 70-72

[2] Library 2 (Ernst et al., 2013).
Diversified residues:
(Region 1) 2, 4, 6, 8-12, 14
(Region 2) 42, 44, 46-49
(Region 3) 62-64, 66, 68, 70-78

Type 'pass' to skip.'''
          )

libraryOptions = {'1': 'Library 1 (Ernst et al., 2013)',
                  '2': 'Library 2 (Ernst et al., 2013)',
                  'pass': 'No library chosen'
                  }
while True:
    libraryInput = input()
    libraryInput = libraryInput.lower()
    # Library 1 (Ernst et al., 2013).
    if libraryInput in libraryOptions:
        # Region 1 formatting.
        if libraryInput == '1':
            worksheet1.write(len(conservedList) + 4, 1, libraryOptions.get('1'), library_format)
            logging.info('Library 1 (Ernst et al., 2013) selected.')
            worksheet1.merge_range(1, 2, 1, 14, 'Region 1', title_format)
            worksheet1.conditional_format(3, 2, len(conservedList) + 2, 2,
                                          {'type': 'no_blanks', 'format': region1_format}
                                          )
            worksheet1.conditional_format(3, 4, len(conservedList) + 2, 4,
                                          {'type': 'no_blanks', 'format': region1_format}
                                          )
            worksheet1.conditional_format(3, 6, len(conservedList) + 2, 6,
                                          {'type': 'no_blanks', 'format': region1_format}
                                          )
            worksheet1.conditional_format(3, 8, len(conservedList) + 2, 12,
                                          {'type': 'no_blanks', 'format': region1_format}
                                          )
            worksheet1.conditional_format(3, 14, len(conservedList) + 2, 14,
                                          {'type': 'no_blanks', 'format': region1_format}
                                          )
            logging.info('Region 1 coloured.')
            # Region 2 formatting.
            worksheet1.merge_range(1, 35, 1, 49, 'Region 2', title_format)
            worksheet1.conditional_format(3, 35, len(conservedList) + 2, 35,
                                          {'type': 'no_blanks', 'format': region2_format}
                                          )
            worksheet1.conditional_format(3, 37, len(conservedList) + 2, 37,
                                          {'type': 'no_blanks', 'format': region2_format}
                                          )
            worksheet1.conditional_format(3, 39, len(conservedList) + 2, 40,
                                          {'type': 'no_blanks', 'format': region2_format}
                                          )
            worksheet1.conditional_format(3, 42, len(conservedList) + 2, 42,
                                          {'type': 'no_blanks', 'format': region2_format}
                                          )
            worksheet1.conditional_format(3, 44, len(conservedList) + 2, 44,
                                          {'type': 'no_blanks', 'format': region2_format}
                                          )
            worksheet1.conditional_format(3, 46, len(conservedList) + 2, 49,
                                          {'type': 'no_blanks', 'format': region2_format}
                                          )
            logging.info('Region 2 coloured.')
            # Region 3 formatting.
            worksheet1.merge_range(1, 62, 1, 72, 'Region 3', title_format)
            worksheet1.conditional_format(3, 62, len(conservedList) + 2, 64,
                                          {'type': 'no_blanks', 'format': region3_format}
                                          )
            worksheet1.conditional_format(3, 66, len(conservedList) + 2, 66,
                                          {'type': 'no_blanks', 'format': region3_format}
                                          )
            worksheet1.conditional_format(3, 68, len(conservedList) + 2, 68,
                                          {'type': 'no_blanks', 'format': region3_format}
                                          )
            worksheet1.conditional_format(3, 70, len(conservedList) + 2, 72,
                                          {'type': 'no_blanks', 'format': region3_format}
                                          )
            logging.info('Region 3 coloured.')
            break

        # Library 2 (Ernst et al., 2013).
        elif libraryInput == '2':
            worksheet1.write(len(conservedList) + 4, 1, libraryOptions.get('2'), library_format)
            logging.info('Library 2 (Ernst et al., 2013) selected.')
            # Region 1 formatting.
            worksheet1.merge_range(1, 2, 1, 14, 'Region 1', title_format)
            worksheet1.conditional_format(3, 2, len(conservedList) + 2, 2,
                                          {'type': 'no_blanks', 'format': region1_format}
                                          )
            worksheet1.conditional_format(3, 4, len(conservedList) + 2, 4,
                                          {'type': 'no_blanks', 'format': region1_format}
                                          )
            worksheet1.conditional_format(3, 6, len(conservedList) + 2, 6,
                                          {'type': 'no_blanks', 'format': region1_format}
                                          )
            worksheet1.conditional_format(3, 8, len(conservedList) + 2, 12,
                                          {'type': 'no_blanks', 'format': region1_format}
                                          )
            worksheet1.conditional_format(3, 14, len(conservedList) + 2, 14,
                                          {'type': 'no_blanks', 'format': region1_format}
                                          )
            logging.info('Region 1 coloured.')
            # Region 2 formatting.
            worksheet1.merge_range(1, 42, 1, 49, 'Region 2', title_format)
            worksheet1.conditional_format(3, 42, len(conservedList) + 2, 42,
                                          {'type': 'no_blanks', 'format': region2_format}
                                          )
            worksheet1.conditional_format(3, 44, len(conservedList) + 2, 44,
                                          {'type': 'no_blanks', 'format': region2_format}
                                          )
            worksheet1.conditional_format(3, 46, len(conservedList) + 2, 49,
                                          {'type': 'no_blanks', 'format': region2_format}
                                          )
            logging.info('Region 2 coloured.')
            # Region 3 formatting.
            worksheet1.merge_range(1, 62, 1, 78, 'Region 3', title_format)
            worksheet1.conditional_format(3, 62, len(conservedList) + 2, 64,
                                          {'type': 'no_blanks', 'format': region3_format}
                                          )
            worksheet1.conditional_format(3, 66, len(conservedList) + 2, 66,
                                          {'type': 'no_blanks', 'format': region3_format}
                                          )
            worksheet1.conditional_format(3, 68, len(conservedList) + 2, 68,
                                          {'type': 'no_blanks', 'format': region3_format}
                                          )
            worksheet1.conditional_format(3, 70, len(conservedList) + 2, 78,
                                          {'type': 'no_blanks', 'format': region3_format}
                                          )
            logging.info('Region 3 coloured.')
            break

        # Pass.
        elif libraryInput == 'pass':
            logging.info('No library design selected.')
            break

    else:
        cyanprint('\nInvalid option.'
                  'Please try again.')

##################
# Analyse biochemistry and diversification patterns in sequences.
##################
# TODO: Add logging.
aaTypes = {'hydrophobic': ['G', 'A', 'V', 'P', 'F', 'M', 'L', 'W', 'I'],
           'polar': ['S', 'N', 'Y', 'H', 'C', 'T', 'G', 'D', 'E', 'R', 'K'],
           'acidic': ['D', 'E'],
           'basic': ['K', 'R', 'H'],
           'aromatic': ['F', 'W', 'Y'],
           'aliphatic': ['A', 'G', 'I', 'L', 'P', 'V'],
           }

if libraryInput == '1':
    resDiversified = []
    resUntargeted = []
    percentUntargeted = []
    resTargeted = []
    percentTargeted = []
    resDiversifiedReg1 = []
    percentDiversifiedReg1 = []
    resDiversifiedReg2 = []
    percentDiversifiedReg2 = []
    resDiversifiedReg3 = []
    percentDiversifiedReg3 = []
    resHydrophobic = []
    percentHydrophobic = []
    resPolar = []
    percentPolar = []
    resAcidic = []
    percentAcidic = []
    resBasic = []
    percentBasic = []
    resAromatic = []
    percentAromatic = []
    resAliphatic = []
    percentAliphatic = []

    region1index = [index - 1 for index in [2, 4, 6, 8, 9, 10, 11, 12, 14]]
    region2index = [index - 1 for index in [35, 37, 39, 40, 42, 44, 46, 47, 48, 49]]
    region3index = [index - 1 for index in [62, 63, 64, 66, 68, 70, 71, 72]]
    allRegionIndex = region1index + region2index + region3index

    for sequence in conservedList:
        # Total diversified residues.
        totalDiversified = 0
        for residue in sequence:
            if residue != '-':
                totalDiversified += 1
        resDiversified.append(totalDiversified)

        # Total diversified residues in targeted regions.
        totalTargeted = 0
        for index in allRegionIndex:
            if sequence[index] != '-':
                totalTargeted += 1
        resTargeted.append(totalTargeted)
        percentTargeted.append(totalTargeted / len(allRegionIndex))

        # Total diversified residues in non-targeted regions.
        totalUntargeted = totalDiversified - totalTargeted
        resUntargeted.append(totalUntargeted)
        try:
            percentUntargeted.append(totalUntargeted / (len(conservedList[0]) - len(allRegionIndex)))
        except ZeroDivisionError:
            percentUntargeted.append('0')

        # Total diversified residues in region 1.
        totalResReg1 = 0
        for index in region1index:
            if sequence[index] != '-':
                totalResReg1 += 1
        resDiversifiedReg1.append(totalResReg1)
        percentDiversifiedReg1.append(totalResReg1 / len(region1index))

        # Total diversified residues in region 2.
        totalResReg2 = 0
        for index in region2index:
            if sequence[index] != '-':
                totalResReg2 += 1
        resDiversifiedReg2.append(totalResReg2)
        percentDiversifiedReg2.append(totalResReg2 / len(region2index))

        # Total diversified residues in region 3.
        totalResReg3 = 0
        for index in region3index:
            if sequence[index] != '-':
                totalResReg3 += 1
        resDiversifiedReg3.append(totalResReg3)
        percentDiversifiedReg3.append(totalResReg3 / len(region3index))

        # Total hydrophobic residues.
        totalHydrophobic = 0
        for residue in sequence:
            if residue in aaTypes['hydrophobic']:
                totalHydrophobic += 1
        resHydrophobic.append(totalHydrophobic)
        try:
            percentHydrophobic.append(totalHydrophobic / totalDiversified)
        except ZeroDivisionError:
            percentHydrophobic.append('0')

        # Total polar residues.
        totalPolar = 0
        for residue in sequence:
            if residue in aaTypes['polar']:
                totalPolar += 1
        resPolar.append(totalPolar)
        try:
            percentPolar.append(totalPolar / totalDiversified)
        except ZeroDivisionError:
            percentPolar.append('0')

        # Total acidic residues.
        totalAcidic = 0
        for residue in sequence:
            if residue in aaTypes['acidic']:
                totalAcidic += 1
        resAcidic.append(totalAcidic)
        try:
            percentAcidic.append(totalAcidic / totalDiversified)
        except ZeroDivisionError:
            percentAcidic.append('0')

        # Total basic residues.
        totalBasic = 0
        for residue in sequence:
            if residue in aaTypes['basic']:
                totalBasic += 1
        resBasic.append(totalBasic)
        try:
            percentBasic.append(totalBasic / totalDiversified)
        except ZeroDivisionError:
            percentBasic.append('0')

        # Total aromatic residues.
        totalAromatic = 0
        for residue in sequence:
            if residue in aaTypes['aromatic']:
                totalAromatic += 1
        resAromatic.append(totalAromatic)
        try:
            percentAromatic.append(totalAromatic / totalDiversified)
        except ZeroDivisionError:
            percentAromatic.append('0')

        # Total aliphatic residues.
        totalAliphatic = 0
        for residue in sequence:
            if residue in aaTypes['aliphatic']:
                totalAliphatic += 1
        resAliphatic.append(totalAliphatic)
        try:
            percentAliphatic.append(totalAliphatic / totalDiversified)
        except ZeroDivisionError:
            percentAliphatic.append('0')

    diversityTable = {'Total Untargeted Diversified': resUntargeted,
                      'Untargeted: Percent of Untargeted Regions': percentUntargeted,
                      'Total Targeted Diversified': resTargeted,
                      'Targeted: Percent of Diversified Regions': percentTargeted,
                      'Region 1 Diversified': resDiversifiedReg1,
                      'Percent of Region 1': percentDiversifiedReg1,
                      'Region 2 Diversified': resDiversifiedReg2,
                      'Percent of Region 2': percentDiversifiedReg2,
                      'Region 3 Diversified': resDiversifiedReg3,
                      'Percent of Region 3': percentDiversifiedReg3,
                      }
    diversityDataframe = pandas.DataFrame(diversityTable)

    biochemicalTable = {'Total Hydrophobic Residues': resHydrophobic,
                        'Hydrophobic: All Diversified Residues': percentHydrophobic,
                        'Total Polar Residues': resPolar,
                        'Polar: All Diversified Residues': percentPolar,
                        'Total Acidic Residues': resAcidic,
                        'Acidic: All Diversified Residues': percentAcidic,
                        'Total Basic Residues': resBasic,
                        'Basic: All Diversified Residues': percentBasic,
                        'Total Aromatic Residues': resAromatic,
                        'Aromatic: All Diversified Residues': percentAromatic,
                        'Total Aliphatic Residues': resAliphatic,
                        'Aliphatic: All Diversified Residues': percentAliphatic
                        }
    biochemicalDataframe = pandas.DataFrame(biochemicalTable)

elif libraryInput == '2':
    resDiversified = []
    resUntargeted = []
    percentUntargeted = []
    resTargeted = []
    percentTargeted = []
    resDiversifiedReg1 = []
    percentDiversifiedReg1 = []
    resDiversifiedReg2 = []
    percentDiversifiedReg2 = []
    resDiversifiedReg3 = []
    percentDiversifiedReg3 = []
    resHydrophobic = []
    percentHydrophobic = []
    resPolar = []
    percentPolar = []
    resAcidic = []
    percentAcidic = []
    resBasic = []
    percentBasic = []
    resAromatic = []
    percentAromatic = []
    resAliphatic = []
    percentAliphatic = []

    region1index = [index - 1 for index in [2, 4, 6, 8, 9, 10, 11, 12, 14]]
    region2index = [index - 1 for index in [42, 44, 46, 47, 48, 49]]
    region3index = [index - 1 for index in [62, 63, 64, 66, 68, 70, 71, 72, 73, 74, 75, 76, 77, 78]]
    allRegionIndex = region1index + region2index + region3index

    for sequence in conservedList:
        # Total diversified residues.
        totalDiversified = 0
        for residue in sequence:
            if residue != '-':
                totalDiversified += 1
        resDiversified.append(totalDiversified)

        # Total diversified residues in targeted regions.
        totalTargeted = 0
        for index in allRegionIndex:
            if sequence[index] != '-':
                totalTargeted += 1
        resTargeted.append(totalTargeted)
        percentTargeted.append(totalTargeted / len(allRegionIndex))

        # Total diversified residues in non-targeted regions.
        totalUntargeted = totalDiversified - totalTargeted
        resUntargeted.append(totalUntargeted)
        try:
            percentUntargeted.append(totalUntargeted / (len(conservedList[0]) - len(allRegionIndex)))
        except ZeroDivisionError:
            percentUntargeted.append('0')

        # Total diversified residues in region 1.
        totalResReg1 = 0
        for index in region1index:
            if sequence[index] != '-':
                totalResReg1 += 1
        resDiversifiedReg1.append(totalResReg1)
        percentDiversifiedReg1.append(totalResReg1 / len(region1index))

        # Total diversified residues in region 2.
        totalResReg2 = 0
        for index in region2index:
            if sequence[index] != '-':
                totalResReg2 += 1
        resDiversifiedReg2.append(totalResReg2)
        percentDiversifiedReg2.append(totalResReg2 / len(region2index))

        # Total diversified residues in region 3.
        totalResReg3 = 0
        for index in region3index:
            if sequence[index] != '-':
                totalResReg3 += 1
        resDiversifiedReg3.append(totalResReg3)
        percentDiversifiedReg3.append(totalResReg3 / len(region3index))

        # Total hydrophobic residues.
        totalHydrophobic = 0
        for residue in sequence:
            if residue in aaTypes['hydrophobic']:
                totalHydrophobic += 1
        resHydrophobic.append(totalHydrophobic)
        try:
            percentHydrophobic.append(totalHydrophobic / totalDiversified)
        except ZeroDivisionError:
            percentHydrophobic.append('0')

        # Total polar residues.
        totalPolar = 0
        for residue in sequence:
            if residue in aaTypes['polar']:
                totalPolar += 1
        resPolar.append(totalPolar)
        try:
            percentPolar.append(totalPolar / totalDiversified)
        except ZeroDivisionError:
            percentPolar.append('0')

        # Total acidic residues.
        totalAcidic = 0
        for residue in sequence:
            if residue in aaTypes['acidic']:
                totalAcidic += 1
        resAcidic.append(totalAcidic)
        try:
            percentAcidic.append(totalAcidic / totalDiversified)
        except ZeroDivisionError:
            percentAcidic.append('0')

        # Total basic residues.
        totalBasic = 0
        for residue in sequence:
            if residue in aaTypes['basic']:
                totalBasic += 1
        resBasic.append(totalBasic)
        try:
            percentBasic.append(totalBasic / totalDiversified)
        except ZeroDivisionError:
            percentBasic.append('0')

        # Total aromatic residues.
        totalAromatic = 0
        for residue in sequence:
            if residue in aaTypes['aromatic']:
                totalAromatic += 1
        resAromatic.append(totalAromatic)
        try:
            percentAromatic.append(totalAromatic / totalDiversified)
        except ZeroDivisionError:
            percentAromatic.append('0')

        # Total aliphatic residues.
        totalAliphatic = 0
        for residue in sequence:
            if residue in aaTypes['aliphatic']:
                totalAliphatic += 1
        resAliphatic.append(totalAliphatic)
        try:
            percentAliphatic.append(totalAliphatic / totalDiversified)
        except ZeroDivisionError:
            percentAliphatic.append('0')

    diversityTable = {'ID': IDlist,
                      'Total Untargeted Diversified': resUntargeted,
                      'Untargeted: Percent of Untargeted Regions': percentUntargeted,
                      'Total Targeted Diversified': resTargeted,
                      'Targeted: Percent of Diversified Regions': percentTargeted,
                      'Region 1 Diversified': resDiversifiedReg1,
                      'Percent of Region 1': percentDiversifiedReg1,
                      'Region 2 Diversified': resDiversifiedReg2,
                      'Percent of Region 2': percentDiversifiedReg2,
                      'Region 3 Diversified': resDiversifiedReg3,
                      'Percent of Region 3': percentDiversifiedReg3,
                      }
    diversityDataframe = pandas.DataFrame(diversityTable)

    biochemicalTable = {'ID': IDlist,
                        'Total Hydrophobic Residues': resHydrophobic,
                        'Hydrophobic: All Diversified Residues': percentHydrophobic,
                        'Total Polar Residues': resPolar,
                        'Polar: All Diversified Residues': percentPolar,
                        'Total Acidic Residues': resAcidic,
                        'Acidic: All Diversified Residues': percentAcidic,
                        'Total Basic Residues': resBasic,
                        'Basic: All Diversified Residues': percentBasic,
                        'Total Aromatic Residues': resAromatic,
                        'Aromatic: All Diversified Residues': percentAromatic,
                        'Total Aliphatic Residues': resAliphatic,
                        'Aliphatic: All Diversified Residues': percentAliphatic
                        }
    biochemicalDataframe = pandas.DataFrame(biochemicalTable)

elif libraryInput == 'pass':
    resDiversified = []
    resHydrophobic = []
    percentHydrophobic = []
    resPolar = []
    percentPolar = []
    resAcidic = []
    percentAcidic = []
    resBasic = []
    percentBasic = []
    resAromatic = []
    percentAromatic = []
    resAliphatic = []
    percentAliphatic = []

    for sequence in conservedList:
        # Total diversified residues.
        totalDiversified = 0
        for residue in sequence:
            if residue != '-':
                totalDiversified += 1
        resDiversified.append(totalDiversified)

        # Total hydrophobic residues.
        totalHydrophobic = 0
        for residue in sequence:
            if residue in aaTypes['hydrophobic']:
                totalHydrophobic += 1
        resHydrophobic.append(totalHydrophobic)
        try:
            percentHydrophobic.append(totalHydrophobic / totalDiversified)
        except ZeroDivisionError:
            percentHydrophobic.append('0')

        # Total polar residues.
        totalPolar = 0
        for residue in sequence:
            if residue in aaTypes['polar']:
                totalPolar += 1
        resPolar.append(totalPolar)
        try:
            percentPolar.append(totalPolar / totalDiversified)
        except ZeroDivisionError:
            percentPolar.append('0')

        # Total acidic residues.
        totalAcidic = 0
        for residue in sequence:
            if residue in aaTypes['acidic']:
                totalAcidic += 1
        resAcidic.append(totalAcidic)
        try:
            percentAcidic.append(totalAcidic / totalDiversified)
        except ZeroDivisionError:
            percentAcidic.append('0')

        # Total basic residues.
        totalBasic = 0
        for residue in sequence:
            if residue in aaTypes['basic']:
                totalBasic += 1
        resBasic.append(totalBasic)
        try:
            percentBasic.append(totalBasic / totalDiversified)
        except ZeroDivisionError:
            percentBasic.append('0')

        # Total aromatic residues.
        totalAromatic = 0
        for residue in sequence:
            if residue in aaTypes['aromatic']:
                totalAromatic += 1
        resAromatic.append(totalAromatic)
        try:
            percentAromatic.append(totalAromatic / totalDiversified)
        except ZeroDivisionError:
            percentAromatic.append('0')

        # Total aliphatic residues.
        totalAliphatic = 0
        for residue in sequence:
            if residue in aaTypes['aliphatic']:
                totalAliphatic += 1
        resAliphatic.append(totalAliphatic)
        try:
            percentAliphatic.append(totalAliphatic / totalDiversified)
        except ZeroDivisionError:
            percentAliphatic.append('0')

    diversityTable = {
    }
    diversityDataframe = pandas.DataFrame(diversityTable)

    biochemicalTable = {'ID': IDlist,
                        'Total Hydrophobic Residues': resHydrophobic,
                        'Hydrophobic: All Diversified Residues': percentHydrophobic,
                        'Total Polar Residues': resPolar,
                        'Polar: All Diversified Residues': percentPolar,
                        'Total Acidic Residues': resAcidic,
                        'Acidic: All Diversified Residues': percentAcidic,
                        'Total Basic Residues': resBasic,
                        'Basic: All Diversified Residues': percentBasic,
                        'Total Aromatic Residues': resAromatic,
                        'Aromatic: All Diversified Residues': percentAromatic,
                        'Total Aliphatic Residues': resAliphatic,
                        'Aliphatic: All Diversified Residues': percentAliphatic
                        }
    biochemicalDataframe = pandas.DataFrame(biochemicalTable)

##################
# Create worksheet for unique amino acid diversity analysis.
##################
# TODO: Finish this worksheet and add biochemical analysis worksheet.
if libraryInput == '1' or '2':
    worksheet2Name = 'Diversity Analysis'
    worksheet2 = workbook.add_worksheet(worksheet2Name)
    worksheet2.hide_gridlines(option=2)
    worksheet2.set_column(0, 0, 10)
    worksheet2.set_column(1, 1, 12)
    worksheet2.set_column(2, 2, 17)
    worksheet2.set_column(3, 3, 12)
    worksheet2.set_column(4, 4, 17)
    worksheet2.freeze_panes(2, 1)
    logging.info('%s worksheet created.' % worksheet2Name)

    # Write unique sequence IDs.
    worksheet2.merge_range(0, 0, 1, 0, 'ID', title_format)
    idRow = 2
    for ID in IDlist:
        worksheet2.write(idRow, 0, ID, general_format)
        idRow += 1
    logging.info('IDs written to %s worksheet.' % worksheet2Name)

    # Write untargeted diversified residues.
    untargetedRow = 2
    untargetedCol = 1
    worksheet2.merge_range(0, 1, 1, 1, 'Untargeted', title_format)
    for untargeted in diversityDataframe['Total Untargeted Diversified']:
        worksheet2.write(untargetedRow, untargetedCol, untargeted, integer_format)
        untargetedRow += 1
    logging.info('Untargeted residues written to %s worksheet.' % worksheet2Name)

    # Write untargeted diversified residues as a percent of untargeted regions.
    untargetedPercentRow = 2
    untargetedPercentCol = 2
    worksheet2.merge_range(0, 2, 1, 2, '% of Untargeted Regions', title_format)
    for untargetedPercent in diversityDataframe['Untargeted: Percent of Untargeted Regions']:
        worksheet2.write(untargetedPercentRow, untargetedPercentCol, untargetedPercent, percent_format)
        untargetedPercentRow += 1
    logging.info('Untargeted residues as a percent of untargeted regions written to %s worksheet.' % worksheet2Name)

    # Write targeted diversified residues.
    targetedRow = 2
    targetedCol = 3
    worksheet2.merge_range(0, 3, 1, 3, 'Targeted', title_format)
    for targeted in diversityDataframe['Total Targeted Diversified']:
        worksheet2.write(targetedRow, targetedCol, targeted, integer_format)
        targetedRow += 1
    logging.info('Untargeted residues as a percent of untargeted regions written to %s worksheet.' % worksheet2Name)

    # Write targeted diversified residues as a percent of targeted regions.
    targetedPercentRow = 2
    targetedPercentCol = 4
    worksheet2.merge_range(0, 4, 1, 4, '% of Targeted Regions', title_format)
    for targetedPercent in diversityDataframe['Targeted: Percent of Diversified Regions']:
        worksheet2.write(targetedPercentRow, targetedPercentCol, targetedPercent, percent_format)
        targetedPercentRow += 1
    logging.info('Untargeted residues as a percent of untargeted regions written to %s worksheet.' % worksheet2Name)

    # Conditional formatting for diversity anlysis.
    worksheet2.conditional_format(2,
                                  1,
                                  len(conservedList) + 2,
                                  1,
                                  {'type': '2_color_scale', 'min_color': '#FAFAFA', 'max_color': '#008000'})
    worksheet2.conditional_format(2,
                                  2,
                                  len(conservedList) + 2,
                                  2,
                                  {'type': '2_color_scale', 'min_color': '#FAFAFA', 'max_color': '#008000'})
    worksheet2.conditional_format(2,
                                  3,
                                  len(conservedList) + 2,
                                  3,
                                  {'type': '2_color_scale', 'min_color': '#FAFAFA', 'max_color': '#008000'})
    worksheet2.conditional_format(2,
                                  4,
                                  len(conservedList) + 2,
                                  4,
                                  {'type': '2_color_scale', 'min_color': '#FAFAFA', 'max_color': '#008000'})
    logging.info('Conditional formatting applied to diversity analysis.')

elif libraryInput == 'pass':
    logging.info('No diversity analysis worksheet created because no library was chosen.')
    pass

##################
# Final workbook formatting.
##################
# TODO: Add tables to worksheet 2.
# Transform data into proper Excel-formatted tables without any design style applied.
if inputFormat == '1':
    worksheet1.add_table(3, 0, len(conservedList) + 2, consensusLen + 7, {'header_row': False,
                                                                          'style': None
                                                                          }
                         )
elif inputFormat == '2':
    worksheet1.add_table(3, 0, len(conservedList) + 2, consensusLen + 2, {'header_row': False,
                                                                          'style': None
                                                                          }
                         )

##################
# Conclusion
##################

workbook.close()
if inputFormat == '1':
    greenprint('\nExcel conserved alignment with ELISA scores saved as %s_conservation.xlsx.' % inFileNameShort)
    logging.info('Excel file exported as %s_conservation.xlsx.' % inFileNameShort)
elif inputFormat == '2':
    greenprint('\nExcel conserved alignment saved as %s_conservation.xlsx.' % inFileNameShort)
    logging.info('Excel file exported as %s_conservation.xlsx.' % inFileNameShort)
greenprint('\n Conservation Analysis program finished running. See log file for details.')
logging.info('Conservation Analysis program finished running.')
logging.shutdown()
