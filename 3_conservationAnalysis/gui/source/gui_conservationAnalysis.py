#! python3
# gui_conservationAnalysis.py - Analyses UbV ELISA and/or sequencing data by showing only non-conserved amino acid
# residues and highlighting regions that were diversified in the phage display library. Final amino acid alignment is
# in xlsx format.

# Usage notes:
# * This code is dependent on the style of the worksheet used as the ELISA data source. This will be entirely based
#   upon the output from the "phageDisplayELISA384well" export format used with the BioTek plate reader.
# * Any assumptions that were made from previous code will be retained.
#   E.g. if the data source is the output from "phageDisplaySeqAnalysis.py" then all alignments will exclude sequences
#   that weren't full length and those that have premature stop codons.
# * This code will assume that there is an extra amino acid residue, K, leftover from the FLAG tag when using the
#   output from the 'phageDisplaySeqAnalysis.py' code. Without this, it will trim incorrectly.

# Compatibility notes:
# * PyCharm is the recommended IDE to use. If using Spyder, avoid version 5 as this version for has conflicts with the
#   xlsxwriter package and will get stuck on importing modules.
# * This code is confirmed to work with the latest version of Python 3 (3.9). Later/earlier versions may work but have
#   not been verified.
# * This code is confirmed to work in Windows and unconfirmed to work in Macs and Linux. It should work in theory
#   but path names may need to be changed to suit Macs and Linux' path formats.


##################
#    Modules
##################

import PySimpleGUI as Sg
from molbiotools import comparestrings
import os
import re
import logging
import xlsxwriter
import pandas
from collections import Counter, OrderedDict


##################
#    Classes
##################

# Ordered list of counts.
class OrderedCounter(Counter, OrderedDict):
    pass


##################
#    GUI
##################

# Choose window theme.
Sg.theme('DarkGrey13')
# Create window layout.
layout = [

    # Title and introduction.
    [Sg.Text('Phage Display - UbV Trimming',
             text_color='#8294cc',
             font=('Segoe UI Semibold', 16),
             expand_x=True)
     ],
    [Sg.Text('''    Analyses ELISA and/or sequencing data for UbVs by showing only non-conserved amino acid
    residues and highlighting regions that were targeted for diversification based on the initial phage
    display library design.
    Final output is in xlsx format.\n''',
             text_color='#8294cc',
             font=('Segoe UI', 12)
             )
     ],

    # Input format prompt.
    [Sg.Text('''\n1. Choose input format:''',
             text_color='white',
             font=('Segoe UI Bold', 10)
             )
     ],
    [Sg.Radio('ELISA and sequencing data',
              'FORMAT',
              default=True,
              key='-FORMAT1-',
              text_color='#bfbfbf',
              font=('Segoe UI Bold', 10),
              enable_events=True
              ),
     Sg.Radio('Sequencing data only',
              'FORMAT',
              default=False,
              key='-FORMAT2-',
              text_color='#bfbfbf',
              font=('Segoe UI Bold', 10),
              enable_events=True
              )
     ],
    [Sg.Text('''    Requires xlsx output from phageDisplayElisaAnalysis.py.
    Requires fasta alignment output from phageDisplaySeqAnalysis.py.\n''',
             text_color='#bfbfbf',
             font=('Segoe UI', 10)
             )
     ],

    # File input prompt.
    [Sg.Text('''\n2. Enter the full path of the ELISA file:''',
             text_color='white',
             font=('Segoe UI Bold', 10),
             key='-FILEBEGINNINGTEXT-'
             )
     ],
    [Sg.Input(key='-FILEINPUT-',
              size=70),
     Sg.FileBrowse()
     ],
    [Sg.Text('''    * Must be in xlsx format.\n''',
             text_color='#bfbfbf',
             font=('Segoe UI', 10),
             key='-FILEENDTEXT-'
             )
     ],

    # Consensus sequence input prompt.
    [Sg.Text('3. Enter the consensus sequence for comparison:',
             text_color='white',
             font=('Segoe UI Bold', 10)
             )
     ],
    [Sg.Input(key='-CONSERVATIONINPUT-',
              size=95,
              default_text='MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGGGG'
              )
     ],
    [Sg.Text('''    * Default is human ubiquitin (PDB: ). In most cases this does not need to be changed.\n''',
             text_color='#bfbfbf',
             font=('Segoe UI', 10)
             )
     ],

    # Library format prompt.
    [Sg.Text('''\n4. Choose library design:''',
             text_color='white',
             font=('Segoe UI Bold', 10)
             )
     ],
    [Sg.Radio('Library 1 (Ernst et al., 2013)',
              'LIBRARY',
              default=True,
              key='-LIBRARY1-',
              text_color='#bfbfbf',
              font=('Segoe UI Bold', 10),
              enable_events=True
              ),
     Sg.Radio('Library 2 (Ernst et al., 2013)',
              'LIBRARY',
              default=False,
              key='-LIBRARY2-',
              text_color='#bfbfbf',
              font=('Segoe UI Bold', 10),
              enable_events=True
              ),
     Sg.Radio('N/A',
              'LIBRARY',
              default=False,
              key='-LIBRARY0-',
              text_color='#bfbfbf',
              font=('Segoe UI Bold', 10),
              enable_events=True
              )
     ],

    [Sg.Text('''    Diversified residues:
        (Region 1) 2, 4, 6, 8-12, 14
        (Region 2) 35, 37, 39-40, 42, 44, 46-49
        (Region 3) 62-64, 66, 68, 70-72\n\n''',
             text_color='#bfbfbf',
             font=('Segoe UI', 10),
             key='-LIBRARYENDTEXT-'
             )
     ],

    # 'Enter' button.
    [Sg.Button('Enter',
               bind_return_key=True,
               font=('Segoe UI Bold', 16),
               size=(10, 0),
               pad=(275, 0)
               )
     ],
]

# Name window, assign layout, and change window behaviour.
window = Sg.Window('Phage Display - UbV Conservation Analysis',
                   layout,
                   alpha_channel=0.95,
                   grab_anywhere=True,
                   resizable=True,
                   size=(725, 775))

# Create a while loop that keeps the window open.
while True:
    event, values = window.read()

    # Break the loop and close the window if 'Exit' or close window buttons are pressed.
    if event == Sg.WIN_CLOSED:
        window.close()
        break

    elif event == '-FORMAT1-':
        window['-FILEBEGINNINGTEXT-'].update('''\n1. Enter the full path of the ELISA file:''')
        window['-FILEENDTEXT-'].update('''    * Must be in xlsx format.\n''')
        continue

    elif event == '-FORMAT2-':
        window['-FILEBEGINNINGTEXT-'].update('''\n2. Enter the full path of the amino acid alignment file:''')
        window['-FILEENDTEXT-'].update('''    * Must be in fasta format.\n''')
        continue

    elif event == '-LIBRARY1-':
        window['-LIBRARYENDTEXT-'].update('''    Diversified residues:
        (Region 1) 2, 4, 6, 8-12, 14
        (Region 2) 35, 37, 39-40, 42, 44, 46-49
        (Region 3) 62-64, 66, 68, 70-72\n\n''')
        continue

    elif event == '-LIBRARY2-':
        window['-LIBRARYENDTEXT-'].update('''    Diversified residues:
        (Region 1) 2, 4, 6, 8-12, 14
        (Region 2) 42, 44, 46-49
        (Region 3) 62-64, 66, 68, 70-78\n\n''')
        continue

    elif event == '-LIBRARY0-':
        window['-LIBRARYENDTEXT-'].update('''    Diversified residues:
        N/A\n\n''')
        continue

    # If 'Enter' is pressed, updates Sg.Text with Sg.Input values.
    elif event == 'Enter':
        inputOptions = {'1': 'ELISA and sequencing',
                        '2': 'sequencing'}
        libraryOptions = {'1': 'Library 1 (Ernst et al., 2013)',
                          '2': 'Library 2 (Ernst et al., 2013)',
                          'pass': 'No library chosen'}
        # TODO: Add consensus sequence settings here.

        if values['-FORMAT1-']:
            inputFormat = '1'
            elisaFilePath = str(values['-FILEINPUT-'])
            elisaFilePath = elisaFilePath.replace('\\', '/')
            # Stops user if no file is found in the working directory.
            if not os.path.exists(elisaFilePath):
                Sg.Popup('The entered ELISA file does not exist in this location.'
                         'Please enter it again.',
                         title='File Not Found',
                         grab_anywhere=True,
                         text_color='#4276ac')
                continue
            path = re.sub(r'[a-zA-Z0-9_]+\.xlsx$', '', elisaFilePath)
            os.chdir(path)
            elisaFileName = re.findall(r'[a-zA-Z0-9_]+\.xlsx$', elisaFilePath)
            elisaFileName = elisaFileName[0]
            logging.info('%s chosen as ELISA absorbance data source.' % elisaFileName)

        if values['-FORMAT2-']:
            inputFormat = '2'
            aaAlignFilePath = str(values['-FILEINPUT-'])
            aaAlignFilePath = aaAlignFilePath.replace('\\', '/')
            # Stops user if no file is found in the working directory.
            if not os.path.exists(aaAlignFilePath):
                Sg.Popup('The entered amino acid alignment file does not exist in this location.'
                         'Please enter it again.',
                         title='File Not Found',
                         grab_anywhere=True,
                         text_color='#4276ac')
                continue
            path = re.sub(r'[a-zA-Z0-9_]+\.fasta$', '', aaAlignFilePath)
            os.chdir(path)
            aaAlignFileName = re.findall(r'[a-zA-Z0-9_]+\.fasta$', aaAlignFilePath)
            aaAlignFileName = aaAlignFileName[0]

        if values['-LIBRARY1-']:
            libraryInput = '1'

        if values['-LIBRARY2-']:
            libraryInput = '2'

        if values['-LIBRARY0-']:
            libraryInput = 'pass'

        else:
            break

##################
# Setup.
##################

# Skip performing any processes if window is closed.
if event == Sg.WIN_CLOSED:
    pass
# Carry on with code otherwise.
else:
    # Logging setup.
    if inputFormat == '1':
        elisaName = elisaFileName.replace('_analysed.xlsx', '')
        logging.basicConfig(filename=path + '/' + elisaName + '.log',
                            level=logging.INFO,
                            format='%(asctime)s - %(message)s',
                            filemode='w'
                            )
        logging.info('%s chosen as data source.' % elisaFileName)

    elif inputFormat == '2':
        aaAlignShortName = re.sub(r'_aaTrimmed_aligned.fasta', '', aaAlignFileName)
        logging.basicConfig(filename=path + '/' + aaAlignShortName + '.log',
                            level=logging.INFO,
                            format='%(asctime)s - %(message)s',
                            filemode='w'
                            )
        logging.info('%s chosen as data source.' % aaAlignFileName)
    logging.info('Working directory changed to %s.' % path)

    ##################
    # Select input format and retrieve/parse data.
    ##################

    # ELISA and sequencing data.
    if inputFormat == '1':
        # Read ELISA file.
        allCells = pandas.read_excel(elisaFileName,
                                     sheet_name=1
                                     )
        logging.info('%s data read.' % elisaFileName)

        # Retrieve statistical data.
        countListFloat = list(allCells['Count'])
        countListFloat = countListFloat[1:]
        countList = []
        for count in countListFloat:
            countList.append(int(count))
        logging.info('Count values extracted from %s.' % elisaFileName)
        maxList = list(allCells['Max.'])
        maxList = maxList[1:]
        logging.info('Maximum values extracted from %s.' % elisaFileName)
        minList = list(allCells['Min.'])
        minList = minList[1:]
        logging.info('Minimum values extracted from %s.' % elisaFileName)
        medianList = list(allCells['Median'])
        medianList = medianList[1:]
        logging.info('Median values extracted from %s.' % elisaFileName)
        meanList = list(allCells['Mean'])
        meanList = meanList[1:]
        logging.info('Mean values extracted from %s.' % elisaFileName)
        devList = list(allCells['St. Dev.'])
        devList = devList[1:]
        logging.info('Standard deviation values extracted from %s.' % elisaFileName)

        # Retrieve well data.
        wellList = list(allCells['Wells'])
        wellList = wellList[1:]
        countID = []
        for wells in wellList:
            trimList = re.findall(r'([A-Z][0-1][0-9])', wells)
            countID.append(trimList)
        logging.info('Wells extracted from %s.' % elisaFileName)

        # Retrieve amino acid sequences from ELISA file.
        trimCells = allCells.iloc[1:, 1:]
        trimCells = trimCells.to_string(index=False)
        aaList = trimCells.replace(' ', '')
        seqRegex = re.compile(r'([ARNDCEQGHILKMFPSTWYVX]{10,})')
        aaList = seqRegex.findall(aaList)

        # Create list of unique nucleotide sequences ordered by frequency.
        uniqueDict = dict(zip(aaList, countList))

    elif inputFormat == '2':
        seqRegex = re.compile(r'([ARNDCEQGHILKMFPSTWYVX]{10,})')
        stopRegex = re.compile(r'([*]+[A-Z]*)')
        with open(aaAlignFileName, 'r') as alignFile:
            allData = alignFile.read()
            # Retrieve amino acid sequences.
            seqClean = allData.replace('\n',
                                       ''
                                       )
            seqClean = stopRegex.sub('',
                                     seqClean
                                     )
            aaList = seqRegex.findall(seqClean)
            logging.info('Amino acid sequences retrieved from %s.' % aaAlignFileName)
            alignFile.close()

        # Retrieve well data.
        countIDregex = re.compile(r'([A-Z][0-1][0-9])')
        wellList = countIDregex.findall(allData)
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
        countID = []
        begin = 0
        for uniqueSeq, count in uniqueDict.items():
            end = int(count) + begin
            countID.append(orderedIndex[begin:end])
            begin += count
        logging.info('List of specific well IDs associated with amino acid sequences created.')

    ##################
    # For each amino acid sequence, replace non-diversified regions with dashes.
    ##################

    # Remove amino acid prior to start codon.
    aaShortList = []
    for seq in aaList:
        aaShortSeq = seq[1:]
        aaShortList.append(aaShortSeq)
    logging.info('Initial non-ubiquitin amino acid residue removed from all sequences.')

    # TODO: Have option to change consensus sequence.
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
    # Export as .xlsx.
    ##################

    # Create workbook.
    if inputFormat == '1':
        workbook = xlsxwriter.Workbook(path + '/' + elisaName + '_conservation.xlsx')
        logging.info('''Excel spreadsheet created as '%s_conservation.xlsx'.''' % elisaName)
    elif inputFormat == '2':
        workbook = xlsxwriter.Workbook(path + '/' + aaAlignShortName + '_conservation.xlsx')
        logging.info('''Excel spreadsheet created as '%s_conservation.xlsx'.''' % aaAlignShortName)

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
    wellTitle_format = workbook.add_format({'bold': True,
                                            'font_size': 12
                                            }
                                           )
    wellTitle_format.set_align('left')
    library_format = workbook.add_format({'font_size': 12})
    library_format.set_align('left')
    library_format.set_align('vcenter')
    # Statistics.
    stats_format = workbook.add_format({'num_format': '#,##0.0'})
    stats_format.set_align('center')
    stats_format.set_align('vcenter')
    # Wells.
    wellList_format = workbook.add_format({'font_size': 11})
    wellID_format = workbook.add_format({'font_size': 12})
    wellID_format.set_align('center')
    # Residue numbers.
    residue_format = workbook.add_format({'font_size': 10})
    residue_format.set_align('center')
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
    worksheet1.freeze_panes(0, 1)
    logging.info('%s worksheet created.' % worksheet1Name)

    # Assign IDs to each unique amino acid sequence.
    worksheet1.write(1, 0, 'ID', title_format)
    numberList = list(range(1,
                            len(uniqueDict) + 1)
                      )
    row1 = 3
    for number in numberList:
        worksheet1.write(row1, 0, number, general_format)
        row1 += 1
    logging.info('IDs written to %s worksheet.' % worksheet1Name)

    # Write amino acid residue numbers above sequences.
    numberList = list(range(1,
                            consensusLen + 1)
                      )
    residueCol = 1
    for number in numberList:
        worksheet1.write(2, residueCol, number, residue_format)
        residueCol += 1

    # Write unique amino acid sequences.
    worksheet1.write(0, 6, 'Amino Acid Sequence', title_format)
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
    worksheet1.write(1, consensusLen + 1, 'Count', title_format)
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
        worksheet1.write(1, consensusLen + 2, 'Max.', title_format)
        for number in maxList:
            worksheet1.write(maxRow, maxCol, number, stats_format)
            maxRow += 1
        logging.info('Maximum values written to %s worksheet.' % worksheet1Name)

        # Min.
        minRow = 3
        minCol = consensusLen + 3
        worksheet1.write(1, consensusLen + 3, 'Min.', title_format)
        for number in minList:
            worksheet1.write(minRow, minCol, number, stats_format)
            minRow += 1
        logging.info('Minimum values written to %s worksheet.' % worksheet1Name)

        # Median.
        medianRow = 3
        medianCol = consensusLen + 4
        worksheet1.write(1, consensusLen + 4, 'Median', title_format)
        for number in medianList:
            worksheet1.write(medianRow, medianCol, number, stats_format)
            medianRow += 1
        logging.info('Median values written to %s worksheet.' % worksheet1Name)

        # Mean.
        meanRow = 3
        meanCol = consensusLen + 5
        worksheet1.write(1, consensusLen + 5, 'Mean', title_format)
        for number in meanList:
            worksheet1.write(meanRow, meanCol, number, stats_format)
            meanRow += 1
        logging.info('Mean values written to %s worksheet.' % worksheet1Name)

        # St. dev.
        stdevRow = 3
        stdevCol = consensusLen + 6
        worksheet1.write(1, consensusLen + 6, 'St. Dev.', title_format)
        for number in devList:
            worksheet1.write(stdevRow, stdevCol, number, stats_format)
            stdevRow += 1
        logging.info('Standard deviation values written to %s worksheet.' % worksheet1Name)

        # Wells.
        wellRow = 3
        wellCol = consensusLen + 7
        worksheet1.write(1, consensusLen + 7, 'Wells', wellTitle_format)
        # Change column width to fit all IDs.
        wellColWidth = len(wellList[0])
        worksheet1.set_column(consensusLen + 7, consensusLen + 7, wellColWidth)
        # Write wells to worksheet.
        for well in wellList:
            worksheet1.write(wellRow, wellCol, well, wellList_format)
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
        # Change column width to fit all IDs.
        wellColWidth = round((len(countID[0]) * 3) * 1.4)
        worksheet1.set_column(consensusLen + 2, consensusLen + 2, wellColWidth)
        # Write specific IDs to worksheet.
        worksheet1.write(1, consensusLen + 2, 'Wells', wellTitle_format)
        wellRow = 3
        wellCol = consensusLen + 2
        countIDregex = re.compile(r'([A-Z][0-1][0-9])')
        sep = ', '
        for wellList in countID:
            wellList = countIDregex.findall(str(wellList))
            wellList = sep.join(wellList)
            worksheet1.write(wellRow, wellCol, wellList, wellList_format)
            wellRow += 1
        logging.info('Amino acid sequence-specific well IDs written to %s worksheet.' % worksheet1Name)
        pass

    ##################
    # Choose phage display library to overlay.
    ##################

    # Library 1 (Ernst et al., 2013).
    # Region 1 formatting.
    if libraryInput == '1':
        worksheet1.write(len(conservedList) + 4, 1, libraryOptions.get('1'), library_format)
        logging.info('Library 1 (Ernst et al., 2013) selected.')
        worksheet1.write(1, 7, 'Region 1', title_format)
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
        worksheet1.write(1, 41, 'Region 2', title_format)
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
        worksheet1.write(1, 70, 'Region 3', title_format)
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

    # Library 2 (Ernst et al., 2013).
    elif libraryInput == '2':
        worksheet1.write(len(conservedList) + 4, 1, libraryOptions.get('2'), library_format)
        logging.info('Library 2 (Ernst et al., 2013) selected.')
        # Region 1 formatting.
        worksheet1.write(1, 7, 'Region 1', title_format)
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
        worksheet1.write(1, 41, 'Region 2', title_format)
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
        worksheet1.write(1, 70, 'Region 3', title_format)
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

    # Pass.
    elif libraryInput == 'pass':
        logging.info('No library design selected.')

    ##################
    # Final workbook formatting.
    ##################

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
        logging.info('Excel file exported as %s_conservation.xlsx.' % elisaFileName)
        Sg.Popup('''Analysis finished. See log file for details.
        \n\nPost-analysis help:
        \nNon-trimmed files are in the 'noTrim' folder and couldn't be trimmed because of one of the following reasons:
        \n      a) Statistical error.
        Statistics will encounter an error if 'overflow' cells are in the raw ELISA data. This is reflected in the'''
                 ''' number of averaged ELISA scores (i.e. the length of cellAve) being less than the total number of'''
                 ''' sequences retrieved from the alignment files. Replace 'OVFLW' wells with '4' to fix.
    \n      b) Indexing error.
    If you encounter an error involving index values being out of range, this is because the sum of all the'''
                 ''' unique sequence counts does not match the total number of ELISA scores. Check the original'''
                 ''' alignment files to make sure sequences aren't repeated.''',
                 title='Analysis Finished',
                 grab_anywhere=True,
                 text_color='#4276ac')
        logging.info('phageDisplayElisaAnalysis.py finished running.')
    elif inputFormat == '2':
        logging.info('Excel file exported as %s_conservation.xlsx.' % aaAlignShortName)
        Sg.Popup('''Analysis finished. See log file for details.
        \n\nPost-analysis help:
        \nNon-trimmed files are in the 'noTrim' folder and couldn't be trimmed because of one of the following reasons:
        \n      a) Statistical error.
        Statistics will encounter an error if 'overflow' cells are in the raw ELISA data. This is reflected in the'''
                 ''' number of averaged ELISA scores (i.e. the length of cellAve) being less than the total number of'''
                 ''' sequences retrieved from the alignment files. Replace 'OVFLW' wells with '4' to fix.
    \n      b) Indexing error.
    If you encounter an error involving index values being out of range, this is because the sum of all the'''
                 ''' unique sequence counts does not match the total number of ELISA scores. Check the original'''
                 ''' alignment files to make sure sequences aren't repeated.''',
                 title='Analysis Finished',
                 grab_anywhere=True,
                 text_color='#4276ac')
        logging.info('phageDisplayElisaAnalysis.py finished running.')
    logging.info('ubvTrims.py finished running.')
    logging.shutdown()
