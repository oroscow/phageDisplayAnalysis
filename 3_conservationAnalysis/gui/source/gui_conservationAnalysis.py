#! python3

##################
#    Modules
##################

import PySimpleGUI as Sg
import os
import re
import logging
import xlsxwriter
import pandas
from collections import Counter, OrderedDict


##################
#    FUNCTIONS
##################

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
    [Sg.Text('Phage Display - Conservation Analysis',
             text_color='#8294cc',
             font=('Segoe UI Semibold', 16),
             expand_x=True)
     ],
    [Sg.Text('''    Analyses UbV ELISA and/or sequencing data by showing only non-conserved amino acid
    residues, highlighting regions that were targeted for diversification in the original 
    phage display library, and calculating diversity and biochemistry metrics.\n''',
             text_color='#8294cc',
             font=('Segoe UI', 12)
             )
     ],

    # Input format prompt.
    [Sg.Text('1. Choose input format:',
             text_color='white',
             font=('Segoe UI Bold', 10)
             )
     ],
    [Sg.Radio('Binding and sequencing data',
              'FORMAT',
              default=True,
              key='-ELISA_SEQ_DATA-',
              text_color='#bfbfbf',
              font=('Segoe UI Bold', 10),
              enable_events=True
              ),
     Sg.Radio('Sequencing data only',
              'FORMAT',
              default=False,
              key='-SEQ_DATA-',
              text_color='#bfbfbf',
              font=('Segoe UI Bold', 10),
              enable_events=True
              )
     ],
    [Sg.Text('''    * Requires xlsx output from Binding Analysis program.\n''',
             text_color='#bfbfbf',
             font=('Segoe UI', 10),
             key='-RADIO_TEXT-'
             )
     ],

    # File input prompt.
    [Sg.Text('2. Enter the full path of the ELISA file:',
             text_color='white',
             font=('Segoe UI Bold', 10),
             key='-INFILE_BEGIN_TEXT-'
             )
     ],
    [Sg.Input(key='-FILEINPUT-',
              size=70,
              pad=(25, 0),
              font=('Segoe UI', 10),
              focus=True
              ),
     Sg.FileBrowse(font=('Segoe UI Bold', 10),
                   size=(10, 0),
                   file_types=(('Excel/Fasta Files', '*.xlsx;*.fasta'), ('All Files', '*.*'),),
                   )
     ],
    [Sg.Text('''    * Must be in xlsx format.''',
             text_color='#bfbfbf',
             font=('Segoe UI', 10),
             key='-INFILE_END_TEXT-'
             )
     ],

    # Trim prompt.
    [Sg.Radio('Trim first amino acid',
              'TRIM',
              default=True,
              key='-TRIM-',
              text_color='#bfbfbf',
              font=('Segoe UI Bold', 10),
              enable_events=True
              ),
     Sg.Radio('Do not trim first amino acid',
              'TRIM',
              default=False,
              key='-NOTRIM-',
              text_color='#bfbfbf',
              font=('Segoe UI Bold', 10),
              enable_events=True
              )
     ],
    [Sg.Text('''    * During sequence analysis, if the second codon wasn't conserved and the codon preceding the first 
      codon had to be used, this option allows that codon to be trimmed and removed from analysis.''',
             text_color='#bfbfbf',
             font=('Segoe UI', 10)
             )
     ],

    # Consensus sequence input prompt.
    [Sg.Text('\n3. Enter the consensus sequence for comparison:',
             text_color='white',
             font=('Segoe UI Bold', 10)
             )
     ],
    [Sg.Input(key='-CONSENSUSINPUT-',
              size=100,
              pad=(25, 0),
              font=('Segoe UI', 10),
              default_text='MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGGGG'
              )
     ],
    [Sg.Text('''    * Not case-sensitive.
    * Default is human ubiquitin (accession: AAA36789, PDB: 1UBQ).\n''',
             text_color='#bfbfbf',
             font=('Segoe UI', 10)
             )
     ],

    # Library format prompt.
    [Sg.Text('4. Choose library design:',
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
        (Region 3) 62-64, 66, 68, 70-72\n''',
             text_color='#bfbfbf',
             font=('Segoe UI', 10),
             key='-LIBRARY_END_TEXT-'
             )
     ],

    # 'Enter' button.
    [Sg.Button('Enter',
               bind_return_key=True,
               font=('Segoe UI Bold', 16),
               size=(10, 0),
               pad=(30, 0),
               use_ttk_buttons=True
               )
     ],
]

# Name window, assign layout, and change window behaviour.
window = Sg.Window('Phage Display - Conservation Analysis',
                   layout,
                   alpha_channel=0.95,
                   grab_anywhere=True,
                   resizable=True,
                   size=(675, 775),
                   ttk_theme='clam'
                   )

# Create a while loop that keeps the window open.
while True:
    event, values = window.read()

    # Break the loop and close the window if 'Exit' or close window buttons are pressed.
    if event == Sg.WIN_CLOSED:
        window.close()
        break

    # Update elements based on radio buttons.
    elif event == '-ELISA_SEQ_DATA-':
        window['-INFILE_BEGIN_TEXT-'].update('2. Enter the full path of the ELISA file:')
        window['-RADIO_TEXT-'].update('    * Requires xlsx output from Binding Analysis program.\n')
        window['-INFILE_END_TEXT-'].update('    * Must be in xlsx format.\n')
        continue
    elif event == '-SEQ_DATA-':
        window['-INFILE_BEGIN_TEXT-'].update('2. Enter the full path of the amino acid alignment file:')
        window['-RADIO_TEXT-'].update('    * Requires amino acid fasta alignment output from Sequence Analysis program'
                                      ' or elsewhere.\n')
        window['-INFILE_END_TEXT-'].update('    * Must be in fasta format.\n')
        continue
    elif event == '-LIBRARY1-':
        window['-LIBRARY_END_TEXT-'].update('''    Diversified residues:
        (Region 1) 2, 4, 6, 8-12, 14
        (Region 2) 35, 37, 39-40, 42, 44, 46-49
        (Region 3) 62-64, 66, 68, 70-72\n''')
        continue
    elif event == '-LIBRARY2-':
        window['-LIBRARY_END_TEXT-'].update('''    Diversified residues:
        (Region 1) 2, 4, 6, 8-12, 14
        (Region 2) 42, 44, 46-49
        (Region 3) 62-64, 66, 68, 70-78\n''')
        continue
    elif event == '-LIBRARY0-':
        window['-LIBRARY_END_TEXT-'].update('''    Diversified residues:
        N/A\n\n\n''')
        continue

    # If 'Enter' is pressed, updates variables with input values.
    elif event == 'Enter':
        inputOptions = {'1': 'ELISA and sequencing',
                        '2': 'sequencing'
                        }
        libraryOptions = {'1': 'Library 1 (Ernst et al., 2013)',
                          '2': 'Library 2 (Ernst et al., 2013)',
                          'pass': 'No library chosen'
                          }
        consensusSeq = str(values['-CONSENSUSINPUT-'])
        consensusSeq = consensusSeq.upper()
        aaSeqRegex = re.compile(r'[ARNDCEQGHILKMFPSTWYVX]{10,}')
        consensusSeqInput = aaSeqRegex.search(consensusSeq, re.IGNORECASE)
        trimInput = 'N'

        if values['-ELISA_SEQ_DATA-']:
            inputFormat = '1'
            inFilePath = str(values['-FILEINPUT-'])
            inFilePath = inFilePath.replace('\\', '/')
            # Stops user if incorrect file type is chosen.
            if not re.search('.xlsx$', inFilePath):
                Sg.Popup('''The entered file is not an xlsx file.
Please choose the file again.''',
                         title='File Not Found',
                         grab_anywhere=True,
                         text_color='#4276ac'
                         )
                continue
            # Stops user if no file is found in the working directory.
            if not os.path.exists(inFilePath):
                Sg.Popup('''The entered ELISA file does not exist in this location.
Please choose the file again.''',
                         title='File Not Found',
                         grab_anywhere=True,
                         text_color='#4276ac'
                         )
                continue
            path = re.sub(r'[a-zA-Z0-9_]+\.xlsx$', '', inFilePath)
            os.chdir(path)
            inFileName = re.findall(r'[a-zA-Z0-9_]+\.xlsx$', inFilePath)
            inFileName = inFileName[0]

            if consensusSeqInput is None:
                Sg.Popup('''Invalid input for consensus sequence.
Please enter a valid IUPAC amino acid sequence at least 10 digits long.''',
                         title='Invalid Consensus Sequence',
                         grab_anywhere=True,
                         text_color='#4276ac',
                         any_key_closes=False
                         )
                continue

        if values['-SEQ_DATA-']:
            inputFormat = '2'
            inFilePath = str(values['-FILEINPUT-'])
            inFilePath = inFilePath.replace('\\', '/')
            # Stops user if incorrect file type is chosen.
            if not re.search('.fasta$', inFilePath):
                Sg.Popup('''The entered file is not a fasta file.
Please choose the file again.''',
                         title='File Not Found',
                         grab_anywhere=True,
                         text_color='#4276ac'
                         )
                continue
            # Stops user if no file is found in the working directory.
            if not os.path.exists(inFilePath):
                Sg.Popup('''The entered amino acid alignment file does not exist in this location.
Please choose the file again.''',
                         title='File Not Found',
                         grab_anywhere=True,
                         text_color='#4276ac'
                         )
                continue
            path = re.sub(r'[a-zA-Z0-9_]+\.fasta$', '', inFilePath)
            os.chdir(path)
            inFileName = re.findall(r'[a-zA-Z0-9_]+\.fasta$', inFilePath)
            inFileName = inFileName[0]

        if values['-TRIM-']:
            trimInput = 'Y'

        if consensusSeqInput is None:
            Sg.Popup('''Invalid input for consensus sequence.
Please enter a valid IUPAC amino acid sequence at least 10 digits long.''',
                     title='Invalid Consensus Sequence',
                     grab_anywhere=True,
                     text_color='#4276ac',
                     any_key_closes=False
                     )
            continue

        if values['-LIBRARY1-']:
            libraryInput = '1'
            break

        if values['-LIBRARY2-']:
            libraryInput = '2'
            break

        if values['-LIBRARY0-']:
            libraryInput = 'pass'
            break

##################
# Setup.
##################

# Skip performing any processes if window is closed.
if event == Sg.WIN_CLOSED:
    pass
# Carry on with code otherwise.
else:
    # Setup logging file.
    if inputFormat == '1':
        outFileNameShort = re.sub(r'[a-zA-Z]{15,}[.]xlsx',
                                  'conservationAnalysis',
                                  inFileName
                                  )
    elif inputFormat == '2':
        outFileNameShort = re.sub(r'[a-zA-Z]{6,}_*[a-zA-Z]+[.]fasta',
                                  'conservationAnalysis',
                                  inFileName
                                  )

    logging.basicConfig(filename=path + '/' + outFileNameShort + '.log',
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

    # ELISA and sequencing data.
    if inputFormat == '1':
        # Read ELISA file.
        allCells = pandas.read_excel(inFileName,
                                     sheet_name=1
                                     )
        # Remove rows with more than eight NaN values.
        allCells = allCells.dropna(axis=0, thresh=8)
        logging.info('%s data read.' % outFileNameShort)

        # Retrieve statistical data.
        countListFloat = list(allCells['Count'])
        countListFloat = countListFloat[1:]
        countList = []
        for count in countListFloat:
            countList.append(int(count))
        logging.info('Count values extracted from %s.' % outFileNameShort)
        maxList = list(allCells['Max.'])
        maxList = maxList[1:]
        logging.info('Maximum values extracted from %s.' % outFileNameShort)
        minList = list(allCells['Min.'])
        minList = minList[1:]
        logging.info('Minimum values extracted from %s.' % outFileNameShort)
        medianList = list(allCells['Median'])
        medianList = medianList[1:]
        logging.info('Median values extracted from %s.' % outFileNameShort)
        meanList = list(allCells['Mean'])
        meanList = meanList[1:]
        logging.info('Mean values extracted from %s.' % outFileNameShort)
        devList = list(allCells['St. Dev.'])
        devList = devList[1:]
        logging.info('Standard deviation values extracted from %s.' % outFileNameShort)

        # Retrieve well data.
        countID = list(allCells['Wells'])
        countID = countID[1:]
        logging.info('Wells extracted from %s.' % outFileNameShort)

        # Retrieve amino acid sequences from ELISA file.
        seqCells = allCells.iloc[1:, 1:]
        seqCells = seqCells.to_string(index=False)
        aaList = seqCells.replace(' ', '')
        aaSeqRegex = re.compile(r'[ARNDCEQGHILKMFPSTWYVX]{10,}')
        aaList = aaSeqRegex.findall(aaList)

        # Create list of unique amino acid sequences ordered by frequency.
        uniqueDict = dict(zip(aaList, countList))

    elif inputFormat == '2':
        # Read alignment file.
        aaSeqRegex = re.compile(r'[ARNDCEQGHILKMFPSTWYVX]{10,}')
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
            aaList = aaSeqRegex.findall(seqClean)
            logging.info('Amino acid sequences retrieved from %s.' % inFileName)

        # Retrieve well data.
        wellList = re.findall(r'([A-H][0-1][0-9])',
                              allData
                              )
        # Create list of unique amino acid sequences ordered by frequency.
        unique = OrderedCounter(aaList)
        unique = unique.most_common()
        uniqueDict = dict(unique)
        # Create ordered list of wells.
        orderedSeq = []
        for seq in uniqueDict.keys():
            orderedSeq.append(seq)
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

    # Trim first amino acid residue.
    aaShortList = []
    if trimInput == 'Y':
        # Remove amino acid prior to start codon.
        for seq in aaList:
            aaShortSeq = seq[1:]
            aaShortList.append(aaShortSeq)
        logging.info('First amino acid residue removed from all sequences.')
    else:
        aaShortList = aaList
        logging.info('First amino acid residue left in all sequences.')
        pass

    # Compare UbV sequences against a consensus sequence and replace conserved amino acids with dashes.
    logging.info('''Consensus sequence set to '%s'.''' % consensusSeq)
    consensusLen = len(consensusSeq)
    conservedList = []
    for querySeq in aaShortList:
        conservedList.append(comparestrings(consensusSeq,
                                            querySeq
                                            )
                             )
    logging.info('List of conserved sequences created.')

    ##################
    # Create lists of arbitrary IDs and residue numbers.
    ##################

    IDlist = [*range(1, len(uniqueDict) + 1)]

    residueList = [*range(1, consensusLen + 1)]

    ##################
    # Export data as a single xlsx file.
    ##################

    # Create workbook.
    workbook = xlsxwriter.Workbook(path + '/' + outFileNameShort + '.xlsx')
    logging.info('''Excel spreadsheet created as '%s.xlsx'.''' % outFileNameShort)

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
    title_format.set_text_wrap()
    wellTitle_format = workbook.add_format({'bold': True,
                                            'font_size': 12
                                            }
                                           )
    wellTitle_format.set_align('left')
    wellTitle_format.set_align('vcenter')
    info_format = workbook.add_format({'font_size': 12})
    info_format.set_align('left')
    info_format.set_align('vcenter')
    # Numbers.
    stats_format = workbook.add_format({'num_format': '#,##0.0'})
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
    worksheet1.set_column(1, consensusLen, 3)
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
        for column in range(consensusLen + 2, consensusLen + 7):
            worksheet1.conditional_format(2,
                                          column,
                                          len(conservedList) + 2,
                                          column,
                                          {'type': '2_color_scale', 'min_color': '#FFFFFF', 'max_color': '#3D85C6'})
        logging.info('Conditional formatting applied to statistics.')

        # Table for statistics.
        worksheet1.add_table(3, 0, len(conservedList) + 2, consensusLen + 7, {'header_row': False,
                                                                              'style': None
                                                                              }
                             )

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

        # Table for counts and wells.
        worksheet1.add_table(3, 0, len(conservedList) + 2, consensusLen + 2, {'header_row': False,
                                                                              'style': None
                                                                              }
                             )

    ##################
    # Choose phage display library to overlay.
    ##################

    # Library 1 (Ernst et al., 2013).
    if libraryInput == '1':
        # Region 1 formatting.
        worksheet1.write(len(conservedList) + 4, 1, libraryOptions.get('1'), info_format)
        logging.info('Library 1 (Ernst et al., 2013) selected.')
        worksheet1.merge_range(1, 2, 1, 14, 'Region 1', title_format)
        for column in [2, 4, 6, 8, 9, 10, 11, 12, 14]:
            worksheet1.conditional_format(3, column, len(conservedList) + 2, column,
                                          {'type': 'no_blanks', 'format': region1_format}
                                          )
        logging.info('Region 1 coloured.')
        # Region 2 formatting.
        worksheet1.merge_range(1, 35, 1, 49, 'Region 2', title_format)
        for column in [35, 37, 39, 40, 42, 44, 46, 47, 48, 49]:
            worksheet1.conditional_format(3, column, len(conservedList) + 2, column,
                                          {'type': 'no_blanks', 'format': region2_format}
                                          )
        logging.info('Region 2 coloured.')
        # Region 3 formatting.
        worksheet1.merge_range(1, 62, 1, 72, 'Region 3', title_format)
        for column in [62, 63, 64, 66, 68, 70, 71, 72]:
            worksheet1.conditional_format(3, 62, len(conservedList) + 2, 64,
                                          {'type': 'no_blanks', 'format': region3_format}
                                          )
        logging.info('Region 3 coloured.')

    # Library 2 (Ernst et al., 2013).
    elif libraryInput == '2':
        worksheet1.write(len(conservedList) + 4, 1, libraryOptions.get('2'), info_format)
        logging.info('Library 2 (Ernst et al., 2013) selected.')
        # Region 1 formatting.
        worksheet1.merge_range(1, 2, 1, 14, 'Region 1', title_format)
        for column in [2, 4, 6, 8, 9, 10, 11, 12, 14]:
            worksheet1.conditional_format(3, column, len(conservedList) + 2, column,
                                          {'type': 'no_blanks', 'format': region1_format}
                                          )
        logging.info('Region 1 coloured.')
        # Region 2 formatting.
        worksheet1.merge_range(1, 42, 1, 49, 'Region 2', title_format)
        for column in [42, 44, 46, 47, 48, 49]:
            worksheet1.conditional_format(3, column, len(conservedList) + 2, column,
                                          {'type': 'no_blanks', 'format': region2_format}
                                          )
        logging.info('Region 2 coloured.')
        # Region 3 formatting.
        worksheet1.merge_range(1, 62, 1, 78, 'Region 3', title_format)
        for column in [62, 63, 64, 66, 68, 70, 71, 72, 73, 74, 75, 76, 77, 78]:
            worksheet1.conditional_format(3, column, len(conservedList) + 2, column,
                                          {'type': 'no_blanks', 'format': region3_format}
                                          )
        logging.info('Region 3 coloured.')

    # Pass.
    elif libraryInput == 'pass':
        logging.info('No library design selected.')

    # Info about how the values were obtained.
    if inputFormat == '1':
        statsInfo = 'The statistics presented above are for binder:control ratios.'
        worksheet1.write(len(conservedList) + 6, 1, statsInfo, info_format)
    else:
        pass

    ##################
    # Analyse diversity and biochemical patterns in sequences.
    ##################

    aaTypes = {'hydrophobic': ['G', 'A', 'V', 'P', 'F', 'M', 'L', 'W', 'I'],
               'polar': ['S', 'N', 'Y', 'H', 'C', 'T', 'G', 'D', 'E', 'R', 'K'],
               'acidic': ['D', 'E'],
               'basic': ['K', 'R', 'H'],
               'aromatic': ['F', 'W', 'Y'],
               'aliphatic': ['A', 'G', 'I', 'L', 'P', 'V'],
               }

    if libraryInput == '1' or libraryInput == '2':
        if libraryInput == '1':
            region1index = [index - 1 for index in [2, 4, 6, 8, 9, 10, 11, 12, 14]]
            region2index = [index - 1 for index in [35, 37, 39, 40, 42, 44, 46, 47, 48, 49]]
            region3index = [index - 1 for index in [62, 63, 64, 66, 68, 70, 71, 72]]
            allRegionIndex = region1index + region2index + region3index
            allRegionLength = len(allRegionIndex)

        elif libraryInput == '2':
            region1index = [index - 1 for index in [2, 4, 6, 8, 9, 10, 11, 12, 14]]
            region2index = [index - 1 for index in [42, 44, 46, 47, 48, 49]]
            region3index = [index - 1 for index in [62, 63, 64, 66, 68, 70, 71, 72, 73, 74, 75, 76, 77, 78]]
            allRegionIndex = region1index + region2index + region3index
            allRegionLength = len(allRegionIndex)

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

            # Total diversified hydrophobic residues.
            totalHydrophobic = 0
            for residue in sequence:
                if residue in aaTypes['hydrophobic']:
                    totalHydrophobic += 1
            resHydrophobic.append(totalHydrophobic)
            try:
                percentHydrophobic.append(totalHydrophobic / allRegionLength)
            except ZeroDivisionError:
                percentHydrophobic.append('0')

            # Total diversified polar residues.
            totalPolar = 0
            for residue in sequence:
                if residue in aaTypes['polar']:
                    totalPolar += 1
            resPolar.append(totalPolar)
            try:
                percentPolar.append(totalPolar / len(sequence))
            except ZeroDivisionError:
                percentPolar.append('0')

            # Total diversified acidic residues.
            totalAcidic = 0
            for residue in sequence:
                if residue in aaTypes['acidic']:
                    totalAcidic += 1
            resAcidic.append(totalAcidic)
            try:
                percentAcidic.append(totalAcidic / len(sequence))
            except ZeroDivisionError:
                percentAcidic.append('0')

            # Total diversified basic residues.
            totalBasic = 0
            for residue in sequence:
                if residue in aaTypes['basic']:
                    totalBasic += 1
            resBasic.append(totalBasic)
            try:
                percentBasic.append(totalBasic / len(sequence))
            except ZeroDivisionError:
                percentBasic.append('0')

            # Total diversified aromatic residues.
            totalAromatic = 0
            for residue in sequence:
                if residue in aaTypes['aromatic']:
                    totalAromatic += 1
            resAromatic.append(totalAromatic)
            try:
                percentAromatic.append(totalAromatic / len(sequence))
            except ZeroDivisionError:
                percentAromatic.append('0')

            # Total diversified aliphatic residues.
            totalAliphatic = 0
            for residue in sequence:
                if residue in aaTypes['aliphatic']:
                    totalAliphatic += 1
            resAliphatic.append(totalAliphatic)
            try:
                percentAliphatic.append(totalAliphatic / len(sequence))
            except ZeroDivisionError:
                percentAliphatic.append('0')
        logging.info('Diversity and biochemical analyses calculated.')

        diversityTable = {'Total Untargeted Diversified': resUntargeted,
                          'Untargeted: Percent of Untargeted Regions': percentUntargeted,
                          'Total Targeted Diversified': resTargeted,
                          'Targeted: Percent of Targeted Regions': percentTargeted,
                          'Region 1 Diversified': resDiversifiedReg1,
                          'Percent of Region 1': percentDiversifiedReg1,
                          'Region 2 Diversified': resDiversifiedReg2,
                          'Percent of Region 2': percentDiversifiedReg2,
                          'Region 3 Diversified': resDiversifiedReg3,
                          'Percent of Region 3': percentDiversifiedReg3,
                          }
        diversityDataframe = pandas.DataFrame(diversityTable)

        biochemicalTable = {'Total Hydrophobic Residues': resHydrophobic,
                            'Hydrophobic: Percent of Entire Sequence': percentHydrophobic,
                            'Total Polar Residues': resPolar,
                            'Polar: Percent of Entire Sequence': percentPolar,
                            'Total Acidic Residues': resAcidic,
                            'Acidic: Percent of Entire Sequence': percentAcidic,
                            'Total Basic Residues': resBasic,
                            'Basic: Percent of Entire Sequence': percentBasic,
                            'Total Aromatic Residues': resAromatic,
                            'Aromatic: Percent of Entire Sequence': percentAromatic,
                            'Total Aliphatic Residues': resAliphatic,
                            'Aliphatic: Percent of Entire Sequence': percentAliphatic
                            }
        biochemicalDataframe = pandas.DataFrame(biochemicalTable)

    elif libraryInput == 'pass':
        resDiversified = []
        resDiversifiedPercent = []
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
            resDiversifiedPercent.append(totalDiversified / len(sequence))

            # Total diversified hydrophobic residues.
            totalHydrophobic = 0
            for residue in sequence:
                if residue in aaTypes['hydrophobic']:
                    totalHydrophobic += 1
            resHydrophobic.append(totalHydrophobic)
            try:
                percentHydrophobic.append(totalHydrophobic / len(sequence))
            except ZeroDivisionError:
                percentHydrophobic.append('0')

            # Total diversified polar residues.
            totalPolar = 0
            for residue in sequence:
                if residue in aaTypes['polar']:
                    totalPolar += 1
            resPolar.append(totalPolar)
            try:
                percentPolar.append(totalPolar / len(sequence))
            except ZeroDivisionError:
                percentPolar.append('0')

            # Total diversified acidic residues.
            totalAcidic = 0
            for residue in sequence:
                if residue in aaTypes['acidic']:
                    totalAcidic += 1
            resAcidic.append(totalAcidic)
            try:
                percentAcidic.append(totalAcidic / len(sequence))
            except ZeroDivisionError:
                percentAcidic.append('0')

            # Total diversified basic residues.
            totalBasic = 0
            for residue in sequence:
                if residue in aaTypes['basic']:
                    totalBasic += 1
            resBasic.append(totalBasic)
            try:
                percentBasic.append(totalBasic / len(sequence))
            except ZeroDivisionError:
                percentBasic.append('0')

            # Total diversified aromatic residues.
            totalAromatic = 0
            for residue in sequence:
                if residue in aaTypes['aromatic']:
                    totalAromatic += 1
            resAromatic.append(totalAromatic)
            try:
                percentAromatic.append(totalAromatic / len(sequence))
            except ZeroDivisionError:
                percentAromatic.append('0')

            # Total diversified aliphatic residues.
            totalAliphatic = 0
            for residue in sequence:
                if residue in aaTypes['aliphatic']:
                    totalAliphatic += 1
            resAliphatic.append(totalAliphatic)
            try:
                percentAliphatic.append(totalAliphatic / len(sequence))
            except ZeroDivisionError:
                percentAliphatic.append('0')
        logging.info('Diversity and biochemical analyses calculated.')

        diversityTable = {'Total Diversified': resDiversified,
                          'Diversified: Percent of Entire Sequence': resDiversifiedPercent
                          }
        diversityDataframe = pandas.DataFrame(diversityTable)

        biochemicalTable = {'Total Hydrophobic Residues': resHydrophobic,
                            'Hydrophobic: Percent of Entire Sequence': percentHydrophobic,
                            'Total Polar Residues': resPolar,
                            'Polar: Percent of Entire Sequence': percentPolar,
                            'Total Acidic Residues': resAcidic,
                            'Acidic: Percent of Entire Sequence': percentAcidic,
                            'Total Basic Residues': resBasic,
                            'Basic: Percent of Entire Sequence': percentBasic,
                            'Total Aromatic Residues': resAromatic,
                            'Aromatic: Percent of Entire Sequence': percentAromatic,
                            'Total Aliphatic Residues': resAliphatic,
                            'Aliphatic: Percent of Entire Sequence': percentAliphatic
                            }
        biochemicalDataframe = pandas.DataFrame(biochemicalTable)

    ##################
    # Create worksheet for unique amino acid diversity analyses.
    ##################

    if libraryInput == '1' or libraryInput == '2':
        worksheet2Name = 'Diversity Analyses'
        worksheet2 = workbook.add_worksheet(worksheet2Name)
        worksheet2.hide_gridlines(option=2)
        worksheet2.set_column(0, 0, 10)
        for column in [1, 3, 5, 7, 9]:
            worksheet2.set_column(column, column, 12)
        for column in [2, 4, 6, 8, 10]:
            worksheet2.set_column(column, column, 17)
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
        logging.info('Targeted residues written to %s worksheet.' % worksheet2Name)

        # Write targeted diversified residues as a percent of targeted regions.
        targetedPercentRow = 2
        targetedPercentCol = 4
        worksheet2.merge_range(0, 4, 1, 4, '% of Targeted Regions', title_format)
        for targetedPercent in diversityDataframe['Targeted: Percent of Targeted Regions']:
            worksheet2.write(targetedPercentRow, targetedPercentCol, targetedPercent, percent_format)
            targetedPercentRow += 1
        logging.info('Targeted residues as a percent of targeted regions written to %s worksheet.' % worksheet2Name)

        # Write region 1 diversified residues.
        reg1Row = 2
        reg1Col = 5
        worksheet2.merge_range(0, 5, 1, 5, 'Region 1 Diversified', title_format)
        for reg1 in diversityDataframe['Region 1 Diversified']:
            worksheet2.write(reg1Row, reg1Col, reg1, integer_format)
            reg1Row += 1
        logging.info('Region 1 diversified residues written to %s worksheet.' % worksheet2Name)

        # Write region 1 diversified residues as a percent of region 1.
        reg1PercentRow = 2
        reg1PercentCol = 6
        worksheet2.merge_range(0, 6, 1, 6, '% of Region 1', title_format)
        for reg1Percent in diversityDataframe['Percent of Region 1']:
            worksheet2.write(reg1PercentRow, reg1PercentCol, reg1Percent, percent_format)
            reg1PercentRow += 1
        logging.info('Region 1 diversified residues as a percent of region 1 written to %s worksheet.' % worksheet2Name)

        # Write region 2 diversified residues.
        reg2Row = 2
        reg2Col = 7
        worksheet2.merge_range(0, 7, 1, 7, 'Region 2 Diversified', title_format)
        for reg2 in diversityDataframe['Region 2 Diversified']:
            worksheet2.write(reg2Row, reg2Col, reg2, integer_format)
            reg2Row += 1
        logging.info('Region 2 diversified residues written to %s worksheet.' % worksheet2Name)

        # Write region 2 diversified residues as a percent of region 2.
        reg2PercentRow = 2
        reg2PercentCol = 8
        worksheet2.merge_range(0, 8, 1, 8, '% of Region 2', title_format)
        for reg2Percent in diversityDataframe['Percent of Region 2']:
            worksheet2.write(reg2PercentRow, reg2PercentCol, reg2Percent, percent_format)
            reg2PercentRow += 1
        logging.info('Region 2 diversified residues as a percent of region 2 written to %s worksheet.' % worksheet2Name)

        # Write region 3 diversified residues.
        reg3Row = 2
        reg3Col = 9
        worksheet2.merge_range(0, 9, 1, 9, 'Region 3 Diversified', title_format)
        for reg3 in diversityDataframe['Region 3 Diversified']:
            worksheet2.write(reg3Row, reg3Col, reg3, integer_format)
            reg3Row += 1
        logging.info('Region 3 diversified residues written to %s worksheet.' % worksheet2Name)

        # Write region 3 diversified residues as a percent of region 3.
        reg3PercentRow = 2
        reg3PercentCol = 10
        worksheet2.merge_range(0, 10, 1, 10, '% of Region 3', title_format)
        for reg3Percent in diversityDataframe['Percent of Region 3']:
            worksheet2.write(reg3PercentRow, reg3PercentCol, reg3Percent, percent_format)
            reg3PercentRow += 1
        logging.info('Region 3 diversified residues as a percent of region 3 written to %s worksheet.' % worksheet2Name)

        # Conditional formatting for diversity analyses.
        for column in range(1, 11):
            worksheet2.conditional_format(2, column, len(conservedList) + 2, column,
                                          {'type': '2_color_scale',
                                           'min_color': '#FFFFFF',
                                           'max_color': '#3D85C6'
                                           }
                                          )

        # Table formatting for diversity analyses.
        worksheet2.add_table(2, 0, len(conservedList) + 1, 10,
                             {'header_row': False,
                              'style': None
                              }
                             )

    elif libraryInput == 'pass':
        # Write total diversified region length.
        worksheet2Name = 'Diversity Analyses'
        worksheet2 = workbook.add_worksheet(worksheet2Name)
        worksheet2.hide_gridlines(option=2)
        worksheet2.set_column(0, 0, 10)
        worksheet2.set_column(1, 1, 12)
        worksheet2.set_column(2, 2, 17)
        worksheet2.freeze_panes(2, 1)
        logging.info('%s worksheet created.' % worksheet2Name)

        # Write unique sequence IDs.
        worksheet2.merge_range(0, 0, 1, 0, 'ID', title_format)
        idRow = 2
        for ID in IDlist:
            worksheet2.write(idRow, 0, ID, general_format)
            idRow += 1
        logging.info('IDs written to %s worksheet.' % worksheet2Name)

        # Write total diversified residues.
        diversifiedRow = 2
        diversifiedCol = 1
        worksheet2.merge_range(0, 1, 1, 1, 'Total Diversified', title_format)
        for diversified in diversityDataframe['Total Diversified']:
            worksheet2.write(diversifiedRow, diversifiedCol, diversified, integer_format)
            diversifiedRow += 1
        logging.info('Diversified residues written to %s worksheet.' % worksheet2Name)

        # Write total diversified residues as a percent of the entire sequence.
        diversifiedPercentRow = 2
        diversifiedPercentCol = 2
        worksheet2.merge_range(0, 2, 1, 2, '% of Entire Sequence', title_format)
        for diversifiedPercent in diversityDataframe['Diversified: Percent of Entire Sequence']:
            worksheet2.write(diversifiedPercentRow, diversifiedPercentCol, diversifiedPercent, percent_format)
            diversifiedPercentRow += 1
        logging.info('Diversified residues as a percent of entire sequence written to %s worksheet.' % worksheet2Name)

        # Conditional formatting for diversity analyses.
        for column in range(1, 3):
            worksheet2.conditional_format(2, column, len(conservedList) + 2, column,
                                          {'type': '2_color_scale',
                                           'min_color': '#FFFFFF',
                                           'max_color': '#3D85C6'
                                           }
                                          )

        # Table formatting for diversity analyses.
        worksheet2.add_table(2, 0, len(conservedList) + 1, 2,
                             {'header_row': False,
                              'style': None
                              }
                             )
    ##################
    # Create worksheet for unique amino acid biochemical analyses.
    ##################

    worksheet3Name = 'Biochemical Analyses'
    worksheet3 = workbook.add_worksheet(worksheet3Name)
    worksheet3.hide_gridlines(option=2)
    worksheet3.set_column(0, 0, 10)
    for column in [1, 3, 5, 7, 9, 11]:
        worksheet3.set_column(column, column, 15)
    for column in [2, 4, 6, 8, 10, 12]:
        worksheet3.set_column(column, column, 17)
    worksheet3.freeze_panes(2, 1)
    logging.info('%s worksheet created.' % worksheet3Name)

    # Write unique sequence IDs.
    worksheet3.merge_range(0, 0, 1, 0, 'ID', title_format)
    idRow = 2
    for ID in IDlist:
        worksheet3.write(idRow, 0, ID, general_format)
        idRow += 1
    logging.info('IDs written to %s worksheet.' % worksheet3Name)

    # Write diversified hydrophobic residues.
    hydrophobicRow = 2
    hydrophobicCol = 1
    worksheet3.merge_range(0, 1, 1, 1, 'Hydrophobic', title_format)
    for hydrophobic in biochemicalDataframe['Total Hydrophobic Residues']:
        worksheet3.write(hydrophobicRow, hydrophobicCol, hydrophobic, integer_format)
        hydrophobicRow += 1
    logging.info('Hydrophobic residues written to %s worksheet.' % worksheet3Name)

    # Write diversified hydrophobic residues as a percent of the entire sequence.
    hydrophobicPercentRow = 2
    hydrophobicPercentCol = 2
    worksheet3.merge_range(0, 2, 1, 2, '% of Entire Sequence', title_format)
    for hydrophobicPercent in biochemicalDataframe['Hydrophobic: Percent of Entire Sequence']:
        worksheet3.write(hydrophobicPercentRow, hydrophobicPercentCol, hydrophobicPercent, percent_format)
        hydrophobicPercentRow += 1
    logging.info('Hydrophobic residues as a percent of the entire sequence written to %s worksheet.' % worksheet3Name)

    # Write diversified polar residues.
    polarRow = 2
    polarCol = 3
    worksheet3.merge_range(0, 3, 1, 3, 'Polar', title_format)
    for polar in biochemicalDataframe['Total Polar Residues']:
        worksheet3.write(polarRow, polarCol, polar, integer_format)
        polarRow += 1
    logging.info('Polar residues written to %s worksheet.' % worksheet3Name)

    # Write diversified polar residues as a percent of the entire sequence.
    polarPercentRow = 2
    polarPercentCol = 4
    worksheet3.merge_range(0, 4, 1, 4, '% of Entire Sequence', title_format)
    for polarPercent in biochemicalDataframe['Polar: Percent of Entire Sequence']:
        worksheet3.write(polarPercentRow, polarPercentCol, polarPercent, percent_format)
        polarPercentRow += 1
    logging.info('Polar residues as a percent of the entire sequence written to %s worksheet.' % worksheet3Name)

    # Write diversified acidic residues.
    acidicRow = 2
    acidicCol = 5
    worksheet3.merge_range(0, 5, 1, 5, 'Acidic', title_format)
    for acidic in biochemicalDataframe['Total Acidic Residues']:
        worksheet3.write(acidicRow, acidicCol, acidic, integer_format)
        acidicRow += 1
    logging.info('Acidic residues written to %s worksheet.' % worksheet3Name)

    # Write diversified acidic residues as a percent of targeted regions.
    acidicPercentRow = 2
    acidicPercentCol = 6
    worksheet3.merge_range(0, 6, 1, 6, '% of Entire Sequence', title_format)
    for acidicPercent in biochemicalDataframe['Acidic: Percent of Entire Sequence']:
        worksheet3.write(acidicPercentRow, acidicPercentCol, acidicPercent, percent_format)
        acidicPercentRow += 1
    logging.info('Acidic residues as a percent of the entire sequence written to %s worksheet.' % worksheet3Name)

    # Write diversified basic residues.
    basicRow = 2
    basicCol = 7
    worksheet3.merge_range(0, 7, 1, 7, 'Basic', title_format)
    for basic in biochemicalDataframe['Total Basic Residues']:
        worksheet3.write(basicRow, basicCol, basic, integer_format)
        basicRow += 1
    logging.info('Basic residues written to %s worksheet.' % worksheet3Name)

    # Write diversified basic residues as a percent of the entire sequence.
    basicPercentRow = 2
    basicPercentCol = 8
    worksheet3.merge_range(0, 8, 1, 8, '% of Entire Sequence', title_format)
    for basicPercent in biochemicalDataframe['Basic: Percent of Entire Sequence']:
        worksheet3.write(basicPercentRow, basicPercentCol, basicPercent, percent_format)
        basicPercentRow += 1
    logging.info('Basic residues as a percent of the entire sequence written to %s worksheet.' % worksheet3Name)

    # Write diversified aromatic residues.
    aromaticRow = 2
    aromaticCol = 9
    worksheet3.merge_range(0, 9, 1, 9, 'Aromatic', title_format)
    for aromatic in biochemicalDataframe['Total Aromatic Residues']:
        worksheet3.write(aromaticRow, aromaticCol, aromatic, integer_format)
        aromaticRow += 1
    logging.info('Aromatic residues written to %s worksheet.' % worksheet3Name)

    # Write diversified aromatic residues as a percent of the entire sequence.
    aromaticPercentRow = 2
    aromaticPercentCol = 10
    worksheet3.merge_range(0, 10, 1, 10, '% of Entire Sequence', title_format)
    for aromaticPercent in biochemicalDataframe['Aromatic: Percent of Entire Sequence']:
        worksheet3.write(aromaticPercentRow, aromaticPercentCol, aromaticPercent, percent_format)
        aromaticPercentRow += 1
    logging.info('Aromatic residues as a percent of the entire sequence written to %s worksheet.' % worksheet3Name)

    # Write diversified aliphatic residues.
    aliphaticRow = 2
    aliphaticCol = 11
    worksheet3.merge_range(0, 11, 1, 11, 'Aliphatic', title_format)
    for aliphatic in biochemicalDataframe['Total Aliphatic Residues']:
        worksheet3.write(aliphaticRow, aliphaticCol, aliphatic, integer_format)
        aliphaticRow += 1
    logging.info('Aliphatic residues written to %s worksheet.' % worksheet3Name)

    # Write diversified aliphatic residues as a percent of the entire sequence.
    aliphaticPercentRow = 2
    aliphaticPercentCol = 12
    worksheet3.merge_range(0, 12, 1, 12, '% of Entire Sequence', title_format)
    for aliphaticPercent in biochemicalDataframe['Aliphatic: Percent of Entire Sequence']:
        worksheet3.write(aliphaticPercentRow, aliphaticPercentCol, aliphaticPercent, percent_format)
        aliphaticPercentRow += 1
    logging.info('Aliphatic residues as a percent of the entire sequence written to %s worksheet.' % worksheet3Name)

    # Conditional formatting for biochemical analyses.
    for column in range(1, 13):
        worksheet3.conditional_format(2, column, len(conservedList) + 2, column,
                                      {'type': '2_color_scale',
                                       'min_color': '#FFFFFF',
                                       'max_color': '#3D85C6'}
                                      )

    # Table formatting for biochemical analyses.
    worksheet3.add_table(2, 0, len(conservedList) + 1, 12,
                         {'header_row': False,
                          'style': None
                          }
                         )

    ##################
    # Conclusion
    ##################

    workbook.close()
    logging.info('Excel file exported as %s_conservationAnalysis.xlsx.' % inFileName)
    logging.info('Conservation Analysis program finished running.')
    # TODO: Add suggestions to final gui popup.
    Sg.Popup('Conservation Analysis program finished running. See log file for details.',
             title='Analysis Finished',
             grab_anywhere=True,
             text_color='#4276ac')
    logging.shutdown()
    window.close()
