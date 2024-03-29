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

# Create window layouts.
infoLayout = [
    # Title and introduction.
    [Sg.Text('\n\n\n\n\n\nPhage Display — Conservation Analysis',
             text_color='#8294cc',
             font=('Segoe UI Bold', 22),
             expand_x=True,
             pad=(50, 0)
             )
     ],
    [Sg.Text('''Compare phage sequence and ELISA (optional) data with the corresponding
wildtype scaffold sequence''',
             text_color='#8294cc',
             font=('Segoe UI', 12),
             pad=(50, 0)
             )
     ],
    [Sg.Text('''Analyses sequencing and ELISA data (optional) by showing non-conserved amino acid residues,
highlighting regions that were targeted for diversification in the original phage display
library, and calculating diversity/biochemical metrics.''',
             text_color='#a0a0a2',
             font=('Segoe UI', 10),
             pad=(70, 40)
             )
     ]
]

inputLayout = [
    # Input format prompt.
    [Sg.Text('\n\n1. Choose input format:',
             text_color='white',
             font=('Segoe UI Bold', 10),
             pad=(20, 0)
             )
     ],
    [Sg.Radio('Binding and sequencing data',
              'FORMAT',
              default=True,
              key='-ELISA_SEQ_DATA-',
              text_color='#bfbfbf',
              font=('Segoe UI Bold', 10),
              pad=(40, 0),
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
    [Sg.Text('• Requires xlsx output from Binding Analysis program.\n',
             text_color='#bfbfbf',
             font=('Segoe UI', 10),
             key='-RADIO_TEXT-',
             pad=(40, 0)
             )
     ],

    # File input prompt.
    [Sg.Text('2. Enter the full path of the ELISA file:',
             text_color='white',
             font=('Segoe UI Bold', 10),
             key='-INFILE_BEGIN_TEXT-',
             pad=(20, 0)
             )
     ],
    [Sg.Input(key='-FILE_INPUT-',
              size=70,
              pad=(40, 10),
              font=('Segoe UI', 10),
              focus=True
              ),
     Sg.FileBrowse(font=('Segoe UI Bold', 10),
                   size=(10, 0),
                   file_types=(('Excel/Fasta Files', '*.xlsx;*.fasta'), ('All Files', '*.*'),),
                   )
     ],
    [Sg.Text('''• Must be in xlsx format.''',
             text_color='#bfbfbf',
             font=('Segoe UI', 10),
             key='-INFILE_END_TEXT-',
             pad=(40, 0)
             )
     ],

    # Trim prompt.
    [Sg.Text('\n3. Choose whether to trim N-termini of aligned sequences:',
             text_color='white',
             font=('Segoe UI Bold', 10),
             pad=(20, 0)
             )
     ],
    [Sg.Checkbox(key='-N_TRIM_CHOICE-',
                 text="If checked, enter the number of residues to trim below",
                 text_color='#bfbfbf',
                 font=('Segoe UI Italic', 10),
                 pad=(40, 0),
                 default=True,
                 enable_events=True
                 ),
     ],
    [Sg.Input(key='-N_TRIM_INPUT-',
              size=3,
              pad=(40, 10),
              font=('Segoe UI', 10),
              default_text='1',
              disabled_readonly_background_color='#161616',
              disabled=False
              )
     ],
    [Sg.Text('''• Trim residues at the beginning and end of the sequence so that only the residues needed
   for alignment are included.
• For UbVs, typically only one residue at the beginning of the amino acid sequences needs to be trimmed.''',
             text_color='#bfbfbf',
             font=('Segoe UI', 10),
             pad=(40, 0)
             )
     ],

    # Reference sequence input prompt.
    [Sg.Text('\n4. Enter the reference sequence for comparison:',
             text_color='white',
             font=('Segoe UI Bold', 10),
             pad=(20, 0)
             )
     ],
    [Sg.Input(key='-REFERENCE_INPUT-',
              size=90,
              pad=(40, 10),
              font=('Segoe UI', 10),
              default_text='MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGGGG'
              )
     ],
    [Sg.Text('''• Not case-sensitive.
• Default is human ubiquitin (accession: AAA36789, PDB: 1UBQ).\n''',
             text_color='#bfbfbf',
             font=('Segoe UI', 10),
             pad=(40, 0)
             )
     ],

    # Library format prompt.
    [Sg.Text('5. Choose library design:',
             text_color='white',
             font=('Segoe UI Bold', 10),
             pad=(20, 0)
             )
     ],
    [Sg.Radio('Library 1 (Ernst et al., 2013)',
              'LIBRARY',
              default=True,
              key='-LIBRARY1-',
              text_color='#bfbfbf',
              font=('Segoe UI Bold', 10),
              pad=(40, 0),
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
              key='-NO_LIBRARY-',
              text_color='#bfbfbf',
              font=('Segoe UI Bold', 10),
              enable_events=True
              )
     ],

    [Sg.Text('''Diversified residues:
    (Region 1) 2, 4, 6, 8-12, 14
    (Region 2) 35, 37, 39-40, 42, 44, 46-49
    (Region 3) 62-64, 66, 68, 70-72\n''',
             text_color='#bfbfbf',
             font=('Segoe UI', 10),
             key='-LIBRARY_END_TEXT-',
             pad=(50, 0)
             )
     ],

    # 'Enter' button.
    [Sg.Button('Enter',
               bind_return_key=True,
               font=('Segoe UI Bold', 16),
               size=(10, 0),
               pad=(40, 30),
               use_ttk_buttons=True
               )
     ],
]

# TODO: Add to this.
troubleshootLayout = [
    [Sg.Text('''Issue #1:
    .''',
             text_color='#f44336',
             font=('Segoe UI', 10),
             pad=((50, 50), (50, 0))
             )
     ],
    [Sg.Text(
        '''Cause:
    .''',
        text_color='#bfbfbf',
        font=('Segoe UI', 10),
        pad=((50, 50), (0, 0))
    )
    ],
    [Sg.Text('''Solutions:
    1. .''',
             text_color='#93c47d',
             font=('Segoe UI', 10),
             pad=((50, 50), (0, 0))
             )
     ],
]

notesLayout = [
    [Sg.Text(
        '''• .''',
        text_color='#bfbfbf',
        font=('Segoe UI', 10),
        pad=((50, 50), (50, 0))
    )
    ]
]

# Create tab layout.
tabGroup = [
    [Sg.TabGroup(
        [
            [Sg.Tab('Info',
                    infoLayout,
                    border_width=40
                    ),
             Sg.Tab('Input Files',
                    inputLayout,
                    border_width=40
                    ),
             Sg.Tab('Troubleshooting',
                    troubleshootLayout,
                    border_width=40
                    ),
             Sg.Tab('Notes',
                    notesLayout,
                    border_width=40
                    )
             ]
        ],
        tab_location='topleft',
        border_width=1,
        font=('Segoe UI Bold', 10),
        title_color='#bfbfbf',
        selected_title_color='#8294cc',
        size=(750, 800)
    )
    ]
]

# Set up window behaviour.
window = Sg.Window('Phage Display - Conservation Analysis',
                   tabGroup,
                   alpha_channel=0.95,
                   size=(750, 825),
                   ttk_theme='clam'
                   )

# Create a while loop that keeps the window open.
while True:
    event, values = window.read()

    # Break the loop and close the window if 'Exit' or close window buttons are pressed.
    if event == Sg.WIN_CLOSED:
        window.close()
        break

    # Update elements based on radio buttons and checkboxes.
    if event == '-N_TRIM_CHOICE-':
        if values['-N_TRIM_CHOICE-']:
            window['-N_TRIM_INPUT-'].update(disabled=False)
            continue
        if not values['-N_TRIM_CHOICE-']:
            window['-N_TRIM_INPUT-'].update(disabled=True)
            continue
    elif event == '-ELISA_SEQ_DATA-':
        window['-INFILE_BEGIN_TEXT-'].update('2. Enter the full path of the ELISA file:')
        window['-RADIO_TEXT-'].update('• Requires xlsx output from Binding Analysis program.\n')
        window['-INFILE_END_TEXT-'].update('• Must be in *.xlsx format.')
        continue
    elif event == '-SEQ_DATA-':
        window['-INFILE_BEGIN_TEXT-'].update('2. Enter the full path of the amino acid alignment file:')
        window['-RADIO_TEXT-'].update('• Requires amino acid fasta alignment output from Sequence Analysis program'
                                      ' or elsewhere.\n')
        window['-INFILE_END_TEXT-'].update('• Must be in *.fasta format.')
        continue
    elif event == '-LIBRARY1-':
        window['-LIBRARY_END_TEXT-'].update('''Diversified residues:
    (Region 1) 2, 4, 6, 8-12, 14
    (Region 2) 35, 37, 39-40, 42, 44, 46-49
    (Region 3) 62-64, 66, 68, 70-72\n''')
        continue
    elif event == '-LIBRARY2-':
        window['-LIBRARY_END_TEXT-'].update('''Diversified residues:
    (Region 1) 2, 4, 6, 8-12, 14
    (Region 2) 42, 44, 46-49
    (Region 3) 62-64, 66, 68, 70-78\n''')
        continue
    elif event == '-NO_LIBRARY-':
        window['-LIBRARY_END_TEXT-'].update('''Diversified residues:
        N/A\n\n\n''')
        continue

    # Updates variables with input values when 'Enter' button is pressed.
    elif event == 'Enter':
        inputOptions = {'1': 'ELISA and sequencing',
                        '2': 'sequencing'
                        }
        libraryOptions = {'1': 'Library 1 (Ernst et al., 2013)',
                          '2': 'Library 2 (Ernst et al., 2013)',
                          'pass': 'No library chosen'
                          }
        referenceSeq = str(values['-REFERENCE_INPUT-'])
        referenceSeq = referenceSeq.upper()
        aaSeqRegex = re.compile(r'[ARNDCEQGHILKMFPSTWYVX]{10,}')
        referenceSeqInput = aaSeqRegex.search(referenceSeq, re.IGNORECASE)

        if values['-ELISA_SEQ_DATA-']:
            inputFormat = '1'
            inFilePath = str(values['-FILE_INPUT-'])
            inFilePath = inFilePath.replace('\\', '/')
            # Stops user if ELISA filetype isn't valid and prompts to retry.
            if not re.search('.xlsx$', inFilePath):
                Sg.Popup('''The entered file is not an xlsx file.
Please choose the file again.''',
                         title='File Not Found',
                         grab_anywhere=True,
                         text_color='#4276ac'
                         )
                continue
            # Stops user if ELISA path isn't valid and prompts to retry.
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

        # Stops user if sequence file path isn't valid and prompts to retry.
        if values['-SEQ_DATA-']:
            inputFormat = '2'
            inFilePath = str(values['-FILE_INPUT-'])
            inFilePath = inFilePath.replace('\\', '/')
            if not re.search('.fasta$', inFilePath):
                Sg.Popup('''The sequence file is not in fasta format.
Please choose the file again.''',
                         title='File Not Found',
                         grab_anywhere=True,
                         text_color='#4276ac'
                         )
                continue
            # Stops user if no file is found in the working directory.
            if not os.path.exists(inFilePath):
                Sg.Popup('''The  amino acid alignment file does not exist in this location.
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

        # Stops user if trim input is invalid and prompts to retry.
        if values['-N_TRIM_CHOICE-']:
            nTrimChoice = 'Y'
            nTrimInput = values['-N_TRIM_INPUT-']
            # Stops user if trim input isn't valid and prompts to retry.
            if re.search(r'[a-zA-Z]+|\s+', nTrimInput):
                Sg.Popup('''Invalid input for trim.
Please enter a number.''',
                         title='Invalid Trim Input',
                         grab_anywhere=True,
                         text_color='#4276ac'
                         )
                continue
            else:
                nTrimInput = int(values['-N_TRIM_INPUT-'])
        else:
            nTrimChoice = 'N'
            nTrimInput = 0

        # Stops user if reference sequence input is invalid and prompts to retry.
        if referenceSeqInput is None:
            Sg.Popup('''Invalid input for reference sequence.
Please enter a valid IUPAC amino acid sequence at least ten residues long.''',
                     title='Invalid Reference Sequence',
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

        if values['-NO_LIBRARY-']:
            libraryInput = 'pass'
            break

##################
#    MAIN
##################

# Skip performing any processes if window is closed.
if event == Sg.WIN_CLOSED:
    pass

##################
# Set up working directory and logging file.
##################

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
    logging.info('%i residues trimmed from sequence N-termini.' % nTrimInput)
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
    if nTrimChoice == 'Y':
        # Remove amino acid prior to start codon.
        for seq in aaList:
            aaShortSeq = seq[nTrimInput:]
            aaShortList.append(aaShortSeq)
        logging.info('First amino acid residue removed from all sequences.')
    elif nTrimChoice == 'N':
        aaShortList = aaList
        logging.info('First amino acid residue left in all sequences.')
        pass

    # Compare binder sequences against a reference sequence and replace conserved amino acids with dashes.
    logging.info('''Reference sequence set to '%s'.''' % referenceSeq)
    referenceLen = len(referenceSeq)
    conservedList = []
    for querySeq in aaShortList:
        conservedList.append(comparestrings(referenceSeq,
                                            querySeq
                                            )
                             )
    logging.info('List of conserved sequences created.')

    ##################
    # Create lists of arbitrary IDs and residue numbers.
    ##################

    IDlist = [*range(1, len(uniqueDict) + 1)]

    residueList = [*range(1, referenceLen + 1)]

    ##################
    # Export data as a single xlsx file.
    ##################

    # Create workbook.
    workbook = xlsxwriter.Workbook(path + '/' + outFileNameShort + '.xlsx')
    logging.info('''Excel spreadsheet created as '%s.xlsx'.''' % outFileNameShort)

    #########
    # Cell formatting rules.
    #########

    # General.
    general_format = workbook.add_format({'font_size': 10,
                                          'font_name': 'Segoe UI'
                                          }
                                         )
    general_format.set_align('center')
    general_format.set_align('vcenter')

    # Titles.
    title_format = workbook.add_format({'bold': True,
                                        'font_size': 10,
                                        'font_name': 'Segoe UI'
                                        }
                                       )
    title_format.set_align('center')
    title_format.set_align('vcenter')
    title_format.set_text_wrap()

    wellTitle_format = workbook.add_format({'bold': True,
                                            'font_size': 10,
                                            'font_name': 'Segoe UI'
                                            }
                                           )
    wellTitle_format.set_align('left')
    wellTitle_format.set_align('vcenter')

    # Information.
    info_format = workbook.add_format({'font_size': 10,
                                       'font_name': 'Segoe UI'
                                       }
                                      )
    info_format.set_align('left')
    info_format.set_align('vcenter')

    # Statistics.
    stats_format = workbook.add_format({'num_format': '#,##0.0',
                                        'font_size': 10,
                                        'font_name': 'Segoe UI'
                                        }
                                       )
    stats_format.set_align('center')
    stats_format.set_align('vcenter')

    integer_format = workbook.add_format({'num_format': '#,##0',
                                          'font_size': 10,
                                          'font_name': 'Segoe UI'
                                          }
                                         )
    integer_format.set_align('center')
    integer_format.set_align('vcenter')

    percent_format = workbook.add_format({'num_format': '#,##0.0%',
                                          'font_size': 10,
                                          'font_name': 'Segoe UI'
                                          }
                                         )
    percent_format.set_align('center')
    percent_format.set_align('vcenter')
    percent_format.set_right(3)

    # Residue numbers.
    residue_format = workbook.add_format({'font_size': 8,
                                          'font_name': 'Segoe UI'
                                          }
                                         )
    residue_format.set_align('center')
    residue_format.set_align('vcenter')

    # Wells.
    wellList_format = workbook.add_format({'font_size': 10,
                                           'font_name': 'Segoe UI'
                                           }
                                          )

    wellID_format = workbook.add_format({'font_size': 10,
                                         'font_name': 'Segoe UI'
                                         }
                                        )
    wellID_format.set_align('center')
    wellID_format.set_align('vcenter')

    # Reference.
    referenceSeq_format = workbook.add_format({'bg_color': '#F2F2F2'
                                               }
                                              )
    referenceSeq_format.set_bottom(3)

    referenceInteger_format = workbook.add_format({'num_format': '#,##0',
                                                   'font_size': 10,
                                                   'font_name': 'Segoe UI'
                                                   }
                                                  )
    referenceInteger_format.set_align('center')
    referenceInteger_format.set_align('vcenter')
    referenceInteger_format.set_bottom(8)

    referencePercent_format = workbook.add_format({'num_format': '#,##0.0%',
                                                   'font_size': 10,
                                                   'font_name': 'Segoe UI'
                                                   }
                                                  )
    referencePercent_format.set_align('center')
    referencePercent_format.set_align('vcenter')
    referencePercent_format.set_bottom(8)
    referencePercent_format.set_right(3)

    # Sequences.
    sequence_format = workbook.add_format({'font_size': 9,
                                           'font_name': 'Lucida Console'
                                           }
                                          )
    sequence_format.set_align('center')
    sequence_format.set_align('vcenter')

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
    worksheet1.set_column(1, referenceLen, 3)
    worksheet1.set_column(referenceLen + 1, referenceLen + 5, 8)
    worksheet1.freeze_panes(3, 1)
    logging.info('%s worksheet created.' % worksheet1Name)

    # Write amino acid residue numbers above sequences.
    residueCol = 1
    for residue in residueList:
        worksheet1.write(2, residueCol, residue, residue_format)
        residueCol += 1

    # Write reference amino acid sequence.
    worksheet1.write(3, 0, "Reference", general_format)
    worksheet1.conditional_format(3, 1, 3, len(referenceSeq) + 1,
                                  {'type': 'no_blanks', 'format': referenceSeq_format}
                                  )
    referenceRow = 3
    referenceCol = 1
    for letter in referenceSeq:
        worksheet1.write(referenceRow, referenceCol, letter, sequence_format)
        referenceCol += 1

    # Assign IDs to each unique amino acid sequence.
    worksheet1.merge_range(0, 0, 2, 0, 'ID', title_format)
    idRow = 4
    for ID in IDlist:
        worksheet1.write(idRow, 0, ID, general_format)
        idRow += 1
    logging.info('IDs written to %s worksheet.' % worksheet1Name)

    # Write unique amino acid sequences.
    worksheet1.merge_range(0, 1, 0, 9, 'Amino Acid Sequence', title_format)
    seqRow = 4
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
    worksheet1.merge_range(0, referenceLen + 1, 2, referenceLen + 1, 'Count', title_format)
    count = list(uniqueDict.values())
    countRow = 4
    countCol = referenceLen + 1
    for number in count:
        worksheet1.write(countRow, countCol, number, general_format)
        countRow += 1
    logging.info('Counts written to %s worksheet.' % worksheet1Name)

    # Write statistics to worksheet.
    if inputFormat == '1':
        # Max.
        maxRow = 4
        maxCol = referenceLen + 2
        worksheet1.merge_range(0, referenceLen + 2, 2, referenceLen + 2, 'Max.', title_format)
        for number in maxList:
            worksheet1.write(maxRow, maxCol, number, stats_format)
            maxRow += 1
        logging.info('Maximum values written to %s worksheet.' % worksheet1Name)

        # Min.
        minRow = 4
        minCol = referenceLen + 3
        worksheet1.merge_range(0, referenceLen + 3, 2, referenceLen + 3, 'Min.', title_format)
        for number in minList:
            worksheet1.write(minRow, minCol, number, stats_format)
            minRow += 1
        logging.info('Minimum values written to %s worksheet.' % worksheet1Name)

        # Median.
        medianRow = 4
        medianCol = referenceLen + 4
        worksheet1.merge_range(0, referenceLen + 4, 2, referenceLen + 4, 'Median', title_format)
        for number in medianList:
            worksheet1.write(medianRow, medianCol, number, stats_format)
            medianRow += 1
        logging.info('Median values written to %s worksheet.' % worksheet1Name)

        # Mean.
        meanRow = 4
        meanCol = referenceLen + 5
        worksheet1.merge_range(0, referenceLen + 5, 2, referenceLen + 5, 'Mean', title_format)
        for number in meanList:
            worksheet1.write(meanRow, meanCol, number, stats_format)
            meanRow += 1
        logging.info('Mean values written to %s worksheet.' % worksheet1Name)

        # St. dev.
        stdevRow = 4
        stdevCol = referenceLen + 6
        worksheet1.merge_range(0, referenceLen + 6, 2, referenceLen + 6, 'St. Dev.', title_format)
        for number in devList:
            worksheet1.write(stdevRow, stdevCol, number, stats_format)
            stdevRow += 1
        logging.info('Standard deviation values written to %s worksheet.' % worksheet1Name)

        # Wells.
        worksheet1.merge_range(0, referenceLen + 7, 2, referenceLen + 7, 'Wells', wellTitle_format)
        wellRow = 4
        wellCol = referenceLen + 7
        # Change column width to fit all IDs.
        wellColWidth = round((len(countID[0]) / 1.16))
        worksheet1.set_column(referenceLen + 7, referenceLen + 7, wellColWidth)
        # Write specific IDs to worksheet.
        for wellList in countID:
            worksheet1.write(wellRow, wellCol, wellList, wellList_format)
            wellRow += 1
        logging.info('Wells written to %s worksheet.' % worksheet1Name)

        # Conditional formatting for statistics.
        for column in range(referenceLen + 2, referenceLen + 7):
            worksheet1.conditional_format(2,
                                          column,
                                          len(conservedList) + 3,
                                          column,
                                          {'type': '2_color_scale',
                                           'min_color': '#FFFFFF',
                                           'max_color': '#3D85C6'})
        logging.info('Conditional formatting applied to statistics.')

        # Table for statistics.
        worksheet1.add_table(4, 0, len(conservedList) + 3, referenceLen + 7, {'header_row': False,
                                                                              'style': None
                                                                              }
                             )

    elif inputFormat == '2':
        # Wells.
        worksheet1.write(1, referenceLen + 2, 'Wells', wellTitle_format)
        wellRow = 4
        wellCol = referenceLen + 2
        # Change column width to fit all IDs.
        wellColWidth = round((len(countID[0]) / 1.16))
        worksheet1.set_column(referenceLen + 2, referenceLen + 2, wellColWidth)

        # Write specific IDs to worksheet.
        for wellList in countID:
            worksheet1.write(wellRow, wellCol, wellList, wellList_format)
            wellRow += 1
        logging.info('Amino acid sequence-specific well IDs written to %s worksheet.' % worksheet1Name)

        # Table for counts and wells.
        worksheet1.add_table(3, 0, len(conservedList) + 3, referenceLen + 2, {'header_row': False,
                                                                              'style': None
                                                                              }
                             )

    ##################
    # Choose phage display library to overlay.
    ##################

    # Library 1 (Ernst et al., 2013).
    if libraryInput == '1':
        # Region 1 formatting.
        worksheet1.write(len(conservedList) + 5, 1, libraryOptions.get('1'), info_format)
        logging.info('Library 1 (Ernst et al., 2013) selected.')
        worksheet1.merge_range(1, 2, 1, 14, 'Region 1', title_format)
        for column in [2, 4, 6, 8, 9, 10, 11, 12, 14]:
            worksheet1.conditional_format(3, column, len(conservedList) + 3, column,
                                          {'type': 'no_blanks', 'format': region1_format}
                                          )
        logging.info('Region 1 coloured.')
        # Region 2 formatting.
        worksheet1.merge_range(1, 35, 1, 49, 'Region 2', title_format)
        for column in [35, 37, 39, 40, 42, 44, 46, 47, 48, 49]:
            worksheet1.conditional_format(3, column, len(conservedList) + 3, column,
                                          {'type': 'no_blanks', 'format': region2_format}
                                          )
        logging.info('Region 2 coloured.')
        # Region 3 formatting.
        worksheet1.merge_range(1, 62, 1, 72, 'Region 3', title_format)
        for column in [62, 63, 64, 66, 68, 70, 71, 72]:
            worksheet1.conditional_format(3, 62, len(conservedList) + 3, 64,
                                          {'type': 'no_blanks', 'format': region3_format}
                                          )
        logging.info('Region 3 coloured.')

    # Library 2 (Ernst et al., 2013).
    elif libraryInput == '2':
        worksheet1.write(len(conservedList) + 5, 1, libraryOptions.get('2'), info_format)
        logging.info('Library 2 (Ernst et al., 2013) selected.')
        # Region 1 formatting.
        worksheet1.merge_range(1, 2, 1, 14, 'Region 1', title_format)
        for column in [2, 4, 6, 8, 9, 10, 11, 12, 14]:
            worksheet1.conditional_format(3, column, len(conservedList) + 3, column,
                                          {'type': 'no_blanks', 'format': region1_format}
                                          )
        logging.info('Region 1 coloured.')
        # Region 2 formatting.
        worksheet1.merge_range(1, 42, 1, 49, 'Region 2', title_format)
        for column in [42, 44, 46, 47, 48, 49]:
            worksheet1.conditional_format(3, column, len(conservedList) + 3, column,
                                          {'type': 'no_blanks', 'format': region2_format}
                                          )
        logging.info('Region 2 coloured.')
        # Region 3 formatting.
        worksheet1.merge_range(1, 62, 1, 78, 'Region 3', title_format)
        for column in [62, 63, 64, 66, 68, 70, 71, 72, 73, 74, 75, 76, 77, 78]:
            worksheet1.conditional_format(3, column, len(conservedList) + 3, column,
                                          {'type': 'no_blanks', 'format': region3_format}
                                          )
        logging.info('Region 3 coloured.')

    # Pass.
    elif libraryInput == 'pass':
        logging.info('No library design selected.')

    # Info about how the values were obtained.
    if inputFormat == '1':
        statsInfo = 'Statistics presented above are for binder:control binding ratios.'
        worksheet1.write(len(conservedList) + 5, referenceLen + 2, statsInfo, info_format)
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
    # Reference sequence analysis.
    resHydrophobicCon = []
    percentHydrophobicCon = []
    resPolarCon = []
    percentPolarCon = []
    resAcidicCon = []
    percentAcidicCon = []
    resBasicCon = []
    percentBasicCon = []
    resAromaticCon = []
    percentAromaticCon = []
    resAliphaticCon = []
    percentAliphaticCon = []

    # Total hydrophobic residues.
    totalHydrophobic = 0
    for residue in referenceSeq:
        if residue in aaTypes['hydrophobic']:
            totalHydrophobic += 1
    resHydrophobicCon.append(totalHydrophobic)
    try:
        percentHydrophobicCon.append(totalHydrophobic / len(referenceSeq))
    except ZeroDivisionError:
        percentHydrophobicCon.append('0')

    # Total polar residues.
    totalPolar = 0
    for residue in referenceSeq:
        if residue in aaTypes['polar']:
            totalPolar += 1
    resPolarCon.append(totalPolar)
    try:
        percentPolarCon.append(totalPolar / len(referenceSeq))
    except ZeroDivisionError:
        percentPolarCon.append('0')

    # Total diversified acidic residues.
    totalAcidic = 0
    for residue in referenceSeq:
        if residue in aaTypes['acidic']:
            totalAcidic += 1
    resAcidicCon.append(totalAcidic)
    try:
        percentAcidicCon.append(totalAcidic / len(referenceSeq))
    except ZeroDivisionError:
        percentAcidicCon.append('0')

    # Total basic residues.
    totalBasic = 0
    for residue in referenceSeq:
        if residue in aaTypes['basic']:
            totalBasic += 1
    resBasicCon.append(totalBasic)
    try:
        percentBasicCon.append(totalBasic / len(referenceSeq))
    except ZeroDivisionError:
        percentBasicCon.append('0')

    # Total aromatic residues.
    totalAromatic = 0
    for residue in referenceSeq:
        if residue in aaTypes['aromatic']:
            totalAromatic += 1
    resAromaticCon.append(totalAromatic)
    try:
        percentAromaticCon.append(totalAromatic / len(referenceSeq))
    except ZeroDivisionError:
        percentAromaticCon.append('0')

    # Total aliphatic residues.
    totalAliphatic = 0
    for residue in referenceSeq:
        if residue in aaTypes['aliphatic']:
            totalAliphatic += 1
    resAliphaticCon.append(totalAliphatic)
    try:
        percentAliphaticCon.append(totalAliphatic / len(referenceSeq))
    except ZeroDivisionError:
        percentAliphaticCon.append('0')

    # Binder sequence analysis.
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
            # Diversity analysis.
            # Total diversified residues.
            totalDiversified = 0
            for residue in sequence:
                if residue != '-':
                    totalDiversified += 1
            resDiversified.append(totalDiversified)

            # Total diversified residues in targeted regions.
            totalTargeted = 0
            for index in allRegionIndex:
                try:
                    if sequence[index] != '-':
                        totalTargeted += 1
                except IndexError:
                    continue
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
                try:
                    if sequence[index] != '-':
                        totalResReg1 += 1
                except IndexError:
                    continue
            resDiversifiedReg1.append(totalResReg1)
            percentDiversifiedReg1.append(totalResReg1 / len(region1index))

            # Total diversified residues in region 2.
            totalResReg2 = 0
            for index in region2index:
                try:
                    if sequence[index] != '-':
                        totalResReg2 += 1
                except IndexError:
                    continue
            resDiversifiedReg2.append(totalResReg2)
            percentDiversifiedReg2.append(totalResReg2 / len(region2index))

            # Total diversified residues in region 3.
            totalResReg3 = 0
            for index in region3index:
                try:
                    if sequence[index] != '-':
                        totalResReg3 += 1
                except IndexError:
                    continue
            resDiversifiedReg3.append(totalResReg3)
            percentDiversifiedReg3.append(totalResReg3 / len(region3index))
        logging.info('Diversity statistics calculated.')

        # Biochemical analysis.
        conBinderSeq = aaShortList
        conBinderSeq.insert(0, referenceSeq)
        for sequence in conBinderSeq:
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
        logging.info('Biochemical statistics calculated.')

        diversityTable = {'Total Untargeted Diversified': resUntargeted,
                          'Untargeted: Percent of Untargeted Residues': percentUntargeted,
                          'Total Targeted Diversified': resTargeted,
                          'Targeted: Percent of Targeted Residues': percentTargeted,
                          'Region 1 Diversified': resDiversifiedReg1,
                          'Percent of Region 1': percentDiversifiedReg1,
                          'Region 2 Diversified': resDiversifiedReg2,
                          'Percent of Region 2': percentDiversifiedReg2,
                          'Region 3 Diversified': resDiversifiedReg3,
                          'Percent of Region 3': percentDiversifiedReg3,
                          }
        diversityDataframe = pandas.DataFrame(diversityTable)

        biochemicalTable = {'Total Hydrophobic Residues': resHydrophobic,
                            'Hydrophobic: Percent of Sequence': percentHydrophobic,
                            'Total Polar Residues': resPolar,
                            'Polar: Percent of Sequence': percentPolar,
                            'Total Acidic Residues': resAcidic,
                            'Acidic: Percent of Sequence': percentAcidic,
                            'Total Basic Residues': resBasic,
                            'Basic: Percent of Sequence': percentBasic,
                            'Total Aromatic Residues': resAromatic,
                            'Aromatic: Percent of Sequence': percentAromatic,
                            'Total Aliphatic Residues': resAliphatic,
                            'Aliphatic: Percent of Sequence': percentAliphatic
                            }
        biochemicalDataframe = pandas.DataFrame(biochemicalTable)
        binderBiochemicalDataframe = biochemicalDataframe.drop([0, 0])
        biochemicalDataframe = pandas.DataFrame(biochemicalTable)
        referenceBiochemicalDataframe = biochemicalDataframe.loc[0]

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
        logging.info('Diversity and biochemical statistics calculated.')

        diversityTable = {'Total Diversified': resDiversified,
                          'Diversified: Percent of Sequence': resDiversifiedPercent
                          }
        diversityDataframe = pandas.DataFrame(diversityTable)

        biochemicalTable = {'Total Hydrophobic Residues': resHydrophobic,
                            'Hydrophobic: Percent of Sequence': percentHydrophobic,
                            'Total Polar Residues': resPolar,
                            'Polar: Percent of Sequence': percentPolar,
                            'Total Acidic Residues': resAcidic,
                            'Acidic: Percent of Sequence': percentAcidic,
                            'Total Basic Residues': resBasic,
                            'Basic: Percent of Sequence': percentBasic,
                            'Total Aromatic Residues': resAromatic,
                            'Aromatic: Percent of Sequence': percentAromatic,
                            'Total Aliphatic Residues': resAliphatic,
                            'Aliphatic: Percent of Sequence': percentAliphatic
                            }
        biochemicalDataframe = pandas.DataFrame(biochemicalTable)
        binderBiochemicalDataframe = biochemicalDataframe.drop([0, 0])
        biochemicalDataframe = pandas.DataFrame(biochemicalTable)
        referenceBiochemicalDataframe = biochemicalDataframe.loc[0]

    ##################
    # Create worksheet for diversity analyses.
    ##################

    if libraryInput == '1' or libraryInput == '2':
        worksheet2Name = 'Diversity Analyses'
        worksheet2 = workbook.add_worksheet(worksheet2Name)
        worksheet2.hide_gridlines(option=2)
        worksheet2.set_column(0, 0, 10)
        for column in [1, 3, 5, 7, 9]:
            worksheet2.set_column(column, column, 15)
        for column in [2, 4, 6, 8, 10]:
            worksheet2.set_column(column, column, 20)
        worksheet2.freeze_panes(3, 1)
        logging.info('%s worksheet created.' % worksheet2Name)

        # Write unique sequence IDs.
        worksheet2.merge_range(0, 0, 2, 0, 'ID', title_format)
        idRow = 3
        for ID in IDlist:
            worksheet2.write(idRow, 0, ID, general_format)
            idRow += 1
        logging.info('IDs written to %s worksheet.' % worksheet2Name)

        # Write untargeted diversified residues.
        untargetedRow = 3
        untargetedCol = 1
        worksheet2.merge_range(0, 1, 2, 1, 'Untargeted Diversified', title_format)
        for untargeted in diversityDataframe['Total Untargeted Diversified']:
            worksheet2.write(untargetedRow, untargetedCol, untargeted, integer_format)
            untargetedRow += 1
        logging.info('Untargeted residues written to %s worksheet.' % worksheet2Name)

        # Write untargeted diversified residues as a percent of untargeted regions.
        untargetedPercentRow = 3
        untargetedPercentCol = 2
        worksheet2.merge_range(0, 2, 2, 2, '% of Untargeted', title_format)
        for untargetedPercent in diversityDataframe['Untargeted: Percent of Untargeted Residues']:
            worksheet2.write(untargetedPercentRow, untargetedPercentCol, untargetedPercent, percent_format)
            untargetedPercentRow += 1
        logging.info('Untargeted residues as a percent of untargeted regions written to %s worksheet.' % worksheet2Name)

        # Write targeted diversified residues.
        targetedRow = 3
        targetedCol = 3
        worksheet2.merge_range(0, 3, 2, 3, 'Targeted Diversified', title_format)
        for targeted in diversityDataframe['Total Targeted Diversified']:
            worksheet2.write(targetedRow, targetedCol, targeted, integer_format)
            targetedRow += 1
        logging.info('Targeted residues written to %s worksheet.' % worksheet2Name)

        # Write targeted diversified residues as a percent of targeted regions.
        targetedPercentRow = 3
        targetedPercentCol = 4
        worksheet2.merge_range(0, 4, 2, 4, '% of Targeted', title_format)
        for targetedPercent in diversityDataframe['Targeted: Percent of Targeted Residues']:
            worksheet2.write(targetedPercentRow, targetedPercentCol, targetedPercent, percent_format)
            targetedPercentRow += 1
        logging.info('Targeted residues as a percent of targeted regions written to %s worksheet.' % worksheet2Name)

        # Write region 1 diversified residues.
        reg1Row = 3
        reg1Col = 5
        worksheet2.merge_range(0, 5, 2, 5, 'Region 1 Diversified', title_format)
        for reg1 in diversityDataframe['Region 1 Diversified']:
            worksheet2.write(reg1Row, reg1Col, reg1, integer_format)
            reg1Row += 1
        logging.info('Region 1 diversified residues written to %s worksheet.' % worksheet2Name)

        # Write region 1 diversified residues as a percent of region 1.
        reg1PercentRow = 3
        reg1PercentCol = 6
        worksheet2.merge_range(0, 6, 2, 6, '% of Region 1', title_format)
        for reg1Percent in diversityDataframe['Percent of Region 1']:
            worksheet2.write(reg1PercentRow, reg1PercentCol, reg1Percent, percent_format)
            reg1PercentRow += 1
        logging.info('Region 1 diversified residues as a percent of region 1 written to %s worksheet.' % worksheet2Name)

        # Write region 2 diversified residues.
        reg2Row = 3
        reg2Col = 7
        worksheet2.merge_range(0, 7, 2, 7, 'Region 2 Diversified', title_format)
        for reg2 in diversityDataframe['Region 2 Diversified']:
            worksheet2.write(reg2Row, reg2Col, reg2, integer_format)
            reg2Row += 1
        logging.info('Region 2 diversified residues written to %s worksheet.' % worksheet2Name)

        # Write region 2 diversified residues as a percent of region 2.
        reg2PercentRow = 3
        reg2PercentCol = 8
        worksheet2.merge_range(0, 8, 2, 8, '% of Region 2', title_format)
        for reg2Percent in diversityDataframe['Percent of Region 2']:
            worksheet2.write(reg2PercentRow, reg2PercentCol, reg2Percent, percent_format)
            reg2PercentRow += 1
        logging.info('Region 2 diversified residues as a percent of region 2 written to %s worksheet.' % worksheet2Name)

        # Write region 3 diversified residues.
        reg3Row = 3
        reg3Col = 9
        worksheet2.merge_range(0, 9, 2, 9, 'Region 3 Diversified', title_format)
        for reg3 in diversityDataframe['Region 3 Diversified']:
            worksheet2.write(reg3Row, reg3Col, reg3, integer_format)
            reg3Row += 1
        logging.info('Region 3 diversified residues written to %s worksheet.' % worksheet2Name)

        # Write region 3 diversified residues as a percent of region 3.
        reg3PercentRow = 3
        reg3PercentCol = 10
        worksheet2.merge_range(0, 10, 2, 10, '% of Region 3', title_format)
        for reg3Percent in diversityDataframe['Percent of Region 3']:
            worksheet2.write(reg3PercentRow, reg3PercentCol, reg3Percent, percent_format)
            reg3PercentRow += 1
        logging.info('Region 3 diversified residues as a percent of region 3 written to %s worksheet.' % worksheet2Name)

        # Conditional formatting for diversity analyses.
        for column in range(1, 11):
            worksheet2.conditional_format(2, column, len(conservedList) + 2, column,
                                          {'type': '2_color_scale',
                                           'min_color': '#FFFFFF',
                                           'max_color': '#8C198C'
                                           }
                                          )

        # Table formatting for diversity analyses.
        worksheet2.add_table(3, 0, len(conservedList) + 2, 10,
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
        worksheet2.merge_range(0, 2, 1, 2, '% of Sequence', title_format)
        for diversifiedPercent in diversityDataframe['Diversified: Percent of Sequence']:
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
    # Create worksheet for biochemical analyses.
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

    # Write reference sequence data.
    worksheet3.write(2, 0, 'Reference', general_format)
    conRow = 2
    conCol = 1
    for integer in referenceBiochemicalDataframe.iloc[range(0, 12, 2)]:
        worksheet3.write(conRow, conCol, integer, referenceInteger_format)
        conCol += 2
    conCol = 2
    for percent in referenceBiochemicalDataframe.iloc[range(1, 12, 2)]:
        worksheet3.write(conRow, conCol, percent, referencePercent_format)
        conCol += 2
    logging.info('Reference sequence data written to %s worksheet.' % worksheet3Name)

    # Write sequence IDs.
    worksheet3.merge_range(0, 0, 1, 0, 'ID', title_format)
    idRow = 3
    for ID in IDlist:
        worksheet3.write(idRow, 0, ID, general_format)
        idRow += 1
    logging.info('IDs written to %s worksheet.' % worksheet3Name)

    # Write diversified hydrophobic residues.
    hydrophobicRow = 3
    hydrophobicCol = 1
    worksheet3.merge_range(0, 1, 1, 1, 'Hydrophobic', title_format)
    for hydrophobic in binderBiochemicalDataframe['Total Hydrophobic Residues']:
        worksheet3.write(hydrophobicRow, hydrophobicCol, hydrophobic, integer_format)
        hydrophobicRow += 1
    logging.info('Hydrophobic residues written to %s worksheet.' % worksheet3Name)

    # Write diversified hydrophobic residues as a percent of the entire sequence.
    hydrophobicPercentRow = 3
    hydrophobicPercentCol = 2
    worksheet3.merge_range(0, 2, 1, 2, '% of Sequence', title_format)
    for hydrophobicPercent in binderBiochemicalDataframe['Hydrophobic: Percent of Sequence']:
        worksheet3.write(hydrophobicPercentRow, hydrophobicPercentCol, hydrophobicPercent, percent_format)
        hydrophobicPercentRow += 1
    logging.info('Hydrophobic residues as a percent of the entire sequence written to %s worksheet.' % worksheet3Name)

    # Write diversified polar residues.
    polarRow = 3
    polarCol = 3
    worksheet3.merge_range(0, 3, 1, 3, 'Polar', title_format)
    for polar in binderBiochemicalDataframe['Total Polar Residues']:
        worksheet3.write(polarRow, polarCol, polar, integer_format)
        polarRow += 1
    logging.info('Polar residues written to %s worksheet.' % worksheet3Name)

    # Write diversified polar residues as a percent of the entire sequence.
    polarPercentRow = 3
    polarPercentCol = 4
    worksheet3.merge_range(0, 4, 1, 4, '% of Sequence', title_format)
    for polarPercent in binderBiochemicalDataframe['Polar: Percent of Sequence']:
        worksheet3.write(polarPercentRow, polarPercentCol, polarPercent, percent_format)
        polarPercentRow += 1
    logging.info('Polar residues as a percent of the entire sequence written to %s worksheet.' % worksheet3Name)

    # Write diversified acidic residues.
    acidicRow = 3
    acidicCol = 5
    worksheet3.merge_range(0, 5, 1, 5, 'Acidic', title_format)
    for acidic in binderBiochemicalDataframe['Total Acidic Residues']:
        worksheet3.write(acidicRow, acidicCol, acidic, integer_format)
        acidicRow += 1
    logging.info('Acidic residues written to %s worksheet.' % worksheet3Name)

    # Write diversified acidic residues as a percent of targeted regions.
    acidicPercentRow = 3
    acidicPercentCol = 6
    worksheet3.merge_range(0, 6, 1, 6, '% of Sequence', title_format)
    for acidicPercent in binderBiochemicalDataframe['Acidic: Percent of Sequence']:
        worksheet3.write(acidicPercentRow, acidicPercentCol, acidicPercent, percent_format)
        acidicPercentRow += 1
    logging.info('Acidic residues as a percent of the entire sequence written to %s worksheet.' % worksheet3Name)

    # Write diversified basic residues.
    basicRow = 3
    basicCol = 7
    worksheet3.merge_range(0, 7, 1, 7, 'Basic', title_format)
    for basic in binderBiochemicalDataframe['Total Basic Residues']:
        worksheet3.write(basicRow, basicCol, basic, integer_format)
        basicRow += 1
    logging.info('Basic residues written to %s worksheet.' % worksheet3Name)

    # Write diversified basic residues as a percent of the entire sequence.
    basicPercentRow = 3
    basicPercentCol = 8
    worksheet3.merge_range(0, 8, 1, 8, '% of Sequence', title_format)
    for basicPercent in binderBiochemicalDataframe['Basic: Percent of Sequence']:
        worksheet3.write(basicPercentRow, basicPercentCol, basicPercent, percent_format)
        basicPercentRow += 1
    logging.info('Basic residues as a percent of the entire sequence written to %s worksheet.' % worksheet3Name)

    # Write diversified aromatic residues.
    aromaticRow = 3
    aromaticCol = 9
    worksheet3.merge_range(0, 9, 1, 9, 'Aromatic', title_format)
    for aromatic in binderBiochemicalDataframe['Total Aromatic Residues']:
        worksheet3.write(aromaticRow, aromaticCol, aromatic, integer_format)
        aromaticRow += 1
    logging.info('Aromatic residues written to %s worksheet.' % worksheet3Name)

    # Write diversified aromatic residues as a percent of the entire sequence.
    aromaticPercentRow = 3
    aromaticPercentCol = 10
    worksheet3.merge_range(0, 10, 1, 10, '% of Sequence', title_format)
    for aromaticPercent in binderBiochemicalDataframe['Aromatic: Percent of Sequence']:
        worksheet3.write(aromaticPercentRow, aromaticPercentCol, aromaticPercent, percent_format)
        aromaticPercentRow += 1
    logging.info('Aromatic residues as a percent of the entire sequence written to %s worksheet.' % worksheet3Name)

    # Write diversified aliphatic residues.
    aliphaticRow = 3
    aliphaticCol = 11
    worksheet3.merge_range(0, 11, 1, 11, 'Aliphatic', title_format)
    for aliphatic in binderBiochemicalDataframe['Total Aliphatic Residues']:
        worksheet3.write(aliphaticRow, aliphaticCol, aliphatic, integer_format)
        aliphaticRow += 1
    logging.info('Aliphatic residues written to %s worksheet.' % worksheet3Name)

    # Write diversified aliphatic residues as a percent of the entire sequence.
    aliphaticPercentRow = 3
    aliphaticPercentCol = 12
    worksheet3.merge_range(0, 12, 1, 12, '% of Sequence', title_format)
    for aliphaticPercent in binderBiochemicalDataframe['Aliphatic: Percent of Sequence']:
        worksheet3.write(aliphaticPercentRow, aliphaticPercentCol, aliphaticPercent, percent_format)
        aliphaticPercentRow += 1
    logging.info('Aliphatic residues as a percent of the entire sequence written to %s worksheet.' % worksheet3Name)

    # Conditional formatting for biochemical analyses.
    for column in range(1, 13):
        worksheet3.conditional_format(2, column, len(conservedList) + 2, column,
                                      {'type': '2_color_scale',
                                       'min_color': '#FFFFFF',
                                       'max_color': '#B77600'}
                                      )

    # Table formatting for biochemical analyses.
    worksheet3.add_table(2, 0, len(conservedList) + 2, 12,
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
    Sg.Popup('''Conservation Analysis program finished running.
\nSee log file for details.''',
             title='Conservation Analysis Completed',
             text_color='#4276ac')
    logging.shutdown()
    window.close()
