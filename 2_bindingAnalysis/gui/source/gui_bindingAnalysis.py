#! python3

##################
#    MODULES
##################

import PySimpleGUI as Sg
import os
import re
import logging
import xlsxwriter
import pandas
import statistics
from Bio import AlignIO
from collections import Counter, OrderedDict


##################
#    CLASSES
##################

# Ordered list of counts.
class OrderedCounter(Counter, OrderedDict):
    pass


##################
#    GUI
##################

# Choose window theme.
Sg.theme('DarkGrey13')

# Establish negative control options.
negControls = ("BSA", "PBS-T", "GST")

# Create window layouts.
infoLayout = [
    # Title and introduction.
    [Sg.Text('\n\n\n\n\n\nPhage Display — Binding Analysis',
             text_color='#8294cc',
             font=('Segoe UI Bold', 22),
             pad=(30, 0)
             )
     ],
    [Sg.Text('Analyse raw phage ELISA with corresponding phage sequence data\n',
             text_color='#8294cc',
             font=('Segoe UI', 12),
             pad=(30, 0)
             )
     ],
    [Sg.Text('''Calculates the average of duplicate absorbances for each protein, normalizes
absorbances against the average of the blanks, and calculates the binder:control
absorbance ratios. ELISA data that don't have corresponding sequencing data are
excluded from the final results.''',
             text_color='#a0a0a2',
             font=('Segoe UI', 10),
             pad=(50, 40)
             )
     ],

    # 'Plate layout' button.
    [Sg.Text('To use this program, ELISA data must correspond to the following plate layout:\n',
             text_color='#a0a0a2',
             font=('Segoe UI', 10),
             pad=(50, 0)
             )
     ],
    [Sg.Button('Plate Layout',
               font=('Segoe UI Bold', 10),
               size=(13, 0),
               pad=(75, 0),
               use_ttk_buttons=True
               )
     ]
]

inputLayout = [
    # ELISA file input prompt.
    [Sg.Text('\n\n1. Enter the full path of the raw ELISA data file:',
             text_color='white',
             font=('Segoe UI Bold', 10),
             pad=(20, 0)
             )
     ],
    [Sg.Input(key='-ELISA_INPUT-',
              size=60,
              pad=(40, 0),
              font=('Segoe UI', 10)
              ),
     Sg.FileBrowse(font=('Segoe UI Bold', 10),
                   size=(10, 0),
                   file_types=(('Excel Files', '*.xlsx'), ('All Files', '*.*'),),
                   )
     ],
    [Sg.Text('''• Must be in *.xlsx format.
• This location will also be the location of the output files.
• Make sure there are no dashes in the name, replace with an underscore if necessary.\n''',
             text_color='#bfbfbf',
             font=('Segoe UI', 10),
             pad=(40, 0)
             )
     ],

    # Input format prompt.
    [Sg.Text('\n2. Choose input format:',
             text_color='white',
             font=('Segoe UI Bold', 10),
             pad=(20, 0)
             )
     ],
    [Sg.Radio('*.xlsx (one file)',
              'FORMAT',
              default=True,
              key='-XLSX_FORMAT-',
              text_color='#bfbfbf',
              font=('Segoe UI Bold', 10),
              pad=(40, 0),
              enable_events=True
              ),
     Sg.Radio('*.fasta (two files)',
              'FORMAT',
              default=False,
              key='-FASTA_FORMAT-',
              text_color='#bfbfbf',
              font=('Segoe UI Bold', 10),
              enable_events=True
              )
     ],
    [Sg.Text('''• Requires xlsx output from Sequence Analysis program.\n''',
             text_color='#bfbfbf',
             font=('Segoe UI', 10),
             key='-RADIO_TEXT-',
             pad=(40, 0)
             )
     ],

    # Xlsx alignment file input prompt (changes to amino acid alignment input prompt depending on radio buttons).
    [Sg.Text('\n3. Enter the full path of the excel alignment file:',
             text_color='white',
             font=('Segoe UI Bold', 10),
             key='-XLSX_INPUT_START_TEXT-',
             pad=(20, 0)
             )
     ],
    [Sg.Input(key='-SWITCH_INPUT-',
              size=60,
              pad=(40, 0),
              font=('Segoe UI', 10),
              do_not_clear=False
              ),
     Sg.FileBrowse(font=('Segoe UI Bold', 10),
                   size=(10, 0),
                   file_types=(('Excel/Fasta Files', '*.xlsx;*.fasta'), ('All Files', '*.*'),)
                   )
     ],
    [Sg.Text('''• Must be in *.fasta format.\n''',
             text_color='#bfbfbf',
             font=('Segoe UI', 10),
             key='-XLSX_INPUT_END_TEXT-',
             pad=(40, 0)
             )
     ],

    # Nucleotide alignment file input prompt.
    [Sg.Text('\n4. Enter the full path of the nucleotide alignment file:',
             text_color='#464646',
             font=('Segoe UI Bold', 10),
             key='-NT_INPUT_START_TEXT-',
             pad=(20, 0)
             )
     ],
    [Sg.Input(key='-NT_INPUT-',
              size=60,
              pad=(40, 0),
              font=('Segoe UI', 10),
              disabled_readonly_background_color='#161616',
              disabled=True,
              do_not_clear=False
              ),
     Sg.FileBrowse(key='-NT_BROWSE-',
                   font=('Segoe UI Bold', 10),
                   size=(10, 0),
                   file_types=(('Fasta Files', '*.fasta'), ('All Files', '*.*'),),
                   disabled=True
                   )
     ],
    [Sg.Text('''• Must be in *.fasta format.\n''',
             text_color='#464646',
             font=('Segoe UI', 10),
             key='-NT_INPUT_END_TEXT-',
             pad=(40, 0)
             )
     ],

    # Negative control selection prompt.
    [Sg.Text('\n5. Select the compound used as the blank:',
             text_color='white',
             font=('Segoe UI Bold', 10),
             pad=(20, 0)
             )
     ],
    [Sg.Combo(key='-BLANK_INPUT-',
              values=negControls,
              default_value=negControls[0],
              size=(15, len(negControls)),
              pad=(40, 5),
              font=('Segoe UI', 10)
              )
     ],

    # Emission absorbance input prompt.
    [Sg.Text('\n6. Enter the wavelength of the emission peak (nm):',
             text_color='white',
             font=('Segoe UI Bold', 10),
             pad=(20, 0)
             )
     ],
    [Sg.Input(key='-EMISSION_INPUT-',
              size=5,
              pad=(40, 5),
              default_text='450',
              font=('Segoe UI', 10)
              )
     ],
    [Sg.Text('''• In most cases it will be 450 nm and does not need to be changed.''',
             text_color='#bfbfbf',
             font=('Segoe UI', 10),
             pad=(40, 0)
             )
     ],

    # Enter button.
    [Sg.Button('Enter',
               bind_return_key=True,
               font=('Segoe UI Bold', 16),
               size=(10, 0),
               pad=(40, 50),
               use_ttk_buttons=True
               )
     ]
]

# TODO: Add to this.
troubleshootLayout = [
    [Sg.Text('''Issue #1:
    Program crashes after pressing "Enter'.''',
             text_color='#f44336',
             font=('Segoe UI', 10),
             pad=((50, 50), (50, 0))
             )
     ],
    [Sg.Text(
        '''Cause:
    The raw ELISA data input is incorrectly formatted.''',
        text_color='#bfbfbf',
        font=('Segoe UI', 10),
        pad=((50, 50), (0, 0))
    )
    ],
    [Sg.Text('''Solutions:
    1. Ensure the correct ELISA file has been chosen.
    2. Ensure the data is in *.xlsx format.
    3. In Gen5, make sure the 'Phage Display' export format has been chosen; this format is specifically made for
       this program, which will not work with any other differently formatted data.''',
             text_color='#93c47d',
             font=('Segoe UI', 10),
             pad=((50, 50), (0, 0))
             )
     ],
]

notesLayout = [
    [Sg.Text(
        '''• Bases that are lowercase in the final output indicate bases that were manually called
   due to being miscalled or not called at all as a result of machine error.
    
• Blanks are not coated with the target protein and contain BSA or another
  appropriate control protein.

• Controls are not coated with the target protein and contain their corresponding
  phage.''',
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
        size=(650, 850)
    )
    ]
]

# Set up window behaviour.
window = Sg.Window('Phage Display - Binding Analysis',
                   tabGroup,
                   alpha_channel=0.95,
                   size=(650, 850),
                   ttk_theme='clam'
                   )

# Create a while loop that keeps the window open.
while True:
    event, values = window.read()

    # Break the loop and close the window if 'Exit' or close window buttons are pressed.
    if event == Sg.WIN_CLOSED:
        window.close()
        break

    # Popup will show the plate layout picture if 'Plate Layout' is pressed.
    elif event == 'Plate Layout':
        # Make sure the path for the image corresponds to the location of the image in the executable's folder.
        Sg.Popup(image='.\\images\\plateLayout.png',
                 title='Plate Layout'
                 )

    # Update elements based on radio buttons.
    elif event == '-XLSX_FORMAT-':
        window['-RADIO_TEXT-'].update('    * Requires *.xlsx output from Sequence Analysis program.\n')
        window['-XLSX_INPUT_START_TEXT-'].update('2. Enter the full path of the excel alignment file:')
        window['-XLSX_INPUT_END_TEXT-'].update('    * Must be in *.xlsx format.\n')
        window['-NT_INPUT_START_TEXT-'].update(text_color='#464646')
        window['-NT_BROWSE-'].update(disabled=True)
        window['-NT_INPUT-'].update(disabled=True)
        window['-NT_INPUT_END_TEXT-'].update(text_color='#464646')
        continue

    elif event == '-FASTA_FORMAT-':
        window['-RADIO_TEXT-'].update('    * Requires *.fasta alignment output from sequence analysis program or'
                                      ' elsewhere.\n')
        window['-XLSX_INPUT_START_TEXT-'].update('2. Enter the full path of the amino acid alignment file:')
        window['-XLSX_INPUT_END_TEXT-'].update('    * Must be in *.fasta format.\n')
        window['-NT_INPUT_START_TEXT-'].update(text_color='white')
        window['-NT_BROWSE-'].update(disabled=False)
        window['-NT_INPUT-'].update(disabled=False)
        window['-NT_INPUT_END_TEXT-'].update(text_color='#bfbfbf')
        continue

    # Updates variables with input values when 'Enter' button is pressed.
    elif event == 'Enter':
        inputOptions = {'1': 'xlsx file',
                        '2': 'fasta files'
                        }
        elisaInFilePath = str(values['-ELISA_INPUT-'].replace('\\',
                                                              '/')
                              )
        blankID = str(values['-BLANK_INPUT-'])
        emissionAbs = str(values['-EMISSION_INPUT-'])

        # Stops user if ELISA filetype isn't valid and prompts to retry.
        if not re.search('.xlsx$', elisaInFilePath):
            Sg.Popup('''The entered file is not an xlsx file.
Please choose the file again.''',
                     title='File Not Found',
                     grab_anywhere=True,
                     text_color='#4276ac'
                     )
            continue
        # Stops user if ELISA file path isn't valid and prompts to retry.
        if not os.path.exists(elisaInFilePath):
            Sg.Popup('''The entered ELISA file does not exist in this location.
Please choose the file again.''',
                     title='File Not Found',
                     grab_anywhere=True,
                     text_color='#8294cc'
                     )
            continue
        else:
            path = re.sub(r'[a-zA-Z0-9_]+\.xlsx$',
                          '',
                          elisaInFilePath
                          )
            os.chdir(path)

        # Stops user if alignment file path isn't valid and prompts to retry.
        if values['-XLSX_FORMAT-']:
            inputFormat = '1'
            xlsxInFilePath = str(values['-SWITCH_INPUT-'].replace('\\',
                                                                  '/')
                                 )
            if not os.path.exists(xlsxInFilePath):
                Sg.Popup('''The alignment file does not exist in this location.
Please choose the file again.''',
                         title='File Not Found',
                         grab_anywhere=True,
                         text_color='#4276ac')
                continue

        # Stops user if amino acid/nucleotide alignment paths aren't valid and prompts to retry.
        if values['-FASTA_FORMAT-']:
            inputFormat = '2'
            aaInFilePath = str(values['-SWITCH_INPUT-'].replace('\\',
                                                                '/')
                               )
            if not os.path.exists(aaInFilePath):
                Sg.Popup('''The amino acid alignment file does not exist in this location.
Please choose the file again.''',
                         title='File Not Found',
                         grab_anywhere=True,
                         text_color='#4276ac')
                continue
            aaInFileName = re.findall(r'[a-zA-Z0-9]+\.fasta$',
                                      aaInFilePath
                                      )

            ntInFilePath = str(values['-NT_INPUT-'].replace('\\',
                                                            '/')
                               )
            if not os.path.exists(ntInFilePath):
                Sg.Popup('''The nucleotide alignment file does not exist in this location.
Please choose the file again.''',
                         title='File Not Found',
                         grab_anywhere=True,
                         text_color='#4276ac')
                continue
            ntInFileName = re.findall(r'[a-zA-Z0-9]+\.fasta$',
                                      ntInFilePath
                                      )

        # Stops user if emission absorbance isn't valid and prompts to retry.
        emissionAbsInput = re.match(r'\d{3}',
                                    emissionAbs
                                    )
        if emissionAbsInput is None:
            Sg.Popup('''Invalid input for emission absorbance.
Please make sure there are three digits and try again.''',
                     title='Invalid Emission Absorbance Input',
                     grab_anywhere=True,
                     text_color='#8294cc'
                     )
            continue

        else:
            elisaInFileName = re.findall(r'[a-zA-Z0-9_]+\.xlsx$',
                                         elisaInFilePath
                                         )
            elisaInFileName = elisaInFileName[0]
            emissionAbs = int(emissionAbs)
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
    # Logging setup.
    elisaInFileNameShort = re.sub(r'[.].*', '', elisaInFileName)
    logging.basicConfig(filename=path + '/' + elisaInFileNameShort + '.log',
                        level=logging.INFO,
                        format='%(asctime)s - %(message)s',
                        filemode='w'
                        )
    logging.info('Working directory changed to %s.' % path)
    logging.info('%s chosen as raw ELISA data source.' % elisaInFileName)

    if inputFormat == '1':
        logging.info('%s chosen as the amino acid and nucleotide sequence source.' % xlsxInFilePath)

    elif inputFormat == '2':
        logging.info('%s chosen as the amino acid sequence source.' % aaInFilePath)
        logging.info('%s chosen as the nucleotide sequence source.' % ntInFilePath)

    logging.info('%s chosen as the blank compound.' % blankID)
    logging.info('%i nm chosen as emission wavelength.' % emissionAbs)

    ##################
    # Retrieve and parse raw ELISA data.
    ##################

    allCells = pandas.read_excel(elisaInFilePath)
    logging.info('Raw data read from ELISA file.')
    # Remove rows where the last column isn't equal to the emission absorbance.
    lastColName = 'Unnamed: ' + str(allCells.shape[1] - 1)
    emissionAbs = 450
    dataCellsRaw = allCells[allCells[lastColName] == emissionAbs]
    # Remove the emission absorbance column.
    dataCellsRaw = dataCellsRaw.iloc[:, :-1]
    dataCellsRaw = dataCellsRaw.dropna(axis=0, thresh=2)
    dataCellsRaw = dataCellsRaw.dropna(axis=1, how='all')
    # Remove rows that contain more than two NaN values.
    dataCellsClean = dataCellsRaw.dropna(axis=0, thresh=2)
    # Replace 'OVRFLW' wells with a float of 4.
    dataCellsClean = dataCellsClean.replace('OVRFLW', float(4))
    # Remove designation column.
    dataCellsClean = dataCellsClean.iloc[:, 1:]

    # Rename rows to make parsing easier.
    dataCellsClean.index = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P']

    # Create a list of all the data absorbances.
    cellValues = [dataCellsClean.loc['A'][value] for value in range(0, 24)] + \
                 [dataCellsClean.loc['C'][value] for value in range(0, 24)] + \
                 [dataCellsClean.loc['E'][value] for value in range(0, 24)] + \
                 [dataCellsClean.loc['G'][value] for value in range(0, 24)] + \
                 [dataCellsClean.loc['I'][value] for value in range(0, 24)] + \
                 [dataCellsClean.loc['K'][value] for value in range(0, 24)] + \
                 [dataCellsClean.loc['M'][value] for value in range(0, 24)] + \
                 [dataCellsClean.loc['O'][value] for value in range(0, 24)]
    logging.info('Data absorbances retrieved from raw data file.')

    # Retrieve control wells.
    controlCells = [dataCellsClean.loc['B'][value] for value in range(1, 24, 2)] + \
                   [dataCellsClean.loc['D'][value] for value in range(1, 24, 2)] + \
                   [dataCellsClean.loc['F'][value] for value in range(1, 24, 2)] + \
                   [dataCellsClean.loc['H'][value] for value in range(1, 24, 2)] + \
                   [dataCellsClean.loc['J'][value] for value in range(1, 24, 2)] + \
                   [dataCellsClean.loc['L'][value] for value in range(1, 24, 2)] + \
                   [dataCellsClean.loc['N'][value] for value in range(1, 24, 2)] + \
                   [dataCellsClean.loc['P'][value] for value in range(1, 24, 2)]
    logging.info('Control absorbances retrieved from raw data file.')

    ##################
    # Extract blank values.
    ##################

    # Retrieve blank wells and average.
    blankCells = [dataCellsClean.loc['B'][value] for value in range(0, 23, 2)]
    # Remove values that contain no data in the Excel cell.
    blankValues = [absorbance for absorbance in blankCells if not (pandas.isna(absorbance))]
    # Average blanks.
    blankAve = statistics.mean(blankValues)
    logging.info('Blank absorbances retrieved from raw data file.')
    logging.info('Blank values averaged.')

    ##################
    # Extract control values.
    ##################

    elisaPlateIDs = ['A01', 'A02', 'A03', 'A04', 'A05', 'A06', 'A07', 'A08', 'A09', 'A10', 'A11', 'A12', 'B01', 'B02',
                     'B03', 'B04', 'B05', 'B06', 'B07', 'B08', 'B09', 'B10', 'B11', 'B12', 'C01', 'C02', 'C03', 'C04',
                     'C05', 'C06', 'C07', 'C08', 'C09', 'C10', 'C11', 'C12', 'D01', 'D02', 'D03', 'D04', 'D05', 'D06',
                     'D07', 'D08', 'D09', 'D10', 'D11', 'D12', 'E01', 'E02', 'E03', 'E04', 'E05', 'E06', 'E07', 'E08',
                     'E09', 'E10', 'E11', 'E12', 'F01', 'F02', 'F03', 'F04', 'F05', 'F06', 'F07', 'F08', 'F09', 'F10',
                     'F11', 'F12', 'G01', 'G02', 'G03', 'G04', 'G05', 'G06', 'G07', 'G08', 'G09', 'G10', 'G11', 'G12',
                     'H01', 'H02', 'H03', 'H04', 'H05', 'H06', 'H07', 'H08', 'H09', 'H10', 'H11', 'H12'
                     ]
    controlDict = {key: value for key, value in zip(elisaPlateIDs, controlCells)}

    ##################
    # Average paired ELISA values for each sequence and normalise to the blank/negative control average.
    ##################

    # Average paired values and append to a new list.
    cellAveList = [statistics.mean([cellValues[value], cellValues[value + 1]]) for value in
                   range(0, len(cellValues), 2)]
    logging.info('Paired ELISA absorbances averaged.')

    # Normalise ELISA scores to the blank/negative control average.
    relAveList = [(value / blankAve) for value in cellAveList]
    logging.info('Averaged absorbances normalised to the blanks/negative control average.')

    ##################
    # Retrieve and parse amino acid/nucleotide sequence data.
    ##################

    if inputFormat == '1':
        aaCells = pandas.read_excel(xlsxInFilePath, sheet_name=0)
        ntCells = pandas.read_excel(xlsxInFilePath, sheet_name=2)

        aaSeqRegex = re.compile(r'[ARNDCEQGHILKMFPSTWYVX]{10,}')
        ntSeqRegex = re.compile('[AGTCNagtc]{10,}')

        # Retrieve amino acid sequences.
        aaSeqCells = aaCells.iloc[1:, 1:]
        aaSeqCells = aaSeqCells.to_string(index=False)
        aaSeqList = aaSeqCells.replace(' ', '')
        aaSeqList = aaSeqRegex.findall(aaSeqList)
        logging.info('Amino acid sequences retrieved from %s.' % xlsxInFilePath)
        aaAlignLen = len(aaSeqList[0])
        logging.info('Amino acid alignment length calculated to be %s.' % aaAlignLen)

        # Retrieve nucleotide sequences.
        ntSeqCells = ntCells.iloc[1:, 1:]
        ntSeqCells = ntSeqCells.to_string(index=False)
        ntSeqList = ntSeqCells.replace(' ', '')
        ntSeqList = ntSeqRegex.findall(ntSeqList)
        logging.info('Nucleotide sequences retrieved from %s.' % xlsxInFilePath)
        ntAlignLen = len(ntSeqList[0])
        logging.info('Nucleotide alignment length calculated to be %s.' % ntAlignLen)

        # Retrieve amino acid sequence names.
        aaNameCells = aaCells.iloc[1:, 0:1]
        aaNameCells = aaNameCells.to_string(index=False)
        aaNameList = re.sub(r'(_M\w*)',
                            '',
                            aaNameCells
                            )
        aaNameList = aaNameList.replace(' ', '')
        aaNameList = aaNameList.split('\n')
        aaNameList = aaNameList[1:]
        logging.info('Amino acid sequence names retrieved from %s.' % xlsxInFilePath)

        # Retrieve nucleotide sequence names.
        ntNameCells = ntCells.iloc[1:, 0:1]
        ntNameCells = ntNameCells.to_string(index=False)
        ntNameList = re.sub(r'(_M\w*)',
                            '',
                            ntNameCells
                            )
        ntNameList = ntNameList.replace(' ', '')
        ntNameList = ntNameList.split('\n')
        ntNameList = ntNameList[1:]
        logging.info('Nucleotide sequence names retrieved from %s.' % xlsxInFilePath)

    elif inputFormat == '2':

        #########
        # Amino Acids
        #########

        # Retrieve amino acid sequence names.
        aaSeqRegex = re.compile(r'([ARNDCEQGHILKMFPSTWYVX]{10,})')
        with open(aaInFilePath, 'r') as file:
            aaAllLines = file.read()
            # Remove primer name and 'aaTrimmed' from fasta name.
            aaAllLinesTrim = re.sub(r'(_M\w*)',
                                    '',
                                    aaAllLines
                                    )
            aaNameList = re.findall(r'>(.*)',
                                    aaAllLinesTrim
                                    )
        logging.info('Amino acid sequence names retrieved from %s.' % aaInFileName)

        # Retrieve amino acid sequences.
        aaAllLinesClean = aaAllLines.replace('\n',
                                             ''
                                             )
        aaSeqList = aaSeqRegex.findall(aaAllLinesClean)
        logging.info('Amino acid sequences retrieved from %s.' % aaInFileName)

        # Retrieve amino acid alignment length.
        aaAlignment = AlignIO.read(aaInFilePath,
                                   'fasta'
                                   )
        aaAlignLen = aaAlignment.get_alignment_length()
        logging.info('Amino acid alignment length calculated to be %s.' % aaAlignLen)

        #########
        # Nucleotides
        #########

        # Retrieve nucleotide sequence names.
        with open(ntInFilePath, 'r') as file:
            ntAllLines = file.read()
            # Remove primer name and 'aaTrimmed' from fasta name.
            ntAllLines = re.sub(r'(_M\w*)',
                                '',
                                ntAllLines
                                )
            ntNameList = re.findall(r'>(.*)',
                                    ntAllLines
                                    )
        logging.info('Nucleotide sequence names retrieved from %s.' % ntInFileName)

        # Retrieve nucleotide sequences.
        ntSeqRegex = re.compile('[AGTCURYNWSMKBHDV]{10,}')
        ntAllLinesClean = ntAllLines.replace('\n',
                                             ''
                                             )
        ntSeqList = ntSeqRegex.findall(ntAllLinesClean)
        logging.info('Nucleotide sequences retrieved from %s.' % ntInFileName)

        # Retrieve nucleotide alignment length.
        ntAlignment = AlignIO.read(ntInFilePath,
                                   'fasta'
                                   )
        ntAlignLen = ntAlignment.get_alignment_length()
        logging.info('Nucleotide alignment length calculated to be %s.' % ntAlignLen)

    ##################
    # Correlate sequence data with ELISA data and remove ELISA data that don't have sequencing counterparts.
    ##################

    # Create a dictionary of all possible ELISA IDs and assign corresponding numbers for indexing.
    wellDict = {key: value for key, value in zip(elisaPlateIDs, range(0, 96))}

    #########
    # Amino acids
    #########

    wellRegex = re.compile(r'[A-H][0-9]{2}')

    # Determine the amino acid IDs present in the sequencing data.
    aaPlateIDs = [wellRegex.findall(name) for name in aaNameList]
    # Turn list of lists into a flat list.
    aaSeqPlateIDs = []
    for sublist in aaPlateIDs:
        for ID in sublist:
            aaSeqPlateIDs.append(ID)

    # Retrieve ELISA absorbances for only the IDs present in the sequencing data.
    aaRawListShort = [cellAveList[wellDict.get(ID)] for ID in aaSeqPlateIDs]
    aaRelAveListShort = [relAveList[wellDict.get(ID)] for ID in aaSeqPlateIDs]
    logging.info('ELISA results without corresponding amino acid sequences removed from analysis.')

    # Relate IDs to their respective ELISA absorbances.
    aaShortNameList = []
    for name in aaNameList:
        if re.findall(r'([A-H])(\d$)', name):
            newName = re.sub(r'([A-H])(\d$)',
                             r'\g<1>0\g<2>',
                             name
                             )
            aaShortNameList.append(newName)
        else:
            aaShortNameList.append(name)
    aaShortNameList.sort()

    # Create a dictionary with the new reorganised names and sequences.
    aaSeqDict = dict(zip(aaShortNameList,
                         aaSeqList)
                     )

    # Create a dictionary with the new reorganised names and ELISA absorbances.
    aaAbsDict = dict(zip(aaShortNameList,
                         aaRelAveListShort)
                     )

    # Create a list of unique amino acid sequences ordered by frequency.
    aaUnique = OrderedCounter(aaSeqList)
    aaUnique = aaUnique.most_common()
    aaUniqueDict = dict(aaUnique)
    logging.info('Dictionary of unique amino acid sequences created.')

    # Remove controls that don't have sequencing counterparts.
    aaControlListRaw = [controlDict[well] for well in aaSeqPlateIDs]
    aaControlListRel = [(absorbance / blankAve) for absorbance in aaControlListRaw]

    #########
    # Nucleotides
    #########

    # Determine the nucleotide IDs present in the sequencing data.
    ntPlateIDs = [wellRegex.findall(name) for name in ntNameList]
    # Turn list of lists into a flat list.
    ntSeqPlateIDs = []
    for sublist in ntPlateIDs:
        for item in sublist:
            ntSeqPlateIDs.append(item)

    # Retrieve ELISA absorbances for only the IDs present in the sequencing data.
    ntRawListShort = [cellAveList[wellDict.get(ID)] for ID in ntSeqPlateIDs]
    ntRelAveListShort = [relAveList[wellDict.get(ID)] for ID in ntSeqPlateIDs]
    logging.info('ELISA results without corresponding nucleotide sequences removed from analysis.')

    # Relate IDs to their respective ELISA absorbances.
    ntShortNameList = []
    for name in ntNameList:
        if re.findall(r'([A-H])(\d$)', name):
            newName = re.sub(r'([A-H])(\d$)',
                             r'\g<1>0\g<2>',
                             name
                             )
            ntShortNameList.append(newName)
        else:
            ntShortNameList.append(name)
    ntShortNameList.sort()

    # Create a dictionary with the new reorganised names and sequences.
    ntSeqDict = dict(zip(ntShortNameList,
                         ntSeqList)
                     )

    # Create a dictionary with the new reorganised names and ELISA absorbances.
    ntAbsDict = dict(zip(ntShortNameList,
                         ntRelAveListShort)
                     )

    # Create a list of unique amino acid sequences ordered by frequency.
    ntUnique = OrderedCounter(ntSeqList)
    ntUnique = ntUnique.most_common()
    ntUniqueDict = dict(ntUnique)
    logging.info('Dictionary of unique nucleotide sequences created.')

    # Remove ELISA data that don't have sequencing counterparts.
    ntControlListRaw = [controlDict[well] for well in ntSeqPlateIDs]
    ntControlListRel = [(absorbance / blankAve) for absorbance in ntControlListRaw]

    ##################
    # Order sequences and wells so they can be attributed to unique sequences. Necessary for subsequent statistics.
    ##################

    #########
    # Amino acids
    #########

    # Create ordered list of wells that correspond to unique sequences.
    aaOrderedSeq = [key for key in aaUniqueDict.keys()]
    aaOrderedNames = []
    for seq in aaOrderedSeq:
        for key, value in aaSeqDict.items():
            if seq in value:
                aaOrderedNames.append(key)
    logging.info('Ordered index of amino acid well IDs created.')

    # Create ordered list of absorbances that correspond to unique sequences.
    aaOrderedAbs = []
    for index in aaOrderedNames:
        for ID, score in aaAbsDict.items():
            if index == ID:
                aaOrderedAbs.append(score)
    logging.info('Ordered list of amino acid absorbances created.')

    # Associate specific well IDs with corresponding unique sequence.
    aaCountID = []
    begin = 0
    for uniqueSeq, count in aaUniqueDict.items():
        end = int(count) + begin
        aaCountID.append(aaOrderedNames[begin:end])
        begin += count
    logging.info('List of specific well IDs associated with amino acid sequences created.')

    #########
    # Nucleotides
    #########

    # Get ordered list of wells that correspond to unique sequences.
    ntOrderedSeq = [key for key in ntUniqueDict.keys()]
    ntOrderedNames = []
    for seq in ntOrderedSeq:
        for key, value in ntSeqDict.items():
            if seq in value:
                ntOrderedNames.append(key)
    logging.info('Ordered index of nucleotide well IDs created.')

    # Get ordered values corresponding to ordered list of wells for unique sequences; necessary for subsequent
    # statistics.
    ntOrderedAbs = []
    for index in ntOrderedNames:
        for ID, score in ntAbsDict.items():
            if index == ID:
                ntOrderedAbs.append(score)
    logging.info('Ordered list of nucleotide absorbances created.')

    # Associate specific well IDs with corresponding unique sequence.
    ntCountID = []
    begin = 0
    for uniqueSeq, count in ntUniqueDict.items():
        end = int(count) + begin
        ntCountID.append(ntOrderedNames[begin:end])
        begin += count
    logging.info('List of specific well IDs associated with nucleotide sequences created.')

    ##################
    # Compare binder and control absorbances and calculate statistics.
    ##################

    #########
    # Amino acids
    #########

    # Normalise absorbances to controls.
    aaBinderControlRatio = [(binder / control) for binder, control in zip(aaRawListShort, aaControlListRaw)]
    aaNameRatioDict = dict(zip(aaShortNameList,
                               aaBinderControlRatio)
                           )
    # Create ordered list of ratios.
    aaOrderedRatios = []
    for name in aaOrderedNames:
        for ID, score in aaNameRatioDict.items():
            if name == ID:
                aaOrderedRatios.append(score)

    # Create ordered list of absorbances that correspond to unique sequences.
    aaCounts = [count for ratio, count in aaUniqueDict.items()]
    aaRatioCountDict = dict(zip(aaBinderControlRatio,
                                aaCounts)
                            )

    # Retrieve max absorbance for ordered values.
    aaRatioMaxList = []
    begin = 0
    for ratio, count in aaRatioCountDict.items():
        end = int(count) + begin
        uniqueMax = max(aaOrderedRatios[begin:end])
        aaRatioMaxList.append(uniqueMax)
        begin += count
    logging.info('List of unique amino acid maximum absorbances created.')

    # Retrieve min absorbance for ordered values.
    aaRatioMinList = []
    begin = 0
    for ratio, count in aaRatioCountDict.items():
        end = int(count) + begin
        uniqueMin = min(aaOrderedRatios[begin:end])
        aaRatioMinList.append(uniqueMin)
        begin += count
    logging.info('List of unique amino acid minimum absorbances created.')

    # Retrieve median absorbance for ordered values.
    aaRatioMedianList = []
    begin = 0
    for ratio, count in aaRatioCountDict.items():
        end = int(count) + begin
        uniqueMedian = statistics.median(aaOrderedRatios[begin:end])
        aaRatioMedianList.append(uniqueMedian)
        begin += count
    logging.info('List of unique amino acid median absorbances created.')

    # Retrieve mean absorbance for ordered values.
    aaRatioMeanList = []
    begin = 0
    for ratio, count in aaRatioCountDict.items():
        end = int(count) + begin
        uniqueMean = statistics.mean(aaOrderedRatios[begin:end])
        aaRatioMeanList.append(uniqueMean)
        begin += count
    logging.info('List of unique amino acid mean absorbances created.')

    # Retrieve standard deviation of absorbances for ordered values.
    aaRatioStdevList = []
    begin = 0
    for ratio, count in aaRatioCountDict.items():
        end = int(count) + begin
        try:
            uniqueStdev = statistics.stdev(aaOrderedRatios[begin:end])
            aaRatioStdevList.append(uniqueStdev)
        # Above statistic won't work if only a single value so append '0.0' value to the list.
        except statistics.StatisticsError:
            aaRatioStdevList.append(0)
        begin += count
    logging.info('List of unique amino acid absorbance standard deviations created.')

    aaRatioStatsTable = {'Ratio Max': aaRatioMaxList,
                         'Ratio Min': aaRatioMinList,
                         'Ratio Median': aaRatioMedianList,
                         'Ratio Mean': aaRatioMeanList,
                         'Ratio St Dev': aaRatioStdevList
                         }

    aaRatioStatsDataframe = pandas.DataFrame(aaRatioStatsTable)

    #########
    # Nucleotides
    #########

    # Normalise absorbances to controls.
    ntBinderControlRatio = [(binder / control) for binder, control in zip(ntRawListShort, ntControlListRaw)]
    ntNameRatioDict = dict(zip(ntShortNameList,
                               ntBinderControlRatio)
                           )
    # Create ordered list of ratios.
    ntOrderedRatios = []
    for name in ntOrderedNames:
        for ID, score in ntNameRatioDict.items():
            if name == ID:
                ntOrderedRatios.append(score)

    # Create ordered list of absorbances that correspond to unique sequences.
    ntCounts = []
    for ratio, count in ntUniqueDict.items():
        ntCounts.append(count)
    ntRatioCountDict = dict(zip(ntBinderControlRatio,
                                ntCounts)
                            )

    # Retrieve max absorbance for ordered values.
    ntRatioMaxList = []
    begin = 0
    for ratio, count in ntRatioCountDict.items():
        end = int(count) + begin
        uniqueMax = max(ntOrderedRatios[begin:end])
        ntRatioMaxList.append(uniqueMax)
        begin += count
    logging.info('List of unique amino acid maximum absorbances created.')

    # Retrieve min absorbance for ordered values.
    ntRatioMinList = []
    begin = 0
    for ratio, count in ntRatioCountDict.items():
        end = int(count) + begin
        uniqueMin = min(ntOrderedRatios[begin:end])
        ntRatioMinList.append(uniqueMin)
        begin += count
    logging.info('List of unique amino acid minimum absorbances created.')

    # Retrieve median absorbance for ordered values.
    ntRatioMedianList = []
    begin = 0
    for ratio, count in ntRatioCountDict.items():
        end = int(count) + begin
        uniqueMedian = statistics.median(ntOrderedRatios[begin:end])
        ntRatioMedianList.append(uniqueMedian)
        begin += count
    logging.info('List of unique amino acid median absorbances created.')

    # Retrieve mean absorbance for ordered values.
    ntRatioMeanList = []
    begin = 0
    for ratio, count in ntRatioCountDict.items():
        end = int(count) + begin
        uniqueMean = statistics.mean(ntOrderedRatios[begin:end])
        ntRatioMeanList.append(uniqueMean)
        begin += count
    logging.info('List of unique amino acid mean absorbances created.')

    # Retrieve standard deviation of absorbances for ordered values.
    ntRatioStdevList = []
    begin = 0
    for ratio, count in ntRatioCountDict.items():
        end = int(count) + begin
        try:
            uniqueStdev = statistics.stdev(ntOrderedRatios[begin:end])
            ntRatioStdevList.append(uniqueStdev)
        # Above statistic won't work if only a single value so append '0.0' value to the list.
        except statistics.StatisticsError:
            ntRatioStdevList.append(0)
        begin += count
    logging.info('List of unique amino acid absorbance standard deviations created.')

    ntRatioStatsTable = {'Ratio Max': ntRatioMaxList,
                         'Ratio Min': ntRatioMinList,
                         'Ratio Median': ntRatioMedianList,
                         'Ratio Mean': ntRatioMeanList,
                         'Ratio St Dev': ntRatioStdevList
                         }

    ntRatioStatsDataframe = pandas.DataFrame(ntRatioStatsTable)

    ##################
    # Export data as a single xlsx file.
    ##################

    # Create workbook.
    workbook = xlsxwriter.Workbook(path + '/' + elisaInFileNameShort + '_bindingAnalysis.xlsx')
    logging.info('Excel spreadsheet created as "%s.xlsx".' % elisaInFileName)

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

    wellTitle_format = workbook.add_format({'bold': True,
                                            'font_size': 10,
                                            'font_name': 'Segoe UI'
                                            }
                                           )
    wellTitle_format.set_align('left')
    wellTitle_format.set_align('vcenter')

    # Absorbance.
    abs_format = workbook.add_format({'num_format': '#,##0.00',
                                      'font_size': 10,
                                      'font_name': 'Segoe UI'
                                      }
                                     )
    abs_format.set_align('center')
    abs_format.set_align('vcenter')

    # Statistics.
    stats_format = workbook.add_format({'num_format': '#,##0.0',
                                        'font_size': 10,
                                        'font_name': 'Segoe UI'
                                        }
                                       )
    stats_format.set_align('center')
    stats_format.set_align('vcenter')

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

    # Residue numbers.
    residue_format = workbook.add_format({'font_size': 8,
                                          'font_name': 'Segoe UI'
                                          }
                                         )
    residue_format.set_align('center')
    residue_format.set_align('vcenter')

    # Sequences.
    sequence_format = workbook.add_format({'font_size': 9,
                                           'font_name': 'Segoe UI'
                                           }
                                          )
    sequence_format.set_align('center')
    sequence_format.set_align('vcenter')

    # Information.
    info_format = workbook.add_format({'font_size': 10,
                                       'font_name': 'Segoe UI'
                                       }
                                      )
    info_format.set_align('left')
    info_format.set_align('vcenter')

    logging.info('Cell formatting rules set.')

    ##################
    # Create worksheet for all amino acid sequences.
    ##################

    # Initialize sheet characteristics.
    worksheet1Name = 'All AA Seq'
    worksheet1 = workbook.add_worksheet(worksheet1Name)
    worksheet1.hide_gridlines(option=2)
    idColWidth = round(len(aaShortNameList[0]) * 1.4)
    worksheet1.set_column(0, 0, idColWidth)
    worksheet1.set_column(1, aaAlignLen, 3)
    worksheet1.set_column(aaAlignLen + 1, aaAlignLen + 4, 12)
    worksheet1.set_column(aaAlignLen + 5, aaAlignLen + 5, 15)
    worksheet1.freeze_panes(2, 1)
    logging.info('%s worksheet created.' % worksheet1Name)

    # Write well IDs.
    worksheet1.merge_range(0, 0, 1, 0, 'ID', title_format)
    idRow = 2
    for name in aaShortNameList:
        worksheet1.write(idRow, 0, name, wellID_format)
        idRow += 1
    logging.info('Amino acid well IDs written to %s worksheet.' % worksheet1Name)

    # Write amino acid sequences.
    worksheet1.merge_range(0, 1, 0, 9, 'Amino Acid Sequence', title_format)
    seqRow = 2
    seqCol = 1
    for aa in aaSeqList:
        letterList = list(aa)
        for letter in letterList:
            worksheet1.write(seqRow, seqCol, letter, sequence_format)
            seqCol += 1
        seqRow += 1
        seqCol = 1
    logging.info('All Amino acid sequences written to %s worksheet.' % worksheet1Name)

    # Write raw binder absorbances.
    absRow = 2
    absCol = aaAlignLen + 1
    worksheet1.merge_range(0, absCol, 0, absCol + 1, 'Raw', title_format)
    worksheet1.merge_range(0, absCol + 2, 0, absCol + 3, 'Normalised', title_format)
    worksheet1.write(1, absCol, 'Binder', title_format)
    for absorbance in aaRawListShort:
        worksheet1.write(absRow, absCol, absorbance, abs_format)
        absRow += 1
    logging.info('All raw binder absorbances written to %s worksheet.' % worksheet1Name)

    # Write raw control absorbances.
    absRow = 2
    absCol = aaAlignLen + 2
    worksheet1.write(1, absCol, 'Control', title_format)
    for absorbance in aaControlListRaw:
        worksheet1.write(absRow, absCol, absorbance, abs_format)
        absRow += 1
    logging.info('All raw control absorbances written to %s worksheet.' % worksheet1Name)

    # Write averaged binder absorbances.
    absRow = 2
    absCol = aaAlignLen + 3
    worksheet1.write(1, absCol, 'Binder', title_format)
    for absorbance in aaRelAveListShort:
        worksheet1.write(absRow, absCol, absorbance, abs_format)
        absRow += 1
    logging.info('All normalised binder absorbances written to %s worksheet.' % worksheet1Name)

    # Write normalised control absorbances.
    absRow = 2
    absCol = aaAlignLen + 4
    worksheet1.write(1, absCol, 'Control', title_format)
    for absorbance in aaControlListRel:
        worksheet1.write(absRow, absCol, absorbance, abs_format)
        absRow += 1
    logging.info('All normalised control absorbances written to %s worksheet.' % worksheet1Name)

    # Write binder:control absorbance ratios.
    absRow = 2
    absCol = aaAlignLen + 5
    worksheet1.merge_range(0, absCol, 1, absCol, 'Binder:Control', title_format)
    for absorbance in aaBinderControlRatio:
        worksheet1.write(absRow, absCol, absorbance, stats_format)
        absRow += 1
    logging.info('All binder:control absorbance ratios written to %s worksheet.' % worksheet1Name)

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
    worksheet2.set_column(aaAlignLen + 1, aaAlignLen + 5, 8)
    worksheet2.freeze_panes(2, 1)
    logging.info('%s worksheet created.' % worksheet2Name)

    # Write unique amino acid sequences.
    worksheet2.merge_range(0, 1, 0, 9, 'Amino Acid Sequence', title_format)
    uniqueRow = 2
    uniqueCol = 1
    for seq in aaUniqueDict.keys():
        letterList = list(seq)
        for letter in letterList:
            worksheet2.write(uniqueRow, uniqueCol, letter, sequence_format)
            uniqueCol += 1
        uniqueRow += 1
        uniqueCol = 1
    logging.info('Unique amino acid sequences written to %s worksheet.' % worksheet2Name)

    # Add counts for each unique amino acid sequence.
    countRow = 2
    countCol = aaAlignLen + 1
    worksheet2.merge_range(0, aaAlignLen + 1, 1, aaAlignLen + 1, 'Count', title_format)
    count = list(aaUniqueDict.values())
    for number in count:
        worksheet2.write_number(countRow, countCol, number, general_format)
        countRow += 1
    logging.info('Amino acid sequence counts written to %s worksheet.' % worksheet2Name)

    # Write amino acid residue numbers above sequences.
    residueCol = 1
    for number in aaResList:
        worksheet2.write(1, residueCol, number, residue_format)
        residueCol += 1
    logging.info('Residue numbers written to %s worksheet.' % worksheet2Name)

    ##################
    # Write statistics to unique amino acid sequences worksheet.
    ##################

    # Max.
    maxRow = 2
    maxCol = aaAlignLen + 2
    worksheet2.merge_range(0, aaAlignLen + 2, 1, aaAlignLen + 2, 'Max.', title_format)
    for seq in aaRatioMaxList:
        worksheet2.write(maxRow, maxCol, seq, stats_format)
        maxRow += 1
    logging.info('List of unique amino acid maximum absorbances written to %s worksheet.' % worksheet2Name)

    # Min.
    minRow = 2
    minCol = aaAlignLen + 3
    worksheet2.merge_range(0, aaAlignLen + 3, 1, aaAlignLen + 3, 'Min.', title_format)
    for seq in aaRatioMinList:
        worksheet2.write(minRow, minCol, seq, stats_format)
        minRow += 1
    logging.info('List of unique amino acid minimum absorbances written to %s worksheet.' % worksheet2Name)

    # Median.
    medianRow = 2
    medianCol = aaAlignLen + 4
    worksheet2.merge_range(0, aaAlignLen + 4, 1, aaAlignLen + 4, 'Median', title_format)
    for seq in aaRatioMedianList:
        worksheet2.write(medianRow, medianCol, seq, stats_format)
        medianRow += 1
    logging.info('List of unique amino acid median absorbances written to %s worksheet.' % worksheet2Name)

    # Mean.
    meanRow = 2
    meanCol = aaAlignLen + 5
    worksheet2.merge_range(0, aaAlignLen + 5, 1, aaAlignLen + 5, 'Mean', title_format)
    for seq in aaRatioMeanList:
        worksheet2.write(meanRow, meanCol, seq, stats_format)
        meanRow += 1
    logging.info('List of unique amino acid mean absorbances written to %s worksheet.' % worksheet2Name)

    # Standard deviation.
    stdevRow = 2
    stdevCol = aaAlignLen + 6
    worksheet2.merge_range(0, aaAlignLen + 6, 1, aaAlignLen + 6, 'St. Dev.', title_format)
    for seq in aaRatioStdevList:
        worksheet2.write(stdevRow, stdevCol, seq, stats_format)
        stdevRow += 1
    logging.info('List of unique amino acid absorbance standard deviations written to %s worksheet.' % worksheet2Name)

    # Change column width to fit all IDs.
    wellColWidth = round((len(aaCountID[0]) * 3) * 1.4)
    worksheet2.set_column(aaAlignLen + 7, aaAlignLen + 7, wellColWidth)
    # Write IDs to worksheet.
    worksheet2.merge_range(0, aaAlignLen + 7, 1, aaAlignLen + 7, 'Wells', wellTitle_format)
    wellRow = 2
    wellCol = aaAlignLen + 7
    sep = ', '
    for wellList in aaCountID:
        wellList = wellRegex.findall(str(wellList))
        wellList = sep.join(wellList)
        worksheet2.write(wellRow, wellCol, wellList, wellList_format)
        wellRow += 1
    logging.info('Amino acid sequence-specific well IDs written to %s worksheet.' % worksheet2Name)

    # Assign arbitrary IDs to each unique amino acid sequence.
    worksheet2.merge_range(0, 0, 1, 0, 'ID', title_format)
    aaIdList = list(range(1,
                          len(aaUniqueDict) + 1)
                    )
    idRow = 2
    for number in aaIdList:
        worksheet2.write_number(idRow, 0, number, general_format)
        idRow += 1
    logging.info('Arbitrary unique amino acid sequence IDs written to %s worksheet.' % worksheet2Name)

    ##################
    # Create worksheet for all nucleotide sequences.
    ##################

    # Initialize sheet characteristics.
    worksheet3Name = 'All NT Seq'
    worksheet3 = workbook.add_worksheet(worksheet3Name)
    worksheet3.hide_gridlines(option=2)
    worksheet3.set_column(0, 0, idColWidth)
    worksheet3.set_column(1, ntAlignLen, 3)
    worksheet3.set_column(ntAlignLen + 1, ntAlignLen + 4, 12)
    worksheet3.set_column(ntAlignLen + 5, ntAlignLen + 5, 15)
    worksheet3.freeze_panes(2, 1)
    logging.info('%s worksheet created.' % worksheet3Name)

    # Write IDs.
    worksheet3.merge_range(0, 0, 1, 0, 'ID', title_format)
    idRow = 2
    for name in ntNameList:
        worksheet3.write(idRow, 0, name, wellID_format)
        idRow += 1
    logging.info('Nucleotide well IDs written to %s worksheet.' % worksheet3Name)

    # Write nucleotide sequences.
    worksheet3.merge_range(0, 1, 0, 7, 'Nucleotide Sequence', title_format)
    seqRow = 2
    seqCol = 1
    for nt in ntSeqList:
        letterList = list(nt)
        for letter in letterList:
            worksheet3.write(seqRow, seqCol, letter, sequence_format)
            seqCol += 1
        seqRow += 1
        seqCol = 1
    logging.info('All nucleotide sequences written to %s worksheet.' % worksheet3Name)

    # Write raw binder absorbances.
    absRow = 2
    absCol = ntAlignLen + 1
    worksheet3.merge_range(0, absCol, 0, absCol + 1, 'Raw', title_format)
    worksheet3.merge_range(0, absCol + 2, 0, absCol + 3, 'Normalised', title_format)
    worksheet3.write(1, absCol, 'Binder', title_format)
    for absorbance in ntRawListShort:
        worksheet3.write(absRow, absCol, absorbance, abs_format)
        absRow += 1
    logging.info('All raw binder absorbances written to %s worksheet.' % worksheet3Name)

    # Write raw control absorbances.
    absRow = 2
    absCol = ntAlignLen + 2
    worksheet3.write(1, absCol, 'Control', title_format)
    for absorbance in ntControlListRaw:
        worksheet3.write(absRow, absCol, absorbance, abs_format)
        absRow += 1
    logging.info('All raw control absorbances written to %s worksheet.' % worksheet3Name)

    # Write averaged binder absorbances.
    absRow = 2
    absCol = ntAlignLen + 3
    worksheet3.write(1, absCol, 'Binder', title_format)
    for absorbance in ntRelAveListShort:
        worksheet3.write(absRow, absCol, absorbance, abs_format)
        absRow += 1
    logging.info('All normalised binder absorbances written to %s worksheet.' % worksheet3Name)

    # Write normalised control absorbances.
    absRow = 2
    absCol = ntAlignLen + 4
    worksheet3.write(1, absCol, 'Control', title_format)
    for absorbance in ntControlListRel:
        worksheet3.write(absRow, absCol, absorbance, abs_format)
        absRow += 1
    logging.info('All normalised control absorbances written to %s worksheet.' % worksheet3Name)

    # Write binder:control absorbance ratios.
    absRow = 2
    absCol = ntAlignLen + 5
    worksheet3.merge_range(0, absCol, 1, absCol, 'Binder:Control', title_format)
    for absorbance in ntBinderControlRatio:
        worksheet3.write(absRow, absCol, absorbance, stats_format)
        absRow += 1
    logging.info('All binder:control absorbance ratios written to %s worksheet.' % worksheet1Name)

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
    worksheet2.set_column(ntAlignLen + 1, ntAlignLen + 5, 8)
    worksheet4.freeze_panes(2, 1)
    logging.info('%s worksheet created.' % worksheet4Name)

    # Write unique nucleotide sequences.
    worksheet4.merge_range(0, 1, 0, 7, 'Nucleotide Sequence', title_format)
    uniqueRow = 2
    uniqueCol = 1
    for seq in ntUniqueDict.keys():
        letterList = list(seq)
        for letter in letterList:
            worksheet4.write(uniqueRow, uniqueCol, letter, sequence_format)
            uniqueCol += 1
        uniqueRow += 1
        uniqueCol = 1
    logging.info('Unique nucleotide sequences written to %s worksheet.' % worksheet4Name)

    # Add counts for each unique sequence.
    countRow = 2
    countCol = ntAlignLen + 1
    worksheet4.merge_range(0, ntAlignLen + 1, 1, ntAlignLen + 1, 'Count', title_format)
    count = list(ntUniqueDict.values())
    for number in count:
        worksheet4.write_number(countRow, countCol, number, general_format)
        countRow += 1
    logging.info('Nucleotide sequence counts written to %s worksheet.' % worksheet4Name)

    # Write nucleotide base pair numbers above sequences.
    bpCol = 1
    for number in ntBpList:
        worksheet4.write(1, bpCol, number, residue_format)
        bpCol += 1
    logging.info('Base pair numbers written to %s worksheet.' % worksheet4Name)

    ##################
    # Write statistics to unique amino acid sequences worksheet.
    ##################

    # Max.
    maxRow = 2
    maxCol = ntAlignLen + 2
    worksheet4.merge_range(0, ntAlignLen + 2, 1, ntAlignLen + 2, 'Max.', title_format)
    for seq in ntRatioMaxList:
        worksheet4.write(maxRow, maxCol, seq, stats_format)
        maxRow += 1
    logging.info('List of unique nucleotide maximum absorbances written to %s worksheet.' % worksheet4Name)

    # Min.
    minRow = 2
    minCol = ntAlignLen + 3
    worksheet4.merge_range(0, ntAlignLen + 3, 1, ntAlignLen + 3, 'Min.', title_format)
    for seq in ntRatioMinList:
        worksheet4.write(minRow, minCol, seq, stats_format)
        minRow += 1
    logging.info('List of unique nucleotide minimum absorbances written to %s worksheet.' % worksheet4Name)

    # Median.
    medianRow = 2
    medianCol = ntAlignLen + 4
    worksheet4.merge_range(0, ntAlignLen + 4, 1, ntAlignLen + 4, 'Median', title_format)
    for seq in ntRatioMedianList:
        worksheet4.write(medianRow, medianCol, seq, stats_format)
        medianRow += 1
    logging.info('List of unique nucleotide median absorbances written to %s worksheet.' % worksheet4Name)

    # Mean.
    meanRow = 2
    meanCol = ntAlignLen + 5
    worksheet4.merge_range(0, ntAlignLen + 5, 1, ntAlignLen + 5, 'Mean', title_format)
    for seq in ntRatioMeanList:
        worksheet4.write(meanRow, meanCol, seq, stats_format)
        meanRow += 1
    logging.info('List of unique nucleotide mean absorbances written to %s worksheet.' % worksheet4Name)

    # Standard deviation.
    stdevRow = 2
    stdevCol = ntAlignLen + 6
    worksheet4.merge_range(0, ntAlignLen + 6, 1, ntAlignLen + 6, 'St. Dev.', title_format)
    for seq in ntRatioStdevList:
        worksheet4.write(stdevRow, stdevCol, seq, stats_format)
        stdevRow += 1
    logging.info('List of unique nucleotide absorbance standard deviations written to %s worksheet.' % worksheet4Name)

    # Change column width to fit all IDs.
    worksheet4.set_column(ntAlignLen + 7, ntAlignLen + 7, wellColWidth)
    # Write IDs to worksheet.
    worksheet4.merge_range(0, ntAlignLen + 7, 1, ntAlignLen + 7, 'Wells', wellTitle_format)
    wellRow = 2
    wellCol = ntAlignLen + 7
    sep = ', '
    for wellList in ntCountID:
        wellList = wellRegex.findall(str(wellList))
        wellList = sep.join(wellList)
        worksheet4.write(wellRow, wellCol, wellList, wellList_format)
        wellRow += 1
    logging.info('Nucleotide sequence-specific well IDs written to %s worksheet.' % worksheet4Name)

    # Assign arbitrary IDs to each unique nucleotide sequence.
    worksheet4.merge_range(0, 0, 1, 0, 'ID', title_format)
    ntIdList = list(range(1,
                          len(ntUniqueDict) + 1)
                    )
    idRow = 2
    for number in ntIdList:
        worksheet4.write_number(idRow, 0, number, general_format)
        idRow += 1
    logging.info('Arbitrary unique nucleotide sequence IDs written to %s worksheet.' % worksheet4Name)

    ##################
    # Final workbook formatting.
    ##################

    # Info about how the values were obtained.
    blankInfo = 'Normalised against the average absorbance of %i %s blanks (%s).' \
                % (len(blankValues), blankID, round(blankAve, 3))
    worksheet1.write(len(aaShortNameList) + 3, aaAlignLen + 1, blankInfo, info_format)
    worksheet3.write(len(ntShortNameList) + 3, 1, blankInfo, info_format)

    calcInfo = 'Manual calculation may produce slightly different results due to rounding.'
    worksheet1.write(len(aaShortNameList) + 5, aaAlignLen + 1, calcInfo, info_format)
    worksheet2.write(len(aaUnique) + 5, aaAlignLen + 2, calcInfo, info_format)
    worksheet3.write(len(ntShortNameList) + 5, ntAlignLen + 1, calcInfo, info_format)
    worksheet4.write(len(ntUnique) + 5, ntAlignLen + 2, calcInfo, info_format)

    statsInfo = 'Statistics presented above are for binder:control binding ratios.'
    worksheet2.write(len(aaUnique) + 3, aaAlignLen + 2, statsInfo, info_format)
    worksheet4.write(len(ntUnique) + 3, ntAlignLen + 2, statsInfo, info_format)

    # Conditionally format columns.
    for column in range(1, 6, 2):
        worksheet1.conditional_format(2, aaAlignLen + column, len(aaShortNameList) + 2, aaAlignLen + column + 1,
                                      {'type': '2_color_scale',
                                       'min_color': '#FFFFFF',
                                       'max_color': '#3D85C6'
                                       }
                                      )

    for column in range(2, 7):
        worksheet2.conditional_format(2, aaAlignLen + column, len(aaUnique) + 2, aaAlignLen + column,
                                      {'type': '2_color_scale',
                                       'min_color': '#FFFFFF',
                                       'max_color': '#3D85C6'
                                       }
                                      )

    for column in range(1, 6, 2):
        worksheet3.conditional_format(2, ntAlignLen + column, len(ntShortNameList) + 2, ntAlignLen + column + 1,
                                      {'type': '2_color_scale',
                                       'min_color': '#FFFFFF',
                                       'max_color': '#3D85C6'
                                       }
                                      )

    for column in range(2, 7):
        worksheet4.conditional_format(2, ntAlignLen + column, len(aaUnique) + 2, ntAlignLen + column,
                                      {'type': '2_color_scale',
                                       'min_color': '#FFFFFF',
                                       'max_color': '#3D85C6'
                                       }
                                      )

    # Transform data into proper Excel-formatted tables without any design style applied.
    worksheet1.add_table(2, 0, len(aaShortNameList) + 1, aaAlignLen + 5,
                         {'header_row': False,
                          'style': None
                          }
                         )
    worksheet2.add_table(2, 0, len(aaUnique) + 1, aaAlignLen + 7,
                         {'header_row': False,
                          'style': None
                          }
                         )
    worksheet3.add_table(2, 0, len(ntNameList) + 1, ntAlignLen + 5,
                         {'header_row': False,
                          'style': None
                          }
                         )
    worksheet4.add_table(2, 0, len(ntUnique) + 1, ntAlignLen + 7,
                         {'header_row': False,
                          'style': None
                          }
                         )

    ##################
    # Conclusion
    ##################

    workbook.close()
    logging.info('Excel file exported as %s_analysed.xlsx.' % elisaInFileNameShort)
    Sg.Popup('''Binding Analysis program finished running.
\nSee log file for details.''',
             title='Binding Analysis Completed',
             grab_anywhere=True,
             text_color='#8294cc'
             )
    logging.info('Binding Analysis program finished running.')
    logging.shutdown()
    window.close()
