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

# Create window layout.
infoLayout = [
    # Title and introduction.
    [Sg.Text('\nPhage Display - Binding Analysis',
             text_color='#8294cc',
             font=('Segoe UI Semibold', 16)
             )
     ],
    [Sg.Text('''        Analyses ELISA binding data along with corresponding sequence data.
        Calculates the average of duplicate absorbances for each protein and
        normalizes them against the average of the negative controls/blanks.
        ELISA data that don't have corresponding sequencing data are excluded
        from the final results.\n\n''',
             text_color='#8294cc',
             font=('Segoe UI', 12)
             )
     ],

    # 'Plate layout' button.
    [Sg.Text('        Click the button below to see the required plate layout for this program.\n',
             text_color='#8294cc',
             font=('Segoe UI', 12)
             )
     ],
    [Sg.Button('Plate Layout',
               font=('Segoe UI Bold', 10),
               size=(13, 0),
               pad=(40, 0),
               use_ttk_buttons=True
               )
     ]
]

inFileLayout = [

    # Input format prompt.
    [Sg.Text('1. Choose input format:',
             text_color='white',
             font=('Segoe UI Bold', 10)
             )
     ],
    [Sg.Radio('Xlsx (one file)',
              'FORMAT',
              default=True,
              key='-XLSX_FORMAT-',
              text_color='#bfbfbf',
              font=('Segoe UI Bold', 10),
              enable_events=True
              ),
     Sg.Radio('Fasta (two files)',
              'FORMAT',
              default=False,
              key='-FASTA_FORMAT-',
              text_color='#bfbfbf',
              font=('Segoe UI Bold', 10),
              enable_events=True
              )
     ],
    [Sg.Text('''    * Requires xlsx output from Sequence Analysis program.\n''',
             text_color='#bfbfbf',
             font=('Segoe UI', 10),
             key='-RADIO_TEXT-'
             )
     ],

    # ELISA file input prompt.
    [Sg.Text('2. Enter the full path of the raw ELISA data file:',
             text_color='white',
             font=('Segoe UI Bold', 10)
             )
     ],
    [Sg.Input(key='-ELISA_INPUT-',
              size=60,
              pad=(25, 0),
              font=('Segoe UI', 10)
              ),
     Sg.FileBrowse(font=('Segoe UI Bold', 10),
                   size=(10, 0),
                   file_types=(('Excel Files', '*.xlsx'), ('All Files', '*.*'),),
                   )
     ],
    [Sg.Text('''    * Must be in xlsx format.
    * This location will also be the location of the output files.
    * Make sure there are no dashes in the name, replace with an underscore if necessary.\n''',
             text_color='#bfbfbf',
             font=('Segoe UI', 10)
             )
     ],

    # Xlsx alignment file input prompt (changes to amino acid alignment input prompt depending on radio buttons).
    [Sg.Text('3. Enter the full path of the excel alignment file:',
             text_color='white',
             font=('Segoe UI Bold', 10),
             key='-XLSX_INPUT_START_TEXT-'
             )
     ],
    [Sg.Input(key='-SWITCH_INPUT-',
              size=60,
              pad=(25, 0),
              font=('Segoe UI', 10),
              do_not_clear=False
              ),
     Sg.FileBrowse(font=('Segoe UI Bold', 10),
                   size=(10, 0),
                   file_types=(('Excel/Fasta Files', '*.xlsx;*.fasta'), ('All Files', '*.*'),)
                   )
     ],
    [Sg.Text('''    * Must be in fasta format.\n''',
             text_color='#bfbfbf',
             font=('Segoe UI', 10),
             key='-XLSX_INPUT_END_TEXT-'
             )
     ],

    # Nucleotide alignment file input prompt.
    [Sg.Text('4. Enter the full path of the nucleotide alignment file:',
             text_color='#464646',
             font=('Segoe UI Bold', 10),
             key='-NT_INPUT_START_TEXT-',
             )
     ],
    [Sg.Input(key='-NT_INPUT-',
              size=60,
              pad=(25, 0),
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
    [Sg.Text('''    * Must be in fasta format.\n''',
             text_color='#464646',
             font=('Segoe UI', 10),
             key='-NT_INPUT_END_TEXT-'
             )
     ]
]

settingsLayout = [

    # Blank well ID input prompt.
    [Sg.Text('1. Enter well IDs that contain ELISA blanks/negative controls:',
             text_color='white',
             font=('Segoe UI Bold', 10)
             )
     ],
    [Sg.Input(key='-BLANK_INPUT-',
              size=20,
              pad=(25, 0),
              font=('Segoe UI', 10)
              )
     ],
    [Sg.Text('''    * Not case-sensitive.
    * Blank wells not associated with any absorbance value will be ignored.
    * For multiple wells, use any kind of separator between them.
      E.g. p22,p23,p24    or    P22-P23-P24\n''',
             text_color='#bfbfbf',
             font=('Segoe UI', 10)
             )
     ],

    # Emission absorbance input prompt.
    [Sg.Text('2. Enter the wavelength of the emission peak (nm):',
             text_color='white',
             font=('Segoe UI Bold', 10)
             )
     ],
    [Sg.Input(key='-EMISSION_INPUT-',
              size=5,
              pad=(25, 0),
              default_text='450',
              font=('Segoe UI', 10)
              )
     ],
    [Sg.Text('''    * In most cases it will be 450 nm and does not need to be changed.\n''',
             text_color='#bfbfbf',
             font=('Segoe UI', 10)
             )
     ],
    [Sg.Button('Enter',
               bind_return_key=True,
               font=('Segoe UI Bold', 16),
               size=(10, 0),
               pad=(30, 25),
               use_ttk_buttons=True
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
                    inFileLayout,
                    border_width=40
                    ),
             Sg.Tab('Settings',
                    settingsLayout,
                    border_width=40
                    )
             ]
        ],
        tab_location='topleft',
        border_width=1,
        font=('Segoe UI Bold', 10),
        title_color='#bfbfbf',
        selected_title_color='#8294cc',
        size=(650, 550)
    )
    ]
]

# Set up window behaviour..
window = Sg.Window('Phage Display - Binding Analysis',
                   tabGroup,
                   alpha_channel=0.95,
                   grab_anywhere=True,
                   ttk_theme='clam'
                   )

# Create a while loop that keeps the window open.
while True:
    event, values = window.read()

    # If 'Exit' or close window buttons are pressed, break the loop and close the window.
    if event == Sg.WIN_CLOSED:
        window.close()
        break

    # If 'Plate Layout' is pressed, a popup will show the plate layout picture.
    elif event == 'Plate Layout':
        # Make sure the path for the image corresponds to the location of the image in the executable's folder.
        Sg.Popup(image='.\\images\\plateLayout.png',
                 title='Plate Layout',
                 grab_anywhere=True
                 )

    elif event == '-XLSX_FORMAT-':
        window['-RADIO_TEXT-'].update('    * Requires xlsx output from Sequence Analysis program.\n')
        window['-XLSX_INPUT_START_TEXT-'].update('2. Enter the full path of the excel alignment file:')
        window['-XLSX_INPUT_END_TEXT-'].update('    * Must be in xlsx format.\n')
        window['-NT_INPUT_START_TEXT-'].update(text_color='#464646')
        window['-NT_BROWSE-'].update(disabled=True)
        window['-NT_INPUT-'].update(disabled=True)
        window['-NT_INPUT_END_TEXT-'].update(text_color='#464646')
        continue

    elif event == '-FASTA_FORMAT-':
        window['-RADIO_TEXT-'].update('    * Requires fasta alignment output from sequence analysis program or'
                                      ' elsewhere.\n')
        window['-XLSX_INPUT_START_TEXT-'].update('2. Enter the full path of the amino acid alignment file:')
        window['-XLSX_INPUT_END_TEXT-'].update('    * Must be in fasta format.\n')
        window['-NT_INPUT_START_TEXT-'].update(text_color='white')
        window['-NT_BROWSE-'].update(disabled=False)
        window['-NT_INPUT-'].update(disabled=False)
        window['-NT_INPUT_END_TEXT-'].update(text_color='#bfbfbf')
        continue

    # If 'Enter' is pressed, updates variables with input values.
    elif event == 'Enter':
        inputOptions = {'1': 'xlsx file',
                        '2': 'fasta files'
                        }
        elisaInFilePath = str(values['-ELISA_INPUT-'])
        blankWells = str(values['-BLANK_INPUT-'])
        emissionAbs = str(values['-EMISSION_INPUT-'])

        # Stops user if no ELISA file is found in the working directory.
        if not os.path.exists(elisaInFilePath):
            Sg.Popup('''The entered raw ELISA xlsx file does not exist in this location.
Please enter it again.''',
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

        if values['-XLSX_FORMAT-']:
            inputFormat = '1'
            xlsxInFilePath = str(values['-SWITCH_INPUT-'])
            if not os.path.exists(xlsxInFilePath):
                Sg.Popup('''The entered xlsx alignment file does not exist in this location.
Please choose the file again.''',
                         title='File Not Found',
                         grab_anywhere=True,
                         text_color='#4276ac')
                continue

        if values['-FASTA_FORMAT-']:
            inputFormat = '2'
            # Stops user if no amino acid alignment file is found in the working directory.
            aaInFilePath = str(values['-SWITCH_INPUT-'])
            if not os.path.exists(aaInFilePath):
                Sg.Popup('''The entered amino acid alignment file does not exist in this location.
Please choose the file again.''',
                         title='File Not Found',
                         grab_anywhere=True,
                         text_color='#4276ac')
                continue
            aaInFileName = re.findall(r'[a-zA-Z0-9]+\.fasta$',
                                      aaInFilePath
                                      )

            # Stops user if no nucleotide alignment file is found in the working directory.
            ntInFilePath = str(values['-NT_INPUT-'])
            if not os.path.exists(ntInFilePath):
                Sg.Popup('''The entered nucleotide alignment file does not exist in this location.
Please choose the file again.''',
                         title='File Not Found',
                         grab_anywhere=True,
                         text_color='#4276ac')
                continue
            ntInFileName = re.findall(r'[a-zA-Z0-9]+\.fasta$',
                                      ntInFilePath
                                      )

        # Stops user if well formatting is incorrect.
        blankWellsInput = re.match(r'[pP]\d{2}[,]*',
                                   blankWells
                                   )
        if blankWellsInput is None:
            Sg.Popup('''Invalid input for blank wells.
Please make sure there are commas (no spaces) between each entry and try again.''',
                     title='Invalid Blank Wells Input',
                     grab_anywhere=True,
                     text_color='#8294cc'
                     )
            continue

        # Stops user if emission absorbance formatting is incorrect.
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

    ##################
    # Retrieve and parse raw ELISA data.
    ##################

    allCells = pandas.read_excel(elisaInFilePath)
    logging.info('Raw data read from ELISA file.')
    # Remove rows where the last column isn't equal to the emission absorbance.
    lastColName = 'Unnamed: ' + str(allCells.shape[1] - 1)
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
    # Rename rows to make parsing easier
    # For dataframes that include empty rows.
    if dataCellsClean.shape[0] == 16:
        dataCellsClean.index = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P']
        # Remove empty rows from the dataframe.
        dataCellsClean = dataCellsClean.drop(['B', 'D', 'F', 'H', 'J', 'L', 'N'], axis=0)
    # For dataframes that don't include empty rows.
    elif dataCellsClean.shape[0] == 9:
        dataCellsClean.index = ['A', 'C', 'E', 'G', 'I', 'K', 'M', 'O', 'P']
    # Create a dictionary of all the blank values and their corresponding well IDs.
    blankCells = {'P1': dataCellsClean.iloc[8][0], 'P2': dataCellsClean.iloc[8][1], 'P3': dataCellsClean.iloc[8][2],
                  'P4': dataCellsClean.iloc[8][3], 'P5': dataCellsClean.iloc[8][4], 'P6': dataCellsClean.iloc[8][5],
                  'P7': dataCellsClean.iloc[8][6], 'P8': dataCellsClean.iloc[8][7], 'P9': dataCellsClean.iloc[8][8],
                  'P10': dataCellsClean.iloc[8][9], 'P11': dataCellsClean.iloc[8][10],
                  'P12': dataCellsClean.iloc[8][11],
                  'P13': dataCellsClean.iloc[8][12], 'P14': dataCellsClean.iloc[8][13],
                  'P15': dataCellsClean.iloc[8][14],
                  'P16': dataCellsClean.iloc[8][15], 'P17': dataCellsClean.iloc[8][16],
                  'P18': dataCellsClean.iloc[8][17],
                  'P19': dataCellsClean.iloc[8][18], 'P20': dataCellsClean.iloc[8][19],
                  'P21': dataCellsClean.iloc[8][20],
                  'P22': dataCellsClean.iloc[8][21], 'P23': dataCellsClean.iloc[8][22],
                  'P24': dataCellsClean.iloc[8][23]
                  }
    logging.info('Blank absorbances retrieved from raw data file.')
    # Create a list of all the data absorbances.
    cellValues = list(dataCellsClean.iloc[0]) + list(dataCellsClean.iloc[1]) + list(dataCellsClean.iloc[2]) + \
                 list(dataCellsClean.iloc[3]) + list(dataCellsClean.iloc[4]) + list(dataCellsClean.iloc[5]) + \
                 list(dataCellsClean.iloc[6]) + list(dataCellsClean.iloc[7])
    logging.info('Data absorbances retrieved from raw data file.')

    ##################
    # Extract blank/negative control values.
    ##################

    # Create a list of ELISA absorbances corresponding to user-inputted blank IDs.
    blankWells = blankWells.upper()
    blankWells = blankWells.split(',')
    blankValues = []
    for well in blankWells:
        if well in blankCells:
            blankValues.append(blankCells.get(well))
    # Remove values that contain no data in the excel cell.
    blankValues = [value for value in blankValues if not (pandas.isna(value))]
    # Average blanks.
    blankAve = statistics.mean(blankValues)
    logging.info('Blank values averaged.')

    ##################
    # Average paired ELISA values for each sequence and normalise to the blank/negative control average.
    ##################

    # Average paired values and append to a new list.
    cellAveList = []
    for value in range(0, len(cellValues), 2):
        cellAveList.append(statistics.mean([cellValues[value], cellValues[value + 1]]))
    logging.info('Paired ELISA absorbances averaged.')

    # Normalise ELISA scores to the blank/negative control average.
    relAveList = []
    for value in cellAveList:
        relAve = value / blankAve
        relAveList.append(relAve)
    logging.info('Averaged absorbances normalised to the blanks/negative control average.')

    ##################
    # Retrieve and parse amino acid/nucleotide sequence data.
    ##################

    if inputFormat == '1':
        aaCells = pandas.read_excel(xlsxInFilePath, sheet_name=0)
        ntCells = pandas.read_excel(xlsxInFilePath, sheet_name=2)

        aaSeqRegex = re.compile(r'[ARNDCEQGHILKMFPSTWYVX]{10,}')
        ntSeqRegex = re.compile('[AGTCURYNWSMKBHDV]{10,}')

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

        wellRegex = re.compile(r'([A-Z][0-1][0-9])')
        # Retrieve amino acid sequence names.
        aaNameCells = aaCells.iloc[1:, 0:1]
        aaNameCells = aaNameCells.to_string(index=False)
        aaNameList = re.sub(r'([_][M][\w]*)',
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
        ntNameList = re.sub(r'([_][M][\w]*)',
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
            aaAllLinesTrim = re.sub(r'([_][M][\w]*)',
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
            ntAllLines = re.sub(r'([_][M][\w]*)',
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

    #########
    # Amino acids
    #########

    # Create a list of all possible ELISA IDs and assign corresponding numbers for indexing.
    elisaPlateIDs = ['A01', 'A02', 'A03', 'A04', 'A05', 'A06', 'A07', 'A08', 'A09', 'A10', 'A11', 'A12', 'B01', 'B02',
                     'B03', 'B04', 'B05', 'B06', 'B07', 'B08', 'B09', 'B10', 'B11', 'B12', 'C01', 'C02', 'C03', 'C04',
                     'C05', 'C06', 'C07', 'C08', 'C09', 'C10', 'C11', 'C12', 'D01', 'D02', 'D03', 'D04', 'D05', 'D06',
                     'D07', 'D08', 'D09', 'D10', 'D11', 'D12', 'E01', 'E02', 'E03', 'E04', 'E05', 'E06', 'E07', 'E08',
                     'E09', 'E10', 'E11', 'E12', 'F01', 'F02', 'F03', 'F04', 'F05', 'F06', 'F07', 'F08', 'F09', 'F10',
                     'F11', 'F12', 'G01', 'G02', 'G03', 'G04', 'G05', 'G06', 'G07', 'G08', 'G09', 'G10', 'G11', 'G12',
                     'H01', 'H02', 'H03', 'H04', 'H05', 'H06', 'H07', 'H08', 'H09', 'H10', 'H11', 'H12'
                     ]
    wellListConversion = {key: value for key, value in zip(elisaPlateIDs, range(0, 96))}

    # Determine the amino acid IDs present in the sequencing data.
    seqPlateRegex = re.compile(r'[A-H][0-9]{2}')
    plateIDs = list()
    for name in aaNameList:
        seqID = seqPlateRegex.findall(name)
        plateIDs.append(seqID)
    # Turn list of lists into a flat list.
    seqPlateIDs = []
    for sublist in plateIDs:
        for ID in sublist:
            seqPlateIDs.append(ID)

    # Retrieve ELISA absorbances for only the IDs present in the sequencing data.
    aaRawListShort = []
    for ID in seqPlateIDs:
        well = wellListConversion.get(ID)
        aaRawListShort.append(cellAveList[well])
    aaRelAveListShort = []
    for ID in seqPlateIDs:
        well = wellListConversion.get(ID)
        aaRelAveListShort.append(relAveList[well])
    logging.info('ELISA results without corresponding amino acid sequences removed from analysis.')

    # Relate IDs to their respective ELISA absorbances.
    aaShortNameList = list()
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

    #########
    # Nucleotides
    #########

    # Determine the nucleotide IDs present in the sequencing data.
    plateIDs = list()
    for name in ntNameList:
        seqID = seqPlateRegex.findall(name)
        plateIDs.append(seqID)
    # Turn list of lists into a flat list.
    seqPlateIDs = []
    for sublist in plateIDs:
        for item in sublist:
            seqPlateIDs.append(item)

    # Retrieve ELISA absorbances for only the IDs present in the sequencing data.
    ntRawListShort = []
    for ID in seqPlateIDs:
        well = wellListConversion.get(ID)
        ntRawListShort.append(cellAveList[well])
    ntRelAveListShort = []
    for ID in seqPlateIDs:
        well = wellListConversion.get(ID)
        ntRelAveListShort.append(relAveList[well])
    logging.info('ELISA results without corresponding nucleotide sequences removed from analysis.')

    # Relate IDs to their respective ELISA absorbances.
    ntNameListShort = list()
    for name in ntNameList:
        if re.findall(r'([A-H])(\d$)', name):
            newName = re.sub(r'([A-H])(\d$)',
                             r'\g<1>0\g<2>',
                             name
                             )
            ntNameListShort.append(newName)
        else:
            ntNameListShort.append(name)
    ntNameListShort.sort()

    # Create dictionaries with the new reorganised names/sequences.
    ntDict = dict(zip(ntNameListShort,
                      ntSeqList)
                  )
    ntAveCellsDict = dict(zip(ntNameListShort,
                              ntRelAveListShort)
                          )

    # Create a list of unique amino acid sequences ordered by frequency.
    uniqueNt = OrderedCounter(ntSeqList)
    uniqueNt = uniqueNt.most_common()
    uniqueNtDict = dict(uniqueNt)
    logging.info('Dictionary of unique nucleotide sequences created.')

    ##################
    # Order sequences and wells so they can be attributed to unique sequences. Necessary for subsequent statistics.
    ##################

    #########
    # Amino acids
    #########

    # Create ordered list of wells that correspond to unique sequences.
    aaOrderedSeq = []
    for key in aaUniqueDict.keys():
        aaOrderedSeq.append(key)
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

    # Get ordered list of wells that correspond to unique sequences; necessary for subsequent statistics.
    ntOrderedSeq = []
    for key in uniqueNtDict.keys():
        ntOrderedSeq.append(key)
    ntOrderedNames = []
    for seq in ntOrderedSeq:
        for key, value in ntDict.items():
            if seq in value:
                ntOrderedNames.append(key)
    logging.info('Ordered index of nucleotide well IDs created.')

    # Get ordered values corresponding to ordered list of wells for unique sequences; necessary for subsequent
    # statistics.
    ntOrderedAbs = []
    for index in ntOrderedNames:
        for ID, score in ntAveCellsDict.items():
            if index == ID:
                ntOrderedAbs.append(score)
    logging.info('Ordered list of nucleotide absorbances created.')

    # Associate specific well IDs with corresponding unique sequence.
    ntCountID = []
    begin = 0
    for uniqueSeq, count in uniqueNtDict.items():
        end = int(count) + begin
        ntCountID.append(ntOrderedNames[begin:end])
        begin += count
    logging.info('List of specific well IDs associated with nucleotide sequences created.')

    ##################
    # Perform statistical analyses for unique sequences.
    ##################

    #########
    # Amino acids
    #########

    # Retrieve max absorbance for ordered values.
    aaUniqueMaxList = []
    begin = 0
    for seq, count in aaUniqueDict.items():
        end = int(count) + begin
        uniqueMax = max(aaOrderedAbs[begin:end])
        aaUniqueMaxList.append(uniqueMax)
        begin += count
    logging.info('List of unique amino acid maximum absorbances created.')

    # Retrieve min absorbance for ordered values.
    aaUniqueMinList = []
    begin = 0
    for seq, count in aaUniqueDict.items():
        end = int(count) + begin
        uniqueMin = min(aaOrderedAbs[begin:end])
        aaUniqueMinList.append(uniqueMin)
        begin += count
    logging.info('List of unique amino acid minimum absorbances created.')

    # Retrieve median absorbance for ordered values.
    aaUniqueMedianList = []
    begin = 0
    for seq, count in aaUniqueDict.items():
        end = int(count) + begin
        uniqueMedian = statistics.median(aaOrderedAbs[begin:end])
        aaUniqueMedianList.append(uniqueMedian)
        begin += count
    logging.info('List of unique amino acid median absorbances created.')

    # Retrieve mean absorbance for ordered values.
    aaUniqueMeanList = []
    begin = 0
    for seq, count in aaUniqueDict.items():
        end = int(count) + begin
        uniqueMean = statistics.mean(aaOrderedAbs[begin:end])
        aaUniqueMeanList.append(uniqueMean)
        begin += count
    logging.info('List of unique amino acid mean absorbances created.')

    # Retrieve standard deviation of absorbances for ordered values.
    aaUniqueStdevList = []
    begin = 0
    for seq, count in aaUniqueDict.items():
        end = int(count) + begin
        try:
            uniqueStdev = statistics.stdev(aaOrderedAbs[begin:end])
            aaUniqueStdevList.append(uniqueStdev)
        # Above statistic won't work if only a single value so append '0.0' value to the list.
        except statistics.StatisticsError:
            aaUniqueStdevList.append(0)
        begin += count
    logging.info('List of unique amino acid absorbance standard deviations created.')

    #########
    # Nucleotides
    #########

    # Retrieve max absorbance for ordered values.
    ntUniqueMaxList = []
    begin = 0
    for seq, count in uniqueNtDict.items():
        end = int(count) + begin
        uniqueMax = max(ntOrderedAbs[begin:end])
        ntUniqueMaxList.append(uniqueMax)
        begin += count
    logging.info('List of unique nucleotide maximum absorbances created.')

    # Retrieve min absorbance for ordered values.
    ntUniqueMinList = []
    begin = 0
    for seq, count in uniqueNtDict.items():
        end = int(count) + begin
        uniqueMin = min(ntOrderedAbs[begin:end])
        ntUniqueMinList.append(uniqueMin)
        begin += count
    logging.info('List of unique nucleotide minimum absorbances created.')

    # Retrieve median absorbance for ordered values.
    ntUniqueMedianList = []
    begin = 0
    for seq, count in uniqueNtDict.items():
        end = int(count) + begin
        uniqueMedian = statistics.median(ntOrderedAbs[begin:end])
        ntUniqueMedianList.append(uniqueMedian)
        begin += count
    logging.info('List of unique nucleotide median absorbances created.')

    # Retrieve mean absorbance for ordered values.
    ntUniqueMeanList = []
    begin = 0
    for seq, count in uniqueNtDict.items():
        end = int(count) + begin
        uniqueMean = statistics.mean(ntOrderedAbs[begin:end])
        ntUniqueMeanList.append(uniqueMean)
        begin += count
    logging.info('List of unique nucleotide mean absorbances created.')

    # Retrieve stdev absorbance for ordered values.
    ntUniqueStdevList = []
    begin = 0
    for seq, count in uniqueNtDict.items():
        end = int(count) + begin
        try:
            uniqueStdev = statistics.stdev(ntOrderedAbs[begin:end])
            ntUniqueStdevList.append(uniqueStdev)
        # Above statistic won't work if only a single value so append '0.0' value to the list.
        except statistics.StatisticsError:
            ntUniqueStdevList.append(0)
        begin += count
    logging.info('List of unique nucleotide absorbance standard deviations created.')

    ##################
    # Export data as a single xlsx file.
    ##################

    # Create workbook.
    workbook = xlsxwriter.Workbook(path + '/' + elisaInFileNameShort + '_analysed.xlsx')
    logging.info('Excel spreadsheet created as "%s.xlsx".' % elisaInFileNameShort)

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
    # Information.
    info_format = workbook.add_format({'font_size': 12})
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
    worksheet1.set_column(1, aaAlignLen, 2)
    worksheet1.set_column(aaAlignLen + 1, aaAlignLen + 2, 12)
    worksheet1.freeze_panes(0, 1)
    logging.info('%s worksheet created.' % worksheet1Name)

    # Write well IDs.
    idRow = 2
    worksheet1.write(0, 0, 'ID', title_format)
    for name in aaShortNameList:
        worksheet1.write(idRow, 0, name, wellID_format)
        idRow += 1
    logging.info('Amino acid well IDs written to %s worksheet.' % worksheet1Name)

    # Write amino acid sequences.
    worksheet1.write(0, 6, 'Amino Acid Sequence', title_format)
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

    # Write raw ELISA absorbances.
    absRow = 2
    absCol = aaAlignLen + 1
    worksheet1.merge_range(0, absCol, 0, absCol + 1, 'Absorbance', title_format)
    worksheet1.write(1, absCol, 'Raw', title_format)
    for absorbance in aaRawListShort:
        worksheet1.write(absRow, absCol, absorbance, stats_format)
        absRow += 1
    logging.info('All raw absorbances written to %s worksheet.' % worksheet1Name)

    # Write averaged ELISA absorbances.
    absRow = 2
    absCol = aaAlignLen + 2
    worksheet1.write(1, absCol, 'Normalised', title_format)
    for absorbance in aaRelAveListShort:
        worksheet1.write(absRow, absCol, absorbance, stats_format)
        absRow += 1
    logging.info('All normalised absorbances written to %s worksheet.' % worksheet1Name)

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
    worksheet2.freeze_panes(0, 1)
    logging.info('%s worksheet created.' % worksheet2Name)

    # Write unique amino acid sequences.
    worksheet2.write(0, 6, 'Amino Acid Sequence', title_format)
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
    worksheet2.write(0, aaAlignLen + 1, 'Count', title_format)
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
    worksheet2.write(0, aaAlignLen + 2, 'Max.', title_format)
    for seq in aaUniqueMaxList:
        worksheet2.write(maxRow, maxCol, seq, stats_format)
        maxRow += 1
    logging.info('List of unique amino acid maximum absorbances written to %s worksheet.' % worksheet2Name)

    # Min.
    minRow = 2
    minCol = aaAlignLen + 3
    worksheet2.write(0, aaAlignLen + 3, 'Min.', title_format)
    for seq in aaUniqueMinList:
        worksheet2.write(minRow, minCol, seq, stats_format)
        minRow += 1
    logging.info('List of unique amino acid minimum absorbances written to %s worksheet.' % worksheet2Name)

    # Median.
    medianRow = 2
    medianCol = aaAlignLen + 4
    worksheet2.write(0, aaAlignLen + 4, 'Median', title_format)
    for seq in aaUniqueMedianList:
        worksheet2.write(medianRow, medianCol, seq, stats_format)
        medianRow += 1
    logging.info('List of unique amino acid median absorbances written to %s worksheet.' % worksheet2Name)

    # Mean.
    meanRow = 2
    meanCol = aaAlignLen + 5
    worksheet2.write(0, aaAlignLen + 5, 'Mean', title_format)
    for seq in aaUniqueMeanList:
        worksheet2.write(meanRow, meanCol, seq, stats_format)
        meanRow += 1
    logging.info('List of unique amino acid mean absorbances written to %s worksheet.' % worksheet2Name)

    # Standard deviation.    
    stdevRow = 2
    stdevCol = aaAlignLen + 6
    worksheet2.write(0, aaAlignLen + 6, 'St. Dev.', title_format)
    for seq in aaUniqueStdevList:
        worksheet2.write(stdevRow, stdevCol, seq, stats_format)
        stdevRow += 1
    logging.info('List of unique amino acid absorbance standard deviations written to %s worksheet.' % worksheet2Name)

    # Change column width to fit all IDs.
    wellColWidth = round((len(aaCountID[0]) * 3) * 1.4)
    worksheet2.set_column(aaAlignLen + 7, aaAlignLen + 7, wellColWidth)
    # Write IDs to worksheet.
    worksheet2.write(0, aaAlignLen + 7, 'Wells', wellTitle_format)
    wellRow = 2
    wellCol = aaAlignLen + 7
    countIDregex = re.compile(r"([A-Z][0-1][0-9])")
    sep = ', '
    for wellList in aaCountID:
        wellList = countIDregex.findall(str(wellList))
        wellList = sep.join(wellList)
        worksheet2.write(wellRow, wellCol, wellList, wellList_format)
        wellRow += 1
    logging.info('Amino acid sequence-specific well IDs written to %s worksheet.' % worksheet2Name)

    # Assign arbitrary IDs to each unique amino acid sequence.
    worksheet2.write(0, 0, 'ID', title_format)
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
    worksheet3.set_column(ntAlignLen + 1, ntAlignLen + 2, 12)
    worksheet3.freeze_panes(0, 1)
    logging.info('%s worksheet created.' % worksheet3Name)

    # Write IDs.
    worksheet3.write(0, 0, 'ID', title_format)
    idRow = 2
    for name in ntNameList:
        worksheet3.write(idRow, 0, name, wellID_format)
        idRow += 1
    logging.info('Nucleotide well IDs written to %s worksheet.' % worksheet3Name)

    # Write nucleotide sequences.
    worksheet3.write(0, 4, 'Nucleotide Sequence', title_format)
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

    # Write raw ELISA absorbances.
    absRow = 2
    absCol = ntAlignLen + 1
    worksheet3.merge_range(0, absCol, 0, absCol + 1, 'Absorbance', title_format)
    worksheet3.write(1, absCol, 'Raw', title_format)
    for absorbance in ntRawListShort:
        worksheet3.write(absRow, absCol, absorbance, stats_format)
        absRow += 1
    logging.info('All raw absorbances written to %s worksheet.' % worksheet1Name)

    # Write averaged ELISA absorbances.
    absRow = 2
    absCol = ntAlignLen + 2
    worksheet3.write(1, absCol, 'Normalised', title_format)
    for result in ntRelAveListShort:
        worksheet3.write(absRow, absCol, result, stats_format)
        absRow += 1
    logging.info('All normalised absorbances written to %s worksheet.' % worksheet3Name)

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
    worksheet4.freeze_panes(0, 1)
    logging.info('%s worksheet created.' % worksheet4Name)

    # Write unique nucleotide sequences.
    worksheet4.write(0, 4, 'Nucleotide Sequence', title_format)
    uniqueRow = 2
    uniqueCol = 1
    for seq in uniqueNtDict.keys():
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
    worksheet4.write(0, ntAlignLen + 1, 'Count', title_format)
    count = list(uniqueNtDict.values())
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
    worksheet4.write(0, ntAlignLen + 2, 'Max.', title_format)
    for seq in ntUniqueMaxList:
        worksheet4.write(maxRow, maxCol, seq, stats_format)
        maxRow += 1
    logging.info('List of unique nucleotide maximum absorbances written to %s worksheet.' % worksheet4Name)

    # Min.
    minRow = 2
    minCol = ntAlignLen + 3
    worksheet4.write(0, ntAlignLen + 3, 'Min.', title_format)
    for seq in ntUniqueMinList:
        worksheet4.write(minRow, minCol, seq, stats_format)
        minRow += 1
    logging.info('List of unique nucleotide minimum absorbances written to %s worksheet.' % worksheet4Name)

    # Median.
    medianRow = 2
    medianCol = ntAlignLen + 4
    worksheet4.write(0, ntAlignLen + 4, 'Median', title_format)
    for seq in ntUniqueMedianList:
        worksheet4.write(medianRow, medianCol, seq, stats_format)
        medianRow += 1
    logging.info('List of unique nucleotide median absorbances written to %s worksheet.' % worksheet4Name)

    # Mean.
    meanRow = 2
    meanCol = ntAlignLen + 5
    worksheet4.write(0, ntAlignLen + 5, 'Mean', title_format)
    for seq in ntUniqueMeanList:
        worksheet4.write(meanRow, meanCol, seq, stats_format)
        meanRow += 1
    logging.info('List of unique nucleotide mean absorbances written to %s worksheet.' % worksheet4Name)

    # Standard deviation.
    stdevRow = 2
    stdevCol = ntAlignLen + 6
    worksheet4.write(0, ntAlignLen + 6, 'St. Dev.', title_format)
    for seq in ntUniqueStdevList:
        worksheet4.write(stdevRow, stdevCol, seq, stats_format)
        stdevRow += 1
    logging.info('List of unique nucleotide absorbance standard deviations written to %s worksheet.' % worksheet4Name)

    # Change column width to fit all IDs.
    worksheet4.set_column(ntAlignLen + 7, ntAlignLen + 7, wellColWidth)
    # Write IDs to worksheet.
    worksheet4.write(0, ntAlignLen + 7, 'Wells', wellTitle_format)
    wellRow = 2
    wellCol = ntAlignLen + 7
    countIDregex = re.compile(r'([A-Z][0-1][0-9])')
    sep = ', '
    for wellList in ntCountID:
        wellList = countIDregex.findall(str(wellList))
        wellList = sep.join(wellList)
        worksheet4.write(wellRow, wellCol, wellList, wellList_format)
        wellRow += 1
    logging.info('Nucleotide sequence-specific well IDs written to %s worksheet.' % worksheet4Name)

    # Assign arbitrary IDs to each unique nucleotide sequence.
    worksheet4.write(0, 0, 'ID', title_format)
    ntIdList = list(range(1,
                          len(uniqueNtDict) + 1)
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
    blankInfo = 'Normalised against the average absorbance of %i blanks (%s).' % (len(blankValues), round(blankAve, 2))
    worksheet1.write(len(aaShortNameList) + 3, 1, blankInfo, info_format)
    worksheet2.write(len(aaUnique) + 3, 1, blankInfo, info_format)
    worksheet3.write(len(ntNameListShort) + 3, 1, blankInfo, info_format)
    worksheet4.write(len(uniqueNt) + 3, 1, blankInfo, info_format)
    calcInfo = 'Values were rounded after calculation, not before. As a result, manual calculation will produce' \
               ' slightly different results.'
    worksheet1.write(len(aaShortNameList) + 5, 1, calcInfo, info_format)
    worksheet2.write(len(aaUnique) + 5, 1, calcInfo, info_format)
    worksheet3.write(len(ntNameListShort) + 5, 1, calcInfo, info_format)
    worksheet4.write(len(uniqueNt) + 5, 1, calcInfo, info_format)

    # Conditionally format columns.
    worksheet1.conditional_format(1, aaAlignLen + 1, len(aaShortNameList) + 1, aaAlignLen + 1,
                                  {'type': '2_color_scale',
                                   'min_color': '#FAFAFA',
                                   'max_color': '#008000'
                                   }
                                  )
    worksheet1.conditional_format(1, aaAlignLen + 2, len(aaShortNameList) + 1, aaAlignLen + 2,
                                  {'type': '2_color_scale',
                                   'min_color': '#FAFAFA',
                                   'max_color': '#008000'
                                   }
                                  )
    worksheet2.conditional_format(1, aaAlignLen + 2, len(aaUnique) + 1, aaAlignLen + 6,
                                  {'type': '2_color_scale',
                                   'min_color': '#FAFAFA',
                                   'max_color': '#008000'
                                   }
                                  )
    worksheet3.conditional_format(1, ntAlignLen + 1, len(ntNameListShort) + 1, ntAlignLen + 1,
                                  {'type': '2_color_scale',
                                   'min_color': '#FAFAFA',
                                   'max_color': '#008000'
                                   }
                                  )
    worksheet3.conditional_format(1, ntAlignLen + 1, len(ntNameListShort) + 1, ntAlignLen + 2,
                                  {'type': '2_color_scale',
                                   'min_color': '#FAFAFA',
                                   'max_color': '#008000'
                                   }
                                  )
    worksheet4.conditional_format(1, ntAlignLen + 2, len(uniqueNt) + 1, ntAlignLen + 6,
                                  {'type': '2_color_scale',
                                   'min_color': '#FAFAFA',
                                   'max_color': '#008000'
                                   }
                                  )

    # Transform data into proper Excel-formatted tables without any design style applied.
    worksheet1.add_table(1, 0, len(aaShortNameList) + 1, aaAlignLen + 2,
                         {'header_row': False,
                          'style': None
                          }
                         )
    worksheet2.add_table(1, 0, len(OrderedCounter(aaSeqList)) + 1, aaAlignLen + 7,
                         {'header_row': False,
                          'style': None
                          }
                         )
    worksheet3.add_table(1, 0, len(ntNameList) + 1, ntAlignLen + 2,
                         {'header_row': False,
                          'style': None
                          }
                         )
    worksheet4.add_table(1, 0, len(OrderedCounter(ntSeqList)) + 1, ntAlignLen + 7,
                         {'header_row': False,
                          'style': None
                          }
                         )

    ##################
    # Conclusion
    ##################

    workbook.close()
    logging.info('Excel file exported as %s_analysed.xlsx.' % elisaInFileNameShort)
    # TODO: Change what the popup says and have earlier popups that address this if the code fails. The code won't even
    #  get to this popup if any of these issues arise.
    Sg.Popup('''Binding Analysis program finished running. See log file for details.
\n\nPost-analysis help:
\nNon-trimmed files are in the 'noTrim' folder and couldn't be trimmed because of one of the following reasons:
\n      a) Statistical error.
Statistics will encounter an error if 'overflow' cells are in the raw ELISA data. This is reflected in the'''
             ''' number of averaged ELISA scores (i.e. the length of cellAveList) being less than the total number of'''
             ''' sequences retrieved from the alignment files. Replace 'OVFLW' wells with '4' to fix.
\n      b) Indexing error.
If you encounter an error involving index values being out of range, this is because the sum of all the'''
             ''' unique sequence counts does not match the total number of ELISA scores. Check the original alignment'''
             ''' files to make sure sequences aren't repeated.''',
             title='Analysis Finished',
             grab_anywhere=True,
             text_color='#8294cc'
             )
    logging.info('Binding Analysis program finished running.')
    logging.shutdown()
    window.close()
