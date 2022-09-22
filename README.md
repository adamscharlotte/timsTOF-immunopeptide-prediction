# timsTOF-immunopeptide



## Identify Spectra
in the **Fasta folder** there is a script to create a fasta file (create_fasta.R) and some code to run this script for each pool to generate all the needed fasta files.

In the **MaxQuant-Queue** folder you can find all the code needed to run MaxQuant on each pool.

We first attempted to identify the spectra with COMET, but we could find a conversion tool that generated the right input. You can find some of the attempts in the **Conversion-Tools** folder.

## Annotate Identified Spectra
To eventually train the model, we need to map the identifications to the spectra, perform some filtering, and annotate the spectra.

| Script | Description |
| --- | --- |
| create_annotation_input.R | Retrieve the raw scans *m/Z* and intensities for the precursors in different frames. This is a time-consuming step, because we need to retrieve this information from the raw Bruker .d folder. The output is written to a folder **precursor-mapped/** and used in further steps. |
| filter_annotation_input.R | The identifications are filter for full-length and N-terminal truncated peptides. Also an info.txt file is generated, containing information on how many full-length and truncated peptides were identified. |
| filter_and_add_ce.R | The identifications are filter for full-length and N-terminal truncated peptides and the collision energies are retrieved from the pasefMsmsScans.txt file (no need). |
| filter_and_sum.R | Filter the identifications, prepare the data to be summed on precursor-level. |
|| The output of the filtering scripts in stored in a **un-annotated/** folder. |
| prosit_annotate.py | This script prepares the data to enable training in the prosit framework. The modified sequence is adjusted, spectra are annotated, onehot the charges, and filter based on the adromeda score. The result is written to an hdf5 file. |
| prosit_annotate_fixed.py | This does the same thing as prosit_annotate.py but creates a fixed ce value. |
| master_spectrum_precursor.py | Creates consensus spectra for each precursor. |
| prosit_sum_annotate.py | Does the same as prosit_annotate |


| Script | Description |
| --- | --- |
| master_spectrum_precursor.py | Generate consensus spectra and annotate the spectra. |
| sum_annotate.py | Seperate annotation script. |
| sum_calibrate.py | Seperate annotation script. |

