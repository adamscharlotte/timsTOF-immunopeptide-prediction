/Users/adams/Downloads/crux-4.1.Darwin.x86_64/bin/crux comet

export PATH $PATH:/Users/adams/Downloads/crux-4.1.Darwin.x86_64/bin
echo $PATH
crux comet
which /Users/adams/Downloads/crux-4.1.Darwin.x86_64/bin/crux


cd /Users/adams/Downloads/Kuster
crux comet --decoy_search 1 --search_enzyme_number 0 --decoy_prefix DECOY_ --overwrite T --output-dir /Users/adams/Downloads/Kuster/crux-output 02446a_GD6-TUM_HLA_138_01_01-DDA-1h-R1.mzML TUM_HLA_138.fasta

/Users/adams/Downloads/crux-4.1.Darwin.x86_64/bin/crux comet --decoy_search 2 --search_enzyme_number 0 --decoy_prefix DECOY_ --overwrite T --output-dir /Users/adams/Downloads/Kuster/crux-output-sep-search 02446a_GD6-TUM_HLA_138_01_01-DDA-1h-R1.mzML TUM_HLA_138.fasta

/Users/adams/Downloads/crux-4.1.Darwin.x86_64/bin/crux comet --decoy_search 2 --search_enzyme_number 0 --decoy_prefix DECOY_ --num_enzyme_termini 1 --allowed_missed_cleavage 5  --output-dir /Users/adams/Downloads/Kuster/crux-output-sep-no-enz 02446a_GD6-TUM_HLA_138_01_01-DDA-1h-R1.mzML TUM_HLA_138.fasta


decoy_search 1
search_enzyme_number 0
# num_enzyme_termini 1
# allowed_missed_cleavage 5


/Users/adams/Downloads/crux-4.1.Darwin.x86_64/bin/crux assign-confidence --decoy-prefix DECOY_ /Users/adams/Downloads/Kuster/crux-output-conc-search/comet.target.txt

cd /Users/adams/Downloads/Kuster/crux-output-sep-search/
cd /Users/adams/Downloads/Kuster/crux-output-sep-no-enz/
/Users/adams/Downloads/crux-4.1.Darwin.x86_64/bin/crux assign-confidence --decoy-prefix DECOY_ --estimation-method mix-max comet.target.txt

cd /Users/adams/Downloads/Kuster/crux-output-conc-search/
/Users/adams/Downloads/crux-4.1.Darwin.x86_64/bin/crux assign-confidence --decoy-prefix DECOY_ --overwrite T comet.target.txt


cd /Users/adams/Downloads/Hela/mzML
/Users/adams/Downloads/crux-4.1.Darwin.x86_64/bin/crux comet --decoy_search 1 --search_enzyme_number 1 --decoy_prefix DECOY_ --overwrite T --output-dir /Users/adams/Downloads/Hela/mzML/10ng/crux-output 10ng/Hela_10ng_1_1_97.mzML uniprot-proteome-reviewed_210622.fasta
cd /Users/adams/Downloads/Hela/mzML/10ng/
/Users/adams/Downloads/crux-4.1.Darwin.x86_64/bin/crux assign-confidence --decoy-prefix DECOY_ --overwrite T crux-output/comet.target.txt

/Users/adams/Downloads/crux-4.1.Darwin.x86_64/bin/crux comet --decoy_search 1 --search_enzyme_number 1 --decoy_prefix DECOY_ --overwrite T --output-dir /Users/adams/Downloads/Hela/mzML/100ng/crux-output 100ng/Hela_100ng_1_1_99.mzML uniprot-proteome-reviewed_210622.fasta
cd /Users/adams/Downloads/Hela/mzML/100ng/
/Users/adams/Downloads/crux-4.1.Darwin.x86_64/bin/crux assign-confidence --decoy-prefix DECOY_ --overwrite T crux-output/comet.target.txt



# params.outputdir = '/Users/adams/Projects/300K/20220613-library-test/comet-output/CFP_aspn_1/'
# params.mzML = '/Users/adams/Projects/300K/20220613-library-test/mzML/E1_50fmol_bruker_S1-E2_1_6339.mzML'
# params.fasta = '/Users/adams/Projects/300K/20220613-library-test/fasta/AspN/CFP_aspn_1.fasta'

/Users/adams/Downloads/crux-4.1.Darwin.x86_64/bin/crux comet --decoy_search 1 --search_enzyme_number 0 --decoy_prefix DECOY_ --overwrite T --output-dir '/Users/adams/Projects/300K/20220613-library-test/comet-output/CFP_aspn_1/' '/Users/adams/Projects/300K/20220613-library-test/mzML/E1_50fmol_bruker_S1-E2_1_6339.mzML' '/Users/adams/Projects/300K/20220613-library-test/fasta/AspN/CFP_aspn_1.fasta'
/Users/adams/Downloads/crux-4.1.Darwin.x86_64/bin/crux comet --decoy_search 1 --search_enzyme_number 0 --decoy_prefix DECOY_ --overwrite T --output-dir '/Users/adams/Projects/300K/20220613-library-test/comet-output/CFP_aspn_37/' '/Users/adams/Projects/300K/20220613-library-test/mzML/H1_50fmol_bruker_S1-H2_1_6346.mzML' '/Users/adams/Projects/300K/20220613-library-test/fasta/AspN/CFP_aspn_37.fasta'

##  MGF
/Users/adams/Downloads/crux-4.1.Darwin.x86_64/bin/crux comet --decoy_search 1 --search_enzyme_number 6 --decoy_prefix DECOY_ --overwrite T --output-dir '/Users/adams/Projects/300K/20220613-library-test/comet-output/CFP_aspn_1/AspN/' '/Users/adams/Projects/300K/20220613-library-test/raw/E1_50fmol_bruker_S1-E2_1_6339.d/E1_50fmol_bruker_S1-E2_1_6339_5.2.216.mgf' '/Users/adams/Projects/300K/20220613-library-test/fasta/AspN/CFP_aspn_1.fasta'
/Users/adams/Downloads/crux-4.1.Darwin.x86_64/bin/crux comet --decoy_search 1 --search_enzyme_number 1 --decoy_prefix DECOY_ --overwrite T --output-dir '/Users/adams/Projects/300K/20220613-library-test/comet-output/CFP_aspn_1/Tryp/' '/Users/adams/Projects/300K/20220613-library-test/raw/E1_50fmol_bruker_S1-E2_1_6339.d/E1_50fmol_bruker_S1-E2_1_6339_5.2.216.mgf' '/Users/adams/Projects/300K/20220613-library-test/fasta/AspN/CFP_aspn_1.fasta'

/Users/adams/Downloads/crux-4.1.Darwin.x86_64/bin/crux assign-confidence --decoy-prefix DECOY_ --overwrite T /Users/adams/Projects/300K/20220613-library-test/comet-output/50fmol/CFP_aspn_1/AspN/comet.target.txt --output-dir '/Users/adams/Projects/300K/20220613-library-test/comet-output/CFP_aspn_1/AspN/'
/Users/adams/Downloads/crux-4.1.Darwin.x86_64/bin/crux assign-confidence --decoy-prefix DECOY_ --overwrite T /Users/adams/Projects/300K/20220613-library-test/comet-output/CFP_aspn_1/Tryp/comet.target.txt --output-dir '/Users/adams/Projects/300K/20220613-library-test/comet-output/CFP_aspn_1/Tryp/'

/Users/adams/Downloads/crux-4.1.Darwin.x86_64/bin/crux comet --decoy_search 1 --search_enzyme_number 6 --decoy_prefix DECOY_ --overwrite T --output-dir '/Users/adams/Projects/300K/20220613-library-test/comet-output/CFP_aspn_37/' '/Users/adams/Projects/300K/20220613-library-test/raw/H1_50fmol_bruker_S1-H2_1_6346.d/H1_50fmol_bruker_S1-H2_1_6346_5.2.216.mgf' '/Users/adams/Projects/300K/20220613-library-test/fasta/AspN/CFP_aspn_37.fasta'

/Users/adams/Downloads/crux-4.1.Darwin.x86_64/bin/crux assign-confidence --decoy-prefix DECOY_ --overwrite T /Users/adams/Projects/300K/20220613-library-test/comet-output/CFP_aspn_37/comet.target.txt --output-dir '/Users/adams/Projects/300K/20220613-library-test/comet-output/CFP_aspn_37/'

/Users/adams/Downloads/crux-4.1.Darwin.x86_64/bin/crux comet --decoy_search 1 --search_enzyme_number 6 --decoy_prefix DECOY_ --overwrite T --output-dir '/Users/adams/Projects/300K/20220613-library-test/comet-output/CFP_aspn_1/AspN/mzML/' '/Users/adams/Projects/300K/20220613-library-test/mzML/E1_50fmol_bruker_S1-E2_1_6339.mzML' '/Users/adams/Projects/300K/20220613-library-test/fasta/AspN/CFP_aspn_1.fasta'
/Users/adams/Downloads/crux-4.1.Darwin.x86_64/bin/crux comet --decoy_search 1 --search_enzyme_number 1 --decoy_prefix DECOY_ --overwrite T --output-dir '/Users/adams/Projects/300K/20220613-library-test/comet-output/CFP_aspn_1/Tryp/' '/Users/adams/Projects/300K/20220613-library-test/raw/E1_50fmol_bruker_S1-E2_1_6339.d/E1_50fmol_bruker_S1-E2_1_6339_5.2.216.mgf' '/Users/adams/Projects/300K/20220613-library-test/fasta/AspN/CFP_aspn_1.fasta'


/Users/adams/Downloads/crux-4.1.Darwin.x86_64/bin/crux comet --decoy_search 1 --enzyme asp-n --decoy_prefix DECOY_ --variable_mod01 15.9949 M 0 3 -1 0 0 --overwrite T --output-dir '/Users/adams/Projects/300K/2022-library-run/comet-output/CFP_aspn_1/' '/Users/adams/Projects/300K/2022-library-run/MGF/AspNLysN-E1_S1-E1_1_6590_5.2.216.mgf' '/Users/adams/Projects/300K/2022-library-run/fasta/AspN/CFP_aspn_1.fasta'
/Users/adams/Downloads/crux-4.1.Darwin.x86_64/bin/crux assign-confidence --decoy-prefix DECOY_ --overwrite T /Users/adams/Projects/300K/2022-library-run/comet-output/CFP_aspn_1/comet.target.txt --output-dir '/Users/adams/Projects/300K/2022-library-run/comet-output/CFP_aspn_1/'

/Users/adams/Downloads/crux-4.1.Darwin.x86_64/bin/crux comet
/Users/adams/Downloads/crux-4.1.Darwin.x86_64/bin/crux comet --decoy_search 1 --enzyme asp-n --decoy_prefix DECOY_ --variable_mod01 15.9949 M 0 3 -1 0 0 --overwrite T --output-dir '/Users/adams/Projects/300K/2022-library-run/comet-output/CFP_aspn_1/' '/Users/adams/Projects/300K/2022-library-run/MGF/AspNLysN-E1_S1-E1_1_6590_5.2.216.mgf' '/Users/adams/Projects/300K/2022-library-run/fasta/AspN/CFP_aspn_1.fasta'

/Users/adams/Downloads/crux-4.1.Darwin.x86_64/bin/crux comet --parameter-file /Users/adams/Projects/300K/2022-library-run/param/asp-comet.params.txt --output-dir '/Users/adams/Projects/300K/2022-library-run/comet-output/CFP_aspn_1/' '/Users/adams/Projects/300K/2022-library-run/MGF/AspNLysN-A1_S1-A1_1_6542_5.2.216.mgf' '/Users/adams/Projects/300K/2022-library-run/fasta/proteomeTools.fasta'
# /Users/adams/Downloads/crux-4.1.Darwin.x86_64/bin/crux assign-confidence --decoy-prefix DECOY_ --overwrite T /Users/adams/Projects/300K/2022-library-run/comet-output/CFP_aspn_1/comet.target.txt --output-dir '/Users/adams/Projects/300K/2022-library-run/comet-output/CFP_aspn_1/'
/Users/adams/Downloads/crux-4.1.Darwin.x86_64/bin/crux assign-confidence --parameter-file /Users/adams/Projects/300K/2022-library-run/param/asp-comet.params.txt /Users/adams/Projects/300K/2022-library-run/comet-output/CFP_aspn_1/comet.target.txt --output-dir '/Users/adams/Projects/300K/2022-library-run/comet-output/CFP_aspn_1/'

/Users/adams/Downloads/crux-4.1.Darwin.x86_64/bin/crux comet --parameter-file /Users/adams/Projects/300K/2022-library-run/param/hlaii-comet.params.txt --output-dir '/Users/adams/Projects/300K/2022-library-run/comet-output/CFP_HLA-II_1/' '/Users/adams/Projects/300K/2022-library-run/MGF/HLAI_1_96_p2-A1_S2-A1_1_6926_5.2.216.mgf' '/Users/adams/Projects/300K/2022-library-run/fasta/HLA-II/TUM_HLA2_1.fasta'
/Users/adams/Downloads/crux-4.1.Darwin.x86_64/bin/crux assign-confidence --decoy-prefix DECOY_ --overwrite T /Users/adams/Projects/300K/2022-library-run/comet-output/CFP_aspn_1/comet.target.txt --output-dir '/Users/adams/Projects/300K/2022-library-run/comet-output/CFP_HLA-II_1/'

/Users/adams/Downloads/crux-4.1.Darwin.x86_64/bin/crux comet --parameter-file /Users/adams/Projects/300K/2022-library-run/param/asp-comet.params.txt --output-dir '/Users/adams/Projects/300K/20220613-library-test/comet-output/CFP_aspn_1/' '/Users/adams/Projects/300K/20220613-library-test/raw/E1_50fmol_bruker_S1-E2_1_6339.d/E1_50fmol_bruker_S1-E2_1_6339_5.2.216.mgf' '/Users/adams/Projects/300K/20220613-library-test/fasta/AspN/proteomeTools.fasta'
# /Users/adams/Downloads/crux-4.1.Darwin.x86_64/bin/crux assign-confidence --decoy-prefix DECOY_ --overwrite T /Users/adams/Projects/300K/20220613-library-test/comet-output/CFP_aspn_1/comet.target.txt --output-dir '/Users/adams/Projects/300K/20220613-library-test/comet-output/CFP_aspn_1/'
/Users/adams/Downloads/crux-4.1.Darwin.x86_64/bin/crux assign-confidence --parameter-file /Users/adams/Projects/300K/2022-library-run/param/asp-comet.params.txt /Users/adams/Projects/300K/20220613-library-test/comet-output/CFP_aspn_1/comet.target.txt --output-dir '/Users/adams/Projects/300K/20220613-library-test/comet-output/CFP_aspn_1/'


# fasta_path='/Users/adams/Projects/300K/20220613-library-test/fasta/AspN/proteomeTools.fasta'
# fasta_path='/Users/adams/Projects/300K/2022-library-run/fasta/AspN/CFP_aspn_1.fasta'
# output_path='/Users/adams/Projects/300K/2022-library-run/comet-output/CFP_aspn_1'
# parameter_path='/Users/adams/Projects/300K/2022-library-run/param/comet.params.txt'
# mgf_path='/Users/adams/Projects/300K/2022-library-run/MGF/AspNLysN-E1_S1-E1_1_6590_5.2.216.mgf'

fasta_path='/Users/adams/Projects/300K/20220613-library-test/fasta/AspN/CFP_aspn_1.fasta'
output_path='/Users/adams/Projects/300K/20220613-library-test/comet-output/CFP_aspn_1'
parameter_path='/Users/adams/Projects/300K/20220613-library-test/param/comet.params.txt'
mgf_path='/Users/adams/Projects/300K/20220613-library-test/raw/E1_50fmol_bruker_S1-E2_1_6339.d/E1_50fmol_bruker_S1-E2_1_6339_5.2.216.mgf'

/Users/adams/Downloads/crux-4.1.Darwin.x86_64/bin/crux comet --parameter-file $parameter_path --output-dir $output_path $mgf_path $fasta_path
/Users/adams/Downloads/crux-4.1.Darwin.x86_64/bin/crux assign-confidence --decoy-prefix DECOY_ --overwrite T  $output_path/comet.target.txt --output-dir $output_path

echo /Users/adams/Downloads/crux-4.1.Darwin.x86_64/bin/crux comet --parameter-file $parameter_path --output-dir $output_path $mgf_path --database_name $fasta_path


fasta_path='/Users/adams/Projects/300K/2022-library-run/fasta/AspN/CFP_aspn_1.fasta'
output_path='/Users/adams/Projects/300K/20220613-library-test/comet-output/CFP_aspn_1'
parameter_path='/Users/adams/Projects/300K/2022-library-run/param/comet.params.txt'
mzml_path='/Users/adams/Projects/300K/20220613-library-test/mzML-test/E1_50fmol_bruker_S1-E2_1_6339-13-07.mzML'

/Users/adams/Downloads/crux-4.1.Darwin.x86_64/bin/crux comet --parameter-file $parameter_path --output-dir $output_path $mzml_path $fasta_path
/Users/adams/Downloads/crux-4.1.Darwin.x86_64/bin/crux assign-confidence --parameter-file $parameter_path $output_path/comet.target.txt --output-dir $output_path


fasta_path='/Users/adams/Projects/300K/2022-library-run/fasta/AspN/CFP_aspn_1.fasta'
output_path='/Users/adams/Projects/300K/20220613-library-test/comet-output/CFP_aspn_1'
parameter_path='/Users/adams/Projects/300K/2022-library-run/param/comet.params.txt'
mzml_path='/Users/adams/Projects/300K/PXD021013/mzml/03210a_BE1-TUM_aspn_1_01_01-DDA-1h-R1.mzML'

/Users/adams/Downloads/crux-4.1.Darwin.x86_64/bin/crux comet --parameter-file $parameter_path --output-dir $output_path $mzml_path $fasta_path
/Users/adams/Downloads/crux-4.1.Darwin.x86_64/bin/crux assign-confidence --parameter-file $parameter_path $output_path/comet.target.txt --output-dir $output_path

fasta_path='/Users/adams/Projects/300K/2022-library-run/fasta/AspN/CFP_aspn_1.fasta'
output_path='/Users/adams/Projects/300K/2022-library-run/comet-output/CFP_aspn_1'
parameter_path='/Users/adams/Projects/300K/2022-library-run/param/asp-comet.params.txt'
mgf_path='/Users/adams/Projects/300K/2022-library-run/MGF/AspNLysN-E1_S1-E1_1_6590_5.2.216.mgf'

/Users/adams/Downloads/crux-4.1.Darwin.x86_64/bin/crux comet --parameter-file $parameter_path --output-dir $output_path $mgf_path $fasta_path
/Users/adams/Downloads/crux-4.1.Darwin.x86_64/bin/crux assign-confidence --decoy-prefix DECOY_ --overwrite T  $output_path/comet.target.txt --output-dir $output_path

fasta_path='/Users/adams/Projects/300K/20220613-library-test/fasta/AspN/proteomeTools.fasta'
# fasta_path='/Users/adams/Projects/300K/2022-library-run/fasta/AspN/CFP_aspn_1.fasta'
output_path='/Users/adams/Projects/300K/2022-library-run/comet-output/CFP_aspn_1'
parameter_path='/Users/adams/Projects/300K/2022-library-run/param/asp-comet.params.txt'
# mgf_path='/Users/adams/Projects/300K/2022-library-run/MGF/AspNLysN-E1_S1-E1_1_6590_5.2.216.mgf'
# mgf_path='/Users/adams/Projects/300K/20220613-library-test/raw/E1_50fmol_bruker_S1-E2_1_6339.d/E1_50fmol_bruker_S1-E2_1_6339_5.2.216.mgf'
# mgf_path='/Users/adams/Projects/300K/20220613-library-test/raw/E1_5fmol_bruker_S1-E4_1_6321.d/E1_5fmol_bruker_S1-E4_1_6321_5.2.216.mgf'
# mgf_path='/Users/adams/Projects/300K/20220613-library-test/raw/E1_5fmol_reg_S1-E3_1_6320.d/E1_5fmol_reg_S1-E3_1_6320_5.2.216.mgf'
mgf_path='/Users/adams/Projects/300K/20220613-library-test/raw/E1_50fmol_reg_S1-E1_1_6338.d/E1_50fmol_reg_S1-E1_1_6338_5.2.216.mgf'

/Users/adams/Downloads/crux-4.1.Darwin.x86_64/bin/crux comet --parameter-file $parameter_path --output-dir $output_path $mgf_path $fasta_path
/Users/adams/Downloads/crux-4.1.Darwin.x86_64/bin/crux assign-confidence --decoy-prefix DECOY_ --overwrite T  $output_path/comet.target.txt --output-dir $output_path

mzml_path='/Users/adams/Projects/300K/20220613-library-test/mzML-test/E1_50fmol_bruker_S1-E2_1_6339-13-07.mzML'

/Users/adams/Downloads/crux-4.1.Darwin.x86_64/bin/crux comet --parameter-file $parameter_path --output-dir $output_path $mzml_path $fasta_path
