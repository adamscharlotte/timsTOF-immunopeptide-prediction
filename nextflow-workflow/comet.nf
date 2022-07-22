nextflow.enable.dsl=1

params.outputdir = '/Users/adams/Projects/300K/20220613-library-test/comet-output/CFP_aspn_1/'
params.mzML = '/Users/adams/Projects/300K/20220613-library-test/mzML/E1_50fmol_bruker_S1-E2_1_6339.mzML'
params.fasta = '/Users/adams/Projects/300K/20220613-library-test/fasta/AspN/CFP_aspn_1.fasta'

process runComet {
  """
  /Users/adams/Downloads/crux-4.1.Darwin.x86_64/bin/crux comet --decoy_search 1 --search_enzyme_number 1 --decoy_prefix DECOY_ --overwrite T --output-dir ${params.outputdir} ${params.mzML} ${params.fasta}
  """
}
