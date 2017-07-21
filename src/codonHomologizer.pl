#!/usr/bin/perl

#USAGE: Run with no options to get usage or with --help for basic details

use CommandLineInterface;
use warnings;
use strict;

##
## Describe the script
##

our $VERSION = '1.7';
setScriptInfo(CREATED => '6/27/2017',
              VERSION => $VERSION,
              AUTHOR  => 'Robert William Leach',
              CONTACT => 'rleach@princeton.edu',
              COMPANY => 'Princeton University',
              LICENSE => 'Copyright 2017',
              HELP    => << 'END_HELP'

This script, given a set of amino acid sequences and a codon usage table, constructs the underlying genetic sequences to optimize them for crossover events.  It takes 2 or more protein sequences, aligns them to optimize codon homology (using a custom amino acid weight matrix as input to a multiple sequence alignment tool (muscle)), and then recodes the amino acid sequences to be more homologous at the DNA level.  It utilizes the supplied codon usage table to prefer codons that are more common in a particular organism.

Note that the alignment step is a multiple sequence alignment, but the sequence construction step is pair-wise.  Each sequence is produced multiple times, optimized for crossover for every other sequence.

To omit codons from incorporation in the resulting sequence, omit them from the codon usage file (or comment them out - refer to the format for the -c file below).

This script produces "Crossover Metrics" intended to be used to evaluate how good each pair of sequences is in terms of their potential to be involved in crossover events.  Here is a listing of the columns in the Crossover Metrics table:

Source             - The file the metrics are based on.
Pair               - Sequential/arbitrary unique pair ID.
Aln Len            - The length of the alignment produced.
Seqs               - Sequence IDs of the pair of seqs optimized for crossover.
%ID                - Percent identity (identicals / all non-gap aln positions).
%Aligned NTs       - All non-gap alignment positions / alignment length.
Contig Ident >= N  - Number of NTs in identical stretches at least N in length.
Rare Cdns(Usg<=N)  - Number of codons incorporated whose usage is <= N.
Rarest Codon       - Codon:usage.  Rarest codon incorporated and its usage.
Frame Shift bases  - Number of NTs aligned out of frame.
Seg Aln Score[0-1] - Segment alignment score, a value between 0 and 1.

The Seg Aln Score is only present is the -d option is used to align for example, similar protein structural domains with one another.

To evaluate a sequence (nt alignment) produced by another tool, you can use the -e option.  All this does is add entries to the Crossover Metrics table at the end of the run.

END_HELP
);

setDefaults(HEADER        => 0,
	    ERRLIMIT      => 3,
	    COLLISIONMODE => 'error', #,merge,rename (when outfile conflict)
	    DEFRUNMODE    => 'usage');

##
## Configure the command line interface
##

my $weighting_methods    = ['crossover1'];
my $weighting_method_def = $weighting_methods->[0];
my $weighting_method     = $weighting_method_def;
addOption(GETOPTKEY   => 'm|method=s',
	  GETOPTVAL   => \$weighting_method,
	  REQUIRED    => 0,
	  DEFAULT     => $weighting_method,
	  HIDDEN      => (scalar(@$weighting_methods) < 2),
	  SMRY_DESC   => 'Weighting method.',
	  DETAIL_DESC => ('The method used to construct the weight matrix ' .
			  'that is used to align the amino acid sequences.'),
	  ACCEPTS     => $weighting_methods);

my $seq_file_type =
  addInfileOption(GETOPTKEY   => 'i|aa-file=s',
		  REQUIRED    => 0,
		  PRIMARY     => 1,
		  DEFAULT     => undef,
		  SMRY_DESC   => 'Amino acid sequence file.',
		  FORMAT_DESC => << 'END_FORMAT'

The amino acid sequence file can be in fasta or fastq format.  There must be 2 or more sequences.  This is an unaligned sequence file.  Any alignment characters present such as gap characters ('-') or spaces will be removed.  It may be upper or lower case.  Codons for any ambiguous nucleotides will be generated as NNN.

END_FORMAT
		 );

my $sequence_formats = ['auto','fasta','fastq'];
my $seq_format = $sequence_formats->[0];
addOption(GETOPTKEY   => 'f|aa-seq-format=s',
	  GETOPTVAL   => \$seq_format,
	  REQUIRED    => 0,
	  DEFAULT     => $seq_format,
	  HIDDEN      => 0,
	  DETAIL_DESC => 'Amino acid sequence file format.  Applies to files submitted using -i and',
	  ACCEPTS     => $sequence_formats);

my $codon_file_type =
  addInfileOption(GETOPTKEY   => 'c|codon-usage-file=s',
		  REQUIRED    => 1,
		  PRIMARY     => 0,
		  DEFAULT     => undef,
		  SMRY_DESC   => 'Codon usage file.',
		  DETAIL_DESC => 'Tab-delimited file [AA,codon,score].',
		  PAIR_WITH   => $seq_file_type,
		  PAIR_RELAT  => '1:1orM',
		  FORMAT_DESC => << 'END_FORMAT'

A tab- or space- delimited file of codon usage scores.  The first column is the single-character amino acid code, the second column is the codon (3 nucleotides), and the third column us the codon's usage score.  The score can be any integer or fraction.  The scores will be normalized, so the score column does not need to sum to any value.  Additional columns are optional and ignored.  Not all codons are required to be present.  All amino acids and the stop character ('*') are required to be present.  To exclude a codon from incorporation into all recoded sequences, you can either remove that row from the file or comment it out by inserting a '#' at the beginning of the line.  Empty lines are permissable.  Example:

G	GGG	0.11	Gly
G	GGA	0.21	Gly
G	GGT	0.49	Gly
G	GGC	0.19	Gly
E	GAG	0.29	Glu
E	GAA	0.71	Glu
D	GAT	0.65	Asp
D	GAC	0.35	Asp
...

END_FORMAT
		 );

my $extension = '.fna';
my $seq_out_type =
  addOutfileSuffixOption(GETOPTKEY     => 'o|dna-suffix|dna-extension=s',
			 FILETYPEID    => $seq_file_type,
			 GETOPTVAL     => \$extension,
			 REQUIRED      => 0,
			 PRIMARY       => 1,
			 DEFAULT       => $extension,
			 HIDDEN        => 0,
			 DETAIL_DESC   => ('The extension of the output DNA ' .
					   'sequence file that is ' .
					   'constructed from and appended ' .
					   'to the input amino acid ' .
					   'sequence file (see ' .
					   '-i).'),
			 FORMAT_DESC   => 'Fasta sequence file.',
			 COLLISIONMODE => 'error');

my $align_codons_def = 1;
my $align_codons     = $align_codons_def;
addOption(GETOPTKEY   => 'a|align-codons!',
	  GETOPTVAL   => \$align_codons,
	  REQUIRED    => 0,
	  DEFAULT     => $align_codons,
	  HIDDEN      => 0,
	  DETAIL_DESC => ('Align the codons in the output DNA sequence file ' .
			  'to match the amino acid alignment (see -o).  ' .
			  'This is the default behavior.  Use --no-align-' .
			  'codons to generate sequences without gap ' .
			  'characters.'));

my($matrix_suffix);
my $mat_out_type =
  addOutfileSuffixOption(GETOPTKEY     => 's|matrix-suffix|matrix-extension=s',
			 FILETYPEID    => $codon_file_type,
			 GETOPTVAL     => \$matrix_suffix,
			 REQUIRED      => 0,
			 PRIMARY       => 0,
			 DEFAULT       => $matrix_suffix,
			 HIDDEN        => 0,
			 DETAIL_DESC   => << 'END_DETAIL'

The extension of the output protein weight matrix file that is appended to and constructed from the codon usage file (see -c).  The generation of this file is not required except for re-use and reproducability.  It can be re-used by supplying it using the -w option.

END_DETAIL
			 ,
			 FORMAT_DESC   => << 'END_FORMAT'

The format of the weight matrix file is the same as is used by sequence alignment tools such as clustalw and muscle.  There is a row and column for every single-letter aminco acid code and each cell contains a score representing how good of a match the pair of amino acids is (i.e. how similar they are).  The last row and column (represented by '*') contains the minimum score.  Examples of these matrices can be found (at the time of this writing) here:

ftp://ftp.ncbi.nih.gov/blast/matrices/

Example (BLOSUM62):

   A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *
A  4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0 -2 -1  0 -4 
R -1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3 -1  0 -1 -4 
N -2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3  3  0 -1 -4 
D -2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3  4  1 -1 -4 
C  0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1 -3 -3 -2 -4 
Q -1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2  0  3 -1 -4 
E -1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4 
G  0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3 -1 -2 -1 -4 
H -2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3  0  0 -1 -4 
I -1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3 -3 -3 -1 -4 
L -1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1 -4 -3 -1 -4 
K -1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2  0  1 -1 -4 
M -1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1 -3 -1 -1 -4 
F -2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1 -3 -3 -1 -4 
P -1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2 -2 -1 -2 -4 
S  1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2  0  0  0 -4 
T  0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0 -1 -1  0 -4 
W -3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3 -4 -3 -2 -4 
Y -2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1 -3 -2 -1 -4 
V  0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4 -3 -2 -1 -4 
B -2 -1  3  4 -3  0  1 -1  0 -3 -4  0 -3 -3 -2  0 -1 -4 -3 -3  4  1 -1 -4 
Z -1  0  0  1 -3  3  4 -2  0 -3 -3  1 -1 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4 
X  0 -1 -1 -1 -2 -1 -1 -1 -1 -1 -1 -1 -1 -1 -2  0  0 -2 -1 -1 -1 -1 -1 -4 
* -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  1 

END_FORMAT
			 ,
			 COLLISIONMODE => 'error');

my($aa_suffix);
my $aa_out_type =
  addOutfileSuffixOption(GETOPTKEY     => 'aa-aln-suffix=s',
			 FILETYPEID    => $seq_file_type,
			 GETOPTVAL     => \$aa_suffix,
			 REQUIRED      => 0,
			 PRIMARY       => 0,
			 DEFAULT       => $aa_suffix,
			 HIDDEN        => 0,
			 DETAIL_DESC   => << 'END_DETAIL'

The aligned amino acid sequence generated by the alignment tool (-x/-y) can be optionally saved to a file.  While this file is not necessary, it can be useful for debugging.  The end goal is a DNA sequence.  Note, the default outfile suffix for the DNA output is '.fna'.  If you want to output the amino acid alignment, it is recommended to use '.faa' to differentiate it from the DNA output file.

END_DETAIL
			 ,
			 FORMAT_DESC   => ('The format output by the tool ' .
					   'selected by -x/-y.'),
			 COLLISIONMODE => 'error');

my $align_tools    = ['muscle','clustalw'];
my $align_tool_def = $align_tools->[0];
my $align_tool     = $align_tool_def;
addOption(GETOPTKEY   => 'x|alignment-tool=s',
	  GETOPTVAL   => \$align_tool,
	  REQUIRED    => 0,
	  DEFAULT     => $align_tool,
	  HIDDEN      => 0,
	  ACCEPTS     => $align_tools,
	  DETAIL_DESC => << 'END_DETAIL'

The system command line tool to be used in a system call to align the amino acid sequences (see -i).  The path to the executable can be provided by -y and the command line options can be tweaked using -z.

END_DETAIL
	 );

my $align_exe_def = $align_tool;
my $align_exe     = $align_exe_def;
addOption(GETOPTKEY   => 'y|alignment-exe-path=s',
	  GETOPTVAL   => \$align_exe,
	  REQUIRED    => 0,
	  DEFAULT     => 'auto', #Changed to match $align_tool
	  HIDDEN      => 0,
	  DETAIL_DESC => << 'END_DETAIL'

The path to the executable selected by -x.  The default is simply the name of the executable selected by -x (e.g. 'muscle' or 'clustalw'), assuming it's in your PATH.

END_DETAIL
	 );

my($align_opts);

my $muscle_def_opts =
  '-gapopen -12 -gapextend -3';
my $muscle_prohib = ['-in','-out','-clw','-fasta','-html','-phys','-pyi',
		     '-msf','-clwstrict','-clwout','-clwstrictout','-msfout',
		     '-htmlout','-physout','-phyiout','-matrix'];

my $clustalw_def_opts = '-pwgapopen -12 -pwgapext -3 -quiet';
my $clustalw_prohib = ['-infile','-pwmatrix','-matrix','-output','-outfile',
		       '-type','-interactive'];
addOption(GETOPTKEY   => 'z|alignment-opts-str=s',
	  GETOPTVAL   => \$align_opts,
	  REQUIRED    => 0,
	  DEFAULT     => 'auto', #Changed to match $align_tool
	  HIDDEN      => 0,
	  DETAIL_DESC => << "END_DETAIL"

Command line options provided to the sequence alignment tool specified by -y.  The default is determined by the selection made with the -x option.  For muscle, the default options are [$muscle_def_opts].  For clustalw, the default options are [$clustalw_def_opts].  Each command has a series of prohibited options that are controlled by the script.

Muscle prohibited options: [@$muscle_prohib]

Clustalw prohibited options: [@$clustalw_prohib]

It is possible that other options other than these that are prohibited could cause problems, but the lists here represent only those specific options added by this script or those that directly conflict with those added by this script.

END_DETAIL
	 );

my $prealigned = 0;
addOption(GETOPTKEY   => 'r|precomputed-aa-alignment=s',
	  GETOPTVAL   => \$prealigned,
	  REQUIRED    => 0,
	  DEFAULT     => $prealigned,
	  HIDDEN      => 0,
	  SMRY_DESC   => ('Treat the input AA sequence file as pre-aligned.'),
	  DETAIL_DESC => << "END_DETAIL"

If you already have an amino acid alignment and want to simply generate nucleotide sequences from them that have as much homology as possible, supply your AA alignment file using -i and supply this option to prevent the script from stripping the alignment characters (e.g. gaps '-') and re-aligning them.

END_DETAIL
	 );

my $matrix_file_type =
  addInfileOption(GETOPTKEY   => 'w|precomputed-weight-matrix=s',
		  REQUIRED    => 0, #Required if -w is not supplied
		  PRIMARY     => 0,
		  DEFAULT     => undef,
		  DETAIL_DESC => << 'END_DETAIL'

If you do not want this script to construct a matrix based on homology and codon usage (using -c) to be used in aligning the amino acid sequences, supply your preferred matrix here.  One reason to do this would be to sacrifice DNA homology to get an alignment that has more biological significance.  The matrix constructed by this script optimizes for homology and ingores biological protein similarity to promote as many crossovers as possible in disparate sequences.  A codon usage table (-c) is still required in order to construct sequences with optimized homology while allowing homologous ties to defer to the more used codons.

END_DETAIL
		  ,
		  PAIR_WITH   => $codon_file_type,
		  PAIR_RELAT  => '1:1',
		  FORMAT_DESC => << 'END_FORMAT'

The format of the weight matrix file is the same as is used by sequence alignment tools such as clustalw and muscle.  There is a row and column for every single-letter aminco acid code and each cell contains a score representing how good of a match the pair of amino acids is (i.e. how similar they are).  The last row and column (represented by '*') contains the minimum score.  Examples of these matrices can be found (at the time of this writing) here:

ftp://ftp.ncbi.nih.gov/blast/matrices/

Example (BLOSUM62):

   A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *
A  4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0 -2 -1  0 -4 
R -1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3 -1  0 -1 -4 
N -2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3  3  0 -1 -4 
D -2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3  4  1 -1 -4 
C  0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1 -3 -3 -2 -4 
Q -1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2  0  3 -1 -4 
E -1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4 
G  0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3 -1 -2 -1 -4 
H -2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3  0  0 -1 -4 
I -1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3 -3 -3 -1 -4 
L -1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1 -4 -3 -1 -4 
K -1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2  0  1 -1 -4 
M -1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1 -3 -1 -1 -4 
F -2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1 -3 -3 -1 -4 
P -1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2 -2 -1 -2 -4 
S  1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2  0  0  0 -4 
T  0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0 -1 -1  0 -4 
W -3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3 -4 -3 -2 -4 
Y -2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1 -3 -2 -1 -4 
V  0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4 -3 -2 -1 -4 
B -2 -1  3  4 -3  0  1 -1  0 -3 -4  0 -3 -3 -2  0 -1 -4 -3 -3  4  1 -1 -4 
Z -1  0  0  1 -3  3  4 -2  0 -3 -3  1 -1 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4 
X  0 -1 -1 -1 -2 -1 -1 -1 -1 -1 -1 -1 -1 -1 -2  0  0 -2 -1 -1 -1 -1 -1 -4 
* -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  1 

END_FORMAT
		 );

my $eval_file_type =
  addInfileOption(GETOPTKEY   => 'e|evaluate-nt-aln-file|nt-file=s',
		  REQUIRED    => 0,
		  PRIMARY     => 0,
		  DEFAULT     => undef,
		  SMRY_DESC   => 'Aligned nucleotide sequence file.',
		  PAIR_RELAT  => '1',
		  FORMAT_DESC => << 'END_FORMAT'

The aligned nucleotide sequence file can be in fasta or fastq format.  There must be 2 or more aligned sequences of the same length and they must encode a protein (i.e. the first character must be in the first 5'->3' frame).  It may be upper or lower case.

END_FORMAT
		 );

my $rpt_stretch_mins = [5,11];
addArrayOption(GETOPTKEY   => 'stretch-reporting-min=s',
	       GETOPTVAL   => $rpt_stretch_mins,
	       REQUIRED    => 0,
	       DEFAULT     => $rpt_stretch_mins,
	       HIDDEN      => 0,
	       INTERPOLATE => 1,
	       DETAIL_DESC => << "END_DETAIL"

After nucleotide sequence construction (see -o) or the submission of a nucleotide alignment (see -e), metrics are reported to gauge the effectiveness of the alignment and codon usage to produce sequences that are likely to be involved in crossover events.  Each possible pair of sequences will report the number of alignment positions involved in a contiguous stretch of identity that is equal to or larger than the number(s) supplied to this option.  For example, if this minimum value is 5 and a pair of sequences each has an aligned segment of ATTGCTA (and no other stretches of identity that are 5 or longer), the number reported for 'Contig Ident >= 5' will be 7.  The cutoff must be an unsigned integer.

END_DETAIL
	 );

my $rpt_usage_maxes = [0.1];
addArrayOption(GETOPTKEY   => 'usage-reporting-max=s',
	       GETOPTVAL   => $rpt_stretch_mins,
	       REQUIRED    => 0,
	       DEFAULT     => $rpt_stretch_mins,
	       HIDDEN      => 0,
	       INTERPOLATE => 1,
	       DETAIL_DESC => << "END_DETAIL"

After nucleotide sequence construction (see -o) or the submission of a nucleotide alignment (see -e), metrics are reported to gauge the effectiveness of the alignment and codon usage to produce sequences that are likely to be involved in crossover events.  Generally, codons more highly used are utilized, but lesser used codons can get used if another codon pair between aligned amino acids has more homology.  Each sequence will have reported the number of rare codons that were included for each pair of sequences.  You can set the max cutoff using this option.  Multiple space-delimited values can be supplied.  Any unsigned number is acceptable.

END_DETAIL
	 );

my $segment_file_type =
  addInfileOption(GETOPTKEY   => 'd|segment-file|domain-file=s',
		  REQUIRED    => 0,
		  PRIMARY     => 0,
		  DEFAULT     => undef,
		  DETAIL_DESC => << 'END_DETAIL'

To align segments of each sequence (e.g. structural or secondary domains) independently and stitch the global alignment back together, you can provide a segment file containing pairs of coordinates that should be aligned sequentially within each sequence.

END_DETAIL
		  ,
		  PAIR_WITH   => $seq_file_type,
		  PAIR_RELAT  => '1:1',
		  FORMAT_DESC => << 'END_FORMAT'

A segment file is a tab-delimited file consisting of a series of amino acid coordinates (starting from the first amino acid as coordinate 1) in each protein.  Each row of the file is for a different sequence in the input sequence file (see -i).  The first column is the sequence ID (which must match all the IDs in the sequence file from the first non-white-space character after the defline character to the next white-space character (or the end of the line)).  E.g. a defline of ">id1 my protein" would match a value in the first column of the segment file of "id1".  See the usage for option -d for details on how the coordinates are treated.

A segment file contains coordinate pairs for each protein that are aligned with one-another independently.  There must be the same number of coordinates on each row.  An even or odd number of coordinates are acceptable, but they are treated in pairs, meaning each odd and even columned coordinate is treated as an inclusive pair such that the AAs corresponding to the pair are aligned with that segment [inclusive].  In the event that there are an odd number of coordinate columns, the coordinate in the last column is paired with the last coordinate in the protein.

Gaps between pairs are also aligned in an exclusive coordinate fashion, but if one protein contains a gap and another has no gap between pairs, e.g. 1-10 and 5-15 (where no characters in the first sequence are aligned with the first 4 characters of the second sequence), gap characters are artificially inserted.  If the first coordinate (1) or last coordinate of a protein are not included, the first and last segments are aligned such that the end coordinate is inclusive and the inner coordinate is exclusive.

Example:

Mj	24	119	160	256	352	482	487	711
Pf	1	100	152	269	360	544	552	769
Sp	1	148	220	342	405	488	500	796
Hs	38	176	230	371	419	514	516	817

Coordinates aligned:
Mj	1-23	24-119	120-159	160-256	257-351	352-482	483-486	487-711	712-end
Pf	*	1-100	101-151	152-269	270-359	360-544	545-551	552-769	770-end
Sp	*	1-148	149-219	220-342	343-404	405-488	489-499	500-796	797-end
Hs	1-37	38-176	177-229	230-371	372-418	419-514	515-515	516-817	818-end

* Sequence excluded from segment alignment and filled in with gap characters.

If a coordinate is larger than the protein provided, an error or warning will be printed.  If the second coordinate in the last pair is larger than the protein, it will be automatically adjusted to the last possible protein coordinate with a warning.  Any other coordinate out of range will generate a fatal error.

Coordinates on each row must be in sequential ascending order.

END_FORMAT
		 );

addOutdirOption(GETOPTKEY   => 'dir|outdir=s',
		REQUIRED    => 0,
		DEFAULT     => undef,
		HIDDEN      => 0,
		DETAIL_DESC => 'All outfiles will be put in this directory.');

##
## Process the command line
##

processCommandLine();

##
## Validate user options (files handled automatically)
##

#If -r is provided, ensure  -w, -x, -y, -z, --aa-aln-suffix, and -m were not
#supplied because they are incompatible.
if($prealigned)
  {
    my @flags     = ('-w','-x','-y','-z','--aa-aln-suffix','-m');
    my @supplied  = (getNumFileGroups($matrix_file_type),
		     $align_tool_def ne $align_tool,
		     $align_exe_def ne $align_exe,
		     defined($align_opts),
		     defined($aa_suffix),
		     $weighting_method_def ne $weighting_method);
    my @incompats = map {$flags[$_]} grep {$supplied[$_]} (0..$#flags);
    if(scalar(@incompats))
      {
	error("[-r] is incompatible with: [",join(',',@incompats),"].",
	      {DETAIL => "When supply a pre-computed alignment, options " .
	       "having to do with aligning the submitted sequences, such as " .
	       "which alignment method to use and what matrix to use, are " .
	       "ambiguous.  Please either eliminate -r or [" .
	       join(',',@incompats) . "] from the script call."});
	quit(7);
      }
  }

if(defined($matrix_suffix))
  {
    my @flags     = ('-w','-r');
    my $num_mat_file_groups = getNumFileGroups($matrix_file_type);
    my @supplied  = ($num_mat_file_groups,$prealigned);
    my @incompats = map {$flags[$_]} grep {$supplied[$_]} (0..$#flags);
    if(scalar(@incompats))
      {
	error("[-s] is incompatible with: [",join(',',@incompats),"].",
	      {DETAIL => "When supplying a matrix outfile extension, " .
	       "supplying a " .
	       (getNumFileGroups($matrix_file_type) ?
		"pre-computed matrix (-w [" .
		getNumFileGroups($matrix_file_type) . " supplied])" : '') .
	       (scalar(@incompats) == 2 ? " or " : '') .
	       ($prealigned ? "pre-computed alignment (-r)" : '') .
	       " is ambiguous.  Please either eliminate -s or [" .
	       join(',',@incompats) . "] from the script call."});
	quit(8);
      }
  }

if(getNumFileGroups($matrix_file_type) &&
   $weighting_method_def ne $weighting_method)
  {
    error("[-w] is incompatible with [-m].",
	  {DETAIL => "When supplying a matrix outfile extension, supplying " .
	   "a pre-computed matrix (-w), selecting a matrix construction " .
	   "method (-m) is ambiguous.  Please either eliminate [-w] or [-m] " .
	   "from the script call."});
    quit(9);
  }

#If align-codons is not the default values (implying it was explicitly
#supplied) and it is true and sequences were not provided
if($align_codons_def != $align_codons &&
   $align_codons && getNumFileGroups($seq_file_type) == 0)
  {
    error("The option to align codons (-a) requires sequences (-i).");
    quit(10);
  }

if($weighting_method_def ne $weighting_method &&
   getNumFileGroups($seq_file_type) == 0)
  {
    error("The option to select a weighting method (-m) requires sequences ",
	  "(-i).");
    quit(11);
  }

#If a pre-computed alignment hasn't been supplied
if(!$prealigned)
  {
    #Set the default executable path
    $align_exe = getExe($align_exe);

    #Make sure that the executable path provided matches the tool selected
    if(($align_exe =~ /muscle[^\/]*$/i  && $align_tool eq 'clustalw') ||
       ($align_exe =~ /clustal[^\/]*$/i && $align_tool eq 'muscle'))
      {
	error("The alignment tool (-x): [$align_tool] does not appear to ",
	      "match the executable (-y): [$align_exe].",
	      {DETAIL => ("The selected tool must correspond to the " .
			  "executable so that the correct (automatically " .
			  "controlled) options can be provided.")});
	quit(1);
      }

    #Make sure the alignment tool is properly installed
    if($align_exe eq '')
      {
	if($align_tool eq 'muscle')
	  {
	    error("Muscle executable not found.  If you have not installed ",
		  "muscle, you can find it here: http://www.drive5.com/",
		  "muscle/downloads.htm");
	    quit(2);
	  }
	elsif($align_tool eq 'clustalw')
	  {
	    error("Clustalw executable not found.  If you have not installed ",
		  "clustalw, you can find it here: http://www.clustal.org/",
		  "download/current/");
	    quit(3);
          }
      }
    if(incompatible($align_exe,$align_tool))
      {
        error("If you have not ",
	      "installed muscle, you can find it here: http://www.drive5.com/",
	      "muscle/downloads.htm");
        quit(4);
      }

    #If alignment options were not explicitly provided
    if(!defined($align_opts))
      {
	#Set the default executable parameters
	if($align_tool eq 'clustalw')
	  {$align_opts = $clustalw_def_opts}
	elsif($align_tool eq 'muscle')
	  {$align_opts = $muscle_def_opts}
      }

    #Make sure they did not supply prohibited opts
    my $prohib_opts =
      ($align_tool eq 'clustalw' ? $clustalw_prohib : $muscle_prohib);
    my @violations = checkAlignOpts($align_opts,$prohib_opts);
    if(scalar(@violations))
      {
	error("These prohibited options were provided via -z: [@violations].",
	      {DETAIL => 'codonHomologizer adds either these options to the ' .
	       'command line itself or adds options that conflict with ' .
	       'these options.  please remove these options from the -z ' .
	       'string.'});
	quit(5);
      }
  }

#If -i is not provided, -s is required
if(getNumFileGroups($seq_file_type) == 0 && !defined($matrix_suffix) &&
   getNumFileGroups($eval_file_type) == 0)
  {
    error("If a sequence file (-i) is not provided, then either a matrix ",
	  "output file extension (-s) or an aligned nucleotide file for ",
	  "evaluation (-e) is required.");
    quit(6);
  }

if((getNumFileGroups($seq_file_type) || getNumFileGroups($eval_file_type)) &&
   scalar(@$rpt_stretch_mins) == 0 || scalar(grep {/\D/} @$rpt_stretch_mins))
  {
    error("Invalid value(s) for --stretch-reporting-min: [",
	  join(' ',@$rpt_stretch_mins),"].",
	  {DETAIL => ('There must be at least 1 value and all values must ' .
		      'be unsigned integers.')});
    quit(7);
  }

if((getNumFileGroups($seq_file_type) || getNumFileGroups($eval_file_type)) &&
   scalar(@$rpt_usage_maxes) == 0 || scalar(grep {/[^\d\.]/}
					    @$rpt_usage_maxes))
  {
    error("Invalid value(s) for --usage-reporting-max: [",
	  join(' ',@$rpt_stretch_mins),"].",
	  {DETAIL => ('There must be at least 1 value and all values must ' .
		      'be unsigned.')});
    quit(8);
  }

#There can only be 1 codon usage file if an evaluation DNA alignment file has
#been provided
if(getNumFileGroups($eval_file_type) &&
   getNumFileGroups($codon_file_type) != 1)
  {
    error("If a file is provided to -e in order to generate metrics for an ",
	  "alignment, exactly 1 codon usage file (-c) is required.  [",
	  getNumFileGroups($codon_file_type),"] codon usage files were ",
	  "provided.  Please provide only 1 codon usage file.",
	  {DETAIL => ('Only 1 codon usage file can be used to generate ' .
		      'metrics for a nucleotide alignment, and if more than ' .
		      '1 is supplied, the selection of which one to use for ' .
		      'the calculation of usage frequencies (such as the ' .
		      'number of rare codons used - see --usage-reporting-' .
		      'max) is ambiguous.')});
    quit(9);
  }

#Segment files are incompatible in pre-computed alignment mode.
if($prealigned && getNumFileGroups($segment_file_type))
  {
    error("Segment files (-d) are inpompatible with pre-computed alignemnts ",
	  "(-r).");
    quit(10);
  }

#Segment files require protien files
if(getNumFileGroups($segment_file_type) &&
   getNumFileGroups($seq_file_type) == 0)
  {
    error("Segment files (-d) require sequence files (-i).");
    quit(11);
  }

#There can only be 1 segment file if an evaluation DNA alignment file has been
#provided.  The same sequences and domains are assumed for each
if(getNumFileGroups($eval_file_type) &&
   getNumFileGroups($segment_file_type) > 1)
  {
    error("If a file is provided to -e in order to generate metrics for an ",
	  "alignment, only 1 segment file (-d) can be supplied.  [",
	  getNumFileGroups($segment_file_type),"] segment files were ",
	  "provided.  Please provide only 1 segment file.",
	  {DETAIL => ('The evaluation option (-e) is intended to allow the ' .
		      'comparison of the performance of this tool with ' .
		      'another tool (or this tool with different ' .
		      'parameters) on the same set of sequences.  Thus, in ' .
		      'order to make a valid comparison, and compare apples ' .
		      'to apples, the same segments must be evaluated in ' .
		      'the evaluation alignment as in the alignment ' .
		      'produced on this run.')});
    quit(12);
  }

#In processing:
#Construct matrix file if extension provided
#Write to AA alignment file if they provided an extension
#Write to NT file using gaps if they selected them
#Print identity and score the identical stretches
#Select codons to extend stretches

my $default_tmpdir   = './';
my $tmpdir           = getTmpDir();
my $tmp_suffix       = '.tmp';

my $file_seen_hash   = {};
my $matrix_buffer    = {};
my $matrix_obj       = {};
my $usage_buffer     = {};
my $usage_obj        = {};
my $aa_aln_buffer    = {};
my $flex_scores      = {};
my $nt_pair_sets     = [];
my($eval_file);
my $source_files     = []; #For metrics reporting
my $segments_objects = [];
my($segments_obj);

#nextFileSet iterates sets of command-line-supplied files processed together
while(nextFileCombo())
  {
    my $codonUsageFile = getInfile($codon_file_type);
    my $matrixFile     = getInfile($matrix_file_type);
    my $aaFile         = getInfile($seq_file_type);
    if(!defined($eval_file))
      {$eval_file = getInfile($eval_file_type)}
    my $segmentFile    = getInfile($segment_file_type);
    my $matrixOutFile  = getOutfile($mat_out_type);
    my $alnOutFile     = getOutfile($aa_out_type);
    my $ntOutFile      = getOutfile($seq_out_type);

    my $tmpmatfile     = (defined($matrixFile) ?
			  $matrixFile :
			  (defined($matrixOutFile) ?
			   $matrixOutFile : getTempFile($codonUsageFile) .
			   '.mat'));
    my $tmpalnoutfile  = (defined($alnOutFile) ?
			  $alnOutFile :
			  ($prealigned ?
			   $aaFile : getTempFile($codonUsageFile) . '.aln'));

    #Read in the usage data if not already read in previously
    if(!exists($file_seen_hash->{$codonUsageFile}))
      {$usage_buffer->{$codonUsageFile} = readUsageFile($codonUsageFile)}
    $usage_obj = $usage_buffer->{$codonUsageFile};
    next unless(defined($usage_obj));

    #If a protein file was supplied, we're going to align it, and a segment
    #file was supplied
    my $aaSeqs = [];
    if(defined($aaFile) && !$prealigned && defined($segmentFile))
      {
	openIn(*AA,$aaFile) || next;
	$aaSeqs = [getNextSeqRec(*AA,0,$aaFile)];
	closeIn(*AA);

	#Read in the segment file - returns a hash keyed on sequence ID and
	#whose values are lists of coordinate pairs.  Each pair in each
	#sequence will be aligned with the respective pair from each other
	#sequence.
	$segments_obj = readSegmentFile($segmentFile,$aaSeqs,$aaFile);
      }
    else
      {undef($segments_obj)}

    #Create/set a matrix object - we need one if any aa sequences were provided
    #or if we will be writing a matrix file
    if(defined($aaFile) || defined($matrixOutFile))
      {
	if(!exists($matrix_buffer->{$codonUsageFile}))
	  {
	    verbose("Creating matrix using codon usage.");
	    $matrix_buffer->{$codonUsageFile} = createMatrix($usage_obj);
	  }
	$matrix_obj = $matrix_buffer->{$codonUsageFile};
	if(!defined($matrix_obj))
	  {
	    error("Unable to create matrix object from codon usage file: ",
		  "[$codonUsageFile].  Skipping.");
	    next;
	  }
      }

    #Create the matrix file if there was not a pre-aligned file supplied, the
    #matrix file wasn't supplied, and we haven't already created one using the
    #usage_obj from the current codonUsageFile
    if(!$prealigned &&
       !defined($matrixFile) &&
       !exists($file_seen_hash->{$codonUsageFile}))
      {
	verbose("Writing matrix file.");
	writeMatrixFile($tmpmatfile,$matrix_obj);
      }

    #Set the file_seen_hash value for the codon usage file
    $file_seen_hash->{$codonUsageFile}++;
    $file_seen_hash->{$tmpmatfile}++;

    #Align the amino acid file if it is not already pre-aligned
    if(!$prealigned && defined($aaFile))
      {
	verbose("Aligning sequences in [$aaFile] using matrix [$tmpmatfile].");

	if(defined($segments_obj))
	  {aaSegmentAlign($aaSeqs,$tmpmatfile,$tmpalnoutfile,$segments_obj)}
	else
	  {aaAlign($aaFile,$tmpmatfile,$tmpalnoutfile)}
      }

    #Construct the NT sequences
    if(defined($aaFile) && defined($tmpalnoutfile))
      {
	openIn(*AAALN,$tmpalnoutfile) || next;
	my $alnSeqs = [getNextSeqRec(*AAALN,0,$tmpalnoutfile)];
	closeIn(*AAALN);

	push(@$nt_pair_sets,constructNTPairs($alnSeqs,$matrix_obj,$usage_obj));
	push(@$source_files,$aaFile);
	push(@$segments_objects,$segments_obj);

	writeNTPairs($nt_pair_sets->[-1],$ntOutFile);
      }
  }

if(defined($eval_file))
  {
    openIn(*NTALN,$eval_file) || next;
    my $alnSeqs = [getNextSeqRec(*NTALN,0,$eval_file)];
    closeIn(*NTALN);
    push(@$nt_pair_sets,getAllNTPairs($alnSeqs,$eval_file));
    push(@$source_files,$eval_file);
    push(@$segments_objects,$segments_obj);
  }

if(scalar(@$nt_pair_sets))
  {
    my $pairnum = 0;
    my $widests = [];
    my @metrics = ();

    push(@metrics,['Source','Pair','Aln Len','Seqs','%ID',
		   '%Aligned NTs',
		   (map {'Contig Ident >= ' . $_} @$rpt_stretch_mins),
		   (map {"Rare Cdns(Usg<=$_)"} @$rpt_usage_maxes),
		   'Rarest Codon','Frame Shift Bases']);

    if(defined($segments_obj))
      {push(@{$metrics[0]},'Seg Aln Score[0-1]')}

    foreach(0..$#{$metrics[0]})
      {if((scalar(@$widests) - 1) < $_ ||
	  length($metrics[0]->[$_]) > $widests->[$_])
	 {$widests->[$_] = length($metrics[0]->[$_])}}

    foreach my $set_of_nt_pairs (@$nt_pair_sets)
      {
	my $source_file = shift(@$source_files); #Assumed same size
	if(length($source_file) > $widests->[0])
	  {$widests->[0] = length($source_file)}

	my $segobj = shift(@$segments_objects);

	#Generate metrics
	foreach my $pair (@$set_of_nt_pairs)
	  {
	    #Source file (alignment or original set of AA sequences)
	    push(@metrics,[$source_file]);

	    #Pair number
	    $pairnum++;
	    push(@{$metrics[-1]},$pairnum);

	    if(length($pairnum) > $widests->[1])
	      {$widests->[1] = length($pairnum)}

	    my $rec1 = $pair->[0];
	    my $rec2 = $pair->[1];
	    my $id1  = $rec1->[0];
	    my $id2  = $rec2->[0];
	    my $seq1 = $rec1->[1];
	    my $seq2 = $rec2->[1];
	    my $aln_len = length($seq1);
	    if(length($seq2) != $aln_len)
	      {error("Alignment length of sequence pair: [$id1,$id2] is not ",
		     "consistent.  Metrics will be inaccurate.")}

	    #Alignment length
	    push(@{$metrics[-1]},$aln_len);

	    if(length($aln_len) > $widests->[2])
	      {$widests->[2] = length($aln_len)}

	    #Sequence IDs
	    if(length($id1) > 12)
	      {$id1 =~ s/(.{9}).*/$1../}
	    if(length($id2) > 12)
	      {$id2 =~ s/(.{9}).*/$1../}
	    push(@{$metrics[-1]},"$id1\n$id2");

	    if(length($id1) > $widests->[3])
	      {$widests->[3] = length($id1)}
	    if(length($id2) > $widests->[3])
	      {$widests->[3] = length($id2)}

	    #Get metrics
	    my($ident,$aligned,$stretch_bases,$rare_occurrences_per_seq,
	       $rarest_codons_per_seq,$frame_shifts) =
		 calculateMetrics($seq1,$seq2,$aln_len);

	    #Percent identity
	    $ident = roundInt($ident * 100) / 100;
	    $ident =~ s/(\.\d{2}).*/$1/;
	    push(@{$metrics[-1]},$ident);

	    if(length($ident) > $widests->[4])
	      {$widests->[4] = length($ident)}

	    #Percent aligned NTs
	    $aligned = roundInt($aligned * 100) / 100;
	    $aligned =~ s/(\.\d{2}).*/$1/;
	    push(@{$metrics[-1]},$aligned);

	    if(length($aligned) > $widests->[5])
	      {$widests->[5] = length($aligned)}

	    #Number of bases involved in stretches of identity above cutoffs
	    push(@{$metrics[-1]},@$stretch_bases);

	    my $wi = 5; #"widests index"

	    foreach(0..$#{$stretch_bases})
	      {
		$wi++;
		if(length($stretch_bases->[$_]) > $widests->[$wi])
		  {$widests->[$wi] = length($stretch_bases->[$_])}
	      }

	    #Number of rare codons per sequence below usage cutoffs
	    push(@{$metrics[-1]},
		 map {join("\n",@$_)} @$rare_occurrences_per_seq);

	    foreach(0..$#{$rare_occurrences_per_seq})
	      {
		$wi++;
		if(length($rare_occurrences_per_seq->[$_]->[0]) >
		   $widests->[$wi])
		  {$widests->[$wi] =
		     length($rare_occurrences_per_seq->[$_]->[0])}
		if(length($rare_occurrences_per_seq->[$_]->[1]) >
		   $widests->[$wi])
		  {$widests->[$wi] =
		     length($rare_occurrences_per_seq->[$_]->[1])}
	      }

	    #Rarest Codon per sequence
	    push(@{$metrics[-1]},join("\n",@$rarest_codons_per_seq));

	    $wi++;

	    if(length($rarest_codons_per_seq->[0]) > $widests->[$wi])
	      {$widests->[$wi] = length($rarest_codons_per_seq->[0])}
	    if(length($rarest_codons_per_seq->[1]) > $widests->[$wi])
	      {$widests->[$wi] = length($rarest_codons_per_seq->[1])}

	    #Number of bases involved in a frame shift
	    push(@{$metrics[-1]},$frame_shifts);

	    $wi++;

	    if(length($frame_shifts) > $widests->[$wi])
	      {$widests->[$wi] = length($frame_shifts)}

	    #Segment alignment score
	    if(defined($segobj))
	      {
		my $segscore =
		  calculateSegmentAlignScore($id1,$id2,$seq1,$seq2,$segobj);
		$segscore = roundInt($segscore * 1000) / 1000;

		#Number of bases involved in a frame shift
		push(@{$metrics[-1]},$segscore);

		$wi++;

		if(length($segscore) > $widests->[$wi])
		  {$widests->[$wi] = length($segscore)}
	      }
	  }
      }

    #Print the metrics evaluating the alignments
    print STDERR ("\nCrossover Metrics\n=================\n",
		  array2DToString(\@metrics,$widests,'  '));
  }












#Globals used: $rpt_stretch_mins, $rpt_usage_maxes, $usage_obj (assumed to only
#be 1)
sub calculateMetrics
  {
    my $seq1 = $_[0];
    my $seq2 = $_[1];
    my $len  = $_[2];

    #Return values
    my $ident                    = 0;  #This is a sum until the end
    my $aligned                  = 0;  #This is a sum until the end
    my $stretch_bases            = [map {0} @$rpt_stretch_mins];
    my $rare_occurrences_per_seq = [map {[0,0]} @$rpt_usage_maxes];
    my $rarest_codons_per_seq    = ['',''];
    my $frame_shifts             = 0;

    if(length($seq1) != length($seq2))
      {
	error("Aligned sequences are not the same length.  Unable to produce ",
	      "alignment metrics.");
	return(undef,undef,undef,undef,undef,undef);
      }

    my $framepos1   = 0;
    my $framepos2   = 0;
    my $codon1      = '';
    my $codon2      = '';
    my $stretch_len = 0;
    my $usage       = undef;
    my $rarest_usg1 = undef;
    my $rarest_usg2 = undef;
    for(my $i = 0;$i < length($seq1);$i++)
      {
	my $aa1 = uc(substr($seq1,$i,1));
	my $aa2 = uc(substr($seq2,$i,1));

	if($aa1 ne '-')
	  {
	    $framepos1++;
	    $codon1 .= $aa1;
	  }
	if($aa2 ne '-')
	  {
	    $framepos2++;
	    $codon2 .= $aa2;
	  }

	if($aa1 ne '-' && $aa2 ne '-')
	  {
	    if($framepos1 != $framepos2)
	      {$frame_shifts++}

	    if($aa1 eq $aa2)
	      {
		$ident++;
		$stretch_len++;
	      }
	    else
	      {
		for(my $mini = 0;$mini < scalar(@$rpt_stretch_mins);$mini++)
		  {
		    my $min = $rpt_stretch_mins->[$mini];
		    if($stretch_len >= $min)
		      {$stretch_bases->[$mini] += $stretch_len}
		  }

		$stretch_len = 0;
	      }
	    $aligned++;
	  }
	else
	  {
	    for(my $mini = 0;$mini < scalar(@$rpt_stretch_mins);$mini++)
	      {
		my $min = $rpt_stretch_mins->[$mini];
		if($stretch_len >= $min)
		  {$stretch_bases->[$mini] += $stretch_len}
	      }

	    $stretch_len = 0;
	  }

	if($framepos1 == 3)
	  {
	    $usage = getUsageScore2($codon1,$usage_obj);

	    if(!defined($rarest_usg1) || $usage < $rarest_usg1)
	      {
		$rarest_usg1 = $usage;
		$rarest_codons_per_seq->[0] = "$codon1:$usage";
	      }

	    for(my $maxi = 0;$maxi < scalar(@$rpt_usage_maxes);$maxi++)
	      {
		my $max = $rpt_usage_maxes->[$maxi];
		if($usage <= $max)
		  {$rare_occurrences_per_seq->[$maxi]->[0]++}
	      }

	    $codon1    = '';
	    $framepos1 = 0;
	  }
	if($framepos2 == 3)
	  {
	    $usage = getUsageScore2($codon2,$usage_obj);

	    if(!defined($rarest_usg2) || $usage < $rarest_usg2)
	      {
		$rarest_usg2 = $usage;
		$rarest_codons_per_seq->[1] = "$codon2:$usage";
	      }

	    for(my $maxi = 0;$maxi < scalar(@$rpt_usage_maxes);$maxi++)
	      {
		my $max = $rpt_usage_maxes->[$maxi];
		if($usage <= $max)
		  {$rare_occurrences_per_seq->[$maxi]->[1]++}
	      }

	    $codon2    = '';
	    $framepos2 = 0;
	  }
      }

    if($aligned == 0)
      {error("No non-gap characters were found to be aligned together.  ",
	     "Unable to calculate percent identity.")}
    else
      {
	$ident /= $aligned;
	$ident *= 100;
      }
    if($len == 0)
      {error("Alignment length submitted is 0.  Unable to calculate percent ",
	     "aligned sequence characters.")}
    else
      {
	$aligned /= $len;
	$aligned *= 100;
      }

    return($ident,$aligned,$stretch_bases,$rare_occurrences_per_seq,
	   $rarest_codons_per_seq,$frame_shifts);
  }

#This calculates a score based on the number of bases the sequence with the
#smaller aligned segment has aligned to the corresponding segment of the other
#sequence over the combined length of all the smaller segments in the pair.
#E.g.: NNNNTTTTNNN
#      NNNNNNTTTTT  This will score a 0.5 because 2 T's out of the smaller 4
#align.  In this case, a 'T' is a segment that's supposed to align in each
#sequence.
sub calculateSegmentAlignScore
  {
    my $id1    = parseIDFromDef($_[0]);
    my $id2    = parseIDFromDef($_[1]);
    my $seq1   = $_[2];
    my $seq2   = $_[3];
    my $segobj = $_[4];

    my($lesserid,$lesserseq,$greaterid,$greaterseq);
    my $score_sum = 0;
    my $score_total = 0;

    my $lesser_pos = 1;
    my $greater_pos = 1;
    my $alnindex = 0;
    foreach my $segindex (0..$#{$segobj->{SEGS}})
      {
	if($segobj->{SEGS}->[$segindex]->{$id1}->{SIZE} <
	   $segobj->{SEGS}->[$segindex]->{$id2}->{SIZE})
	  {
	    $lesserid  = $id1;
	    $lesserseq = $seq1;
	    $greaterid = $id2;
	    $greaterseq = $seq2;
	  }
	else
	  {
	    $lesserid  = $id2;
	    $lesserseq = $seq2;
	    $greaterid = $id1;
	    $greaterseq = $seq1;
	  }

	$score_total += $segobj->{SEGS}->[$segindex]->{$lesserid}->{SIZE};

	my $start = $alnindex;
	foreach($start..(length($seq1) - 1))
	  {
	    $alnindex = $_;

	    my $lesser_char  = substr($lesserseq,$alnindex,1);
	    my $greater_char = substr($greaterseq,$alnindex,1);

	    #If there's no gap in this position and we're in both
	    #segments' positions for both sequences, increment the sum
	    if($lesser_char ne '-' && $greater_char ne '-' &&
	       $lesser_pos >=
	       $segobj->{SEGS}->[$segindex]->{$lesserid}->{COORDS}->[0] &&
	       $lesser_pos <=
	       $segobj->{SEGS}->[$segindex]->{$lesserid}->{COORDS}->[1] &&
	       $greater_pos >=
	       $segobj->{SEGS}->[$segindex]->{$greaterid}->{COORDS}->[0] &&
	       $greater_pos <=
	       $segobj->{SEGS}->[$segindex]->{$greaterid}->{COORDS}->[1])
	      {$score_sum++}

	    if($lesser_char ne '-')
	      {$lesser_pos++}
	    if($greater_char ne '-')
	      {$greater_pos++}

	    #If we are past the current segment for either sequence, they will
	    #not possibly align, so move on to the next segment
	    if($lesser_pos >
	       $segobj->{SEGS}->[$segindex]->{$lesserid}->{COORDS}->[1] ||
	       $greater_pos >
	       $segobj->{SEGS}->[$segindex]->{$greaterid}->{COORDS}->[1])
	      {
		#We don't want to redo this alignment index.  We already just
		#looked at it and we also already incremented the real
		#positions in each of the sequences
		$alnindex++;
		last;
	      }
	  }
      }

    if($score_total == 0)
      {
	error("Unexpected error.  The combined size of all the lesser sized ",
	      "segments between sequences [$id1] and [$id2] is 0.");
	return(0);
      }

    return($score_sum / $score_total);
  }

sub getAllNTPairs
  {
    my $alnSeqs = $_[0]; #array of 2 (or 3) element arrays with defline,
                         #sequence, and quality
    my $file    = $_[1];
    my $pairs   = [];

    if(scalar(@$alnSeqs) < 2)
      {
	error("Multiple aligned sequences are required to get all pairs of ",
	      "sequences.");
	return($pairs);
      }

    my $got = 0;

    foreach my $i (0..($#{$alnSeqs} - 1))
      {
	my $rec1 = $alnSeqs->[$i];
	if(!defined($rec1))
	  {
	    error("Undefined nucleotide sequence record encountered.");
	    next;
	  }
	foreach my $j (($i+1)..$#{$alnSeqs})
	  {
	    my $rec2 = $alnSeqs->[$j];
	    if(!defined($rec2))
	      {
		error("Undefined nucleotide sequence record encountered.");
		next;
	      }
	    if(length($rec1->[1]) ne length($rec2->[1]))
	      {
		warning("Alignment sequences [$rec1->[0]] and [$rec2->[0]] ",
			"are not equal lengths.  Skipping.",
			{DETAIL => ('You can ignore this warning if ' .
				    'multiple separate nucleotide sequence ' .
				    'alignments were entered into a single ' .
				    'file.')});
		next;
	      }
	    $got++;
	    my $pair = [$rec1,$rec2];
	    push(@$pairs,$pair);
	  }
      }

    if($got == 0)
      {error("No equal length alignment sequences were found",
	     (defined($file) ? " in file [$file]." : '.'))}

    return($pairs);
  }

#Creates a hash from a tab delimited file like this:
#G	GGG	0.11	Gly
#G	GGA	0.21	Gly
#where:
#hash->{G}->{GGG} = {SCORE => 0.11,DATA => ['Gly']}
#hash->{G}->{GGG} = {SCORE => 0.21,DATA => ['Gly']}
sub readUsageFile
  {
    my $file = $_[0];
    my $usg  = {};

    openIn(*USG,$file) || return(undef);

    my $line    = 0;
    my $got_one = 0;
    my $cdnchk  = {};

    while(getLine(*USG))
      {
	verboseOverMe({FREQUENCY => 10,LEVEL => 2},"Reading line [",++$line,
		      "] of codon usage file [$file].");

	next if(/^#/ || /^\s*$/);

	chomp;
	s/^ +//g;
	s/ +$//g;

	my @data = split(/ *\t */,$_);

	if(scalar(@data) < 3)
	  {
	    error("Unable to parse line [$line] of codon usage file ",
		  "[$file].  Skipping.",
		  {DETAIL => ("Not enough tab-delimited columns on line: [" .
			      join(@data,"\t") . "].  A minimum of 3 on " .
			      "each line are required.")});
	    next;
	  }
	elsif(length($data[0]) != 1)
	  {
	    error("The first column on line [$line] of codon usage file ",
		  "[$file] is not a single character.  Skipping.",
		  {DETAIL => ("The first column should contain the single " .
			      "character amino acid codes.  Any single " .
			      "characters are acceptable as long as they " .
			      "match what's in the sequence file.")});
	    next;
	  }
	elsif($data[2] !~ /\d/)
	  {
	    error("The codon usage score on line [$line]: [$data[2]] of ",
		  "codon usage file [$file] must be a number.  Skipping.",
		  {DETAIL => ("The scores are used to break ties between " .
			      "different codon pairs.  Any number is " .
			      "acceptable.")});
	    next;
	  }
	elsif(exists($cdnchk->{$data[1]}))
	  {
	    error("Redundant codon encountered: [$data[1]] on line [$line] ",
		  "of codon usage file [$file].  Skipping.",
		  {DETAIL => "A codon may only be present once in the file."});
	    next;
	  }

	if($data[1] !~ /^[ATGCBDHVRYKMSWN]{3}$/i)
	  {warning("The codon on line [$line]: [$data[1]] of codon usage ",
		   "file [$file] appears mal-formed.",
		   {DETAIL => ("The second column should contain 3 " .
			       "nucleotide characters for construction of " .
			       "the nucleotide sequence.")})}

	my $rest = [];
	if(scalar(@data) > 3)
	  {$rest = [@data[3..$#data]]}

	$cdnchk->{$data[1]} = 1;
	$got_one = 1;

	$usg->{uc($data[0])}->{uc($data[1])} = {SCORE => $data[2],
						DATA  => $rest};
      }

    closeIn(*USG);

    if($got_one == 0)
      {
	error("No acceptable usage records were found in codon usage file: ",
	      "[$file].",
	      {DETAIL => ("Please consult the codon usage file format in " .
			  "the --help output.")});
	return(undef);
      }

    return($usg);
  }

#Given an AA and codon, returns the usage score
sub getUsageScore
  {
    my $aa     = uc($_[0]);
    my $codon  = uc($_[1]);
    my $usgobj = $_[2];
    if(scalar(@_) != 3)
      {
	error("3 parameters are required.",
	      {DETAIL => ("A single character amino acid code, a codon, and " .
			  "the usage object.")});
	return(undef);
      }
    elsif(!defined($usgobj))
      {
	error("Usage object (3rd parameter) not defined.");
	return(undef);
      }
    elsif(!defined($aa))
      {return(getUsageScore2($codon,$usgobj))}
    elsif(exists($usgobj->{$aa}) && exists($usgobj->{$aa}->{$codon}))
      {return($usgobj->{$aa}->{$codon}->{SCORE})}
    error("Score not found for amino acid [",(defined($aa) ? $aa : 'undef'),
	  "] and codon [",(defined($codon) ? $codon : 'undef'),"].");
    return(undef);
  }

sub getUsageScore2
  {
    my $codon  = uc($_[0]);
    my $usgobj = $_[1];

    #Assumes that there's only 1 instance of a codon in the entire object
    foreach my $aa (keys(%$usgobj))
      {if(exists($usgobj->{$aa}->{$codon}))
	 {return($usgobj->{$aa}->{$codon}->{SCORE})}}

    return(undef);
  }

sub getMostUsedCodon
  {
    my $aa     = uc($_[0]);
    my $usgobj = $_[1];
    if(scalar(@_) != 2)
      {
	error("2 parameters are required.",
	      {DETAIL => ("A single character amino acid code and the usage " .
			  "object.")});
	return('');
      }
    elsif(!defined($usgobj))
      {
	error("Usage object (2nd parameter) not defined.");
	return('');
      }
    elsif(exists($usgobj->{$aa}))
      {
	if(scalar(keys(%{$usgobj->{$aa}})) == 0)
	  {
	    error("Amino acid [$aa] does not have any codons.");
	    return('');
	  }

	return((sort {$usgobj->{$aa}->{$b}->{SCORE} <=>
			$usgobj->{$aa}->{$a}->{SCORE}}
		keys(%{$usgobj->{$aa}}))[0]);
      }

    error("Score/codon not found for amino acid [",
	  (defined($aa) ? $aa : 'undef'),"].");

    return('NNN');
  }

#Returns a list of codons in order of descending score (then ascending alpha)
sub getCodons
  {
    my $aa     = uc($_[0]);
    my $usgobj = $_[1];
    my $codons = [];
    if(!defined($usgobj))
      {
	error("Usage object (2nd parameter) not defined.");
	return(undef);
      }
    elsif(exists($usgobj->{$aa}))
      {$codons = [sort {$usgobj->{$aa}->{$b} <=> $usgobj->{$aa}->{$a} ||
			  $a cmp $b} keys(%{$usgobj->{$aa}})]}
    else
      {error("Amino acid: [",(defined($aa) ? $aa : 'undef'),
	     "] not found in usage object.")}
    return(wantarray ? @$codons : $codons);
  }

sub getAminoAcids
  {
    my $usgobj = $_[0];
    my $aas    = [];
    if(!defined($usgobj) || ref($usgobj) ne 'HASH')
      {
	error("Usage object (1st parameter) not defined or invalid.");
	return($aas);
      }
    $aas = [keys(%$usgobj)];
    return(wantarray ? @$aas : $aas);
  }

#Takes a usage object and creates a homology matrix object containing the
#scores for the AA combos (which may or may not be used for writing a matrix
#file) and the list of codon pairs by score that are used to construct the NT
#sequence.  Note, there's no reason to ever read a matrix file in this script.
#It is only used by the alignment program.  It doesn't contain the codon pairs
#used to construct an NT sequence from an aligned AA sequence.  So if a matrix
#file is supplied, it is only used by the aligner and sequence construction
#happens using the codon usage table.
sub createMatrix
  {
    my $usgobj = $_[0];
    if(!defined($usgobj) || ref($usgobj) ne 'HASH')
      {
	error("Invalid usage object of type [",
	      (defined($usgobj) ? (ref($usgobj) eq '' ? 'SCALAR' :
				   ref($usgobj) . ' reference') : 'undef'),
	      "] passed in.");
	return(undef);
      }

    if($weighting_method eq 'crossover1')
      {return(createMatrixCrossover1($usgobj))}
    else
      {error('Unrecognized weighting method: [',
	     (defined($weighting_method) ? $weighting_method : 'undef'),'].')}

    return(undef);
  }

#Should only be called from createMatrix
#Global used: $flex_scores
sub createMatrixCrossover1
  {
    my $usgobj = $_[0]; #Assume valid
    my $matrix = {};
    my @aalist = getAminoAcids($usgobj);

    #This is deprecated
    my $max_homo_scores = {3 => 0.0,
			   2 => 0.0,
			   1 => 0.0,
			   0 => 0.0};

    #The sum of all the scores in 1 section (plus the max homo score) should be
    #less than a single score from the section above it (plus its homo score)
    $flex_scores     = {3 => {LMR => 100.0},  #perfect match
			2 => {LM  => 20.0,    #left and middle match
			      LR  => 10.0,    #left and right match
			      MR  => 20.0},   #middle and right match
			1 => {L   => 2.0,    #left matches
			      M   => 1.0,    #middle matches
			      R   => 2.0},   #right matches
			0 => {''  => 0.0}};  #no matches

    #The lesser usage score of any pair of codons
    my $max_usg_score   = 0.09;
    my $max_global_usage = getMaxUsage($usgobj);

    foreach my $i (0..$#aalist)
      {
	my $aa1 = $aalist[$i];
	foreach my $j ($i..$#aalist)
	  {
	    my $aa2 = $aalist[$j];
	    #The lookup will always be indexed in alphabetical order
	    my($lesser,$greater) = sort {$a cmp $b} ($aa1,$aa2);

	    my $max_matches = 0;
	    if($i == $j)
	      {$max_matches = 3}
	    else
	      {$max_matches = getMaxMatches([getCodons($lesser,$usgobj)],
					    [getCodons($greater,$usgobj)])}

	    $matrix->{$lesser}->{$greater} =
	      getCrossover1ScoreHash($lesser,
				     $greater,
				     $max_matches,
				     $flex_scores->{$max_matches},
				     $usgobj,
				     $max_usg_score,
				     $max_homo_scores->{$max_matches},
				     $max_global_usage);
	  }
      }

    verbose({LEVEL => 2},matrixToString($matrix));

    return($matrix);
  }

sub matrixToString
  {
    my $matrix = $_[0];

    #Include the protein weight matrix
    my $string = "Alignment Weight Matrix:\n\n" . pwMatrixToString($matrix) .
      "\n";

    #Include the flex codes used
    $string .= "Codon Match Flexibility scores:\n\n";
    foreach my $matches (sort {$b <=> $a} keys(%$flex_scores))
      {
	foreach my $flex_code (sort {length($b) <=> length($a) || $a cmp $b}
			       keys(%{$flex_scores->{$matches}}))
	  {$string .= "$flex_code\t$flex_scores->{$matches}->{$flex_code}\n"}
      }
    $string .= "\n";

    #Include the sequence construction matrix
    $string .= "Codon Sequence Construction Matrix:\n\n";

    #Keep track of widest cell width
    my $widest = 0;

    my $data = [['',sort {$a eq '*' || $b eq '*' ?
			    $b cmp $a : $a cmp $b} keys(%$matrix)]];
    foreach my $aa1 (sort {$a eq '*' || $b eq '*' ? $b cmp $a : $a cmp $b}
		     keys(%$matrix))
      {
	if(length($aa1) > $widest)
	  {$widest = length($aa1)}
	push(@$data,[$aa1]);
	foreach my $aa2 (sort {$a eq '*' || $b eq '*' ? $b cmp $a : $a cmp $b}
			 keys(%$matrix))
	  {
	    if(length($aa2) > $widest)
	      {$widest = length($aa2)}
	    my($lesser,$greater) = sort {$a cmp $b} ($aa1,$aa2);
	    my $cell = "$matrix->{$lesser}->{$greater}->{SCORE}\n" .
	      "$lesser/$greater\n" .
		"$matrix->{$lesser}->{$greater}->{BESTPAIR}->[0]/" .
		  "$matrix->{$lesser}->{$greater}->{BESTPAIR}->[1]\n" .
		    "$matrix->{$lesser}->{$greater}->{HOMO}\n";
	    if(length($matrix->{$lesser}->{$greater}->{SCORE}) > $widest)
	      {$widest = length($matrix->{$lesser}->{$greater}->{SCORE})}
	    if(length("$lesser/$greater") > $widest)
	      {$widest = length("$lesser/$greater")}
	    if(length("$matrix->{$lesser}->{$greater}->{BESTPAIR}->[0]/" .
		      "$matrix->{$lesser}->{$greater}->{BESTPAIR}->[1]") >
	       $widest)
	      {$widest = length($matrix->{$lesser}->{$greater}->{BESTPAIR}->[0]
				. "/" .
				$matrix->{$lesser}->{$greater}->{BESTPAIR}
				->[1])}
	    if(length($matrix->{$lesser}->{$greater}->{HOMO}) > $widest)
	      {$widest = length($matrix->{$lesser}->{$greater}->{HOMO})}
	    foreach my $flex_code (keys(%{$matrix->{$lesser}->{$greater}
					    ->{PAIRS}}))
	      {
		my $best = "$flex_code:";
		if(length("$flex_code:") > $widest)
		  {$widest = length("$flex_code:")}
		$best .= "\n";
		$best .= join(",\n",
			      map {my $s = "  $_->[0]/$_->[1]=$_->[2]";
				   if((length($s) + 1) > $widest)
				     {$widest = length($s) + 1}   #1 for comma
				   $s}
			      @{$matrix->{$lesser}->{$greater}->{PAIRS}
				  ->{$flex_code}});
		$cell .= $best . "\n";
	      }
	    push(@{$data->[-1]},$cell);
	  }
      }

    $string .= array2DToString($data,$widest,' ') . "\n";

    return($string);
  }

sub array2DToString
  {
    my $array  = $_[0];
    my $width  = $_[1];
    my $spacer = $_[2];

    my $string = '';

    foreach my $row (@$array)
      {
	my $height = (sort {$b <=> $a} map {scalar(split(/\n/,$_))} @$row)[0];
	if($height == 0)
	  {$string .= "\n"}
	else
	  {
	    foreach my $i (0..($height - 1))
	      {
		my $subrow = [map {my $cell=[split(/\n/,$_)];
				   if(scalar(@$cell) <= $i){''}
				   else{$cell->[$i]}} @$row];
		$string .=
		  join($spacer,
		       map {my $l = ((ref($width) eq 'ARRAY' ?
				      $width->[$_] : $width) -
				     length($subrow->[$_]));
			    $subrow->[$_] . ($l < 1 ? '' : (' ' x $l))}
		       (0..$#{$subrow})) . "\n";
	      }
	  }
      }

    return($string);
  }

sub pwMatrixToString
  {
    my $matrix = $_[0];
    my $colwid = 2;

    my $string =
      '   ' . join('  ',sort {$a eq '*' || $b eq '*' ?
				$b cmp $a : $a cmp $b} keys(%$matrix)) . "\n";
    foreach my $aa1 (sort {$a eq '*' || $b eq '*' ? $b cmp $a : $a cmp $b}
		     keys(%$matrix))
      {
	$string .= $aa1;
	foreach my $aa2 (sort {$a eq '*' || $b eq '*' ? $b cmp $a : $a cmp $b}
			 keys(%$matrix))
	  {
	    my($lesser,$greater) = sort {$a cmp $b} ($aa1,$aa2);
	    my $score = roundInt($matrix->{$lesser}->{$greater}->{SCORE});
	    my $len = length($score);
	    my $l = $colwid - $len;
	    $string .= (' ' .                         #Column spacer
			($l < 1 ? '' : (' ' x $l)) .  #Right-align
			$score);                      #Score
	  }
	$string .= "\n";
      }

    return($string);
  }

#Returns a hash with keys SCORE and PAIRS.  The value for SCORE is a number.
#The value for PAIRS is a hash keyed on the category of flex_score keys that
#correspond to max_matches
sub getCrossover1ScoreHash
  {
    my $lesser           = $_[0];
    my $greater          = $_[1];
    my $matches          = $_[2];
    my $flex_scores      = $_[3];
    my $usgobj           = $_[4];
    my $max_usg_score    = $_[5];
    my $homo_score       = $_[6];
    my $max_global_usage = $_[7];
    my $hash             = {SCORE    => 0.0,
			    HOMO     => $homo_score,
			    BESTPAIR => [],
			    PAIRS    => {}}; #{PAIRS => {$flex_code =>
                                             #                [[cdn1,cdn2,$sc],
                                             #                 [cdn1,cdn2,$sc],
                                             #                 ...]}}
                                             #  The keys are normalized lesser
                                             #  usage scores. Only codon pairs
                                             #  at the $matches level are
                                             #  included

    my $testcase = $lesser eq 'A' && $greater eq 'C';

    my $codons1  = [getCodons($lesser,$usgobj)];  #Assumed to be sorted
    my $codons2  = [getCodons($greater,$usgobj)];

    debug("Best homology among codons: [$matches].") if($testcase);
    debug("codons for $lesser: [",join(',',@$codons1),"].") if($testcase);
    debug("codons for $greater: [",join(',',@$codons2),"].") if($testcase);

    #This hash will track whether a codon pair with a matching pattern defined
    #by the flex scores exists or not, e.g. ATT/GTT is an MR match
    my $best_flex_min_usage = {};
    my $best_pair           = [];

    foreach my $codon1 (@$codons1)
      {
	my $min_usage_score = getUsageScore($lesser,$codon1,$usgobj);

	foreach my $codon2 (grep {getNumCodonMatches($codon1,$_) == $matches}
			    @$codons2)
	  {
	    #Determine the lesser of the 2 usage scores for this codon pair
	    if(getUsageScore($greater,$codon2,$usgobj) < $min_usage_score)
	      {$min_usage_score = getUsageScore($greater,$codon2,$usgobj)}

	    #Save the greater of the min usage scores among all pairs
	    my $flex_code = getCrossover1FlexCode($codon1,$codon2);
	    if(!exists($best_flex_min_usage->{$flex_code}) ||
	       $min_usage_score > $best_flex_min_usage->{$flex_code})
	      {
		$best_flex_min_usage->{$flex_code} = $min_usage_score;
		$best_pair = [$codon1,$codon2];
	      }

	    my $indiv_score =
	      (($min_usage_score / $max_global_usage) * $max_usg_score);
	    push(@{$hash->{PAIRS}->{$flex_code}},
		 [$codon1,$codon2,$indiv_score]);
	  }
      }

    $hash->{BESTPAIR} = $best_pair;

    #Initialize the score with the static homology score
    my $score = $homo_score;

    debug("Homology score: [$homo_score].") if($testcase);

    #Add the incremental flexibility scores to the total score
    foreach my $flex_code (keys(%$best_flex_min_usage))
      {
	debug("Flex score for code [$flex_code]: [$flex_scores->{$flex_code}].") if($testcase);

	debug("Usage score: [($best_flex_min_usage->{$flex_code} / $max_global_usage) * $max_usg_score = ",(($best_flex_min_usage->{$flex_code} / $max_global_usage) * $max_usg_score),"].") if($testcase);

	#Add the flexibility score
	$score += $flex_scores->{$flex_code} +

	  #And add the usage score, normalized to make the max available usage
	  #score a max of the available allotted max score (max_usg_score)
	  (($best_flex_min_usage->{$flex_code} / $max_global_usage) *
	   $max_usg_score);
      }

    $hash->{SCORE} = $score;

    return($hash);
  }

sub getMaxUsage
  {
    my $usgobj = $_[0];
    my $max    = 0;

    foreach my $aa (keys(%$usgobj))
      {
	foreach my $codon (keys(%{$usgobj->{$aa}}))
	  {
	    my $score = getUsageScore($aa,$codon,$usgobj);
	    if($score > $max)
	      {$max = $score}
	  }
      }

    return($max);
  }

sub getMaxMatches
  {
    my $codons1 = $_[0]; #Assumed to be sorted
    my $codons2 = $_[1];
    my $max     = 0;
    if(!defined($codons1)     || !defined($codons2) ||
       scalar(@$codons1) == 0 || scalar($codons2) == 0)
      {return($max)}
    foreach my $codon1 (@$codons1)
      {
	foreach my $codon2 (@$codons2)
	  {
	    my $ms = getNumCodonMatches($codon1,$codon2);
	    #Assuming the codons are sorted, any match greater than 1 can be
	    #immediately returned because a match of 3 is the best possible and
	    #if a match of 3 is not immediately obtained, then a match of 2 is
	    #the next best possible
	    if($ms > 1)
	      {return($ms)}
	    elsif($ms > $max)
	      {$max = $ms}
	  }
      }

    return($max);
  }

sub getNumCodonMatches
  {
    my $codon1 = $_[0];
    my $codon2 = $_[1];
    my $num    = 0;
    $num++ if(substr($codon1,0,1) eq substr($codon2,0,1));
    $num++ if(substr($codon1,1,1) eq substr($codon2,1,1));
    $num++ if(substr($codon1,2,1) eq substr($codon2,2,1));
    return($num);
  }

sub getNumEndCodonMatches
  {
    my $codon1 = $_[0];
    my $codon2 = $_[1];
    my $num    = 0;
    $num++ if(             substr($codon1,2,1) eq substr($codon2,2,1));
    $num++ if($num == 1 && substr($codon1,1,1) eq substr($codon2,1,1));
    $num++ if($num == 2 && substr($codon1,0,1) eq substr($codon2,0,1));
    return($num);
  }

sub getCrossover1FlexCode
  {
    my $codon1 = $_[0];
    my $codon2 = $_[1];
    my $code = '';
    $code .= 'L' if(substr($codon1,0,1) eq substr($codon2,0,1));
    $code .= 'M' if(substr($codon1,1,1) eq substr($codon2,1,1));
    $code .= 'R' if(substr($codon1,2,1) eq substr($codon2,2,1));
    return($code);
  }

sub writeMatrixFile
  {
    my $file   = $_[0];
    my $matrix = $_[1];
    my $colwid = 2;

    if(!defined($file))
      {return(0)}

    openOut(*MATO,$file,0) || return(0);

    print MATO << "END_HEAD"
#This matrix was produced by codonHomologizer.pl version 1.0 [8-Jul-17]
#Crossover Homology and Usage Matrix
#Max possible score (though may not be present based on usage omissions) = 25
#Min possible score = 0
#
END_HEAD
      ;

    print MATO ('   ',join('  ',sort {$a eq '*' || $b eq '*' ?
					$b cmp $a : $a cmp $b} keys(%$matrix)),
		"\n");

    foreach my $aa1 (sort {$a eq '*' || $b eq '*' ? $b cmp $a : $a cmp $b}
		     keys(%$matrix))
      {
	print MATO $aa1;
	foreach my $aa2 (sort {$a eq '*' || $b eq '*' ? $b cmp $a : $a cmp $b}
			 keys(%$matrix))
	  {
	    my($lesser,$greater) = sort {$a cmp $b} ($aa1,$aa2);
	    my $score = roundInt($matrix->{$lesser}->{$greater}->{SCORE});
	    my $len = length($score);
	    my $l = $colwid - $len;
	    print MATO (' ',                          #Column spacer
			($l < 1 ? '' : (' ' x $l)),   #Right-align
			$score);                      #Score
	  }
	print MATO "\n";
      }

    closeOut(*MATO);

    return(1);
  }

sub roundInt
  {
    my $num = $_[0];
    if($num !~ /^[0-9\.\-\+]*$/)
      {
	error("Invalid number: [$num].",
	      {DETAIL => "Only decimal and integer values are allowed."});
	return(undef);
      }
    if($num !~ /\./)
      {return($num)}
    elsif($num =~ /\.[01234]/)
      {return(int($num))}
    return(int($num) + 1);
  }

sub getTempFile
  {
    my $seed = $_[0];
    $seed =~ s%.*/%%;
    my $file = getTmpDir();
    $file .= ($file =~ m%/$% ? '' : '/') . $seed . $tmp_suffix;
    return($file);
  }

sub constructNTPairs
  {
    my $alnSeqs    = $_[0]; #array of 2 (or 3) element arrays with defline,
                            #sequence, and quality
    my $matrix_obj = $_[1]; #complex hash object
    my $usage_obj  = $_[2]; #simple usage hash object

    my $pairs      = [];

    if(scalar(@$alnSeqs) < 2)
      {
	error("Multiple aligned sequences are required to contruct all pairs ",
	      "of sequences.");
	return($pairs);
      }

    foreach my $i (0..($#{$alnSeqs} - 1))
      {
	my $rec1 = $alnSeqs->[$i];
	foreach my $j (($i+1)..$#{$alnSeqs})
	  {
	    my $rec2 = $alnSeqs->[$j];
	    my $pair = constructNTPair($rec1,$rec2,$matrix_obj,$usage_obj);
	    push(@$pairs,$pair) if(defined($pair));
	  }
      }

    return($pairs);
  }

#Creates pairs of NT sequence records from an array of aligned AA sequences
#Globals used: $align_codons
sub constructNTPair
  {
    my $rec1 = $_[0];
    my $rec2 = $_[1];
    my $mat  = $_[2];
    my $usg  = $_[3];

    my $aaseq1 = $rec1->[1];
    my $aaseq2 = $rec2->[1];

    if(length($aaseq1) != length($aaseq2))
      {
	error("The length of the aligned sequences [$rec1->[0]] and ",
	      "[$rec2->[0]] must be the same.");
	return(undef);
      }

    my $ntseq1       = '';
    my $ntseq2       = '';
    my $stretch_size = 0;

    foreach my $ptn (0..(length($aaseq1) - 1))
      {
	my $aa1 = substr($aaseq1,$ptn,1);
	my $aa2 = substr($aaseq2,$ptn,1);
	my $codon1 = '';
	my $codon2 = '';
	my $flex_code = '';

	if($aa1 eq '-' && $aa2 ne '-')
	  {
	    $codon1 = '---';
	    $codon2 = getMostUsedCodon($aa2,$usg);
	    $stretch_size = 0;
	  }
	elsif($aa1 ne '-' && $aa2 eq '-')
	  {
	    $codon1 = getMostUsedCodon($aa1,$usg);
	    $codon2 = '---';
	    $stretch_size = 0;
	  }
	elsif($aa1 eq '-' && $aa2 eq '-')
	  {
	    $codon1 = '---';
	    $codon2 = '---';
	  }
	#When equal, it's a no-brainer - best usage score
	elsif($aa1 eq $aa2)
	  {
	    $codon1 = $codon2 = getMostUsedCodon($aa1,$usg);

	    debug({LEVEL => 2},"Picked Best Used codon: [$codon1] from among ",
		  "[(",join('),(',map {"$_,$usg->{$aa1}->{$_}->{SCORE}"}
			    keys(%{$usg->{$aa1}})),")].");

	    $stretch_size += 3;
	  }
	else
	  {
	    my $can_extend = getExtendLength($mat,$aa1,$aa2);
	    my $can_start  = getStartLength($mat,$aa1,$aa2);

	    #The following numbers are based on 2 references from a book
	    #(chapter 3: Shuffle Optimizer: A Program to Optimize DNA Shuffling
	    #for Protein Engineering), claiming that crossovers won't happen
	    #without at least 5 bases of contiguous identity and occur best
	    #when contiguous identity is > 10nts:

	    #Moore GL, Maranas CD (2000) Modeling DNA mutation and
	    #recombination for directed evolution experiments. J Theor Biol
	    #205(3):483-503. doi:10.1006/jtbi.2000.2082

	    #He L, Friedman AM, Bailey-Kellogg C (2012) Algorithms for
	    #optimizing cross-overs in DNA shuffling. BMC Bioinformatics
	    #13(Suppl 3):S3. doi:10.1186/1471-2105-13-s3-s3

	    if(($stretch_size    == 3  && $can_extend   == 2)  ||
		  ($stretch_size == 4  && $can_extend   == 1)  ||
		  ($stretch_size >  4  && $stretch_size <  11) ||
		  ($stretch_size >= 11 && $can_extend   >= $can_start))
	      {($codon1,$codon2) = getBestExtendPair($mat,$aa1,$aa2)}
	    else
	      {($codon1,$codon2) = getBestStartPair($mat,$aa1,$aa2)}

	    $stretch_size = getNumEndCodonMatches($codon1,$codon2);
	  }

	$ntseq1 .= $codon1;
	$ntseq2 .= $codon2;
      }

    my $def1 = $rec1->[0];
    my $def2 = $rec2->[0];
    $def1 =~ s/^[>\@\+]//;
    $def2 =~ s/^[>\@\+]//;

    my $newrec1 = ["$def1 optimized for $def2",$ntseq1];
    my $newrec2 = ["$def2 optimized for $def1",$ntseq2];

    return(wantarray ? ($newrec1,$newrec2) : [$newrec1,$newrec2]);
  }

sub getExtendLength
  {
    my $mat = $_[0];
    my $aa1 = $_[1];
    my $aa2 = $_[2];

    my($lesser,$greater) = sort {$a cmp $b} ($aa1,$aa2);

    if(!exists($mat->{$lesser}) || !exists($mat->{$lesser}->{$greater}))
      {
	error("The amino acid pair: [$aa1/$aa2] does not exist in the ",
	      "protein weight matrix.");
	return(0);
      }

    if(exists($mat->{$lesser}->{$greater}->{PAIRS}->{LMR}))
      {return(3)}
    elsif(exists($mat->{$lesser}->{$greater}->{PAIRS}->{LM}))
      {return(2)}
    elsif(exists($mat->{$lesser}->{$greater}->{PAIRS}->{L}) ||
	  exists($mat->{$lesser}->{$greater}->{PAIRS}->{LR}))
      {return(1)}
    else
      {return(0)}
  }

sub getBestExtendPair
  {
    my $mat = $_[0];
    my $aa1 = $_[1];
    my $aa2 = $_[2];

    my($lesser,$greater) = sort {$a cmp $b} ($aa1,$aa2);
    my $swap = ($aa2 lt $aa1);

    if(!exists($mat->{$lesser}) || !exists($mat->{$lesser}->{$greater}))
      {
	error("The amino acid pair: [$aa1/$aa2] does not exist in the ",
	      "protein weight matrix.");
	return(wantarray ? ('NNN','NNN') : ['NNN','NNN']);
      }

    my $flex_code = '';
    my($best_codon_pair);
    if(exists($mat->{$lesser}->{$greater}->{PAIRS}->{LMR}))
      {$flex_code = 'LMR'}
    elsif(exists($mat->{$lesser}->{$greater}->{PAIRS}->{LM}))
      {$flex_code = 'LM'}
    elsif(exists($mat->{$lesser}->{$greater}->{PAIRS}->{LR}))
      {$flex_code = 'LR'}
    elsif(exists($mat->{$lesser}->{$greater}->{PAIRS}->{L}))
      {$flex_code = 'L'}
    else
      {return(wantarray ? getBestStartPair($mat,$aa1,$aa2,1) :
	      scalar(getBestStartPair($mat,$aa1,$aa2,1)))}

    if(!exists($mat->{$lesser}->{$greater}->{PAIRS}) ||
       !exists($mat->{$lesser}->{$greater}->{PAIRS}->{$flex_code}) ||
       ref($mat->{$lesser}->{$greater}->{PAIRS}->{$flex_code}) ne 'ARRAY' ||
       scalar(@{$mat->{$lesser}->{$greater}->{PAIRS}->{$flex_code}}) == 0)
      {
	error("The amino acid pair: [$aa1/$aa2] does not have any codon ",
	      "pairs protein weight matrix.");

	return(wantarray ? ('NNN','NNN') : ['NNN','NNN']);
      }

    $best_codon_pair =
      (map {$swap ? [$_->[1],$_->[0]] : [$_->[0],$_->[1]]}
       sort {$b->[2] <=> $a->[2]}
       @{$mat->{$lesser}->{$greater}->{PAIRS}->{$flex_code}})[0];

    debug({LEVEL => 2},"Picked Best Extend pair: [",
	  join(',',@$best_codon_pair),"] from among [(",
	  join('),(', map {join(',',@$_)}
	       map {$swap ?
		      [$_->[1],$_->[0],$_->[2]] :
			[$_->[0],$_->[1],$_->[2]]} sort {$b->[2] <=> $a->[2]}
	       @{$mat->{$lesser}->{$greater}->{PAIRS}->{$flex_code}}),")].");

    return(wantarray ? @$best_codon_pair : $best_codon_pair);
  }

sub getStartLength
  {
    my $mat = $_[0];
    my $aa1 = $_[1];
    my $aa2 = $_[2];

    my($lesser,$greater) = sort {$a cmp $b} ($aa1,$aa2);

    if(!exists($mat->{$lesser}) || !exists($mat->{$lesser}->{$greater}))
      {
	error("The amino acid pair: [$aa1/$aa2] does not exist in the ",
	      "protein weight matrix.");
	return(0);
      }

    if(exists($mat->{$lesser}->{$greater}->{PAIRS}->{LMR}))
      {return(3)}
    elsif(exists($mat->{$lesser}->{$greater}->{PAIRS}->{MR}))
      {return(2)}
    elsif(exists($mat->{$lesser}->{$greater}->{PAIRS}->{R}) ||
	  exists($mat->{$lesser}->{$greater}->{PAIRS}->{LR}))
      {return(1)}
    else
      {return(0)}
  }

sub getBestStartPair
  {
    my $mat = $_[0];
    my $aa1 = $_[1];
    my $aa2 = $_[2];
    my $end = defined($_[3]) ? $_[3] : 0;

    my($lesser,$greater) = sort {$a cmp $b} ($aa1,$aa2);
    my $swap = ($aa2 lt $aa1);

    if(!exists($mat->{$lesser}) || !exists($mat->{$lesser}->{$greater}))
      {
	error("The amino acid pair: [$aa1/$aa2] does not exist in the ",
	      "protein weight matrix.");
	return(wantarray ? ('NNN','NNN') : ['NNN','NNN']);
      }

    my $flex_code = '';
    my($best_codon_pair);
    if(exists($mat->{$lesser}->{$greater}->{PAIRS}->{LMR}))
      {$flex_code = 'LMR'}
    elsif(exists($mat->{$lesser}->{$greater}->{PAIRS}->{MR}))
      {$flex_code = 'MR'}
    elsif(exists($mat->{$lesser}->{$greater}->{PAIRS}->{LR}))
      {$flex_code = 'LR'}
    elsif(exists($mat->{$lesser}->{$greater}->{PAIRS}->{R}))
      {$flex_code = 'R'}
    else
      {
	if(!exists($mat->{$lesser}->{$greater}->{BESTPAIR}) ||
	   !defined($mat->{$lesser}->{$greater}->{BESTPAIR}) ||
	   ref($mat->{$lesser}->{$greater}->{BESTPAIR}) ne 'ARRAY' ||
	   scalar(@{$mat->{$lesser}->{$greater}->{BESTPAIR}}) < 2)
	  {
	    error("The amino acid pair: [$aa1/$aa2] does not have codons in ",
		  "the protein weight matrix.");
	    return(wantarray ? ('NNN','NNN') : ['NNN','NNN']);
	  }

	$best_codon_pair = $mat->{$lesser}->{$greater}->{BESTPAIR};
	if($swap)
	  {$best_codon_pair =
	     [reverse(@{$mat->{$lesser}->{$greater}->{BESTPAIR}})]}

	debug({LEVEL => 2},"Picked Best Overall(start) pair",
	      ($end ? ' (because we could not extend or start)' : ''),": [",
	      join(',',@$best_codon_pair),"].");

	return(wantarray ? @$best_codon_pair : $best_codon_pair);
      }

    if(!exists($mat->{$lesser}->{$greater}->{PAIRS}) ||
       !exists($mat->{$lesser}->{$greater}->{PAIRS}->{$flex_code}) ||
       ref($mat->{$lesser}->{$greater}->{PAIRS}->{$flex_code}) ne 'ARRAY' ||
       scalar(@{$mat->{$lesser}->{$greater}->{PAIRS}->{$flex_code}}) == 0)
      {
	error("The amino acid pair: [$aa1/$aa2] does not have any codon ",
	      "pairs protein weight matrix.");

	return(wantarray ? ('NNN','NNN') : ['NNN','NNN']);
      }

    $best_codon_pair =
      (map {$swap ? [$_->[1],$_->[0]] : [$_->[0],$_->[1]]}
       sort {$b->[2] <=> $a->[2]}
       @{$mat->{$lesser}->{$greater}->{PAIRS}->{$flex_code}})[0];

    debug({LEVEL => 2},"Picked Best Start pair",
	  ($end ? ' (because we could not extend)' : ''),": [",
	  join(',',@$best_codon_pair),"] from among [(",
	  join('),(', map {join(',',@$_)}
	       map {$swap ?
		      [$_->[1],$_->[0],$_->[2]] :
			[$_->[0],$_->[1],$_->[2]]} sort {$b->[2] <=> $a->[2]}
	       @{$mat->{$lesser}->{$greater}->{PAIRS}->{$flex_code}}),")].");

    return(wantarray ? @$best_codon_pair : $best_codon_pair);
  }

sub writeNTPairs
  {
    my $nt_pairs = $_[0];
    my $outfile  = $_[1];

    openOut(*NT,$outfile);

    foreach my $pair (@$nt_pairs)
      {
	foreach my $rec (@$pair)
	  {
	    my $defline = $rec->[0];
	    my $sequence = $rec->[1];

	    #Not sure if beginning defline character is included in the defline
	    #in the record, but we'll remove it in case it's there
	    $defline =~ s/^[>\@\+]//;
	    chomp($defline);
	    chomp($sequence);

	    unless($align_codons)
	      {$sequence =~ s/-+//g}

	    print(">$defline\n");
	    print(formatSequence($sequence));
	  }
      }

    closeOut(*NT);
  }






#Must be sure to output alignments in fasta format
#Use globals for alignment executable
#Globals used: $align_tool
sub aaAlign
  {
    my $aafile     = $_[0];
    my $matfile    = $_[1];
    my $alnoutfile = $_[2];

    if($align_tool eq 'muscle')
      {getMuscleMultipleAlignment($aafile,$matfile,$alnoutfile)}
    elsif($align_tool eq 'clustalw')
      {
	error("Clustalw has not yet been configured to run in ",
	      "codonHomologizer.  Please use muscle.");
      }
    else
      {
	error("Unidentified alignment tool: [$align_tool].",
	      {DETAIL => ("Only these tools are supported: [" .
			  join(',',@$align_tools) . "].")});
      }

###FINISH
  }

sub aaSegmentAlign
  {
    my $aaseqs       = $_[0];
    my $matfile      = $_[1];
    my $alnoutfile   = $_[2];
    my $segments_obj = $_[3];

    my $segfilebase   = getTempFile('segment');
    my $segment_files = [];
    my $first         = 1;

    #For each segment in the sequences, create a sequence file to feed to the
    #alignment executable
    my $alnseqs = {};
    my $cnt     = 0;
    my $counts  = {};
    my $seen   = {};
    foreach my $rec (@$aaseqs)
      {
	my $def = $rec->[0];
	my $id  = parseIDFromDef($def);
	my $seq = $rec->[1];
	$alnseqs->{$id} = {DEF   => $def,
			   SEQ   => '',
			   ORDER => $cnt++};

	if(!exists($segments_obj->{ALL}) ||
	   !exists($segments_obj->{ALL}->{$id}))
	  {
	    error("Sequence ID [$id] not found (or mal-formed) in the ",
		  "segments object.",
		  {DETAIL => ("Please check your segments file (-d) and " .
			      "make sure the IDs in the first column match " .
			      "the IDs in your amino acid sequence files " .
			      "(-i).")});
	    next;
	  }

	my $segnum = 0;
	foreach my $coordpair (@{$segments_obj->{ALL}->{$id}})
	  {
	    $segnum++;
	    my $segfile = $segfilebase . '.' . $segnum . '.fna';

	    #This assumes the same number of coordinate pairs for each sequence
	    #which is a safe assumption, because that's all checked when the
	    #segments file is checked when it's read in
	    if($first)
	      {push(@$segment_files,$segfile);}

	    #The segment doesn't exist for a sequence if it is undefined, which
	    #means the user's coordinates abutt one another in at least one
	    #sequence and not in another
	    if(defined($coordpair))
	      {
		openOut(HANDLE => *SEGOUT,
			FILE   => $segfile,
			SELECT => 1,
			QUIET  => isDebug(),
			HEADER => 0,
			APPEND => exists($seen->{$segfile}));

		$seen->{$segfile}++;

		#This counts the number of sequences written to each segment
		#file
		$counts->{$segfile}++;

		print(">$id\n",
		      substr($seq,
			     $coordpair->[0] - 1,
			     $coordpair->[1] - $coordpair->[0] + 1),"\n");

		closeOut(*SEGOUT);
	      }
	  }

	$first = 0;
      }

    #Now we can align each of the segment files.  This assumes any completely
    #empty segments have been removed (i.e. when all the user's pairs of
    #coordinates abutt one another).  The removal is done by the last loop in
    #insertGapCoords()
    foreach my $segfile (@$segment_files)
      {
	my $segoutfile = $segfile;

	#If there is more than 1 sequence in the current segment fasta file
	if($counts->{$segfile} > 1)
	  {
	    $segoutfile .= '.aln';

	    #Perform the alignment
	    aaAlign($segfile,$matfile,$segoutfile);
	  }
	#Else there's only 1 sequence - so we'll be reading in the segment
	#file's 1 sequence and fill the rest of the sequences with gaps
	#(inefficient, I know, but easier to code)

	#Read in and append the sequences
	openIn(*AAALNSEG,$segoutfile) || next;

	#Turn the records into a hash for easy lookup, and grab the length of
	#the alignment while doing it (assuming the length of each aligned
	#sequence is consistent
	my $len = 0;
	my $segrecs = {map {$len = length($_->[1]) unless($len > 0);
			    parseIDFromDef($_->[0]) => $_}
		       getNextSeqRec(*AAALNSEG,0,$segoutfile)};

	#Close the file
	closeIn(*AAALNSEG);

	#Append each sequence (or gap characters)
	foreach my $id (keys(%$alnseqs))
	  {
	    if(exists($segrecs->{$id}))
	      {
		chomp($segrecs->{$id}->[1]);
		$alnseqs->{$id}->{SEQ} .= $segrecs->{$id}->[1];
	      }
	    else
	      {$alnseqs->{$id}->{SEQ} .= '-' x $len}
	  }
      }

    #Finally, write out the aligned and stitched together segments to one
    #output alignment file
    openOut(*ALNOUT,$alnoutfile) || return(0);
    foreach my $id (sort {$alnseqs->{$a}->{ORDER} <=>
			    $alnseqs->{$b}->{ORDER}} keys(%$alnseqs))
      {print("$alnseqs->{$id}->{DEF}\n$alnseqs->{$id}->{SEQ}\n")}
    closeOut(*ALNOUT);

    return(1);
  }


#Returns alignment in clustalw format of 2 very similar sequences using muscle
#Uses global: $muscle, $verbose
sub getMuscleMultipleAlignment
  {
    my $fasta      = $_[0];
    my $matrix     = $_[1];
    my $outfile    = $_[2];
    my $muscle_exe = (defined($align_exe) && $align_exe ne '' ?
		      $align_exe : 'muscle');

    my $muscle_command = "$muscle_exe -in '$fasta' -out '$outfile' " .
      "-matrix '$matrix' -fasta " .
	(isVerbose() ? '' : ' -quiet') .
	  ($align_opts ? ' ' . $align_opts : '');

    verbose("Running muscle.");

    debug($muscle_command);

    my $output = `$muscle_command`;

    debug("Muscle alignment done:\n$output\n");

    return($output);
  }






sub getExe
  {
    my $cmd = $_[0];
    my $exe = '';

    if(eval("use File::Which;1;") ||
       eval("use local::lib;use File::Which;1;"))
      {$exe = which($cmd)}
    else
      {
	verbose("File::Which not found, switching to backup method.");
	$exe = `which $cmd`;
	chomp($exe);
	if($exe =~ /which: Command not found./ || $exe !~ /\S/)
	  {
	    error("System command 'which' does not appear to exist.  Please ",
		  "install the perl module File::Which.");
	    $exe = '';
	  }
	elsif($exe =~ /not found/i)
	  {$exe = ''}
      }

    return($exe);
  }

sub incompatible
  {
    my $exe  = $_[0];
    my $tool = $_[1];

    if(!defined($exe) || $exe eq '' || !-e $exe || !-x $exe)
      {
	error("The executable [$exe] appears to either not be in ",
	      "your path, not exist, not have execute permissions, or you ",
	      "have not created a symbolic link to the full ",
	      "name of the executable (i.e. with version number).");
	return(1);
      }

    if($tool eq 'muscle')
      {
	my $version = `$exe -version`;
	chomp($version);
	if($version =~ /MUSCLE v(\S+) by Robert C. Edgar/)
	  {
	    my $vnum = $1;
	    my $confirmed = [3,8,31];
	    my $vnums = [split(/\./,$vnum,-1)];
	    my $ok = 1;
	    my $i = 0;
	    for($i = 0;$i < scalar(@$vnums) && $i < scalar(@$confirmed);$i++)
	      {
		if($vnums->[$i] != $confirmed->[$i])
		  {
		    if($vnums->[$i] < $confirmed->[$i])
		      {$ok = 0}
		    else
		      {$ok = 1}
		    last;
		  }
	      }
	    warning("This script was tested with Muscle version 3.8.31.  ",
		    "Your version appears to be [$vnum], thus it may not ",
		    "work properly.") unless($ok);
	  }
	else
	  {warning("This script was tested with Muscle version 3.8.31.  It ",
		   "may not work properly with the version you are using.")}
      }
    elsif($tool eq 'clustalw')
      {
	my $version = `$exe -help`;
	chomp($version);
	if($version =~ /CLUSTAL\S*\s+(\S+)\s*Multiple Sequence Alignments/)
	  {
	    my $vnum = $1;
	    my $confirmed = [2,1];
	    my $vnums = [split(/\./,$vnum,-1)];
	    my $ok = 1;
	    my $i = 0;
	    for($i = 0;$i < scalar(@$vnums) && $i < scalar(@$confirmed);$i++)
	      {
		if($vnums->[$i] != $confirmed->[$i])
		  {
		    if($vnums->[$i] < $confirmed->[$i])
		      {$ok = 0}
		    else
		      {$ok = 1}
		    last;
		  }
	      }
	    warning("This script was tested with clustalw version 2.1.  ",
		    "Your version appears to be [$vnum], thus it may not ",
		    "work properly.") unless($ok);
	  }
	else
	  {warning("This script was tested with clustalw version 2.1.  It ",
		   "may not work properly with the version you are using.")}
      }
    else
      {error("Unrecognized alignment tool selection: [$tool].  Unable to ",
	     "test compatibility.")}

    return(0);
  }

sub checkAlignOpts
  {
    my $align_opts  = $_[0];
    my $prohib_opts = $_[1];
    my @violations  = ();

    foreach my $prohib_opt (@$prohib_opts)
      {
	if($align_opts =~ /\Q$prohib_opt\E(=| |\Z)/)
	  {push(@violations,$prohib_opt)}
      }

    return(@violations);
  }

#Uses global variables: lastfiletype & seq_format
#Edited to use $seq_format instead of global $filetype variable 7/7/17 -Rob
sub getNextSeqRec
  {
    #The first 2 parameters are passed along to the sub that reads the specific
    #file format detected (input handle and the no_format flag).  The third
    #parameter used by this sub is only for constructing errors and warnings
    #that mention the affected file name
    my $input_file = $_[2];

    my $filetype = $seq_format;

    debug("Determining previous type");

    if(!defined($main::lastfiletype) || $filetype ne 'auto')
      {
	if($filetype eq 'fasta')
	  {$main::getnextsub = \&getNextFastaRec}
	elsif($filetype eq 'fastq')
	  {$main::getnextsub = \&getNextFastqRec}
      }
    elsif(defined($main::lastfiletype) &&
	  exists($main::lastfiletype->{$input_file}))
      {
	if($main::lastfiletype->{$input_file} eq 'fasta')
	  {$main::getnextsub = \&getNextFastaRec}
	elsif($main::lastfiletype->{$input_file} eq 'fastq')
	  {$main::getnextsub = \&getNextFastqRec}
      }

    if($filetype eq 'auto' &&
       (!defined($main::lastfiletype) ||
	!exists($main::lastfiletype->{$input_file})))
      {
	debug("Determining type");

	if($input_file eq '-')
	  {
	    error("`-t auto` cannot be used when the input file is supplied ",
		  "on standard input.  Please supply the exact file type.");
	    quit(2);
	  }

	if(!-e $input_file)
	  {
	    error("`-t auto` cannot be used when the input file does not ",
		  "exist.  Please supply the exact file type.");
	    quit(8);
	  }

	my($num_fastq_defs);
	if(-e $input_file)
	  {
	    $num_fastq_defs =
	      `head -n 50 "$input_file" | grep -c -E '^[\@\+]'`;
	    debug({LEVEL => 2},"System output from: [",
		  qq(head -n 50 "$input_file" | grep -c -E '^[\@\+]'),
		  "]:\n$num_fastq_defs");
	    $num_fastq_defs =~ s/^\D+//;
	    $num_fastq_defs =~ s/\D.*//;
	  }
	else
	  {$num_fastq_defs = 0}

	if($num_fastq_defs > 0)
	  {
	    $main::getnextsub = \&getNextFastqRec;
	    $main::lastfiletype->{$input_file} = 'fastq';
	  }
	else
	  {
	    my($num_fasta_defs);
	    if(-e $input_file)
	      {
		$num_fasta_defs = `head -n 50 "$input_file" | grep -c -E '^>'`;

		debug({LEVEL => 2},"System output from: [",
		      qq(head -n 50 "$input_file" | grep -c -E '^>'),
		      "]:\n$num_fasta_defs");

		$num_fasta_defs =~ s/^\D+//;
		$num_fasta_defs =~ s/\D.*//;
	      }
	    else
	      {$num_fasta_defs = 0}

	    if($num_fasta_defs > 0)
	      {
		$main::getnextsub = \&getNextFastaRec;
		$main::lastfiletype->{$input_file} = 'fasta';
	      }
	    else
	      {
		if(!defined($main::lastfiletype) ||
		   !exists($main::lastfiletype->{$input_file}))
		  {
		    debug("Num fasta deflines: [$num_fasta_defs].");
		    error("Unable to determine file type.  Skipping file ",
			  "[$input_file].");
		    return(undef);
		  }
		warning("Unable to determine file type.  Defaulting to ",
			"[$main::lastfiletype->{$input_file}].");
		if($main::lastfiletype->{$input_file} eq 'fasta')
		  {$main::getnextsub = \&getNextFastaRec}
		else
		  {$main::getnextsub = \&getNextFastqRec}
	      }
	  }
      }

    debug("Returning record");

    return($main::getnextsub->(@_));
  }

#Copied from DNAstiffness.pl on 2/12/2014 so as to be independent -Rob
sub getNextFastaRec
  {
    my $handle    = $_[0];      #File handle or file name
    my $no_format = $_[1];

    if(exists($main::{FASTABUFFER}) && exists($main::{FASTABUFFER}->{$handle}))
      {
	if(scalar(@{$main::{FASTABUFFER}->{$handle}}) > 0)
	  {
	    if(wantarray)
	      {
		my @array = (@{$main::{FASTABUFFER}->{$handle}});
		@{$main::{FASTABUFFER}->{$handle}} = ();
		return(@array);
	      }
	    return(shift(@{$main::{FASTABUFFER}->{$handle}}));
	  }
	elsif(eof($handle))
	  {return(undef)}
      }

    my $parent_id_check = {};
    my $first_loop      = 0;
    my $line_num        = 0;
    my $verbose_freq    = 1000;
    my $line            = '';
    my $defline         = '';
    my $seq_lines       = 0;
    my($seq);

    #For each line in the current input file
    while(getLine($handle))
      {
	$line_num++;

	verboseOverMe("Reading line [$line_num].")
	  unless($line_num % $verbose_freq);

	$line = $_;

	next if($line !~ /\S/ || $line =~ /^\s*#/);
	if($line =~ />/)
	  {
	    if($defline)
	      {
		my $solidseq =
		  ($seq_lines == 1 || $no_format ?
		   $seq : formatSequence($seq));
		chomp($solidseq);
		chomp($defline);

		push(@{$main::{FASTABUFFER}->{$handle}},[$defline,$solidseq]);
	      }
	    $defline   = $line;
	    $seq_lines = 0;

	    my $tmp_id = $defline;
	    $tmp_id =~ s/^\s*>\s*//;
	    $tmp_id =~ s/\s.*//;
	    if($tmp_id eq '')
	      {warning("No Defline ID on line: [$line_num] of current file.  ",
		       " Universal coordinates will be used if some were ",
		       "supplied either via command line arguments of via ",
		       "coordinate file with no parent sequence ID.")}
	    elsif(exists($parent_id_check->{$tmp_id}))
	      {
		error("Two sequences found with the same ID on the ",
		      "defline: [$tmp_id] in current fasta file.  The same ",
		      "pairs of coordinates will be used for each sequence.");
	      }

	    undef($seq);
	  }
	elsif($line =~ /^([^\t]+?) *\t\s*(.*)/)
	  {
	    $defline = $1;
	    $seq     = $2;

	    my $solidseq =
	      ($seq_lines == 1 || $no_format ? $seq : formatSequence($seq));
	    chomp($solidseq);
	    chomp($defline);

	    push(@{$main::{FASTABUFFER}->{$handle}},[$defline,$solidseq]);

	    undef($seq);
	  }
	else
	  {
	    $seq .= $line;
	    $seq_lines++;
	  }
      }

    #Handle the last sequence (if there were any sequences)
    if(defined($seq))
      {
	my $solidseq =
	  ($no_format ? $seq :
	   formatSequence($seq));
	chomp($solidseq);
	chomp($defline);

	push(@{$main::{FASTABUFFER}->{$handle}},[$defline,$solidseq]);
      }

    #Return the first sequence (if sequence was parsed)
    if(exists($main::{FASTABUFFER}) && exists($main::{FASTABUFFER}->{$handle}))
      {
	if(scalar(@{$main::{FASTABUFFER}->{$handle}}) > 0)
	  {
	    if(wantarray)
	      {
		my @array = (@{$main::{FASTABUFFER}->{$handle}});
		@{$main::{FASTABUFFER}->{$handle}} = ();
		return(@array);
	      }
	    return(shift(@{$main::{FASTABUFFER}->{$handle}}));
	  }
	else
	  {return(undef)}
      }
    else
      {return(undef)}
  }

#Copied from DNAstiffness.pl on 2/12/2014 so as to be independent -Rob
sub formatSequence
  {
    #1. Read in the parameters.
    my $sequence          = $_[0];
    my $chars_per_line    = $_[1];
    my $coords_left_flag  = $_[2];
    my $coords_right_flag = $_[3];
    my $start_coord       = $_[4];
    my $coords_asc_flag   = $_[5];
    my $coord_upr_bound   = $_[6];
    my $uppercase_flag    = $_[7];
    my $print_flag        = $_[8];
    my $nucleotide_flag   = $_[9];

    my($formatted_sequence,
       $sub_string,
       $sub_sequence,
       $coord,
       $max_num_coord_digits,
       $line_size_left,
       $lead_spaces,
       $line);
    my $coord_separator = '  ';
    my $tmp_sequence = $sequence;
    $tmp_sequence =~ s/\s+//g;
    $tmp_sequence =~ s/<[^>]*>//g;
    my $seq_len = length($tmp_sequence);

    #2. Error check the parameters and set default values if unsupplied.
    my $default_chars_per_line    = ''; #Infinity
    my $default_coords_left_flag  = 0;
    my $default_coords_right_flag = 0;
    my $default_start_coord       = (!defined($coords_asc_flag) ||
				     $coords_asc_flag ? 1 : $seq_len);
    my $default_coords_asc_flag   = 1;
    my $default_coord_upr_bound   = undef();  #infinity (going past 1 produces
    my $default_uppercase_flag    = undef();  #          negative numbers)
    my $default_print_flag        = 0;

    if(!defined($chars_per_line) || $chars_per_line !~ /^\d+$/)
      {
        if(defined($chars_per_line) &&
	   $chars_per_line !~ /^\d+$/ && $chars_per_line =~ /./)
	  {print("WARNING:seq-lib.pl:formatSequence: Invalid ",
	         "chars_per_line: [$chars_per_line] - using default: ",
		 "[$default_chars_per_line]<BR>\n")}
        #end if(chars_per_line !~ /^\d+$/)
	$chars_per_line = $default_chars_per_line;
      }
    elsif(!$chars_per_line)
      {$chars_per_line = ''}
    #end if(!defined($chars_per_line) || $chars_per_line !~ /^\d+$/)
    if(!defined($coords_left_flag))
      {$coords_left_flag = $default_coords_left_flag}
    #end if(!defined(coords_left_flag))
    if(!defined($coords_right_flag))
      {$coords_right_flag = $default_coords_right_flag}
    #end if(!defined(coords_right_flag))
    if(!defined($start_coord) || $start_coord !~ /^\-?\d+$/)
      {
        if(defined($start_coord) &&
           ($coords_left_flag || $coords_right_flag))
          {print("WARNING:formatSequence.pl:formatSequence: Invalid ",
                 "start_coord: [$start_coord] - using default: ",
                 "[$default_start_coord]\n")}
        #end if($start_coord !~ /^\d+$/)
        $start_coord = $default_start_coord;
      }
    #end if(!defined($start_coord) || $start_coord !~ /^\d+$/)
    if(!defined($coords_asc_flag))
      {$coords_asc_flag = $default_coords_asc_flag}
    #end if(!defined(coords_right_flag))
    if(defined($coord_upr_bound) && $coord_upr_bound !~ /^\d+$/)
      {undef($coord_upr_bound)}
    if(!defined($print_flag))
      {$print_flag = $default_print_flag}
    #end if(!defined($print_flag))

    if(defined($coord_upr_bound) && $start_coord < 1)
      {$start_coord = $coord_upr_bound + $start_coord}
    elsif($start_coord < 1)
      {$start_coord--}
    elsif(defined($coord_upr_bound) && $start_coord > $coord_upr_bound)
      {$start_coord -= $coord_upr_bound}

    #3. Initialize the variables used for formatting.  (See the DATASTRUCTURES
    #   section.)
    if($coords_asc_flag)
      {
        if(defined($coord_upr_bound) &&
           ($seq_len + $start_coord) > $coord_upr_bound)
          {$max_num_coord_digits = length($coord_upr_bound)}
        else
          {$max_num_coord_digits = length($seq_len + $start_coord - 1)}

        $coord = $start_coord - 1;
      }
    else
      {
        if(defined($coord_upr_bound) && ($start_coord - $seq_len + 1) < 1)
          {$max_num_coord_digits = length($coord_upr_bound)}
        elsif(!defined($coord_upr_bound) &&
              length($start_coord - $seq_len - 1) > length($start_coord))
          {$max_num_coord_digits = length($start_coord - $seq_len - 1)}
        else
          {$max_num_coord_digits = length($start_coord)}

        $coord = $start_coord + 1;
      }
    $line_size_left = $chars_per_line;
    $lead_spaces    = $max_num_coord_digits - length($start_coord);

    #5. Add the first coordinate with spacing if coords_left_flag is true.
    $line = ' ' x $lead_spaces . $start_coord . $coord_separator
      if($coords_left_flag);

    #6. Foreach sub_string in the sequence where sub_string is either a
    #   sub_sequence or an HTML tag.
    foreach $sub_string (split(/(?=<)|(?<=>)/,$sequence))
      {
        #6.1 If the substring is an HTML tag
        if($sub_string =~ /^</)
          #6.1.1 Add it to the current line of the formatted_sequence
          {$line .= $sub_string}
        #end if(sub_string =~ /^</)
        #6.2 Else
        else
          {
            $sub_string =~ s/\s+//g;

	    if($nucleotide_flag)
	      {
		my(@errors);
		(@errors) = ($sub_string =~ /([^ATGCBDHVRYKMSWNX])/ig);
		$sub_string =~ s/([^ATGCBDHVRYKMSWNX])//ig;
		if(scalar(@errors))
		  {print STDERR ("WARNING:formatSequence.pl:formatSequence:",
				 scalar(@errors),
				 " bad nucleotide characters were ",
				 "filtered out of your sequence: [",
				 join('',@errors),
				 "].\n")}
	      }

            #6.2.1 If the sequence is to be uppercased
            if(defined($uppercase_flag) && $uppercase_flag)
              #6.2.1.1 Uppercase the sub-string
              {$sub_string = uc($sub_string)}
            #end if(defined($uppercase_flag) && $uppercase_flag)
            #6.2.2 Else if the sequence is to be lowercased
            elsif(defined($uppercase_flag) && !$uppercase_flag)
              #6.2.2.1 Lowercase the sub-string
              {$sub_string = lc($sub_string)}
            #end elsif(defined($uppercase_flag) && !$uppercase_flag)

            #6.2.3 While we can grab enough sequence to fill the rest of a line
            while($sub_string =~ /(.{1,$line_size_left})/g)
              {
                $sub_sequence = $1;
                #6.2.3.1 Add the grabbed sequence to the current line of the
                #        formatted sequence
                $line .= $sub_sequence;
                #6.2.3.2 Increment the current coord by the amount of sequence
                #        grabbed
                my $prev_coord = $coord;
                if($coords_asc_flag)
                  {
                    $coord += length($sub_sequence);
                    if(defined($coord_upr_bound)      &&
                       $prev_coord <= $coord_upr_bound &&
                       $coord > $coord_upr_bound)
                      {$coord -= $coord_upr_bound}
                  }
                else
                  {
                    $coord -= length($sub_sequence);
                    if(defined($coord_upr_bound) &&
                       $prev_coord >= 1 && $coord < 1)
                      {$coord = $coord_upr_bound + $coord - 1}
                    elsif($prev_coord >= 1 && $coord < 1)
                      {$coord--}
                  }
                #6.2.3.3 If the length of the current sequence grabbed
                #        completes a line
                if($line_size_left eq '' ||
		   length($sub_sequence) == $line_size_left)
                  {
                    $lead_spaces = $max_num_coord_digits - length($coord);
                    #6.2.3.3.1 Conditionally add coordinates based on the
                    #          coords flags
                    $line .= $coord_separator . ' ' x $lead_spaces . $coord
                      if($coords_right_flag);

                    #6.2.3.3.2 Add a hard return to the current line of the
                    #          formatted sequence
                    $line .= "\n";

                    #6.2.3.3.3 Add the current line to the formatted_sequence
                    $formatted_sequence .= $line;
                    #6.2.3.3.4 Print the current line if the print_flag is true
                    print $line if($print_flag);

                    #6.2.3.3.5 Start the next line
                    $lead_spaces = $max_num_coord_digits - length($coord+1);
                    $line = '';
                    $line = ' ' x $lead_spaces
                          . ($coords_asc_flag ? ($coord+1) : ($coord-1))
                          . $coord_separator
                      if($coords_left_flag);

                    #6.2.3.3.6 Reset the line_size_left (length of remaining
                    #          sequence per line) to chars_per_line
                    $line_size_left = $chars_per_line;
                  }
                #end if(length($sub_sequence) == $line_size_left)
                #6.2.3.4 Else
                else
                  #6.2.3.4.1 Decrement line_size_left (length of remaining
                  #          sequence per line) by the amount of sequence
                  #          grabbed
                  {$line_size_left -= length($sub_sequence)}
                #end 6.2.3.4 Else
              }
            #end while($sub_string =~ /(.{1,$line_size_left})/g)
          }
        #end 6.2 Else
      }
    #end foreach $sub_string (split(/(?=<)|(?<=>)/,$sequence))
    #7. Add the last coodinate with enough leadin white-space to be lined up
    #   with the rest coordinates if the coords_right_flag is true
    $lead_spaces = $max_num_coord_digits - length($coord);
    $line .= ' ' x $line_size_left . $coord_separator . ' ' x $lead_spaces
          . $coord
      if($coords_right_flag && $line_size_left != $chars_per_line);
    $line =~ s/^\s*\d+$coord_separator\s*$// if($coords_left_flag);

    #8. Add the ending PRE tag to the last line of the formatted sequence
    $line =~ s/\n+$/\n/s;

    #9. Add the last line to the formatted_sequence
    $formatted_sequence .= $line;
    #10. Print the last line if the print_flag is true
    print "$line\n" if($print_flag);

    if($coord < 1 && ($coords_left_flag || $coords_right_flag))
      {print("WARNING: The sequence straddles the origin.  Coordinates are ",
             "inaccurate.")}

    #11. Return the formatted_sequence
    return $formatted_sequence;
  }

#Merged the above getNextFastaRec subroutine with the code from convertSeq.pl
#on 2/12/2014
sub getNextFastqRec
  {
    my $handle    = $_[0];      #File handle or file name
    my $no_format = $_[1];

    if(exists($main::{FASTQBUFFER}) && exists($main::{FASTQBUFFER}->{$handle}))
      {
	if(scalar(@{$main::{FASTQBUFFER}->{$handle}}) > 0)
	  {
	    if(wantarray)
	      {
		my @array = (@{$main::{FASTQBUFFER}->{$handle}});
		@{$main::{FASTQBUFFER}->{$handle}} = ();
		return(@array);
	      }
	    return(shift(@{$main::{FASTQBUFFER}->{$handle}}));
	  }
	elsif(eof($handle))
	  {return(undef)}
      }

    my $parent_id_check  = {};
    my $first_loop       = 1;
    my $line_num         = 0;
    my $line             = '';
    my $defline          = '';
    my $seq              = '';
    my $qual             = '';
    my $getting_sequence = 0;
    my $comment_buffer   = '';
    my $seq_lines        = 0;
    my $qual_lines       = 0;
    my $verbose_freq     = 1000;

    #For each line in the current input file
    while(getLine($handle))
      {
	$line_num++;

	verboseOverMe("Reading line [$line_num].")
	  unless($line_num % $verbose_freq);

	$line = $_;

	next if($line !~ /\S/ || ($first_loop && $line =~ /^\s*#/));

	$first_loop = 0;

	#If this is the defline, or the quality length is the same as the seq
	if(length($qual) >= length($seq) && /^\s*\@[^\n\r]*/)
	  {
	    if($defline ne '' || $seq ne '' || $qual ne '')
	      {
		my $solidseq =
		  ($seq_lines == 1 || $no_format ? $seq :
		   formatSequence($seq));
		$qual =~ s/[\s\r\n\t]+//g if(!$no_format && $qual_lines > 1);
		chomp($solidseq);
		chomp($qual);
		chomp($defline);

		push(@{$main::{FASTQBUFFER}->{$handle}},
		     [$defline,$solidseq,$qual]);
	      }
	    $defline    = $line;
	    $seq_lines  = 0;
	    $qual_lines = 0;

	    my $tmp_id = $defline;
	    $tmp_id =~ s/^\s*\@\s*//;
	    $tmp_id =~ s/\s.*//;

	    if($tmp_id eq '')
	      {warning("No Defline ID on line: [$line_num] of current file.")}
	    elsif(exists($parent_id_check->{$tmp_id}))
	      {error("Two sequences found with the same ID on the ",
		     "defline: [$tmp_id] in current fastq file.")}

	    $seq  = '';
	    $qual = '';
	    $getting_sequence = 1;
	  }
	elsif($getting_sequence && /^\s*\+[^\n\r]*/)
	  {$getting_sequence = 0}
	elsif($getting_sequence)
	  {
	    s/\s+//g;
	    $seq_lines++;
	    if(/^[A-Za-z\n\.~]*$/)
	      {$seq .= $_}
	    else
	      {
		error("Expected a sequence character string, but ",
		      "got: [$_].  Appending anyway.");
		$seq .= $_;
	      }
	  }
	elsif($seq =~ /./)
	  {
	    s/\s+//g;
	    $qual_lines++;
	    if(/^[\!-\~]*$/)
	      {$qual .= $_}
	    else
	      {
		error("Expected a quality character string, but ",
		      "got: [$_].  Appending anyway.");
		$qual .= $_;
	      }
	  }
	#else must be a comment, ignore it
      }

    #Handle the last sequence (if there were any sequences)
    if($seq ne '')
      {
	my $solidseq =
	  ($no_format ? $seq :
	   formatSequence($seq));
	$qual =~ s/[\s\r\n\t]+//g unless($no_format);
	chomp($solidseq);
	chomp($defline);
	chomp($qual);

	push(@{$main::{FASTQBUFFER}->{$handle}},[$defline,$solidseq,$qual]);
      }

    #Return the first sequence (if sequence was parsed)
    if(exists($main::{FASTQBUFFER}) && exists($main::{FASTQBUFFER}->{$handle}))
      {
	if(scalar(@{$main::{FASTQBUFFER}->{$handle}}) > 0)
	  {
	    if(wantarray)
	      {
		my @array = (@{$main::{FASTQBUFFER}->{$handle}});
		@{$main::{FASTQBUFFER}->{$handle}} = ();
		return(@array);
	      }
	    return(shift(@{$main::{FASTQBUFFER}->{$handle}}));
	  }
	else
	  {return(undef)}
      }
    else
      {return(undef)}
  }

#Globals used: $default_tmpdir
sub getTmpDir
  {
    if(exists($ENV{TMPDIR}) &&
       -e $ENV{TMPDIR} && -d $ENV{TMPDIR} && -w $ENV{TMPDIR})
      {return($ENV{TMPDIR})}
    elsif(exists($ENV{TMP}) && -e $ENV{TMP} && -d $ENV{TMP} && -w $ENV{TMP})
      {return($ENV{TMP})}
    elsif(exists($ENV{TEMP}) &&
	  -e $ENV{TEMP} && -d $ENV{TEMP} && -w $ENV{TEMP})
      {return($ENV{TEMP})}
    elsif(-e '/tmp' && -d '/tmp' && -w '/tmp')
      {return('/tmp')}
    else
      {return($default_tmpdir)}
  }

#Returns a segments object.  A segment's object has 3 main hash keys:
#ALL, which is a hash keyed on sequence ID, whose values are lists of
#coordinate pairs.  Each coordinate pair will be aligned with the respective
#pair of every other sequence.  Gaps are filled in.  The pair will be undefined
#if it does not contain a gap that another sequence has (i.e. SeqA: 1-10,11-20
#SeqB: 1-12, 20-29.  In this scenario, SeqA doesn't have a gap, but SeqB does).
#TARGET, which is a similar hash, but it only contains the coordinate pairs
#that were read from the file.  It is currently unused, but is intended to be
#used later to allow unequal alignments to consume neighboring gap sequence.
#SEGS, which is an array of each segment containing a hash keyed on the
#sequence ID and containing a hash with keys for segment SIZE and segment
#COORDS of only target segments.  This is used by the segment alignment scoring
#code to determine how good alignments of the desired aligned segments are,
#i.e. did the segments actually align well together?
sub readSegmentFile
  {
    my $file   = $_[0];
    my $aaseqs = $_[1];
    my $aafile = $_[2];

    my $seqid_lookup = getSeqLookup($aaseqs);

    my $segobj = {};

    if(!defined($file))
      {return($segobj)}

    unless(openIn(*SEGS,$file))
      {
	error("Error reading segments file (-d): [$file].");
	return(undef);
      }

    my $raw_coords   = [];
    my $raw_rows     = [];
    my $num_cols     = 0;
    my $line_num     = 0;
    my $errors       = 0;
    my $ucheck       = {};
    while(getLine(*SEGS))
      {
	$line_num++;
	next if(/^\s*$/ || /^\s*#/);
	chomp($_);
	debug("Reading segment line: [$_].");
	$raw_coords = [split(/ *\t */,$_)];
	if($num_cols == 0)
	  {
	    $num_cols = scalar(@$raw_coords);
	    if(scalar(@$raw_coords) < 2)
	      {
		error("No coordinates could be parsed from line [$line_num] ",
		      "of segment file: [$file].");
		$errors++;
	      }
	  }
	if(scalar(@$raw_coords) != $num_cols)
	  {
	    error("Columns inconsistent in segment file: [$file].  Previous ",
		  "rows had [$num_cols] columns, but the row on line ",
		  "[$line_num] as [",scalar(@$raw_coords),"] columns.");
	    $errors++;
	  }
	if(!exists($seqid_lookup->{$raw_coords->[0]}))
	  {
	    error("Sequence ID: [$raw_coords->[0]] on line [$line_num] of ",
		  "segment file [$file] was not found as an ID on a defline ",
		  "in the corresponding sequence file: [$aafile].");
	    $errors++;
	  }
	if(scalar(grep {$raw_coords->[$_] =~ /\D/} (1..$#{$raw_coords})))
	  {
	    error("Non-integer charcters found among coordinate columns on ",
		  "line [$line_num] of segment file [$file].");
	    $errors++;
	  }
	if(scalar(grep {$raw_coords->[$_] !~ /\d/} (1..$#{$raw_coords})))
	  {
	    error("Integer charcters found to be absent among at least 1 ",
		  "coordinate column on line [$line_num] of segment file ",
		  "[$file].");
	    $errors++;
	  }
	if(exists($ucheck->{$raw_coords->[0]}))
	  {
	    error("Duplicate sequence ID found on line [$line_num] of ",
		  "segment file [$file].");
	    $errors++;
	  }
	my $c = 0;
	#If first in a pair is less than or equal to a previous column or the
	#second of a pair is less than a previous column, it's an error
	if(scalar(grep {my $res =
			  ($_ % 2 ? $raw_coords->[$_] <= $c :
			   $raw_coords->[$_] < $c);
			$c = $raw_coords->[$_];$res} (1..$#{$raw_coords})))
	  {
	    error("Each pair of coordinate columns must have coordinates in ",
		  "ascending order.  Coordinates on line [$line_num] in file ",
		  "[$file] are not properly ordered.",
		  {DETAIL => ("The first coordinate of each pair must be " .
			      "greater than the coordinate in the previous " .
			      "column.  The second of each pair must be " .
			      "greater than or equal to the coordinate in " .
			      "the previous column.")});
	    $errors++;
	  }
	#If anything other than the last column has a coordinate that is out of
	#range
	if(exists($seqid_lookup->{$raw_coords->[0]}) &&
	   scalar(grep {$raw_coords->[$_] >
			  length($seqid_lookup->{$raw_coords->[0]}->[1])}
		  (1..($#{$raw_coords} - 1))))
	  {
	    error("Coordinates out of range on line [$line_num] of segment ",
		  "file [$file].",
		  {DETAIL => ("The size of the amino acid " .
			      "[$raw_coords->[0]] is [",
			      length($seqid_lookup->{$raw_coords->[0]}->[1]),
			      "].  Only the coordinate in the last column " .
			      "is allowed to be outside the coordinate " .
			      "range (because it's automatically " .
			      "adjusted).")});
	    $errors++;
	  }

	$ucheck->{$raw_coords->[0]}++;

	push(@$raw_rows,$raw_coords);
      }

    closeIn(*SEGS);

    if($errors)
      {return(undef)}

    debug("Raw rows from segment file [$file]:\n",
	  join("\n",map {join(',',@$_)} @$raw_rows));

    $segobj = insertGapCoords($raw_rows,$seqid_lookup);

    if(scalar(keys(%{$segobj->{ALL}})) != scalar(keys(%$seqid_lookup)))
      {
	error("The number of sequences: [",scalar(keys(%$seqid_lookup)),
	      "] in the protein file: [$aafile] does not match the number of ",
	      "sequence segment rows/IDs [",scalar(keys(%{$segobj->{ALL}})),
	      "] in the segment file: [$file].");
	return(undef);
      }

    return($segobj);
  }

#Takes a 2D array and returns a segments object.
sub insertGapCoords
  {
    my $array = $_[0];
    my $hash  = $_[1];
    my $obj = {};

    my $real_indexes = [];

    foreach my $row (@$array)
      {
	my $id   = $row->[0];
	my $size = length($hash->{$id}->[1]);
	my $c    = 0;
	for(my $i = 1;$i < scalar(@$row);$i += 2)
	  {
	    #Insert a pair of gap coords
	    if($row->[$i] == ($c + 1))
	      {
		push(@{$obj->{ALL}->{$id}},undef);
		#There's a real index value for each pair of coordinates which
		#tracks whether the entire column is all undef values or not
		#(because all the pairs of coordinates in 1 column may abutt
		#all the pairs of coordinates in another column.  Any column
		#without a "real" value will be filtered out at the bottom of
		#this sub.  The 0 below is set whenever the segments object has
		#grown larger than the real_indexes array.
		if(scalar(@$real_indexes) < scalar(@{$obj->{ALL}->{$id}}))
		  {$real_indexes->[$#{$obj->{ALL}->{$id}}] = 0}
	      }
	    else
	      {
		push(@{$obj->{ALL}->{$id}}, [($c + 1),($row->[$i] - 1)]);
		$real_indexes->[$#{$obj->{ALL}->{$id}}] = 1;
	      }

	    #Determine the next/end coordinate in the pair
	    my $end = 0;
	    if(scalar(@$row) > ($i + 1))
	      {$end = $row->[$i + 1]}
	    else
	      {$end = $size}
	    if($end > $size)
	      {$end = $size}

	    #The target pairs of coordinates are the ones the user has
	    #specified in their file that they want to have aligned together
	    push(@{$obj->{TARGET}->{$id}},[$row->[$i],$end]);

	    #Keep track of which sequence has the smaller segment for each
	    #column of segments.  This will be used in the segment alignment
	    #score calculation.  We will need to know which sequence has the
	    #smallest segment for any pair of seqences and each column of
	    #segments and what coordinates it corresponds to.  And so that any
	    #DNA alignment can be scored by domain alignment, we will convert
	    #the coordinates to nucleotide coordinates.  These values are
	    #indexed essentially be the "segment index" (i.e. The current size
	    #of the current TARGET array of coordinate pars is the current
	    #segment number we're on).
	    my $nt_size = ($end - $row->[$i] + 1) * 3;
	    $obj->{SEGS}->[$#{$obj->{TARGET}->{$id}}]->{$id}->{SIZE} =
	      $nt_size;
	    $obj->{SEGS}->[$#{$obj->{TARGET}->{$id}}]->{$id}->{COORDS} =
	      [($row->[$i] - 1) * 3 + 1,$end * 3];

	    #Push the current pair on
	    push(@{$obj->{ALL}->{$id}},[$row->[$i],$end]);
	    $real_indexes->[$#{$obj->{ALL}->{$id}}] = 1;
	    $c = $end;
	  }

	#Append a gap at the end
	if($obj->{ALL}->{$id}->[-1]->[1] == $size)
	  {
	    push(@{$obj->{ALL}->{$id}},undef);
	    if(scalar(@$real_indexes) < scalar(@{$obj->{ALL}->{$id}}))
	      {$real_indexes->[$#{$obj->{ALL}->{$id}}] = 0}
	  }
	else
	  {
	    push(@{$obj->{ALL}->{$id}},[($c + 1),$size]);
		$real_indexes->[$#{$obj->{ALL}->{$id}}] = 1;
	  }
      }

    debug("Segment coordinate pairs before missing gap filtering:\n[",
	  join("\n",
	       map {my $tid = $_;$tid . ':(' .
		      join('),(',
			   map {defined($_) ?
				  join(',',@$_) : 'undef'}
			   @{$obj->{ALL}->{$tid}}) . ')'}
	       keys(%{$obj->{ALL}})),'].');

    #Eliminate non-existent gaps.  If all rows from the segment file had no
    #gaps between pairs of columns, grep out those indexes
    foreach my $id (keys(%{$obj->{ALL}}))
      {$obj->{ALL}->{$id} = [map {$obj->{ALL}->{$id}->[$_]}
			     grep {$real_indexes->[$_]}
			     (0..$#{$obj->{ALL}->{$id}})]}

    debug("Segment coordinate pairs after missing gap filtering:\n[",
	  join("\n",
	       map {my $tid = $_;$tid . ':(' .
		      join('),(',
			   map {defined($_) ?
				  join(',',@$_) : 'undef'}
			   @{$obj->{ALL}->{$tid}}) . ')'}
	       keys(%{$obj->{ALL}})),'].');

    return($obj);
  }

sub getSeqLookup
  {
    my $recs = $_[0];
    my $hash = {};

    if(scalar(@$recs) == 0)
      {return({})}

    foreach my $rec (@$recs)
      {
	my $id = parseIDFromDef($rec->[0]);
	if(exists($hash->{$id}))
	  {
	    error("Sequence IDs are not unique.  ID: [$id] found multiple ",
		  "times.  Skipping duplicates.");
	    next;
	  }
	elsif($id !~ /\S/)
	  {
	    error("Sequence IDs must have non-white-space characters.  ID: ",
		  "[$id] invalid.  Skipping.");
	    next;
	  }
	$hash->{$id} = $rec;
      }

    return($hash);
  }

sub parseIDFromDef
  {
    my $def = $_[0];

    $def =~ s/^[>\@\+]\s*//;
    $def =~ s/^[>\@\+]\s*//;
    $def =~ s/ .*//;

    return($def);
  }
