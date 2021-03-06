#!/usr/bin/perl

#USAGE: Run with no options to get usage or with --help for basic details

use CommandLineInterface;
use warnings;
use strict;
use CodonHomologizer;

##
## Describe the script
##

our $VERSION = '1.036';

setScriptInfo(CREATED => '10/5/2017',
              VERSION => $VERSION,
              AUTHOR  => 'Robert William Leach',
              CONTACT => 'rleach@princeton.edu',
              COMPANY => 'Princeton University',
              LICENSE => 'Copyright 2018',
              HELP    => << 'END_HELP'

This script, given a series of pairwise aligned nucleotide fasta files produced by codonHomologizer.pl, will identify identical homologous subsequences of length N (N=11 is the default) in each submitted alignment pair and select stretches from each alignment (that do not conflict*) and produce 1 sequence for each sequence involved in the submitted alignment pairs that is a hybrid of all the sequences it was paired with, without changing the amino acids encoded.  The goal is to maximize homologous recombination opportunities between all pairs submitted.  This is based on the premise that crossovers occur best when contiguous identity is > 10nts^.

END_HELP
	      ,
	      DETAILED_HELP => << 'END_DETAIL'

This script is written specifically for very divergent proteins with different functions and very little homology.  If your sequences are very similar, you would be better off simply performing a multiple sequence alignment and selecting codons based purely on usage.

It is possible to submit all possible pairs of a set of sequences that have been pairwise aligned with all other sequences, but it will also accept a subset of pairwise aligned sequences, as long as at least 1 sequence is present in more than one alignment file.  For example, if you want to modify 1 sequence to generate a hybrid that has stretches that match in the best spots with a handful of other sequences, you can submit all alignments that include the target sequence and omit all others.  Each sequence that occurs in multiple alignments must be the same length, but they do not have to be the same sequence...

The sequences with the same ID in each alignment must be composed of different codons in order for this script to be of any use.  codonHomologizer, in addition to aligning sequences in a pairwise fashion, selects codons to create the best homology between a pair of sequences to encourage homologous recombinations.  So the same sequence ID in different alignment pairs encode the same protein, but is composed of different codons.  (Note, ad future version of codonHomologizer may allow limited amino acid changes to optimize homology without affecting function.)

Here's how it works.  Each 2-sequence alignment is scanned for homologous stretches.  Each sequence with the same ID (composed of different codons optimized for different partners) is then reconstructed as a mix of codons from its optimized versions with other sequences, containing as many stretches of homology to match as many of the sequences present among the alignments as possible.  It also attempts to spread out the stretches to create homologous stretches in a set of disparate locations.  It outputs a single sequence for each sequence ID involved in an alignment.  Each sequence contains a mix of segments that are optimized to match homologous stretches of its various partners.

* A stretch is defined as a string of contiguous identical bases between 2 or more sequences.  2 stretches conflict if they overlap and the complement of sequences that compose each stretch includes a common sequence and the overlapping portion of any one of those common sequences differs.

^ Numbers of contiguous homologous bases required for crossover events are based on 2 references from a book (chapter 3: Shuffle Optimizer: A Program to Optimize DNA Shuffling for Protein Engineering), claiming that crossovers won't happen without at least 5 bases of contiguous identity and occur best when contiguous identity is > 10nts:

Moore GL, Maranas CD (2000) Modeling DNA mutation and recombination for directed evolution experiments. J Theor Biol 205(3):483-503. doi:10.1006/jtbi.2000.2082

He L, Friedman AM, Bailey-Kellogg C (2012) Algorithms for optimizing cross-overs in DNA shuffling. BMC Bioinformatics 13(Suppl 3):S3. doi:10.1186/1471-2105-13-s3-s3

END_DETAIL
	     );

setDefaults(HEADER        => 0,
	    ERRLIMIT      => 3,
	    COLLISIONMODE => 'error',
	    DEFRUNMODE    => 'usage');

my $seq_file_type =
  addInfileOption(GETOPTKEY   => 'i|aln-files=s',
		  REQUIRED    => 1,
		  PRIMARY     => 1,
		  DEFAULT     => undef,
		  SMRY_DESC   => 'Nucleotide alignment files.',
		  DETAIL_DESC => 'Nucleotide alignment files.',
		  FORMAT_DESC => << 'END_FORMAT'

The nucleotide alignment files must be in aligned fasta format, as output by codonHomologizer.pl.  At least 2 files are required.  Each file must contain exactly 2 sequences.  At least 1 sequence ID must be present in more than 1 file, their unaligned length must be the same, and their codons should differ (otherwise, this script will simply output the same sequence that was input, which would be pointless).  2 alignment files must not contain the same pair of sequence IDs.  The case of the sequences does not matter.

Example file 1:

>hs optimized for mj
---------------TAT---AGTGGTGCT---------GGACCAGCACTTGCTCCA---
CCAGCTCCTCCTCCTCCTATTCAAGGTTAT---------GCTTTTAAACCACCTCCTAGA
---CCTGATTTT---GGTACTTCTGGTAGGACTATT------------AAA---------
TTGCAA---------GCTAAT---TTTTTTGAA---ATGGAT------------ATTCCA
AAAATTGATATTTACCATTAT------------GAATTGGAT---ATTAAACCAGAAAAA
TGTCCAAGACGTGTTAAT
>mj optimized for hs
TTGAATAAAGTTACTTATAAAATTAATGCTTATAAAATTAAAGAAGAATTTATTCCAAAA
GAAGTTCATTTTTATCGTATTAAAAGTTTTGTTAATGAAGCTTTTAAT---TTTTATAGA
TTTGTTAATTTTTATGGT------GGTATGATTATTAATAAAAAAGATAAATCTTTTGTT
TTGCCATATAAAGTTGATAATAAAGTTCTTAAATATAAGGATGGTAATAATGAAATT---
CCAATTGATATTGAA---TATATTAAATCTTTGAAATTGGAGTATGTTAAACCAGAAATA
GCTGAAAAACTTGTT---

Example file 2:

>hs optimized for mp WP_014295921 Marinitoga piezophila Argonaute (MpAgo)
TATTCTGGTGCTGGTCCAGCTTTGGCTCCACCTGCTCCTCCACCACCAATTCAAGGT---
TATGCTTTT---AAACCACCTCCGAGGCCA---GATTTTGGTACA------TCTGGTAGA
------ACTATTAAATTGCAA---GCTAAT---TTTTTTGAGATGGATATTCCAAAGATT
GAT------ATTTATCAT---TATGAACTTGATATTAAACCTGAAAAATGTCCAAGAAGA
GTTAAT
>mp WP_014295921 Marinitoga piezophila Argonaute (MpAgo) optimized for hs
TAT------------CTAAATTTGTATAAAATTGATATTCCAAAAAAAATTAAACGTTTG
TAT---TTTTATAATCCAGATATGGAGCCAAAACTTTTTGCTAGAAATTTGTCT---AGA
GTTAATAATTTTAAATTTCAAGATTCTAATGATCTTGTTTGGATTGAGATTCCAGATATT
GATTTTCAAATTACTCCTAAAAATGTATTTCAGTATAAAGTTGAAAAA------------
------

Note in the examples, which are submitted together in 1 run, each have the sequence ID "hs" (the first string after the ">" up to the first space or new line character).  The hs sequences are different, though they each encode the same protein.  One was optimized for homology to sequence mj and the other for sequence mp.

END_FORMAT
		 );

my $stretch_mins     = [];
my $stretch_mins_def = [5,11];
addArrayOption(GETOPTKEY   => 's|stretch-min=s',
	       GETOPTVAL   => $stretch_mins,
	       REQUIRED    => 0,
	       DEFAULT     => $stretch_mins_def,
	       HIDDEN      => 0,
	       INTERPOLATE => 1,
	       ACCEPTS     => ['>=3'],
	       SMRY_DESC   => 'Length of contiguous identity to search & mix.',
	       DETAIL_DESC => << "END_DETAIL"

This is the number of contiguous identical nucleotides between the sequences in the supplied alignments (see -i) for searching and mixing into the output sequences.  Every identical stretch of this length will be a candidate for weaving into the final output sequence for each input sequence present among the alignments.  Multiple lengths can be supplied.  The algorithm will attempt to insert them in descending size order.  The default values of 11 and 5 are based on the characteristics of recombination.  See --help --extended for more information.

END_DETAIL
	 );

my $codon_file_type =
  addInfileOption(GETOPTKEY   => 'c|codon-file|codon-usage-file=s',
		  REQUIRED    => 0,
		  PRIMARY     => 0,
		  DEFAULT     => undef,
		  SMRY_DESC   => 'Codon usage input file.',
		  DETAIL_DESC => << 'END_DETAIL'

Tab-delimited file containing the single letter amino acid (AA) code, codon, and [optional] usage score (see -n).  This script will use the codons in this file to recode any sequence sections that have changed to include an identity optimize homology after a hybrid sequence is mixed together from multiple alignments.  In other words, if sequence 1 is changed to include a homologous stretch matching sequence 2, that portion of sequence 1 no longer has the optimal homology possible with sequence 3, so sequence 3 is recoded to best match sequence 1's new set of codons.

All AAs must be present.  If usage is unimportant, make the usage scores all the same (e.g. 1.0) or supply -n to ignore the usage column.  Commented out lines in the file will be ignored.  See --help for more details.

All sequences are assumed to be produced in the same organism with a single usage profile.

END_DETAIL
		  ,
		  PAIR_WITH   => $seq_file_type,
		  PAIR_RELAT  => 'ONETOMANY',
		  FORMAT_DESC => << 'END_FORMAT'

A tab- or space- delimited file of codon usage scores.  The first column is the single-character amino acid code, the second column is the codon (3 nucleotides), and the third column us the codon's usage score.  The score can be any integer or fraction.  The scores will be normalized, so the score column does not need to sum to any value.  Additional columns are optional and ignored.  Not all codons are required to be present.  All amino acids and the stop character ('*') are required to be present.  To exclude a codon from incorporation into all recoded sequences, you can either remove that row from the file or comment it out by inserting a '#' at the beginning of the line.  Empty lines are permissible.  Example:

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

my $seq_out_type =
  addOutfileOption(GETOPTKEY     => 'o|outfile|output-file=s',
		   REQUIRED      => 0,
		   PRIMARY       => 1,
		   HIDDEN        => 0,
		   SMRY_DESC     => 'DNA sequence output file.',
		   DETAIL_DESC   => 'DNA sequence output file.',
		   FORMAT_DESC   => 'Nucleotide fasta file.',
		   COLLISIONMODE => 'merge',
		   PAIR_WITH     => $seq_file_type,
		   PAIR_RELAT    => 'ONETOMANY');

debug("Seq out type ID: [$seq_out_type].");

my $map_out_type =
  addOutfileOption(GETOPTKEY   => 'p|map-outfile=s',
		   REQUIRED    => 0,
		   PRIMARY     => 1,
		   HIDDEN      => 0,
		   SMRY_DESC   => 'Sequence map output file.',
		   DETAIL_DESC => 'Sequence map output file.',
		   FORMAT_DESC => << 'END_FORMAT'

Coded fasta file/map.  The purpose of this map is to convey which partner sequence a section of sequence is based on (using the alignments).  The sequences are mixed to try to include stretches of identity from multiple pairwise alignments.  Each sequence ID is assigned an arbitrary single letter code (only 26 sequences supported).  A legend with the character codes and which sequences they represent is printed at the top of the file as a comment.  There will be 1 fasta record for each sequence, consisting of the coded characters from the legend.  The character code at each position of the sequence indicates which sequence the nucleotide in the output sequence file (see --outfile) is based on.  The capitalization of the coded characters indicate whether the sequence in that position came directly from the pairwise alignment with the sequence represented by the coded character or whether the sequence was recoded to best match a version of the sequence indicated by the coded character.  In other words, a lower case coded character means that the sequence which with this portion was aligned was optimized with a third sequence.

END_FORMAT
		   ,
		   COLLISIONMODE => 'merge',
		   PAIR_WITH     => $seq_file_type,
		   PAIR_RELAT    => 'ONETOMANY'
		  );

my $report_out_type =
  addOutfileOption(GETOPTKEY     => 'r|report-outfile=s',
		   REQUIRED      => 0,
		   PRIMARY       => 1,
		   HIDDEN        => 0,
		   SMRY_DESC     => 'Report output file.',
		   DETAIL_DESC   => 'Report output file.',
		   FORMAT_DESC   => 'Free text file.',
		   COLLISIONMODE => 'merge',
		   PAIR_WITH     => $seq_file_type,
		   PAIR_RELAT    => 'ONETOMANY');

my $no_usage = 0;
addOption(GETOPTKEY   => 'n|ignore-codon-usage!',
	  GETOPTVAL   => \$no_usage,
	  REQUIRED    => 0,
	  DEFAULT     => $no_usage,
	  HIDDEN      => 0,
	  DETAIL_DESC => << 'END_DETAIL'

If codon usage is not important, supply this option to ignore the usage values.  Any usage values in the codon usage file (see -c) will be ignored and sequence codon editing will only be based on homology (and homology placement/flexibility).  Note, the hierarchy of factors is: homology, homology flexibility, and usage.

END_DETAIL
	 );

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
			  'that is used to re-code portions of sequence ' .
			  'that overlap an identity segment between 2 other ' .
			  "pairs of sequence.\n\nFor example, if sequences " .
			  'A and B are identical from coordinates 20 to 35 ' .
			  'and 21 to 36 respectively, and sequence A was ' .
			  'also aligned with sequence C, but is not ' .
			  'identical in that region, sequence C will be ' .
			  "recoded in that region to match sequence A's " .
			  'codons that were optimized for sequence B.'),
	  ACCEPTS     => $weighting_methods);

addOutdirOption(GETOPTKEY   => 'dir|outdir=s',
		REQUIRED    => 0,
		DEFAULT     => undef,
		HIDDEN      => 0,
		DETAIL_DESC => 'All outfiles will be put in this directory.');

##
## Process the command line
##

processCommandLine();

#Require at least 2 alignment files
if(scalar(grep {$_ < 2} getFileGroupSizes($seq_file_type)))
  {
    error("Too few alignment files submitted.",
	  {DETAIL => ("2 or more are required.  Refer to the extended usage " .
		      "output for details (i.e. run this script with only " .
		      "the --extended option).")});
    quit(1);
  }

if(scalar(@$stretch_mins) == 0 || scalar(grep {/\D/} @$stretch_mins))
  {
    if(scalar(@$stretch_mins) == 0)
      {@$stretch_mins = @$stretch_mins_def}
    else
      {
	error("Invalid value(s) for --stretch-min: [",
	      join(' ',@$stretch_mins),"].",
	      {DETAIL => ('There must be at least 1 value and all values ' .
			  'must be unsigned integers.')});
	quit(2);
      }
  }
elsif(scalar(grep {$_ < 3} @$stretch_mins))
  {
    error("Value(s) for --stretch-min below allowed minimum: [",
	  join(' ',grep {$_ < 3} @$stretch_mins),"].",
	  {DETAIL => ('The minimum value for a stretch size is 3.')});
    quit(3);
  }

my $codon_data_hash  = {};
my $matrix_data_hash = {};
my $input_data_hash  = {}; #->{$ntOutFile} = {CDNHASH => $cuh, DATA => $hash}

verbose("Loading input data.");

#nextFileSet iterates sets of command-line-supplied files processed together
while(nextFileCombo())
  {
    my $codonUsageFile = getInfile($codon_file_type);
    my $alnFile        = getInfile($seq_file_type);
    my $ntOutFile      = getOutfile($seq_out_type);
    my $mapOutFile     = getOutfile($map_out_type);
    my $reportOutFile  = getOutfile($report_out_type);

    #Read in this usage file if we haven't already
    if(defined($codonUsageFile) &&
       (!exists($codon_data_hash->{$codonUsageFile}) ||
        !exists($matrix_data_hash->{$codonUsageFile}) ||
        !defined($matrix_data_hash->{$codonUsageFile})))
      {
	$codon_data_hash->{$codonUsageFile} = readUsageFile($codonUsageFile,
							    $no_usage);

	verbose("Creating matrix using codon usage.");
	$matrix_data_hash->{$codonUsageFile} =
	  createMatrix($codon_data_hash->{$codonUsageFile});
	if(!defined($matrix_data_hash->{$codonUsageFile}))
	  {
	    error("Unable to create matrix object from codon usage file: ",
		  "[$codonUsageFile].  Skipping.");
	    next;
	  }

	$input_data_hash->{$ntOutFile}->{CDNHASH} =
	  {FWD => $codon_data_hash->{$codonUsageFile},
	   REV => makeReverseCodonHash($codon_data_hash->{$codonUsageFile})};
	$input_data_hash->{$ntOutFile}->{DATA}    = {};
	$input_data_hash->{$ntOutFile}->{MTXHASH} =
	  $matrix_data_hash->{$codonUsageFile};
	$input_data_hash->{$ntOutFile}->{MAPFILE} = $mapOutFile;
	$input_data_hash->{$ntOutFile}->{RPTFILE} = $reportOutFile;
      }

    loadAlignment($alnFile,
		  $stretch_mins,
		  $input_data_hash->{$ntOutFile}->{DATA});
  }

my $max_size = (sort {$b <=> $a} @$stretch_mins)[0];
my($divider_warnings);

foreach my $outfile (sort {$a cmp $b} keys(%$input_data_hash))
  {
    my $solution = {};
    $divider_warnings = {};
    #Add segments to the solution starting with the largest segments and going
    #down in size (to squeeze less optimally sized segments in)
    foreach my $size (sort {$b <=> $a} @$stretch_mins)
      {
	verbose("Doing stretch size: [$size].");

	my $sorted_pair_ids = sortPairIDs($input_data_hash->{$outfile}->{DATA},
					  $size);

	verbose("Searching [",scalar(@$sorted_pair_ids),"] alignments for ",
		"the best assortment of identity stretches.");

	$solution = getMaxWeavedSegments($input_data_hash->{$outfile}->{DATA},
					 $size,
					 $solution,
					 $sorted_pair_ids,
					 $max_size,
					 $input_data_hash->{$outfile}
					 ->{CDNHASH},
					 $stretch_mins,
					 1);
      }

    #Keep the last verboseOverMe message on the screen by printing a "\n"
    if(isVerbose())
      {print STDERR ("\n")}

    my $score = scoreSolution($solution,
			      $input_data_hash->{$outfile}->{DATA},
			      $stretch_mins);

    debug({LEVEL => 2},"Raw unreduced solution with unmerged overlapping ",
	  "identity segments:\n",solutionToString($solution),"\n");

    my($soln_stats,$validation) =
      outputHybrids(generateHybrids($solution,
				    $input_data_hash->{$outfile}->{DATA},
				    $input_data_hash->{$outfile}->{CDNHASH},
				    $input_data_hash->{$outfile}->{MTXHASH},
				    $stretch_mins),
		    $outfile,
		    $input_data_hash->{$outfile}->{MAPFILE});

    outputReport($score,
		 $input_data_hash->{$outfile}->{DATA},
		 $max_size,
		 $soln_stats,
		 $stretch_mins,
		 $validation,
		 $input_data_hash->{$outfile}->{RPTFILE});
  }










##
## Subroutines...
##

#This sub takes the input_data_hash (which I will refer to as the segments
#object) and tries to select the most non-overlapping segments for inclusion in
#the final sequence output.
#Definitions of terms used in the variables:
#  pair - a pair of aligned sequences
#  region - a location in a sequence that represents "optimal" spacing of
#           desired crossover events.  Crossover events happen in regions of
#           sequence identity (see 'segment' below).  There is one number of
#           identical segments of sequence per pair, but the number of regions
#           where they reside may be different per sequence if the segments of
#           identity are not well spaced.  However a pair has 1 number of
#           regions which is set by the sequence with the largest number of
#           populated regions.  For example, imagine an alignment of 2
#           sequences with 2 segments of identity, but there's a large
#           alignment gap between the segments.  In sequence 1, the segments
#           are right next to eachother and in sequence 2, they are far apart.
#           The optimal spacing for 2 segments of identity is 1/3 and 2/3 into
#           each sequence.  Likely, sequence 1 will have 1 populated region
#           where its segments are found while sequence 2 has a segment in each
#           of the 2 regions defined by optimal spacing.
#  segment - An N-long series of identical bases (where N is "size")
#  size - The length of the segments being added.  All segments added are this
#         size.  This method can be called multiple times with different sizes,
#         but the size must exist in the hash.  Each call to this method will
#         try and fit segements of the supplied size in with segments added of
#         the previous size.  Thus, it is recommended that larger segments be
#         added first and smallest last.
#  least or first - A pair or region array index representing the smallest
#                   index containing unprocessed segments
#  start - the pair or region index at which to start processing segments after
#          a recursive call
#  counts - These are the saved locations of where segment processing currently
#           is.  Each pair contains a series of segments which have been sorted
#           into regions.  Thus to maintain state, we need to know, for each
#           pair/region combo, which segments have been processed and which
#           haven't.  Counts maintain that.  But since we alternate pair and
#           region, we also need to know where to pick up those alternations
#           and when we alternate back to the beginning, where the beginning is
#           (excluding those with no more segments left to process).
#  soln - a solution is a collection of segments that have been selected to not
#         conflict with every other contained segment and is ideally optimal.
#         An optimal solution maximizes the number of segments from each pair,
#         maximizes the spacings of those segments in their respective
#         sequences, and includes as many segments and unique bases involved in
#         segments.
#  seqid1 & 2 - The IDs of the sequences in the alignments, sorted in
#               alphanumeric order.
#  pos - The start position of the identity segment (starting from 1)
#  first codon start & last codon end - The first & last base positions of the
#                                       left and right edge codons that the
#                                       segment overlaps.
sub getMaxWeavedSegments
  {
    #Might rewrite later because recursion is a resource hog.
    #http://www.perlmonks.org/?node_id=238938
    no warnings 'recursion'; #This is recursive, by design.

    my $hash         = $_[0];
    my $size         = $_[1];
    my $soln         = $_[2];
    my $pairs        = $_[3];
    my $max_size     = $_[4];
    my $codon_hash   = $_[5];
    my $stretch_mins = $_[6];
    my $greedy       = (defined($_[7])  ? $_[7]  : 1);
    my $pair_least   = (defined($_[8])  ? $_[8]  : 0);           #DO NOT SUPPLY
    my $pair_start   = (defined($_[9])  ? $_[9]  : $pair_least); #BELOW THIS
    my $region_first = (defined($_[10]) ? $_[10]  : [(0) x scalar(@$pairs)]);
    my $region_start = (defined($_[11]) ? $_[11] : [@$region_first]);
    my $counts       = (defined($_[12]) ? $_[12] :
			[map {[(0) x $hash->{$pairs->[$_]}->{NUMUNIQSEGS}
			       ->{$size}]}
			 0..$#{$pairs}]);
    my $depth        = (defined($_[13]) ? $_[13] : 0);
    my $best         = $soln;
    my $best_score   = scoreSolution($soln,$hash,$stretch_mins);
    my $greedy_done  = 0;
    my $progress_msg = '';

    #While the first pair index still containing segments to try and add to the
    #solution is less than the number of pairs
    while($pair_least < scalar(@$pairs))
      {
	foreach my $curpair ($pair_start..$#{$pairs})
	  {
	    my $pairid = $pairs->[$curpair];

	    #Arbitrarily order the sequence IDs
	    my($seqid1,$seqid2) =
	      sort {$a cmp $b} keys(%{$hash->{$pairid}->{SEQS}});

	    my $num_regions =
	      $hash->{$pairid}->{NUMUNIQSEGS}->{$size};

	    while($region_first->[$curpair] < $num_regions)
	      {
		foreach my $curregion ($region_start->[$curpair]..
				       ($num_regions - 1))
		  {
		    my $num_segs = scalar(@{$hash->{$pairid}->{REGS}->{$size}
					      ->[$curregion]});
		    #Current position in the search space = Size of identity,
		    #the current pair (alternates), the current region of the
		    #current pair (alternates), and the current segment of the
		    #current region of the current pair (sequential)
		    my $sol_size = 0;
		    foreach my $sk (keys(%$best))
		      {$sol_size += scalar(grep {$_->{TYPE} eq 'seg'}
					   @{$best->{$sk}})}
		    my $num_seq_keys = scalar(keys(%$best));
		    $progress_msg =
		      join('',("Solution search: Depth:$depth size$size:",
			       "pair$curpair($pair_least-",
			       (scalar(@$pairs) - 1),"):reg",
			       "$curregion($region_first->[$curpair]-",
			       ($num_regions - 1),
			       "):seg$counts->[$curpair]->[$curregion]/",
			       (scalar(@{$hash->{$pairs->[$curpair]}->{REGS}
					   ->{$size}->[$curregion]}) - 1),
			       " SolnScore:",scoreToString($best_score),
			       " NumSolSegs:$sol_size"));
		    verboseOverMe($progress_msg);

		    if($counts->[$curpair]->[$curregion] >= $num_segs)
		      {
			#While the first region index (with segments yet to be
			#evaluated) is less than the number of regions and the
			#segment index we are on for that pair/region is larger
			#than the last segment index for that region
			while($region_first->[$curpair] < $num_regions &&
			      $counts->[$curpair]
			      ->[$region_first->[$curpair]] >
			      $#{$hash->{$pairid}->{REGS}->{$size}
				 ->[$region_first->[$curpair]]})
			  {$region_first->[$curpair]++}
			next;
		      }

		    debug("Region [$curregion] of [$num_regions] for pair ",
			  "[$curpair:$pairid], sequences [$seqid1] and ",
			  "[$seqid2] defined for size [$size]: [",
			  (defined($hash->{$pairid}->{REGS}->{$size}
				   ->[$curregion]) ?
			   'yes' : 'no'),"].\nCounts for this pair/region ",
			  "are [",(defined($counts->[$curpair]->[$curregion]) ?
				   'defined' : 'undef'),"].\nThe number of ",
			  "optimal regions for this alignment at size ",
			  "[$size] are: ",
			  "[$hash->{$pairid}->{NUMUNIQSEGS}->{$size}].",
			  {LEVEL => 10});
		    #If there are still more segments in this region
		    if($counts->[$curpair]->[$curregion] < $num_segs)
		      {
			#The counts keep track of indexes into the segment
			#indexes contained in each region in the REGS hash.
			#They are sorted by each segment's distance from the
			#region center.
			my $reg_seg_index = $counts->[$curpair]->[$curregion];
			#Each segmentUsing the region segment index, we can
			#obtain the actual segment index so we can get the
			#start position of the segment and its spacing score.
			my $index = $hash->{$pairid}->{REGS}->{$size}
			  ->[$curregion]->[$reg_seg_index];

			my $pos1 = $hash->{$pairid}->{SEGS}->{$size}->{$seqid1}
			  ->[$index]->[0];
			my $pos2 = $hash->{$pairid}->{SEGS}->{$size}->{$seqid2}
			  ->[$index]->[0];

			#Spacing scores
			my $spc1 = $hash->{$pairid}->{SEGS}->{$size}->{$seqid1}
			  ->[$index]->[1];
			my $spc2 = $hash->{$pairid}->{SEGS}->{$size}->{$seqid2}
			  ->[$index]->[1];

			#Determine the coordinates after expansion to include
			#every touched codon
			my $first_codon_start1 = $pos1 - (($pos1 - 1) % 3);
			my $last_iden_pos1     = $pos1 + $size - 1;
			my $last_codon_end1    = ($last_iden_pos1 +
						  ($last_iden_pos1 % 3 ? 3 -
						   ($last_iden_pos1 % 3) : 0));
			my $first_codon_start2 = $pos2 - (($pos2 - 1) % 3);
			my $last_iden_pos2     = $pos2 + $size - 1;
			my $last_codon_end2    = ($last_iden_pos2 +
						  ($last_iden_pos2 % 3 ? 3 -
						   ($last_iden_pos2 % 3) : 0));

			my $new_soln = addToSolution($soln,
						     $pairid,
						     $seqid1,$seqid2,
						     $pos1,($pos1 + $size - 1),
						     $pos2,($pos2 + $size - 1),
						     $first_codon_start1,
						     $last_codon_end1,
						     $first_codon_start2,
						     $last_codon_end2,
						     $max_size,
						     $hash,
						     $spc1,$spc2,
						     $codon_hash,
						     $size);

			$counts->[$curpair]->[$curregion]++;

			while($region_first->[$curpair] < $num_regions &&
			      $counts->[$curpair]
			      ->[$region_first->[$curpair]] >
			      $#{$hash->{$pairid}->{REGS}->{$size}
				 ->[$region_first->[$curpair]]})
			  {$region_first->[$curpair]++}

			#If the solution will accept this segment
			if(defined($new_soln))
			  {
			    ##
			    ## Determine where to start in the next recursion
			    ##

			    #Cycle the pair to start looking in the next
			    #alignment for a new segment to add
			    my $next_pair = $curpair + 1;
			    #If we've gotten to the end of the pairs we're
			    #cycling through and we can start over at a place
			    #that's valid
			    if($next_pair >= scalar(@$pairs) &&
			       $pair_least < scalar(@$pairs))
			      {$next_pair = $pair_least}
			    #Keep track of where we left off in this pair so
			    #when we come back to it, we can continue by
			    #searching the next region in the rotation.  This
			    #is separate from which segment to look at in each
			    #region which is tracked in the counts variable.
			    my $next_reg  = [@$region_start];
			    $next_reg->[$curpair] = $curregion + 1;
			    if($next_reg->[$curpair] >= $num_regions)
			      {$next_reg->[$curpair] =
				 $region_first->[$curpair]}

			    #Obtain a candidate solution
			    my($cand);
			    #If we're out of pairs, set the new solution as the
			    #one we were going to send down the recursion (and
			    #don't do the recursion)
			    if($next_pair >= scalar(@$pairs))
			      {$cand = $new_soln}
			    else
			      {
				#We copy the counts that keep track of where we
				#are and the first regions (for each pair)
				#because the positions in the recursive calls
				#will be different from the positions in this
				#call.  This is what allows us to search all
				#possibilities.
				$cand = getMaxWeavedSegments($hash,
							     $size,
							     $new_soln,
							     $pairs,
							     $max_size,
							     $codon_hash,
							     $stretch_mins,
							     $greedy,
							     $pair_least,
							     $next_pair,
							     [@$region_first],
							     $next_reg,
							     [map {[@$_]}
							      @$counts],
							     $depth + 1);
			      }

			    #See if this new candidate solution is better
			    my $cand_score = scoreSolution($cand,
							   $hash,
							   $stretch_mins);
			    if(candidateBetter($cand_score,$best_score))
			      {
				$best_score = $cand_score;
				$best       = $cand;

				$progress_msg =
				  join('',
				       ("Solution search: Depth:$depth size",
					"$size:pair$curpair($pair_least-",
					(scalar(@$pairs) - 1),"):reg",
					"$curregion($region_first->[$curpair]",
					"-",($num_regions - 1),"):seg",
					"$counts->[$curpair]->[$curregion]/",
					(scalar(@{$hash->{$pairs->[$curpair]}
						    ->{REGS}->{$size}
						      ->[$curregion]}) - 1),
					" SolnScore:",
					scoreToString($best_score),
					" NumSolSegs:$sol_size"));
			      }

			    $greedy_done = 1;
			    last if($greedy && $greedy_done);
			  }
		      }
		  }

		#While the segment number we are on in the first region is
		#greater than the segment array's last index (i.e. we are out
		#of segments to try and add from the first region), increment
		#the first region index
		while($region_first->[$curpair] < $num_regions &&
		      $counts->[$curpair]->[$region_first->[$curpair]] >=
		      scalar(@{$hash->{$pairid}->{REGS}->{$size}
				 ->[$region_first->[$curpair]]}))
		  {$region_first->[$curpair]++}

		#This exists to accommodate recursive calls that set
		#$start_region (where the search should continue from) but we
		#want to go back to the first region with segments in it once
		#we reach the end of this set of loop iterations
		$region_start->[$curpair] = $region_first->[$curpair];

		last if($greedy && $greedy_done);
	      }

	    last if($greedy && $greedy_done);
	  }

	#While there are more pairs to start from and the first region index
	#(the one to start from because previous regions are out of segments)
	#is larger than or equal to the number of regions in the starting pair
	#(pair_least), increment pair_least.  Regions are continually
	#alternated through until they have run out of identity segments.  When
	#ever a first region runs out, the region_first value for that pair is
	#incremented.  Once that starting region index goes over the number of
	#regions, it's time to increment the starting pair, because all of the
	#segments from that pair have been added or have tried to have been
	#added to a solution.
	while($pair_least < scalar(@$pairs) &&
	      $region_first->[$pair_least] >=
	      $hash->{$pairs->[$pair_least]}->{NUMUNIQSEGS}->{$size})
	  {$pair_least++}

	#This exists to accommodate recursive calls that set $pair_least (where
	#the search should continue from) but we want to go back to the first
	#pair with segments in it once we reach the end of this set of loop
	#iterations
	$pair_start = $pair_least;

	last if($greedy && $greedy_done);
      }

    verboseOverMe($progress_msg," Returning best solution at depth $depth");

    return($best);
  }

sub sortPairIDs
  {
    my $hash = $_[0];
    my $size = $_[1];

    #Sort the pairs by the ascending number of non-overlapping identical
    #segments it contains (so the ones with the fewest are placed first), then
    #by how poorly spaced the segments are (so the ones with the poorest spaced
    #segments are placed first), then by the number of identical bases (so the
    #ones with the fewest identical bases are placed first), then by the key
    #(alphanumerically) so that ordering can be done consistently
    return([sort {

      $hash->{$a}->{NUMUNIQSEGS}->{$size} <=>
	$hash->{$b}->{NUMUNIQSEGS}->{$size} ||

	  $hash->{$b}->{SPACINGSCORE}->{$size} <=>
	    $hash->{$a}->{SPACINGSCORE}->{$size} ||

	      $hash->{$a}->{NUMBASES}->{$size} <=>
		$hash->{$b}->{NUMBASES}->{$size} ||

		  $a cmp $b}

	    keys(%$hash)]);
  }

#Adds a subsequence from a segment to a solution object while maintaining the
#order of the subsequences by the identity start coordinate.
sub addSegment
  {
    my $soln              = $_[0];
    my $pairid            = $_[1];

    my $seqid             = $_[2];
    my $iden_start        = $_[3];
    my $iden_stop         = $_[4];
    my $first_codon_start = $_[5];
    my $last_codon_end    = $_[6];

    my $type              = $_[7];
    my $alt_codon         = defined($_[8]) ? $_[8] : '';

    my $final_start       = $_[9];
    my $final_stop        = $_[10];

    my $spacing_score     = $_[11];

    my $add_size          = $_[12]; #Search size used to find this segment

    my $insert_index =
      findSegmentInsertPos($soln->{$seqid},$iden_start,$iden_stop);
    if(!defined($insert_index))
      {
	error("Insert of segment into [$seqid] failed.");
	return(1);
      }

    splice(@{$soln->{$seqid}},$insert_index,0,
	   {TYPE              => $type,
	    IDEN_START        => $iden_start,
	    IDEN_STOP         => $iden_stop,
	    FIRST_CODON_START => $first_codon_start,
	    LAST_CODON_STOP   => $last_codon_end,
	    FINAL_START       => $final_start,
	    FINAL_STOP        => $final_stop,
	    PAIR_ID           => $pairid,
	    ALT_CODON         => $alt_codon,
	    SPACING_SCORE     => $spacing_score,
	    ADD_SIZE          => $add_size});

    return(0);
  }

sub solutionToString
  {
    my $soln = $_[0];
    my $str  = '';

    foreach my $seqid (sort {$a cmp $b} keys(%$soln))
      {
	$str .= "$seqid\n";
	my $cnt = 0;
	foreach my $id_segment (@{$soln->{$seqid}})
	  {
	    $cnt++;
	    $str .= "\tSegment $cnt\n";
	    $str .= "\t\t" .
	      join("\n\t\t",
		   map {my $spcr = ' ' x (17 - length($_));
			"$_$spcr => " .
			  (defined($id_segment->{$_}) ? $id_segment->{$_} :
			   'undef')} sort {$a cmp $b}
		   keys(%$id_segment)) . "\n";
	  }
      }

    return($str);
  }

#This returns the index where the record should go.  The ordering is by
#ascending identity start coordinate or descending identity stop coordinate
sub findSegmentInsertPos
  {
    my $seg_array  = $_[0];
    my $iden_start = $_[1];
    my $iden_stop  = $_[2];
    my $index = 0;

    if(!defined($seg_array) || scalar(@$seg_array) == 0)
      {return($index)}

    my $mini = 0;
    my $maxi = $#{$seg_array};
    my $i    = int($maxi / 2);

    #There could be multiple records with the same start (from different
    #pairs), so let's close in on where it is
    while(($maxi - $mini) > 0)
      {
	if($i == $mini && $maxi > $mini)
	  {$i++}
	elsif($i == $maxi && $maxi > $mini)
	  {$i--}

	if(($i == $mini || $i == $maxi))
	  {last}

	#If $i's segment starts to the right of the start coordinate, move the
	#right bound left
	if($iden_start < $seg_array->[$i]->{IDEN_START})
	  {$maxi = $i}
	#Else if $i's segment starts to the left of the start coordinate, move
	#the left bound right
	elsif($iden_start > $seg_array->[$i]->{IDEN_START})
	  {$mini = $i}
	#Else $i's segment starts in the same place as the current segment
	else
	  {
	    while($i > $mini &&
		  $iden_start == $seg_array->[$i-1]->{IDEN_START} &&
		  $iden_stop > $seg_array->[$i]->{IDEN_STOP})
	      {$i--}
	    while($i < $maxi &&
		  $iden_start == $seg_array->[$i+1]->{IDEN_START} &&
		  $iden_stop < $seg_array->[$i]->{IDEN_STOP})
	      {$i++}
	    return($i + 1);
	  }

	$i = int(($maxi - $mini) / 2) + $mini;
      }

    if($iden_start > $seg_array->[$mini]->{IDEN_START} &&
       $iden_start < $seg_array->[$maxi]->{IDEN_START})
      {return($maxi)}
    elsif($iden_start < $seg_array->[$mini]->{IDEN_START})
      {return($mini)}
    elsif($iden_start >= $seg_array->[$maxi]->{IDEN_START})
      {return($maxi + 1)}
    elsif($iden_start == $seg_array->[$mini]->{IDEN_START})
      {return($mini + 1)}
    else
      {
	error("Binary search of a [",scalar(@$seg_array),"] member list ",
	      "failed for identity start coordinate [$iden_start].",
	      {DETAIL => ("Min/Max segment position thus far: " .
			  "[$seg_array->[$mini]->{IDEN_START}/" .
			  "$seg_array->[$maxi]->{IDEN_START}].")});
	return(undef);
      }
  }

sub copySolution
  {
    my $soln = $_[0];
    my $copy = {};

    foreach my $seqid (keys(%$soln))
      {for(my $i = 0;$i <= $#{$soln->{$seqid}};$i++)
	 {$copy->{$seqid}->[$i] = {%{$soln->{$seqid}->[$i]}}}}

    return($copy);
  }

#Adds a pair of subsequences from an identity segment of a pair to a solution
#if it fits.  Returns a copy of the solution with the added segment or undef.
#It also adds alternate codons to the solution if another codon for the same
#amino acid exists to satisfy the identity regions from different segments that
#overlap the same (edge) codon
sub addToSolution
  {
    my $soln               = $_[0];
    my $pairid             = $_[1];
    my $seqid1             = $_[2];
    my $seqid2             = $_[3];
    my $iden_start1        = $_[4];
    my $iden_stop1         = $_[5];
    my $iden_start2        = $_[6];
    my $iden_stop2         = $_[7];
    my $first_codon_start1 = $_[8];
    my $last_codon_stop1   = $_[9];
    my $first_codon_start2 = $_[10];
    my $last_codon_stop2   = $_[11];
    my $max_size           = $_[12]; #Needed to optimize overlap search
    my $hash               = $_[13];
    my $spacing_score1     = $_[14];
    my $spacing_score2     = $_[15];
    my $codon_hash         = $_[16];
    my $add_size           = $_[17]; #Used to know how to score the solution

    #See if the first sequence fits (could generate alt codons on both sides)
    my $fit_data1 = segmentFits($soln,
				$seqid1,
				$first_codon_start1,
				$last_codon_stop1,
				$max_size,
				$iden_start1,
				$iden_stop1,
				$hash,
				$pairid,
			        $codon_hash);
    if(!defined($fit_data1))
      {return(undef)}

    #See if the second sequence fits (could generate alt codons on both sides)
    my $fit_data2 = segmentFits($soln,
				$seqid2,
				$first_codon_start2,
				$last_codon_stop2,
				$max_size,
				$iden_start2,
				$iden_stop2,
				$hash,
				$pairid,
			        $codon_hash);
    if(!defined($fit_data2))
      {return(undef)}

    my $new_soln = copySolution($soln);

    addSegment($new_soln,
	       $pairid,

	       $seqid1,
	       $iden_start1,
	       $iden_stop1,
	       $first_codon_start1,
	       $last_codon_stop1,

	       'seg',
	       undef,

	       $first_codon_start1,
	       $last_codon_stop1,

	       $spacing_score1,
	       $add_size);

    addSegment($new_soln,
	       $pairid,

	       $seqid2,
	       $iden_start2,
	       $iden_stop2,
	       $first_codon_start2,
	       $last_codon_stop2,

	       'seg',
	       undef,

	       $first_codon_start2,
	       $last_codon_stop2,

	       $spacing_score2,
	       $add_size);

    if(defined($fit_data1->{LEFT_CODON}))
      {
	addSegment($new_soln,
		   $pairid,

		   $seqid1,
		   $iden_start1,
		   $first_codon_start1 + 2,
		   $first_codon_start1,
		   $first_codon_start1 + 2,

		   'alt',
		   $fit_data1->{LEFT_CODON},

		   $first_codon_start1,
		   $first_codon_start1 + 2,

		   undef,
		   undef);
      }

    if(defined($fit_data1->{RIGHT_CODON}))
      {
	addSegment($new_soln,
		   $pairid,

		   $seqid1,
		   $last_codon_stop1 - 2,
		   $iden_stop1,
		   $last_codon_stop1 - 2,
		   $last_codon_stop1,

		   'alt',
		   $fit_data1->{RIGHT_CODON},

		   $last_codon_stop1 - 2,
		   $last_codon_stop1,

		   undef,
		   undef);
      }

    if(defined($fit_data2->{LEFT_CODON}))
      {
	addSegment($new_soln,
		   $pairid,

		   $seqid2,
		   $iden_start2,
		   $first_codon_start2 + 2,
		   $first_codon_start2,
		   $first_codon_start2 + 2,

		   'alt',
		   $fit_data2->{LEFT_CODON},

		   $first_codon_start2,
		   $first_codon_start2 + 2,

		   undef,
		   undef);
      }

    if(defined($fit_data2->{RIGHT_CODON}))
      {
	addSegment($new_soln,
		   $pairid,

		   $seqid2,
		   $last_codon_stop2 - 2,
		   $iden_stop2,
		   $last_codon_stop2 - 2,
		   $last_codon_stop2,

		   'alt',
		   $fit_data2->{RIGHT_CODON},

		   $last_codon_stop2 - 2,
		   $last_codon_stop2,

		   undef,
		   undef);
      }

    #I do not need to adjust any "final" coordinates to be used when
    #constructing the sequences.  All I need to do is follow some simple rules.
    #Identity trumps edge codons' non-identical portions, so use the edge
    #codons with larger amounts of identity.  Alternate codons trump
    #everything.

    return($new_soln);
  }

sub codonToHashOfCodons
  {
    my $codon      = $_[0];
    my $codon_hash = $_[1];

    return($codon_hash->{FWD}->{$codon_hash->{REV}->{$codon}});
  }

sub getAA
  {
    my $codon      = $_[0];
    my $codon_hash = $_[1];

    return($codon_hash->{REV}->{$codon});
  }

#Given 2 codons that do not match along with the positions in that codon that
#must not change, find a third alternate codon that meets the criteria of both
#codons (i.e. their fixed positions do not change).  It is assumed that the
#codons submitted do not match but that if they share any fixed relative
#poistions, they do not differ in this positions (Note, the subroutine would
#still return the correct result (no alternative exists), but it does not check
#to make sure that shared fixed positions are the same and instead just checks
#that every other codon alternative meets the fixed position criteria of each
#codon).  This sub requires that both codons encode the same amino acid.  It
#will error out and return undef if they do not.  If there was no error, but no
#alternative codon exists to satisfy the requirements, returns an empty string.
#Globals used: codon hash
sub reCodeMerge
  {
    my $codon_hash      = $_[0];
    my $new_codon       = uc($_[1]);
    my $new_fixed_poses = $_[2];
    my $cur_codon       = uc($_[3]);
    my $cur_fixed_poses = $_[4];
    my $opposite_ends   = defined($_[5]) ? $_[5] : 0;
    #'opposite ends' means that the codons are edge codons of opposite sides.
    #Knowing this allows us to select one of the 2 codons as satisfying both's
    #identity (because one may not, e.g. GAA and GAG - If the left segment's
    #identity ends in GAA and only the G is involved in identity, but the
    #right's segment starts with GAG and only the second G is involved in
    #identity, then the 2 codons are compatible and the alternate is GAG.

    #Do a quick error-check to make sure the codons each encode the same AA
    if(getAA($new_codon,$codon_hash) ne getAA($cur_codon,$codon_hash))
      {
	my $new_aa = getAA($new_codon,$codon_hash);
	my $cur_aa = getAA($cur_codon,$codon_hash);
	error("The codons [$new_codon,$cur_codon] from the same sequence ",
	      "optimized by different alignments do not encode the same ",
	      "amino acid [$new_aa,$cur_aa].  Cannot recode.");
	return(undef);
      }

    my $new_codon_hash = codonToHashOfCodons($new_codon,$codon_hash);

    #Create an array of all possible alternate codons to choose from
    my $alternates = [sort {$a cmp $b}
		      grep {$opposite_ends ?
			      ($_ ne $new_codon || $_ ne $cur_codon) :
			      ($_ ne $new_codon && $_ ne $cur_codon)}
		      keys(%$new_codon_hash)];

    debug({LEVEL => 10},"Alt codons to: $new_codon (@$new_fixed_poses) & ",
	  "$cur_codon (@$cur_fixed_poses): [@$alternates].");

    my $choices = [];
    foreach my $alt (@$alternates)
      {
	my $matches = 1;
	foreach my $new_pos (map {($_ % 3) ? $_ % 3 : 3} @$new_fixed_poses)
	  {
	    if($new_pos > length($new_codon))
	      {
		error("Position: [$new_pos] is larger than the sequence ",
		      "length: [",length($new_codon),"].  Cannot recode.");
		return(undef);
	      }
	    elsif($new_pos > length($alt))
	      {
		error("Position: [$new_pos] is larger than the sequence ",
		      "length: [",length($alt),"].  Cannot recode.");
		return(undef);
	      }
	    if(substr($new_codon,$new_pos - 1,1) ne
	       substr($alt,$new_pos - 1,1))
	      {
		$matches = 0;
		last;
	      }
	  }
	next if(!$matches);
	foreach my $cur_pos (map {($_ % 3) ? $_ % 3 : 3} @$cur_fixed_poses)
	  {
	    if(substr($cur_codon,$cur_pos - 1,1) ne
	       substr($alt,$cur_pos - 1,1))
	      {
		$matches = 0;
		last;
	      }
	  }
	if($matches)
	  {push(@$choices,$alt)}
      }

    if(scalar(@$choices) == 0)
      {return('')}

    #We're no longer concerned with degree of homology for this particular
    #codon.  All we want is to satisfy the homology stretches, so the only
    #deciding factor between multiple matching codons is the usage score
    return((sort {$new_codon_hash->{$b}->{SCORE} <=>
		    $new_codon_hash->{$a}->{SCORE}}
	    @$choices)[0])
  }

sub reCodeOneCodon
  {
    my $fixed_codon  = uc($_[0]);
    my $change_codon = uc($_[1]);
    my $matrix       = $_[2];
    my $codon_hash   = $_[3];

    my $fixed_aa  = getAA($fixed_codon,$codon_hash);
    my $change_aa = getAA($change_codon,$codon_hash);

    my($lesser_aa,$greater_aa,$change_codon_index,$fixed_codon_index);
    if($change_aa le $fixed_aa)
      {
	$change_codon_index = 0;
	$fixed_codon_index  = 1;
	$lesser_aa          = $change_aa;
	$greater_aa         = $fixed_aa;
      }
    else
      {
	$fixed_codon_index  = 0;
	$change_codon_index = 1;
	$lesser_aa          = $fixed_aa;
	$greater_aa         = $change_aa;
      }

    my $flex_sizes = {};
    foreach my $flex_code (keys(%{$matrix->{$lesser_aa}->{$greater_aa}
				    ->{PAIRS}}))
      {$flex_sizes->{length($flex_code)} = 0}

    my $best_change_codon = '';
    my($best_score);
    foreach my $flex_size (sort {$b <=> $a} keys(%$flex_sizes))
      {
	my $flex_codes =
	  [sort {$a cmp $b}
	   grep {length($_) == $flex_size}
	   keys(%{$matrix->{$lesser_aa}->{$greater_aa}->{PAIRS}})];

	foreach my $flex_code (@$flex_codes)
	  {
	    foreach my $codon_pair (sort {$b->[2] <=> $a->[2]}
				    grep {$_->[$fixed_codon_index] eq
					    $fixed_codon}
				    @{$matrix->{$lesser_aa}->{$greater_aa}
					->{PAIRS}->{$flex_code}})
	      {
		if(!defined($best_score) || $codon_pair->[2] > $best_score)
		  {
		    $best_score = $codon_pair->[2];
		    $best_change_codon = $codon_pair->[$change_codon_index];
		  }
		last;
	      }
	  }
	if($best_change_codon ne '')
	  {last}
      }

    return($best_change_codon);
  }

sub reCodeOneSubseq
  {
    my $fixed_partner_seq = $_[0];
    my $fixed_start       = $_[1];
    my $change_this_seq   = $_[2];
    my $change_start      = $_[3];
    my $matrix            = $_[4];
    my $codon_hash        = $_[5];

    #Make sure the sequences both start in the same frame
    if(!framesSynced(@_))
      {
	error("Sequnce alignment does not contain aligned frames.  Unable to ",
	      "recode.");
	return();
      }

    my $fixed_pos  = $fixed_start  - 1;
    my $change_pos = $change_start - 1;
    my $aln_pos    = 0;
    my $min_len    = (length($fixed_partner_seq) < length($change_this_seq) ?
		      length($fixed_partner_seq) : length($change_this_seq));
    my $cur_fixed_codon  = '';
    my $cur_change_codon = '';
    my $new_change_seq = '';

    while($aln_pos < $min_len)
      {
	my $fixed_char  = substr($fixed_partner_seq,$aln_pos,1);
	my $change_char = substr($change_this_seq,  $aln_pos,1);

	if($fixed_char ne '-')
	  {$fixed_pos++}
	if($change_char ne '-')
	  {$change_pos++}

	my $fixed_frame  = $fixed_pos  % 3;
	my $change_frame = $change_pos % 3;

	if($fixed_char ne '-' && $fixed_frame == 1)
	  {$cur_fixed_codon = $fixed_char}
	elsif($fixed_char ne '-')
	  {$cur_fixed_codon .= $fixed_char}

	if($change_char ne '-' && $change_frame == 1)
	  {$cur_change_codon = $change_char}
	elsif($change_char ne '-')
	  {$cur_change_codon .= $change_char}

	if($fixed_char ne '-' && $change_char ne '-' &&
	   $fixed_frame == 3 && $change_frame == 3 &&
	   length($cur_fixed_codon) == 3 && length($cur_change_codon) == 3)
	  {
	    my $tmp_change_codon = reCodeOneCodon($cur_fixed_codon,
						  $cur_change_codon,
						  $matrix,
						  $codon_hash);
	    if($tmp_change_codon eq '')
	      {
		warning("Unable to recode codon [$cur_change_codon].",
		        {DETAIL =>
			 join('',("This should not have happened.  There ",
				  "should always be a best codon match, even ",
				  "if the result is the same codon as you ",
				  "started with.  Make sure your codon usage ",
				  "file has a codon for every amino acid."))});
		$new_change_seq .= $cur_change_codon;
	      }
	    else
	      {$new_change_seq .= $tmp_change_codon}
	  }
	elsif(length($cur_fixed_codon) != length($cur_change_codon) ||
	      length($cur_change_codon) >= 3 ||
	      length($cur_fixed_codon) != $fixed_frame ||
	      length($cur_change_codon) != $change_frame)
	  {
	    $cur_fixed_codon = '';
	    if(length($cur_change_codon))
	      {$new_change_seq .= $cur_change_codon}
	    $cur_change_codon = '';
	  }

	$aln_pos++;
      }

    return($new_change_seq);
  }

#Looks in the aligned sequences for the first position that contains no gap and
#Determines whether the frames are aligned in that position for each of the
#sequences
sub framesSynced
  {
    my $fixed_partner_seq = $_[0];
    my $fixed_start       = $_[1];
    my $change_this_seq   = $_[2];
    my $change_start      = $_[3];

    if(length($fixed_partner_seq) != length($change_this_seq))
      {error("Aligned sequence lengths are not the same.",
	     {DETAIL => ("$fixed_start: $fixed_partner_seq\n" .
			 "$change_start: $change_this_seq")})}

    my $fixed_pos  = $fixed_start  - 1;
    my $change_pos = $change_start - 1;
    my $aln_pos    = 0;
    my $min_len    = (length($fixed_partner_seq) < length($change_this_seq) ?
		      length($fixed_partner_seq) : length($change_this_seq));

    while($aln_pos < $min_len)
      {
	my $fixed_char  = substr($fixed_partner_seq,$aln_pos,1);
	my $change_char = substr($change_this_seq,  $aln_pos,1);

	if($fixed_char ne '-')
	  {$fixed_pos++}
	if($change_char ne '-')
	  {$change_pos++}

	if($fixed_char ne '-' && $change_char ne '-')
	  {
	    my $fixed_frame = $fixed_pos % 3;
	    my $change_frame = $change_pos % 3;

	    return($fixed_frame == $change_frame);
	  }

	$aln_pos++;
      }

    #If we got here, return true if the lengths are the same, false otherwise
    return(length($fixed_partner_seq) == length($change_this_seq));
  }

sub getOverlapSegList
  {
    my $soln     = $_[0];
    my $seqid    = $_[1];
    my $start    = $_[2];
    my $stop     = $_[3]; #Assumed to be larger than $start
    my $max_size = $_[4] + 4;
    #We're adding 4 because the size is just the identity portion,
    #which can be anything, but overlapped edge codons can expand the
    #outer coords by at most 2 on each side.  This is a worst-case
    #scenario.  The max could be smaller, but we'll just keep it
    #simple.

    my $list = [];

    if(scalar(keys(%$soln)) == 0 || !exists($soln->{$seqid}))
      {return($list)}

    my $mini = 0;
    my $maxi = $#{$soln->{$seqid}};
    my $i    = int($maxi / 2);

    #I need to find the first record that overlaps.  All the records are
    #(/should be) sorted by their identity start coordinate.  I will not assume
    #that multiple records cannot have the same identity start, even though
    #that's how I will code this, just to be on the safe side.
    while(($maxi - $mini) > 0)
      {
	#If $i's segment starts to the left of the stop coordinate, move the
	#right bound left
	if($soln->{$seqid}->[$i]->{FIRST_CODON_START} > $stop)
	  {$maxi = $i}
	elsif(($soln->{$seqid}->[$i]->{FIRST_CODON_START} + $max_size - 1) <
	      $start)
	  {$mini = $i}
	else
	  {last}

	$i = int(($maxi - $mini) / 2) + $mini;

	if($i == $mini || $i == $maxi)
	  {last}
      }

    $i = $mini;

    #Now $i is where we will start.  We will move right through the previously
    #added segments until the start + max size is greater than the stop
    #coordinate
    do
      {
	#If there's overlap, add the inner hash to the list
	if(coordsOverlap($start,$stop,
			 $soln->{$seqid}->[$i]->{FIRST_CODON_START},
			 $soln->{$seqid}->[$i]->{LAST_CODON_STOP}))
	  {push(@$list,$soln->{$seqid}->[$i])}

	$i++;
      }
	while($i <= $#{$soln->{$seqid}} &&
	      $soln->{$seqid}->[$i]->{FIRST_CODON_START} <= $stop);

    return(wantarray ? @$list : $list);
  }

sub segmentFits
  {
    my $soln              = $_[0];
    my $seqid             = $_[1];
    my $first_codon_start = $_[2];
    my $last_codon_stop   = $_[3];
    my $max_size          = $_[4];
    my $iden_start        = $_[5];
    my $iden_stop         = $_[6];
    my $hash              = $_[7];
    my $pairid            = $_[8];
    my $codon_hash        = $_[9];

    #Variables to hold alternative codon sequences, named right and left
    #relative to the new segmet side the codon is on.
    my($alt_left_cdn,$alt_right_cdn);

    #These track whether the alternate codon has already been added earlier
    my $no_left_cdn  = 1;
    my $no_right_cdn = 1;

    #This returns a series of references to inner hashes from the solution
    my $seqlist = getOverlapSegList($soln,
				    $seqid,
				    $first_codon_start,
				    $last_codon_stop,
				    $max_size);

    if(isDebug() == 1)
      {debug("Checking overlapping $seqid record starts: [",
	     join(',',map {$_->{FIRST_CODON_START}} @$seqlist),"] out of [",
	     (exists($soln->{$seqid}) ? scalar(@{$soln->{$seqid}}) : '0'),
	     "] records to see if [$seqid:$first_codon_start] fits")}

    foreach my $rec (@$seqlist)
      {
	if($rec->{TYPE} eq 'seg')
	  {
	    next if($rec->{PAIR_ID} eq $pairid);

	    my $common_iden_start = 0;
	    if($iden_start >= $rec->{IDEN_START} &&
	       $iden_start <= $rec->{IDEN_STOP})
	      {$common_iden_start = $iden_start}
	    elsif(($iden_stop >= $rec->{IDEN_START} &&
		   $iden_stop <= $rec->{IDEN_STOP}) ||
		  ($iden_start < $rec->{IDEN_START} &&
		   $iden_stop  > $rec->{IDEN_STOP}))
	      {$common_iden_start = $rec->{IDEN_START}}

	    my $common_iden_stop  = 0;
	    if($iden_stop >= $rec->{IDEN_START} &&
	       $iden_stop <= $rec->{IDEN_STOP})
	      {$common_iden_stop = $iden_stop}
	    elsif(($iden_start >= $rec->{IDEN_START} &&
		   $iden_start <= $rec->{IDEN_STOP}) ||
		  ($iden_start < $rec->{IDEN_START} &&
		   $iden_stop  > $rec->{IDEN_STOP}))
	      {$common_iden_stop = $rec->{IDEN_STOP}}

	    #If the identities overlap, the sequence from each pair must be
	    #identical
	    if($common_iden_start > 0 && $common_iden_stop > 0)
	      {
		#Obtain each sequence segment
		my $subseq1 =
		  substr($hash->{$pairid}->{SEQS}->{$seqid},
			 $common_iden_start - 1,
			 $common_iden_stop - $common_iden_start + 1);

		my $subseq2 =
		  substr($hash->{$rec->{PAIR_ID}}->{SEQS}->{$seqid},
			 $common_iden_start - 1,
			 $common_iden_stop - $common_iden_start + 1);

		if($subseq1 ne $subseq2)
		  {return(undef)}
	      }

	    #Things we do not care about: 1. no overlap at all.  2. overlap
	    #of an edge codon with an identity region.  As long as
	    #overlapping identity regions match and we're comparing the
	    #sequences composed of just different codons that each encode
	    #the same AA sequence.  All we have to do is use sequence from
	    #one pair or another in whole codon increments

	    #If both left codons are not 100% covered by the identity
	    #region and they overlap (assuming the identity is at least 3
	    #bases long) and their identities do not start in the same
	    #place.  We're not going to worry about whether these records
	    #have alternate codons that were previously selected because
	    #either the same conclusion will be reached and the overlapping
	    #alternate codon record will agree and it'll be added or it
	    #won't.
	    if(#both left codons are not 100% covered by the identity reg
	       $iden_start != $first_codon_start &&
	       $rec->{IDEN_START} != $rec->{FIRST_CODON_START} &&

	       #the codons overlap (i.e. are in the same [frame] location)
	       $first_codon_start == $rec->{FIRST_CODON_START} &&

	       #their identities do not stop in the same place
	       $rec->{IDEN_START} != $iden_start)
	      {
		#The codon in the sequence being added
		my $left_codon_new =
		  substr($hash->{$pairid}->{SEQS}->{$seqid},
			 $first_codon_start - 1,
			 3);
		#The codon in the sequence that was already added
		my $left_codon_cur =
		  substr($hash->{$rec->{PAIR_ID}}->{SEQS}->{$seqid},
			 $rec->{FIRST_CODON_START} - 1,
			 3);

		if($left_codon_new ne $left_codon_cur)
		  {
		    #The relative positions in the codon being added that
		    #must stay the same
		    my $fixed_poses_new =
		      [$iden_start..($first_codon_start + 2)];

		    #The relative positions in the codon already added that
		    #must stay the same
		    my $fixed_poses_cur =
		      [$rec->{IDEN_START}..
		       ($rec->{FIRST_CODON_START} + 2)];

		    my $tmp_codon =
		      reCodeMerge($codon_hash,
				  $left_codon_new,$fixed_poses_new,
				  $left_codon_cur,$fixed_poses_cur,
				  ($iden_start != $first_codon_start &&
				   $rec->{IDEN_STOP} !=
				   $rec->{LAST_CODON_STOP}));

		    #If we couldn't select an alternative codon that
		    #satisfies both sequences
		    if(!defined($tmp_codon) || $tmp_codon eq '')
		      {
			if(!defined($tmp_codon))
			  {error("Could not recode codons for sequence ",
				 "[$seqid] at (sequence-relative, i.e. ",
				 "not alignment) positions: ",
				 "[$first_codon_start,",
				 "$rec->{FIRST_CODON_START}] (which ",
				 "should be the same position) in ",
				 "alignments [$pairid,$rec->{PAIR_ID}].")}
			return(undef);
		      }
		    #Else if an alternative codon was already selected to
		    #match another segment and it differs from what we got
		    #here
		    elsif(defined($alt_left_cdn) &&
			  $alt_left_cdn ne $tmp_codon)
		      {return(undef)}
		    #Else if this is a new alternate codon
		    elsif(!defined($alt_left_cdn))
		      {$alt_left_cdn = $tmp_codon}
		  }
	      }

	    #If both right codons are not 100% covered by the identity
	    #region and they overlap (assuming the identity is at least 3
	    #bases long) and their identities do not stop in the same place
	    #We're not going to worry about whether these records have
	    #alternate codons that were previously selected because either
	    #the same conclusion will be reached and the overlapping
	    #alternate codon record will agree and it'll be added or it
	    #won't.
	    if(#both right codons are not 100% covered by the identity reg
	       $iden_stop != $last_codon_stop &&
	       $rec->{IDEN_STOP} != $rec->{LAST_CODON_STOP} &&

	       #the codons overlap (i.e. are in the same [frame] location)
	       $last_codon_stop == $rec->{LAST_CODON_STOP} &&

	       #their identities do not stop in the same place
	       $rec->{IDEN_STOP} != $iden_stop)
	      {
		#The codon in the sequence being added
		my $right_codon_new =
		  substr($hash->{$pairid}->{SEQS}->{$seqid},
			 #-2 to get the first base in the codon and -1 so
			 #the coordinate system starts with 0
			 $last_codon_stop - 2 - 1,
			 3);
		#The codon in the sequence that was already added
		my $right_codon_cur =
		  substr($hash->{$rec->{PAIR_ID}}->{SEQS}->{$seqid},
			 #-2 to get the first base in the codon and -1 so
			 #the coordinate system starts with 0
			 $rec->{LAST_CODON_STOP} - 2 - 1,
			 3);

		if($right_codon_new ne $right_codon_cur)
		  {
		    #The relative positions in the codon being added that
		    #must stay the same
		    my $fixed_poses_new =
		      [($last_codon_stop - 2)..$iden_stop];

		    #The relative positions in the codon already added that
		    #must stay the same
		    my $fixed_poses_cur =
		      [($rec->{LAST_CODON_STOP} - 2)..
		       $rec->{IDEN_STOP}];

		    my $tmp_codon =
		      reCodeMerge($codon_hash,
				  $right_codon_new,$fixed_poses_new,
				  $right_codon_cur,$fixed_poses_cur,
				  ($iden_start != $first_codon_start &&
				   $rec->{IDEN_STOP} !=
				   $rec->{LAST_CODON_STOP}));

		    #If we couldn't select an alternative codon that
		    #satisfies both sequences
		    if(!defined($tmp_codon) || $tmp_codon eq '')
		      {
			if(!defined($tmp_codon))
			  {error("Could not recode codons for sequence ",
				 "[$seqid] at (sequence-relative, i.e. ",
				 "not alignment) positions: ",
				 "[",($last_codon_stop - 2),",",
				 ($rec->{LAST_CODON_STOP} - 2),"] (which ",
				 "should be the same position) in ",
				 "alignments [$pairid,$rec->{PAIR_ID}].")}
			return(undef);
		      }
		    #Else if an alternative codon was already selected to
		    #match another segment and it differs from what we got
		    #here
		    elsif(defined($alt_right_cdn) &&
			  $alt_right_cdn ne $tmp_codon)
		      {return(undef)}
		    #Else if this is a new alternate codon
		    elsif(!defined($alt_right_cdn))
		      {$alt_right_cdn = $tmp_codon}
		  }
	      }

	    #If the left codon of the candidate segment being added and the
	    #right codon of the current segment already in the solution are
	    #not 100% covered by the identity region and they overlap
	    #(assuming the identity is at least 3 bases long) and we do not
	    #need to check if their identity boundaries are the same
	    #because they come in from different sides
	    #We're not going to worry about whether these records have
	    #alternate codons that were previously selected because either
	    #the same conclusion will be reached and the overlapping
	    #alternate codon record will agree and it'll be added or it
	    #won't.
	    if(#both codons are not 100% covered by the identity regs
	       $iden_start != $first_codon_start &&
	       $rec->{IDEN_STOP} != $rec->{LAST_CODON_STOP} &&

	       #the new segment's left edge codon and the already added
	       #segment's right edge codon overlap (assumed to be in the same
	       #frame)
	       $first_codon_start == ($rec->{LAST_CODON_STOP} - 2))
	      {
		#The codon in the sequence being added
		my $left_codon_new =
		  substr($hash->{$pairid}->{SEQS}->{$seqid},
			 $first_codon_start - 1,
			 3);
		#The codon in the sequence that was already added
		my $right_codon_cur =
		  substr($hash->{$rec->{PAIR_ID}}->{SEQS}->{$seqid},
			 #-2 to get the first base in the codon and -1 so
			 #the coordinate system starts with 0
			 $rec->{LAST_CODON_STOP} - 2 - 1,
			 3);

		if($left_codon_new ne $right_codon_cur)
		  {
		    #The relative positions in the codon being added that
		    #must stay the same
		    my $fixed_poses_new =
		      [$iden_start..($first_codon_start + 2)];

		    #The relative positions in the codon already added that
		    #must stay the same
		    my $fixed_poses_cur =
		      [($rec->{LAST_CODON_STOP} - 2)..
		       $rec->{IDEN_STOP}];

		    my $tmp_codon =
		      reCodeMerge($codon_hash,
				  $left_codon_new, $fixed_poses_new,
				  $right_codon_cur,$fixed_poses_cur,
				  1);

		    if(isDebug() == 1)
		      {debug("Alternate codon selected 1: [",
			     (defined($tmp_codon) ? $tmp_codon : 'undef'),
			     "] for position [$seqid:$first_codon_start].")}

		    #If we couldn't select an alternative codon that
		    #satisfies both sequences
		    if(!defined($tmp_codon) || $tmp_codon eq '')
		      {
			if(!defined($tmp_codon))
			  {error("Could not recode codons for sequence ",
				 "[$seqid] at (sequence-relative, i.e. ",
				 "not alignment) positions: ",
				 "[$first_codon_start,",
				 ($rec->{LAST_CODON_STOP} - 2),"] (which ",
				 "should be the same position) in ",
				 "alignments [$pairid,$rec->{PAIR_ID}].")}
			return(undef);
		      }
		    #Else if an alternative codon was already selected to
		    #match another segment and it differs from what we got
		    #here
		    elsif(defined($alt_left_cdn) &&
			  $alt_left_cdn ne $tmp_codon)
		      {return(undef)}
		    #Else if this is a new alternate codon
		    elsif(!defined($alt_left_cdn))
		      {$alt_left_cdn = $tmp_codon}
		  }
	      }

	    #If the right codon of the candidate segment being added and
	    #the left codon of the current segment already in the solution
	    #are not 100% covered by the identity region and they overlap
	    #(assuming the identity is at least 3 bases long) and we do not
	    #need to check if their identity boundaries are the same
	    #because they come in from different sides
	    #We're not going to worry about whether these records have
	    #alternate codons that were previously selected because either
	    #the same conclusion will be reached and the overlapping
	    #alternate codon record will agree and it'll be added or it
	    #won't.
	    if(#both codons are not 100% covered by the identity regs
	       $iden_stop != $last_codon_stop &&
	       $rec->{IDEN_START} != $rec->{FIRST_CODON_START} &&

	       #the new segment's right edge codon and the already added
	       #segment's left edge codon overlap (assumed to be in the same
	       #frame)
	       ($last_codon_stop - 2) == $rec->{FIRST_CODON_START})
	      {
		#The codon in the sequence being added
		my $right_codon_new =
		  substr($hash->{$pairid}->{SEQS}->{$seqid},
			 #-2 to get the first base in the codon and -1 so
			 #the coordinate system starts with 0
			 $last_codon_stop - 2 - 1,
			 3);
		#The codon in the sequence that was already added
		my $left_codon_cur =
		  substr($hash->{$rec->{PAIR_ID}}->{SEQS}->{$seqid},
			 $rec->{FIRST_CODON_START} - 1,
			 3);

		if($right_codon_new ne $left_codon_cur)
		  {
		    #The relative positions in the codon being added that
		    #must stay the same
		    my $fixed_poses_new =
		      [$iden_start..($first_codon_start + 2)];

		    #The relative positions in the codon already added that
		    #must stay the same
		    my $fixed_poses_cur =
		      [($rec->{LAST_CODON_STOP} - 2)..
		       $rec->{IDEN_STOP}];

		    my $tmp_codon =
		      reCodeMerge($codon_hash,
				  $right_codon_new,$fixed_poses_new,
				  $left_codon_cur, $fixed_poses_cur,
				  1);

		    if(isDebug() == 1)
		      {debug("Alternate codon selected 2: [",
			     (defined($tmp_codon) ? $tmp_codon : 'undef'),
			     "] for position [$seqid:$first_codon_start].")}

		    #If we couldn't select an alternative codon that
		    #satisfies both sequences
		    if(!defined($tmp_codon) || $tmp_codon eq '')
		      {
			if(!defined($tmp_codon))
			  {error("Could not recode codons for sequence ",
				 "[$seqid] at (sequence-relative, i.e. ",
				 "not alignment) positions: ",
				 "[",($last_codon_stop - 2),",",
				 "$rec->{FIRST_CODON_START}] (which ",
				 "should be the same position) in ",
				 "alignments [$pairid,$rec->{PAIR_ID}].")}
			return(undef);
		      }
		    #Else if an alternative codon was already selected to
		    #match another segment and it differs from what we got
		    #here
		    elsif(defined($alt_right_cdn) &&
			  $alt_right_cdn ne $tmp_codon)
		      {return(undef)}
		    #Else if this is a new alternate codon
		    elsif(!defined($alt_right_cdn))
		      {$alt_right_cdn = $tmp_codon}
		  }
	      }

	    #At this point, if we've gotten here, the segment can still be
	    #added to the solution (as far as up to this number of
	    #overlapping segments looked at thus far goes).  We will keep
	    #looping to check the other overlaps.
	  }
	else
	  {
	    #This doesn't treat edge codons completely covered by the
	    #identity region any differently because the distinction is
	    #unimportant - what is done is the same test in each case

	    #If the already added alternate codon overlaps the candidate's
	    #left edge codon
	    if($rec->{FIRST_CODON_START} == $first_codon_start)
	      {
		#If the candidate already has an alternate codon & differs
		if(defined($alt_left_cdn) &&
		   $alt_left_cdn ne $rec->{ALT_CODON})
		  {return(undef)}
		#If this alternate codon already exists in the solution as
		#an alternate codon
		elsif(defined($alt_left_cdn))
		  {$no_left_cdn = 0}
		else
		  {
		    #Make sure the identity portion of the sequence being
		    #added is the same as the corresponding portion of the
		    #alternate codon that was already added

		    #Obtain the subseq for this segment from the identity
		    #start to the end of the left edge codon
		    my $subseq1 =
		      substr($hash->{$pairid}->{SEQS}->{$seqid},
			     $iden_start - 1,
			     #+2 to get the last codon base coord +1 to get
			     #the size after getting the diff
			     $first_codon_start + 2 - $iden_start + 1);

		    #Obtain the subseq of the alternate codon from the
		    #relative identity start to the end of the codon
		    #Don't need to add 1 because start = 0 in substr()
		    my $relstart = $first_codon_start - $iden_start;
		    my $relsize  = 3 - $relstart;
		    my $subseq2 =
		      substr($rec->{ALT_CODON},$relstart,$relsize);

		    if($subseq1 ne $subseq2)
		      {return(undef)}

		    #There is an existing compatible alt codon, so no need to
		    #create another
		    $no_left_cdn = 0;
		  }
	      }

	    #If the already added alternate codon overlaps the candidate's
	    #right edge codon
	    if($rec->{LAST_CODON_STOP} == $last_codon_stop)
	      {
		#If the candidate already has an alternate codon & differs
		if(defined($alt_right_cdn) &&
		   $alt_right_cdn ne $rec->{ALT_CODON})
		  {return(undef)}
		#If this alternate codon already exists in the solution as
		#an alternate codon
		elsif(defined($alt_right_cdn))
		  {$no_right_cdn = 0}
		else
		  {
		    #Make sure the identity portion of the sequence being
		    #added is the same as the corresponding portion of the
		    #alternate codon that was already added

		    #-2 to get the last codon first base coord, +1
		    #to get the size after getting the diff
		    my $relidsize =
		      $iden_stop - ($last_codon_stop - 2) + 1;

		    #Obtain the subseq for this segment from the identity
		    #start to the end of the left edge codon
		    my $subseq1 =
		      substr($hash->{$pairid}->{SEQS}->{$seqid},
			     $last_codon_stop - 2 - 1,
			     $relidsize);

		    #Obtain the subseq of the alternate codon from the
		    #codon start to the relative identity stop
		    my $subseq2 =
		      substr($rec->{ALT_CODON},0,$relidsize);

		    if($subseq1 ne $subseq2)
		      {return(undef)}

		    #There is an existing compatible alt codon, so no need to
		    #create another
		    $no_right_cdn = 0;
		  }
	      }

	    #If the already added alternate codon overlaps the candidate's
	    #identity region body (not including the edge codons' identity),
	    #Make sure that that portion of the segment being added matches any
	    #pre-existing alternate codons
	    if($rec->{IDEN_START} >= ($first_codon_start + 3) &&
	       $rec->{IDEN_START} <= ($last_codon_stop - 3))
	      {
		#If the candidate has any alternate codons for its edges,
		#they do not matter in this case because they cannot
		#overlap.

		my $subseq =
		  substr($hash->{$pairid}->{SEQS}->{$seqid},
			 $rec->{FINAL_START} - 1,
			 3);

		#If the middle of the segment doesn't match an overlapping alt
		#codon, return undef
		if($subseq ne $rec->{ALT_CODON})
		  {return(undef)}
	      }
	  }
      }

    #If no alternate codons already in the solution, return them
    return({LEFT_CODON  => ($no_left_cdn  ? $alt_left_cdn  : undef),
	    RIGHT_CODON => ($no_right_cdn ? $alt_right_cdn : undef)});
  }

#Creates an array of hierarchical score values in the following order:
#1 Follow this procedure to get a "weakest link" inclusion score - maximize
#  where S = largest stretch size...
#  Sum the number of unique base positions included in stretches of size S from
#  each pair/sequence combo and int(divide that by S).  If the result is == to
#  the number of unique stretches possible, skip it.  Otherwise, keep track of
#  which pair/sequence's result is the smallest among those whose value is not
#  maxed out (i.e. all possible unique pairs from the alignment have been
#  included).  That will be the score at this level.
#  Note that each pair should have a single score, so I either need to
#  arbitrarily select a requence to represent a sequence to represent the pair/
#  alignment, or double the max possible segments number to decide when a pair
#  has reached its max.
#2 Raw number of unique segments (per sequence) included - not just the weakest
#  link - maximize
#3 Int(sqrt()) of the sum of deviation distances (the shortest distance per
#  region) from the ideal spacing per pair/sequence, normalized by the number
#  of regions.  Each sequence in each pair has a static number of non-
#  overlapping identity segments.  Optimal spacing per pair/sequence is even
#  distribution throughout the sequence length.  For each of the ideal
#  locations per pair/sequence, the closest identity segment (in that region)
#  has its distance from the region center used as a score.  All the scores for
#  each pair/sequence/region are summed.  If no segment is in a region, the
#  region size is added.  The truncated (e.g. int value) of the square root of
#  this average distance from the ideal spacing is used for 2 reasons:
#  succinct display in verbose more and more importantly: So that the next base
#  score can have more influence over the resulting solution. - minimize
#4 Number of bases (times number of sources) included in stretches of any size
#  - maximize
sub scoreSolution
  {
    my $soln      = $_[0];
    my $hash      = $_[1]; #Needed to get: ->{$pairid}->{NUMUNIQSEGS}->{$sz}
    my $sizes     = $_[2]; #Stretch minimums
    my $score_obj = {};

    my $ubase_sums  = {}; #Hash used to sum of unique identity base positions
                          #included per pair/sequence
    my $dist_scores = {}; #A hash to keep track of the segment closest to a
                          #region center
    my $occupied    = {}; #A hash to count the number of segments included

    my $lowest_inclusion_score = {map {$_ => undef} @$sizes};
    my $total_inclusion_score  = {map {$_ => 0} @$sizes};
    my $spacing_score          = {};
    my $base_score             = 0;
    my $num_regs               = {};

    #For each sequence we are constructing (present in the solution)
    foreach my $seqid (sort {$a cmp $b} keys(%$soln))
      {
	#For each identity segment record in the solution
	foreach my $subseqrec (grep {$_->{TYPE} eq 'seg'} @{$soln->{$seqid}})
	  {
	    # Spacing score is required to be present
	    if(!defined($subseqrec->{SPACING_SCORE}))
	      {
		error("Cannot call scoreSolution after solution has been ",
		      "reduced by reduceSolutionForIndirectConflicts.");
		return(undef);
	      }

	    my $pairid   = $subseqrec->{PAIR_ID};
	    my $add_size = $subseqrec->{ADD_SIZE};

	    ## Inclusion score calculations... Record the unique base locations

	    #For each position of identity in this segment, record the position
	    #as the key in the hash (then we'll count the number of keys)
	    foreach my $pos ($subseqrec->{IDEN_START}..
			     $subseqrec->{IDEN_STOP})
	      {$ubase_sums->{$add_size}->{$pairid}->{$seqid}->{$pos} = 0}

	    ## Num score (number of non-identity-overlapped segments, any size)

	    foreach my $pos ($subseqrec->{IDEN_START}..$subseqrec->{IDEN_STOP})
	      {$occupied->{$pairid}->{$seqid}->{$pos} = 0}

	    ## Spacing score calculations for each region for each pair/seq...
	    ## Record the distance of the closest identity sequence from the
	    ## region center

	    my $middle = int($subseqrec->{IDEN_START} + $add_size / 2);
	    my $period_len =
	      int(length($hash->{$pairid}->{SEQS}->{$seqid}) /
		  ($hash->{$pairid}->{NUMUNIQSEGS}->{$add_size} + 1));
	    my $region = int($middle / $period_len);
	    if(!exists($dist_scores->{$add_size}) ||
	       !exists($dist_scores->{$add_size}->{$pairid}) ||
	       !exists($dist_scores->{$add_size}->{$pairid}->{$seqid}) ||
	       !exists($dist_scores->{$add_size}->{$pairid}->{$seqid}
		       ->{$region}) ||
	       $subseqrec->{SPACING_SCORE} <
	       $dist_scores->{$add_size}->{$pairid}->{$seqid}->{$region})
	      {$dist_scores->{$add_size}->{$pairid}->{$seqid}->{$region} =
		 $subseqrec->{SPACING_SCORE}}
	  }
      }

    foreach my $size (@$sizes)
      {
	#For all possible pairs available
	foreach my $pairid (sort {$a cmp $b} keys(%$hash))
	  {
	    my $inclusion_score = 0;

	    ## Spacing score calculations

	    #For all possible sequences
	    foreach my $seqid (keys(%{$hash->{$pairid}->{SEQS}}))
	      {
		my $period_len =
		  int(length($hash->{$pairid}->{SEQS}->{$seqid}) /
		      ($hash->{$pairid}->{NUMUNIQSEGS}->{$size} + 1));
		#For all possible regions
		foreach my $reg (0..($hash->{$pairid}->{NUMUNIQSEGS}->{$size} -
				     1))
		  {
		    $num_regs->{$size}++;
		    if(exists($dist_scores->{$size}) &&
		       exists($dist_scores->{$size}->{$pairid}) &&
		       exists($dist_scores->{$size}->{$pairid}->{$seqid}) &&
		       exists($dist_scores->{$size}->{$pairid}->{$seqid}
			      ->{$reg}))
		      {$spacing_score->{$size} +=
			 $dist_scores->{$size}->{$pairid}->{$seqid}->{$reg}}
		    else
		      {$spacing_score->{$size} += $period_len}
		  }
	      }

	    ## Inclusion score calculations...

	    #For each sequence in the pair (Note: the score for each should
	    #technically be the same, so doing it for both here is useless,
	    #though it would be good to do a sanity check here in the future)
	    foreach my $seqid (sort {$a cmp $b}
			       keys(%{$ubase_sums->{$size}->{$pairid}}))
	      {
		my $cur_size            = 1;
		my $tmp_inclusion_score = 0;
		my($last_pos);
		foreach my $pos (sort {$a <=> $b}
				 keys(%{$ubase_sums->{$size}->{$pairid}
					  ->{$seqid}}))
		  {
		    if(defined($last_pos) && ($last_pos + 1) == $pos)
		      {$cur_size++}
		    elsif(defined($last_pos))
		      {
			#Calculate the number of unique/independent segments
			my $usegs = int($cur_size / $size);
			if($usegs > 0)
			  {$tmp_inclusion_score += $usegs}
			$cur_size = 1;
		      }
		    $last_pos = $pos;
		  }

		#Calculate the last number of unique/independent segments
		my $usegs = int($cur_size / $size);
		if($usegs > 0)
		  {$tmp_inclusion_score += $usegs}

		if(!defined($lowest_inclusion_score->{$size}) ||
		   $tmp_inclusion_score < $lowest_inclusion_score->{$size})
		  {
		    if($tmp_inclusion_score >
		       $hash->{$pairid}->{NUMUNIQSEGS}->{$size})
		      {error("The number of unique segments included in the ",
			     "solution: [$tmp_inclusion_score] is larger ",
			     "than the max possible calculated earlier: ",
			     "[$hash->{$pairid}->{NUMUNIQSEGS}->{$size}].  ",
			     "This should not have happened.  There must ",
			     "be a bug somewhere in the code.")}
		    $lowest_inclusion_score->{$size} = $tmp_inclusion_score;
		  }
		$total_inclusion_score->{$size} += $tmp_inclusion_score;
		$inclusion_score += $tmp_inclusion_score;
	      }

	    if(!defined($lowest_inclusion_score->{$size}) ||
	       $inclusion_score < $lowest_inclusion_score->{$size})
	      {$lowest_inclusion_score->{$size} = $inclusion_score}
	  }
      }

    ## Bases score calculation

    foreach my $pairid (keys(%$occupied))
      {foreach my $seqid (keys(%{$occupied->{$pairid}}))
	 {$base_score += scalar(keys(%{$occupied->{$pairid}->{$seqid}}))}}
    $score_obj->{BASESCORE} = $base_score;

    foreach my $size (@$sizes)
      {
	## Normalize the spacing score by the number of total regions in all
	## pairs/sequences
	$spacing_score->{$size} /= $num_regs->{$size};
	#Take the square root and runcate both for display and so that base
	#score has more influence
	$spacing_score->{$size} = int(sqrt($spacing_score->{$size}));

	$score_obj->{STRETCHES}->{$size} = [$lowest_inclusion_score->{$size},
					    $total_inclusion_score->{$size},
					    $spacing_score->{$size}];
      }

    return($score_obj);
  }

sub candidateBetter
  {
    my $cand_obj = $_[0];
    my $best_obj = $_[1];

    #If one of the solutions is not fully defined (defaulting to cand = better
    if(!defined($best_obj) || !exists($best_obj->{STRETCHES}))
      {return(1)}
    elsif(!defined($cand_obj) || !exists($cand_obj->{STRETCHES}))
      {return(0)}

    my $cand_sizes = [sort {$b <=> $a} keys(%{$cand_obj->{STRETCHES}})];
    my $best_sizes = [sort {$b <=> $a} keys(%{$best_obj->{STRETCHES}})];

    #If one of the solutions is empty in terms of stretch sizes
    if(scalar(@$best_sizes) && scalar(@$cand_sizes) == 0)
      {return(0)}
    elsif(scalar(@$best_sizes) == 0 && scalar(@$cand_sizes))
      {return(1)}

    #If the solutions consist of different complements of stretch sizes
    if(scalar(@$best_sizes) != scalar(@$cand_sizes) ||
       join(',',@$best_sizes) ne join(',',@$cand_sizes))
      {
	for(my $i = 0;$i <= $#{$best_sizes};$i++)
	  {
	    #If the candidate has fewer size segments
	    if($#{$cand_sizes} < $i)
	      {return(0)}

	    if($best_sizes->[$i] > $cand_sizes->[$i])
	      {return(0)}
	    elsif($best_sizes->[$i] < $cand_sizes->[$i])
	      {return(1)}
	    else
	      {
		my $size = $best_sizes->[$i];
		my $best = $best_obj->{STRETCHES}->{$size};
		my $cand = $cand_obj->{STRETCHES}->{$size};
		my $cmp  = candidateSizeCmp($cand,$best);

		if($cmp == -1)
		  {return(0)}
		elsif($cmp == 1)
		  {return(1)}
	      }
	  }
	#The candidate must have more segment sizes (and everything else is the
	#same)
	return(1);
      }

    #If we get here, each solution has the same complement of stretch size data
    #(as a type)
    foreach my $size (@$best_sizes)
      {
	my $best = $best_obj->{STRETCHES}->{$size};
	my $cand = $cand_obj->{STRETCHES}->{$size};
	my $cmp  = candidateSizeCmp($cand,$best);

	if($cmp == -1)
	  {return(0)}
	elsif($cmp == 1)
	  {return(1)}
      }

    #The stretches are all equivalent, so let's use the number of uniquely
    #covered bases
    if($cand_obj->{BASESCORE} > $best_obj->{BASESCORE})
      {return(1)}

    #The candidate is not better
    return(0);
  }

sub candidateSizeCmp
  {
    my $cand = $_[0];
    my $best = $_[1];

    if(!defined($best) || scalar(@$best) < 3 ||
       $cand->[0] > $best->[0] ||
       ($cand->[0] == $best->[0] && $cand->[1] > $best->[1]) ||
       ($cand->[0] == $best->[0] && $cand->[1] == $best->[1] &&
	$cand->[2] < $best->[2]))
      {return(1)}
    elsif($cand->[0] == $best->[0] && $cand->[1] == $best->[1] &&
	  $cand->[2] == $best->[2])
      {return(0)}

    return(-1);
  }

#Creates a hash of hybrid sequences
sub generateHybrids
  {
    my $soln   = $_[0];
    my $data   = $_[1];
    my $usage  = $_[2];
    my $matrix = $_[3];
    my $sizes  = $_[4];

    verbose("Reducing solution...\n",
	    "Merging overlapping identity segments from common alignments...");

    $soln = reduceSolutionMergeOverlaps($soln);

    verbose("Solution:\n",solutionToString($soln));

    #This has to be called on the merged overlapping segments version of the
    #solution, which contains all the duplicate segments from various pairwise
    #alignments, before the solution is further reduced to address conflicts
    #and duplicates
    my $soln_stats = getSolutionStats($soln,$data,$sizes);

    verbose("Adjusting overlapping identity segments from different ",
	    "alignments...");

    $soln = reduceSolutionForDirectConflicts($soln);

    debug("Reduced solution before splitting segments for the recoding map:",
	  "\n",solutionToString($soln));

    #Get the alignments each sequence is involved in, to be used as a reference
    #when constructing the map of how each sequence should be composed (i.e.
    #defining which alignment each portion of the sequence should come from, or
    #which partner sequence it should be recoded to best match).
    my $partners_hash = getPartnersHash($soln,$data);

    verbose("Splitting identity segments based on indirect overlap for the ",
	    "recoding map...");

    $soln = reduceSolutionForIndirectConflicts($soln,$partners_hash,$data);

    debug("Reduced solution with merged overlapping identity segments, ",
	  "adjusted overlapping identity segments from different sources, ",
	  "and split identity segments for determining recoding:\n",
	  solutionToString($soln),{LEVEL => 2});

    verbose("Updating solution final coordinates...");

    $soln = updateFinalCoords($soln,$data);

    debug("Reduced solution with updated final coordinates and basic ",
	  "filtering for duplicate regions:\n",solutionToString($soln),
	  {LEVEL => 2});

    verbose("Expanding solution coordinates...");

    #Now generate the map of where each portion of the sequences should come
    #from, or be based on (for recoding to match a portion of a sequence that
    #was optimized for identity with another sequence)
    my $map = getFilledMap($soln,$data,$partners_hash);

    issueDividerWarnings($map);

    verbose("Building solution sequences...");

    my($woven_seqs,$source_seqs,$validation_message) =
      weaveSeqs($map,$data,$partners_hash,$matrix,$usage,$soln);

    return([$woven_seqs,$source_seqs,$soln_stats,$validation_message]);
  }

#Merges overlapping segments derived from the same pair by extending a single
#record to represent an entire area of overlap
sub reduceSolutionMergeOverlaps
  {
    my $soln     = $_[0];
    my $new_soln = {};

    foreach my $seqid (sort {$a cmp $b} keys(%$soln))
      {
	my @ovlp_buffer = ();
	foreach my $currec (@{$soln->{$seqid}})
	  {
	    #Add any previous records that no longer overlap the current record
	    #to the new solution or else add them to the tmp_ovlp_buffer.  It
	    #is assumed that the records are in order of IDEN_START, and that
	    #that's always less than (or equal to) IDEN_STOP (and
	    #LAST_CODON_STOP.  We're also assuming that FIRST_CODON_START
	    #overlaps IDEN_START.  So once the IDEN_START is past
	    #LAST_CODON_STOP, then there is no more overlap with that record,
	    #thus we can add it to the reduced solution.
	    #Note: Nothing in the ovlp_buffer overlaps (except alternate codons
	    #which we are ignoring), because it has been explicitly edited to
	    #not do so.  It represents segments that did overlap, but they are
	    #kept here to see if they will overlap the next record
	    my @tmp_ovlp_buffer = ();
	    foreach my $lastrec (@ovlp_buffer)
	      {
		#If we have gone past this "last" record, it's good to add.
		#DO NOT edit to merge abutting segments. There should always be
		#overlap because the segments are added exhaustively and exist
		#in a staggered fashion.  If they abutt and there's no
		#overlapping piece that connects them, that means that the
		#partner's segments do not abutt.
		if($lastrec->{IDEN_STOP} < $currec->{IDEN_START})
		  {push(@{$new_soln->{$seqid}},$lastrec)}
		else
		  {push(@tmp_ovlp_buffer,$lastrec)}
	      }
	    @ovlp_buffer = @tmp_ovlp_buffer;

	    #If there are no remaining records in the ovlp_buffer or if this
	    #is an alternate codon (which does not get merged with other
	    #records), add this record to the ovlp_buffer
	    if(scalar(@ovlp_buffer) == 0 || $currec->{TYPE} eq 'alt')
	      {push(@ovlp_buffer,{%$currec})}
	    else
	      {
		my $add = 1;
		#The remaining records all overlap the current record being
		#added.  Search through them to find the one from the same pair
		#source and edit its record to extend the stops (unless one
		#from the same pair doesn't exist - then "$add" it to the
		#ovlp_buffer as its own thing to be searched for overlap on the
		#next iteration)
		foreach my $lastrec (grep {$_->{TYPE} eq 'seg' &&
					     $_->{PAIR_ID} eq
					       $currec->{PAIR_ID}}
				     @ovlp_buffer)
		  {
		    if($currec->{IDEN_STOP} > $lastrec->{IDEN_STOP})
		      {
			$lastrec->{IDEN_STOP}       = $currec->{IDEN_STOP};
			$lastrec->{LAST_CODON_STOP} = $currec
			  ->{LAST_CODON_STOP};
			$lastrec->{FINAL_STOP}      = $currec->{FINAL_STOP};
			if($lastrec->{ADD_SIZE} < $currec->{ADD_SIZE})
			  {error("Smaller stretch size ",
				 "[$lastrec->{ADD_SIZE}] merging with a ",
				 "larger stretch size identity segment ",
				 "[$currec->{ADD_SIZE}].",
				 {DETAIL =>
				  join('',('This should not happen because ',
					   'of the ordering of the search ',
					   'space and will corrupt the ',
					   'solution score calculations.'))})}
			#Keep the score of the better placed piece of identity
			if($currec->{SPACING_SCORE} <
			   $lastrec->{SPACING_SCORE})
			  {$lastrec->{SPACING_SCORE} =
			     $currec->{SPACING_SCORE}}
		      }
		    $add = 0;
		    last;
		  }
		if($add)
		  {push(@ovlp_buffer,{%$currec})}
	      }
	  }
	#All the records have been traversed.  Anything left in @ovlp_buffer
	#can be added to the reduced solution
	push(@{$new_soln->{$seqid}},@ovlp_buffer);
      }

    return($new_soln);
  }

#This must only be called after merging overlaps and before calling any
#reduceSolution* methods
sub getSolutionStats
  {
    my $soln  = $_[0];
    my $data  = $_[1];
    my $sizes = $_[2];
    my $stats = {};

    foreach my $seqid (sort {$a cmp $b} keys(%$soln))
      {
	#Smaller stretch sizes can merge when included to represent a
	#duplicate segment record as the larger ones.
	my $seen = {};

	foreach my $rec (grep {$_->{TYPE} eq 'seg'} @{$soln->{$seqid}})
	  {
	    my $pair_id = $rec->{PAIR_ID};
	    my($seqid1,$seqid2) =
	      sort {$a cmp $b} keys(%{$data->{$pair_id}->{SEQS}});

	    #We will cycle through all the seqids in the loop above and set the
	    #reciprocal values by making sure that seqid1 is always the seqid
	    #from the above loop.  Then when we get to its partner in the above
	    #loop, the reciprocal value will be set
	    if($seqid1 ne $seqid)
	      {
		$seqid2 = $seqid1;
		$seqid1 = $seqid;
	      }

	    foreach my $i (sort {$sizes->[$a] <=> $sizes->[$b]} 0..$#{$sizes})
	      {
		my $min1 = $sizes->[$i];
		my $min2 = $i == $#{$sizes} ? 0 : $sizes->[$i + 1];
		my $size = $rec->{IDEN_STOP} - $rec->{IDEN_START} + 1;
		if(!exists($seen
			   ->{"$seqid:$rec->{IDEN_START}-$rec->{IDEN_STOP}:" .
			      "$seqid2"})
		   && $size >= $min1 && ($i == $#{$sizes} || $size < $min2))
		  {
		    my $times = int($size / $min1);
		    if($i == 0 && !$times)
		      {error("Segment size in the solution for [$seqid] is ",
			     "below the minimum: [$size / $min1].")}
		    $stats->{$seqid1}->{$seqid2}->{$min1} += $times;

		    $seen->{"$seqid:$rec->{IDEN_START}-$rec->{IDEN_STOP}:" .
			    "$seqid2"} = 0;
		  }
	      }
	  }
      }

    return($stats);
  }

#Eliminates overlap between segments of the same sequence, but derived from
#different pairs by either removing whole records that are completely
#overlapped or by adjusting their coordinates by shrinking them inward.
#Reductions are done by codon boundaries.
sub reduceSolutionForDirectConflicts
  {
    my $soln = $_[0];
    my $new_soln = {};
    foreach my $seqid (sort {$a cmp $b} keys(%$soln))
      {
	#This should really be fixed elsewhere.  FINAL_START and FINAL_STOP
	#should be defined from the beginning
	foreach my $rec (@{$soln->{$seqid}})
	  {
	    if(!defined($rec->{FINAL_START}))
	      {$rec->{FINAL_START} = $rec->{FIRST_CODON_START}}
	    if(!defined($rec->{FINAL_STOP}))
	      {$rec->{FINAL_STOP} = $rec->{LAST_CODON_STOP}}
	  }
	my @ovlp_buffer = ();
	foreach my $oldrec (sort {$a->{FINAL_START} <=> $b->{FINAL_START}}
			    @{$soln->{$seqid}})
	  {
	    #Copy the soln record
	    my $currec = {%$oldrec};
	    if(!exists($currec->{FINAL_START}) ||
	       !defined($currec->{FINAL_START}))
	      {$currec->{FINAL_START} = $currec->{FIRST_CODON_START}}
	    if(!exists($currec->{FINAL_STOP}) ||
	       !defined($currec->{FINAL_STOP}))
	      {$currec->{FINAL_STOP} = $currec->{LAST_CODON_STOP}}

	    #Add any previous records that no longer overlap the current record
	    #to the new solution or else add them to the tmp_ovlp_buffer.  It
	    #is assumed that the records are in order of IDEN_START, and that
	    #that's always less than (or equal to) IDEN_STOP (and
	    #LAST_CODON_STOP.  We're also assuming that FIRST_CODON_START
	    #overlaps IDEN_START.  So once the IDEN_START is past
	    #LAST_CODON_STOP, then there is no more overlap with that record,
	    #thus we can add it to the reduced solution.
	    #Note: Nothing in the ovlp_buffer overlaps (except alternate codons
	    #which we are ignoring), because it has been explicitly edited to
	    #not do so.  It represents segments that did overlap, but they are
	    #kept here to see if they will overlap the next record
	    my @tmp_ovlp_buffer = ();
	    foreach my $lastrec (@ovlp_buffer)
	      {
		#If we have gone past this "last" record, it's good to add
		if($lastrec->{FINAL_START} < $currec->{FINAL_STOP})
		  {push(@{$new_soln->{$seqid}},$lastrec)}
		else
		  {push(@tmp_ovlp_buffer,$lastrec)}
	      }
	    @ovlp_buffer = @tmp_ovlp_buffer;

	    #If there are no remaining records in the ovlp_buffer or if this
	    #is an alternate codon (which does not get merged with other
	    #records), add this record to the ovlp_buffer
	    if(scalar(@ovlp_buffer) == 0 || $currec->{TYPE} eq 'alt')
	      {push(@ovlp_buffer,$currec)}
	    else
	      {
		my @lastrecs =
		  sort {$a->{FINAL_STOP} <=> $b->{FINAL_STOP}}
		    grep {$_->{TYPE} eq 'seg' &&
			    $currec->{FINAL_STOP} > $_->{FINAL_STOP}}
		      @ovlp_buffer;
		my $add = 1;

		if(scalar(@lastrecs) == 0)
		  {
		    push(@ovlp_buffer,$currec);
		    next;
		  }
		#The remaining records all overlap the current record being
		#added.  Search through them to find the one from the same pair
		#source and edit its record to adjust the starts and stops
		#(unless an overlap of type seg doesn't exist - then "$add" it
		#to the #ovlp_buffer as its own thing to be searched for
		#overlap on the next iteration)
		my $lastrec = $lastrecs[-1];

		##
		## The point here is to figure out:
		## 1. Whether we will add or skip the current record to the
		##    overlap buffer,
		## 2. whether we should trim back the end of the last
		##    lastrec (which could be necessary so that the
		##    identity of what's being added will be used instead
		##    of any extra non-identity codon bases in the last
		##    lastrec) (it's assumed that the last lastrec is the
		##    one that will be adjusted), and
		## 3. how much we need to increase the start of the record
		##    being added to not have any final overlap
		##

		#If the end of the segment being added is overlapping
		#anything previously added (complete overlap is assumed
		#because we're guaranteed to have a start overlap and
		#everything should be contiguous), don't add it and stop
		if($currec->{FINAL_STOP} <= $lastrec->{FINAL_STOP})
		  {
		    $add = 0;
		    next;
		  }

		#Determine the last base of the previously added record
		#that we are definitely going to keep, which is the last
		#frame 3 position of identity.  Everything else will either
		#be identical to the overlapping thing being added or will
		#get replaced by an alt codon
		my $last_keep_to = $lastrec->{FINAL_STOP};
		if($last_keep_to > $lastrec->{IDEN_STOP})
		  {$last_keep_to -= 3}
		while($currec->{FINAL_START} <= $last_keep_to)
		  {
		    $currec->{FINAL_START} += 3;
		    if($currec->{FINAL_START} > $currec->{IDEN_STOP} &&
		       defined($currec->{TYPE}) && $currec->{TYPE} eq 'SOURCE')
		      {$currec->{TYPE} = ''}
		  }

		#If the current record is no longer valid or no longer
		#contains any identity
		if($currec->{FINAL_START} > $currec->{FINAL_STOP} ||
		   $currec->{FINAL_START} > $currec->{IDEN_STOP})
		  {
		    $add = 0;
		    next;
		  }
		else
		  {
		    #If we adjusted the the end of the lastrec temporarily,
		    #we adjusted the start of the current record being
		    #added, make the adjustment to the last rec permanent.
		    #We will assume that this will only happen to the last
		    #lastrec, ordered by final start
		    if(($last_keep_to + 3) ==
		       $lastrec->{LAST_CODON_STOP} &&
		       $currec->{FINAL_START} !=
		       $currec->{FIRST_CODON_START} &&
		       ($last_keep_to + 1) == $currec->{FINAL_START})
		      {
			if($last_keep_to + 1 != $currec->{FINAL_START} ||
			   $currec->{FINAL_START} % 3 != 1)
			  {error("Internal error: Bad frame boundary.  This ",
				 "should not have happened.")}
			$lastrec->{FINAL_STOP} = $last_keep_to;
			if($lastrec->{FINAL_STOP} < $lastrec->{IDEN_START} &&
			   defined($lastrec->{TYPE}) &&
			   $lastrec->{TYPE} eq 'SOURCE')
			  {$lastrec->{TYPE} = ''}
		      }
		  }
		if($add)
		  {push(@ovlp_buffer,$currec)}
	      }
	  }
	#All the records have been traversed.  Anything left in @ovlp_buffer
	#can be added to the reduced solution
	push(@{$new_soln->{$seqid}},@ovlp_buffer);
      }

    return($new_soln);
  }

sub calcFrame
  {
    my $coord = $_[0];
    if(!defined($coord))
      {
	error("Coordinate not defined.");
	return(0);
      }
    return(($coord % 3) ? ($coord % 3) : 3);
  }

#Eliminates overlap between segments of different sequences, derived from
#unrelated pairs by splitting records that are overlapping.  For example, say
#we have alignments between sequences A&B, A&C, and C&D.  A segment in A&B has
#corresponding coordinates in C (where A aligns with it) and those coordinates
#correspond to coordinates in D.  Those coordinates can overlap.  The
#overlapping portions don't change, but the non-overlapping portions adjacent
#to each segment involved are to be recode to best match the other segment.
#And the way the assigned recoding works, by copying dividers from one sequence
#to another, those dividers need to be split up the the proper sourcing of the
#sequence (whether it's obtained directly from the alignment where it's
#identity segment is or whether it's to be recoded, needs to not get copied
#into the middle of a segment.  By splitting it up and keeping the TYPEs and
#pair IDs up to date, those dividers can be properly sourced.
#Reductions are done by codon boundaries.
sub reduceSolutionForIndirectConflicts
  {
    my $soln          = $_[0];
    my $partners_hash = $_[1];
    my $data          = $_[2];

    #We're not going to be removing any records due to this type of overlap -
    #only dividing them up into separate records, so we'll copy the solution to
    #be able to make modifications to the copy without altering the one sent in
    my $new_soln      = {};
    foreach my $seqid (sort {$a cmp $b} keys(%$soln))
      {
	foreach my $oldrec (sort {$a->{FINAL_START} <=> $b->{FINAL_START}}
			    @{$soln->{$seqid}})
	  {push(@{$new_soln->{$seqid}},{%$oldrec})}
      }

    my $still_finding_overlap = 0;
    my $warn_limit = 4;
    my $iter_num = 0;

    do
      {
	$iter_num++;

	if($iter_num >= $warn_limit)
	  {verbose("The search for indirect overlap is going longer than ",
		   "expected.  Found $still_finding_overlap splits on ",
		   "iteration ",($iter_num - 1)," (the number of new splits ",
		   "found each iteration should be (generally) trending ",
		   "down).")}

	$still_finding_overlap = 0;

	#For each independent sequence representation in the solution
	foreach my $seqid (sort {$a cmp $b} keys(%$new_soln))
	  {
	    #For each segment record in an individual sequence
	    foreach my $rec (sort {$a->{FINAL_START} <=> $b->{FINAL_START}}
			     @{$new_soln->{$seqid}})
	      {
		#Note, this sequence has a direct partner, but it's possible
		#that a shrunken segment could create a divider that interrupts
		#a direct partner sequence because the 2 partners are
		#indirectly related, so I will include it in the partner check
		#below instead of skipping it.  E.g. hs20-30 overlaps hs27-37
		#(could be an alt codon in the middle, but we'll ignore that
		#for clarity).  One of the 2 gets shrunk by
		#reduceSolutionForDirectConflicts.  The partners in each case
		#respectively may be mj:15-25 and pf:30-40.  If the first hs
		#segment gets shrunk to hs:20-26, a divider would get copied to
		#the middle of mj:15-25, so mj should be split to anticipate
		#this.  Look at each sequence this sequence is involved in an
		#alignment with (note, only ones where identity segments were
		#found will be included)
		foreach my $indirect_partner (@{$partners_hash->{$seqid}})
		  {
		    my $indirect_partner_pairid = $indirect_partner->[0];
		    my $indirect_partner_seqid  = $indirect_partner->[1];
		    my $converted_final_start =
		      convertDivider($seqid,
				     $rec->{FINAL_START},
				     $indirect_partner_pairid,
				     $indirect_partner_seqid,
				     $data);
		    my $converted_final_stop =
		      convertAlnSeqCoord($seqid,
					 $rec->{FINAL_STOP},
					 $indirect_partner_pairid,
					 $indirect_partner_seqid,
					 $data);

		    if(calcFrame($converted_final_start) != 1)
		      {error("Converted alignment start coordinate ",
			     "[$indirect_partner_seqid:",
			     "$converted_final_start] (converted from final ",
			     "start: [$seqid:$rec->{FINAL_START}] in ",
			     "alignment [$indirect_partner_pairid]) not in ",
			     "frame 1.")}
		    if(calcFrame($converted_final_stop) != 3)
		      {error("Converted alignment stop coordinate ",
			     "[$indirect_partner_seqid:",
			     "$converted_final_stop] (converted from final ",
			     "stop: [$seqid:$rec->{FINAL_STOP}] in alignment ",
			     "[$indirect_partner_pairid]) not in frame 3.")}

		    #Now we need to search to see if there's overlap with any
		    #segments in the partner
		    foreach my $partner_rec (sort {$a->{FINAL_START} <=>
						     $b->{FINAL_START}}
					     @{$new_soln
						 ->{$indirect_partner_seqid}})
		      {
			my $partner_final_start = $partner_rec->{FINAL_START};
			my $partner_final_stop  = $partner_rec->{FINAL_STOP};

			my $msg =
			  join('',("Indirect overlap search: Iter: $iter_num ",
				   "Splits: $still_finding_overlap ",
				   "$seqid:$rec->{FINAL_START}-",
				   "$rec->{FINAL_STOP} vs ",
				   "($indirect_partner_seqid:",
				   "$partner_final_start-$partner_final_stop ",
				   "-> $converted_final_start-",
				   "$converted_final_stop)"));
			verboseOverMe($msg);

			#We're going to change the partners.  The partner will
			#or already did change us.

			#If a converted edge coord is inside the partner record
			if(($converted_final_start > $partner_final_start &&
			    $converted_final_start <= $partner_final_stop) ||
			   ($converted_final_stop >= $partner_final_start &&
			    $converted_final_stop < $partner_final_stop))
			  {
			    verboseOverMe($msg," Splitting overlap.");
			    splitSegment($new_soln,
					 $indirect_partner_seqid,
					 $partner_rec,
					 $converted_final_start,
					 $converted_final_stop);
			    $still_finding_overlap++;
			  }

			#If we're past the possibility of finding overlap
			last if($partner_final_start > $converted_final_stop);
		      }
		  }
	      }
	  }
      } while($still_finding_overlap);
    #The above while is there in case I'm creating new potential dividing
    #points in a sequence that was already processed in a previous foreach loop
    #iteration

    #Some of the coordinates changed, so let's make sure things are still
    #ordered as intended
    foreach my $seqid (keys(%$new_soln))
      {@{$new_soln->{$seqid}} =
	 sort {$a->{IDEN_START} <=> $b->{IDEN_START} ||
		 $b->{IDEN_STOP} <=> $a->{IDEN_STOP}} @{$new_soln->{$seqid}}}

    #Keep the last verboseOverMe message on the screen
    if(isVerbose())
      {print STDERR ("\n")}

    return($new_soln);
  }

#Assumes that everything will happen on codon boundaries because codons are
#aligned and all segment boundaries should be on codon boundaries.  Only call
#this from reduceSolutionForIndirectConflicts (because it invalidates the
#solution for other methods such as scoreSolution and any other method that
#uses SPACING_SCOREs or depends on identity segments being unmanipulated/split
#up).
sub splitSegment
  {
    my $soln  = $_[0];
    my $seqid = $_[1];
    my $rec   = $_[2];
    my $start = $_[3];
    my $stop  = $_[4];

    #Not to be used on alt codons
    if($rec->{TYPE} eq 'alt')
      {return()}

    #This won't change
    my $orig_size = $rec->{ADD_SIZE};

    #If the start overlaps anything but the start and the stop overlaps
    #anything but the stop and the segment is larger than 2 codons
    if($start > $rec->{FINAL_START} && $start <= $rec->{FINAL_STOP} &&
       $stop >= $rec->{FINAL_START} && $stop  <  $rec->{FINAL_STOP} &&
       ($rec->{FIRST_CODON_START} + 5) > $rec->{LAST_CODON_STOP})
      {
	#The record must be broken into 3 pieces.  Edit the record to represent
	#the left side and add the other 2

	#The kept piece on the left side:
	my $saved_final_stop    = $rec->{FINAL_STOP};
	$rec->{FINAL_STOP}      = $start - 1;
	my $saved_iden_stop     = $rec->{IDEN_STOP};
	$rec->{IDEN_STOP}       = $start - 1;
	$rec->{LAST_CODON_STOP} = $start - 1;

	#The completely overlapping piece:
	addSegment($soln,
		   $rec->{PAIR_ID},
		   $seqid,
		   $start,
		   $stop,
		   $start,
		   $stop,
		   $rec->{TYPE},
		   $rec->{ALT_CODON},
		   $start,
		   $stop,
		   undef,                #Don't replicate the spacing score
		   $orig_size);


	#Add a new piece on the right side
	addSegment($soln,
		   $rec->{PAIR_ID},
		   $seqid,
		   $stop + 1,
		   $saved_iden_stop,
		   $stop + 1,
		   $saved_final_stop,
		   $rec->{TYPE},
		   $rec->{ALT_CODON},
		   $stop + 1,
		   $saved_final_stop,
		   undef,                #Don't replicate the spacing score
		   $orig_size);
      }
    #Else if the start overlaps anything but the start and the stop overlaps
    #anything but the stop and the segment is 2 codons in size
    elsif($start > $rec->{FINAL_START} && $start <= $rec->{FINAL_STOP} &&
	  $stop >= $rec->{FINAL_START} && $stop  <  $rec->{FINAL_STOP} &&
	  ($rec->{FIRST_CODON_START} + 5) == $rec->{LAST_CODON_STOP})
      {
	#The record must be broken into 3 pieces.  Edit the record to represent
	#the left side and add the other 2

	#The kept piece on the left side:
	my $saved_final_stop    = $rec->{FINAL_STOP};
	$rec->{FINAL_STOP}      = $start - 1;
	my $saved_iden_stop     = $rec->{IDEN_STOP};
	$rec->{IDEN_STOP}       = $start - 1;
	$rec->{LAST_CODON_STOP} = $start - 1;

	#Add a new piece on the right side
	addSegment($soln,
		   $rec->{PAIR_ID},
		   $seqid,
		   $stop + 1,
		   $saved_iden_stop,
		   $stop + 1,
		   $saved_final_stop,
		   $rec->{TYPE},
		   $rec->{ALT_CODON},
		   $stop + 1,
		   $saved_final_stop,
		   undef,                #Don't replicate the spacing score
		   $orig_size);
      }
    #Elsif the start overlaps anything but the start
    elsif($start > $rec->{FINAL_START} && $start <= $rec->{FINAL_STOP})
      {
	#Edit the record to represent the left side and add the other

	#The kept piece on the left side:
	my $saved_final_stop    = $rec->{FINAL_STOP};
	$rec->{FINAL_STOP}      = $start - 1;
	my $saved_iden_stop     = $rec->{IDEN_STOP};
	$rec->{IDEN_STOP}       = $start - 1;
	$rec->{LAST_CODON_STOP} = $start - 1;

	#Add a new piece on the right side
	addSegment($soln,
		   $rec->{PAIR_ID},
		   $seqid,
		   $start,
		   $saved_iden_stop,
		   $start,
		   $saved_final_stop,
		   $rec->{TYPE},
		   $rec->{ALT_CODON},
		   $start,
		   $saved_final_stop,
		   undef,                #Don't replicate the spacing score
		   $orig_size);
      }
    #Elsif the stop overlaps anything but the stop
    elsif($stop >= $rec->{FINAL_START} && $stop < $rec->{FINAL_STOP})
      {
	#Edit the record to represent the right side and add the other

	#The kept piece on the right side:
	my $saved_final_start     = $rec->{FINAL_START};
	$rec->{FINAL_START}       = $stop + 1;
	my $saved_iden_start      = $rec->{IDEN_START};
	$rec->{IDEN_START}        = $stop + 1;
	$rec->{FIRST_CODON_START} = $stop + 1;

	#Add a new piece on the left side
	addSegment($soln,
		   $rec->{PAIR_ID},
		   $seqid,
		   $saved_iden_start,
		   $stop,
		   $saved_final_start,
		   $stop,
		   $rec->{TYPE},
		   $rec->{ALT_CODON},
		   $saved_final_start,
		   $stop,
		   undef,                #Don't replicate the spacing score
		   $orig_size);
      }
    else
      {error("Unexpected case.  No overlap was found for splitting the ",
	     "submitted record using the start and stop coordinates ",
	     "provided.")}
  }

#This edits the final starts and stops of portions of sequence, originating
#from different pairs, in each record in the solution, to not overlap.  It does
#not handle overlap of portions of sequence from the same pair.  That's handled
#by reduceSolution.  It has a different purpose - to eliminate overlapping
#records by expanding those from the same pair.  This subroutine however
#shrinks records to be able to tell where portions of sequence should be
#derived from in the final mixed sequence output.  Note that all overlap of
#identity coords is guaranteed to be the same sequence - even though they come
#from different pair sources.  Overlap of edge codons is NOT the same and must
#be dealt with as an entire unit so that the encoded AA does not change.
sub updateFinalCoords
  {
    my $soln     = $_[0];
    my $data     = $_[1];
    my $new_soln = {};

    foreach my $seqid (sort {$a cmp $b} keys(%$soln))
      {
	my @current = ();
	#For each segment record (in order of ascending IDEN_START or
	#descending IDEN_STOP).
	foreach my $currec (@{$soln->{$seqid}})
	  {
	    #Error-check to make sure that the FINAL values are fresh
	    if(!defined($currec->{FINAL_START}) ||
	       $currec->{FIRST_CODON_START} != $currec->{FINAL_START})
	      {
		warning("Final start coordinate has an unexpected value.  ",
			"Resetting.",
			{DETAIL =>
			 join(',',("updateFinalCoords should only be called ",
				   "once on a solution.  It expects the ",
				   "final start and stop to be the same as ",
				   "the first codon start and the last codon ",
				   "stop."))})
		  if(defined($currec->{FINAL_STOP}));
		$currec->{FINAL_START} = $currec->{FIRST_CODON_START};
	      }
	    if(!defined($currec->{FINAL_STOP}) ||
	       $currec->{LAST_CODON_STOP} != $currec->{FINAL_STOP})
	      {
		warning("Final stop coordinate has an unexpected value.  ",
			"Resetting.",
			{DETAIL =>
			 join(',',("updateFinalCoords should only be called ",
				   "once on a solution.  It expects the ",
				   "final start and stop to be the same as ",
				   "the first codon start and the last codon ",
				   "stop."))})
		  if(defined($currec->{FINAL_STOP}));
		$currec->{FINAL_STOP} = $currec->{LAST_CODON_STOP};
	      }

	    #Add any previous records that no longer overlap the current record
	    #to the new solution or else add them to the tmpcurrent array.  It
	    #is assumed that the records are in order of IDEN_START, and that
	    #that's always less than (or equal to) IDEN_STOP (and
	    #LAST_CODON_STOP.  We're also assuming that FIRST_CODON_START
	    #overlaps IDEN_START.  So once the IDEN_START is past
	    #LAST_CODON_STOP, then there is no more overlap with that record,
	    #thus we can add it to the reduced solution.
	    my @tmpcurrent = ();
	    foreach my $lastrec (@current)
	      {
		#If we have gone past this "last" record, it's good to add
		#Assumes FINAL_START has not changed yet
		if($lastrec->{FINAL_STOP} < $currec->{FINAL_START})
		  {push(@{$new_soln->{$seqid}},$lastrec)}
		else
		  {push(@tmpcurrent,$lastrec)}
	      }
	    @current = @tmpcurrent;

	    #If there are no remaining records in the current array, add this
	    #record to the current array
	    if(scalar(@current) == 0)
	      {push(@current,{%$currec})}
	    else
	      {
		#The remaining records all overlap.  Search through them to
		#find any whose final coordinates need to be adjusted
		foreach my $lastrec (@current)
		  {
		    #Only look at records that overlap.  Some may have been
		    #adjusted and can be skipped
		    unless(coordsOverlap($currec->{FINAL_START},
					 $currec->{FINAL_STOP},
					 $lastrec->{FINAL_START},
					 $lastrec->{FINAL_STOP}))
		      {next}

		    #If lastrec is an alternate codon (and currec is not)
		    if($lastrec->{TYPE} eq 'alt' && $currec->{TYPE} eq 'seg')
		      {clipForAltCodon($lastrec,$currec)}
		    #Else if currec is an alternate codon (and lastrec is not)
		    elsif($lastrec->{TYPE} eq 'seg' &&
			  $currec->{TYPE}  eq 'alt')
		      {clipForAltCodon($currec,$lastrec)}
		    #Else if neither are alternate codon
		    elsif($lastrec->{TYPE} eq 'seg' &&
			  $currec->{TYPE}  eq 'seg')
		      {clipSegments($lastrec,$currec)}
		    #Else error
		    else
		      {
			error("Overlapping alternate codons encountered at ",
			      "position: [$seqid:$lastrec->{FINAL_START}]: ",
			      "[$lastrec->{ALT_CODON}] versus ",
			      "[$currec->{ALT_CODON}].  Skipping ",
			      "[$currec->{ALT_CODON}] & arbitrarily ",
			      "selecting: [$lastrec->{ALT_CODON}].",
			      {DETAIL => ("Proposed alt position: " .
					  "[$seqid:$currec->{FINAL_START}]")})
			  unless($lastrec->{ALT_CODON} eq
				 $currec->{ALT_CODON});
		      }
		  }

		#When segments are completely overlapped by other segments,
		#their FINAL_STOP is set to 0 (by clipSegments above).  Other
		#math in adjusting for overlap can also result in the
		#FINAL_START being larger than the FINAL_STOP when total
		#overlap exists.  These segments can be skipped.  However, when
		#a FIRST_START_CODON and LAST_STOP_CODON have the same coords,
		#but their identity segments DO NOT overlap, this is the only
		#case where overlap is permitted - all this is handled in
		#clipSegments by splitting the records at overlap bounds and
		#manipulating the coordinates.
		if($currec->{FINAL_START} <= $currec->{FINAL_STOP})
		  {push(@current,{%$currec})}
	      }
	  }
	#All the records have been traversed.  Anything left in @current can be
	#added to the updated solution
	push(@{$new_soln->{$seqid}},@current);
      }

    return($new_soln);
  }

#This edits the final coordinates of a segment record that overlaps an
#alternate codon record
sub clipForAltCodon
  {
    my $altrec = $_[0];
    my $segrec = $_[1];

    if(($segrec->{FIRST_CODON_START} + 2) == $segrec->{LAST_CODON_STOP})
      {
	debug({LEVEL => 6},"Alternate codon overlaps a split segment edge ",
	      "codon.  Alt codons were originally implemented before segment ",
	      "splitting was implemented.  This is expected to happen now, ",
	      "but a new way of dealing with this has not yet been ",
	      "implemented.  These edge codons either need to be ignored or ",
	      "the identity coordinates they contain must be stored ",
	      "elsewhere (possibly inside the alt codon record itself.");
	return();
      }

    #If the alternate codon overlaps the first codon of the segment and
    #does not overlap the second codon (assumes everything's in
    #the same frame).  Note FINAL coords of alternate codon
    #do not change.
    if($altrec->{FINAL_START} == $segrec->{FIRST_CODON_START} &&
       $altrec->{FINAL_STOP}  != $segrec->{LAST_CODON_STOP})
      {
	#If the FINAL start of the current record is still
	#inside the first edge codon
	if($segrec->{FINAL_START} >= $segrec->{FIRST_CODON_START} &&
	   $segrec->{FINAL_START} <= ($segrec->{FIRST_CODON_START} + 2))
	  {
	    $segrec->{FINAL_START} = $segrec->{FIRST_CODON_START} + 3;

	    if($segrec->{FINAL_STOP} < $segrec->{FINAL_START})
	      {error("Stop less than start.")}
	  }
      }
    #Else if the alternate codon overlaps the last edge
    #codon and does not overlap the first edge codon
    #(assumes all are in the same frame).  Note FINAL
    #coords of alternate codon do not change.
    elsif($altrec->{FINAL_STOP}  == $segrec->{LAST_CODON_STOP} &&
	  $altrec->{FINAL_START} != $segrec->{FIRST_CODON_START})
      {
	#If the FINAL stop of the current record is still
	#inside the last edge codon
	if($segrec->{FINAL_STOP} >= $segrec->{LAST_CODON_STOP} &&
	   $segrec->{FINAL_STOP} <= ($segrec->{LAST_CODON_STOP} - 2))
	  {
	    $segrec->{FINAL_STOP} =
	      $segrec->{LAST_CODON_STOP} - 3;

	    if($segrec->{FINAL_STOP} < $segrec->{FINAL_START})
	      {error("Stop less than start.")}
	  }
      }
    #Else if the alternate codon overlaps both (assumes all
    #are in the same frame).  Note FINAL coords of
    #alternate codon do not change.
    elsif($altrec->{FINAL_START} == $segrec->{FIRST_CODON_START} &&
	  $altrec->{FINAL_STOP}  == $segrec->{LAST_CODON_STOP})
      {error("Alternate codon overlaps both edge codons.")}
    #Else - overlap is inside identity
    else
      {warning("Alternate codon overlaps non-edge codon.",
	       {DETAIL =>
		join(',',
		     ("This is not expected to happen, though it could ",
		      "theoretically happen if the compromise between 2 ",
		      "versions of the same sequence matches a 3rd version ",
		      "of the sequence.  If this warning occurs, the output ",
		      "should be checked for correctness, as code was not ",
		      "written to handle this specific case."))})}
  }

#Shifts the final start & stop of either segment inward to eliminate overlap.
#If one segment is contained inside another, its final stop is changed to 0 to
#indicate that it should be skipped/eliminated.
sub clipSegments
  {
    my $seg1 = $_[0]; #Should start before seg2*
    my $seg2 = $_[1];

    #If seg1 starts after seg2, then it is assumed to have had its start
    #adjusted, which means that seg2 will overlap *something else* and it just
    #got here first before being compared with the other segord which caused
    #seg1's start to have been adjusted.  Thus, the start of seg2 will always
    #be moved past the stop of seg1.

    if($seg2->{FINAL_STOP} == 0)
      {return()}

    #If the stop of segment 2 occurs at or before the stop of segment 1 (and
    #overlaps), set the stop of segment 2 to 0 so that it will get eliminated
    if($seg2->{IDEN_STOP} <= $seg1->{IDEN_STOP} &&
       $seg2->{IDEN_STOP} >= $seg1->{IDEN_START})
      {
	#If the stop of segment 2 is inside segment 1, every other portion of
	#segment 2 to the left is assumed to overlap because of the way the
	#calling method keeps a buffer of overlapping segments.  However, if
	#the overlapping segments are both end codons, there's a chance that
	#their identities do not (fully) overlap (if each was split off of
	#segments to the left and right).  They would have been split because
	#their final coords (not their identity coords) overlap.  This is the
	#only case where we will keep overlapping segments as they are, because
	#we need to make sure we keep any identity coordinate that does not
	#have overlap.
	if(!(#These are both single codon segments
	     ($seg1->{FINAL_START} + 2) == $seg1->{FINAL_STOP} &&
	     ($seg2->{FINAL_START} + 2) == $seg2->{FINAL_STOP}))
	   {
	     if(#Both have no identity on opposite sides, i.e. both are
		#required to preserve all identity coordinates
		($seg1->{IDEN_START} != $seg1->{FINAL_START} &&
		 $seg2->{FINAL_STOP} != $seg2->{IDEN_STOP}) ||
		($seg1->{IDEN_STOP} != $seg1->{FINAL_STOP} &&
		 $seg2->{FINAL_START} != $seg2->{IDEN_START}))
	       {$seg2->{FINAL_STOP} = 0}

	     #Only the second segment is up for possibly skipping.  Must rely
	     #on the ordering of the records, so we don't need to test for seg2
	     #containing seg1.
	   }
        #Else - they're both single codon segments
        else
          {
            if(#Segment 1 completely contains the identity of segment 2
	       $seg1->{IDEN_START} <= $seg2->{IDEN_START} &&
	       $seg1->{IDEN_STOP}  >= $seg2->{IDEN_STOP})
	       {$seg2->{FINAL_STOP} = 0}
          }

	#This should be the only place I have to do this.  There should not
	#exist a case where the start codon of segment 2 overlaps the stop
	#codon of segment 1, where one or both are longer than a singkle codon,
	#because we're assuming that they were split earlier.
      }
    #Else if there's not really overlap, it might be one that was adjusted
    #before, so leave it alone
    elsif($seg2->{FINAL_STOP} < $seg1->{FINAL_START})
      {return()}
    #Else if the start of the segment being added is less than the stop...
    #NOTE: We can assume that overlap exists because segments are added in
    #order of ascending start (FIRST_CODON_START).  Even if the added segment
    #was adjusted to not technically overlap the segment being added, its start
    #could only have been moved to the right to accommodate other overlap which
    #is guaranteed to overlap this segment - so whether it's the current
    #already-added segment or one that caused the already-added segment to
    #adjust, we know there is overlap with something for seg2
    elsif($seg2->{FINAL_START} <= $seg1->{FINAL_STOP})
      {
	#If seg2/seg1 first codon starts & last codon stops overlap each other
	#Example: seg1: ***---------***
	#         seg2: ***---------***
	if($seg2->{FIRST_CODON_START} == $seg1->{FIRST_CODON_START} &&
	   $seg2->{LAST_CODON_STOP}   == $seg1->{LAST_CODON_STOP})
	  {
	    ##
	    ##Split if the identity starts or stops differ
	    ##

	    #If seg2 IDEN_START < seg1 IDEN_START - this should not happen, so
	    #throw error
	    if($seg2->{IDEN_START} < $seg1->{IDEN_START})
	      {
		error("Identity starts not being processed in ascending ",
		      "order.  Unable to add segment.");
		$seg2->{FINAL_STOP} = 0;
	      }

	    #Else if seg2 IDEN_STOP > seg1 IDEN_STOP - keep left side of seg1
	    #and right side of seg2 starting at the position after seg1's
	    #identity stop
	    elsif($seg2->{IDEN_STOP} > $seg1->{IDEN_STOP})
	      {
		#If there is current overlap
		if($seg2->{IDEN_START} <= $seg1->{IDEN_STOP})
		  {
		    $seg1->{FINAL_STOP}        = $seg1->{LAST_CODON_STOP} - 3;
		    $seg1->{IDEN_STOP}         = $seg1->{FINAL_STOP};
		    $seg1->{LAST_CODON_STOP}   = $seg1->{FINAL_STOP};
		    $seg2->{FINAL_START}       = $seg2->{LAST_CODON_STOP} - 2;
		    $seg2->{IDEN_START}        = $seg2->{FINAL_START};
		    $seg2->{FIRST_CODON_START} = $seg2->{FINAL_START};
		  }
	      }

	    #Else (total identity overlap) set seg2 FINAL_STOP to 0 to
	    #eliminate it
	    else
	      {$seg2->{FINAL_STOP} = 0}
	  }
	#else if seg2 first codon start overlaps seg1 first codon start
	#Example: seg1: ***---------***  OR: ***--------***
	#         seg2: ***---***            ***---------------***
	elsif($seg2->{FIRST_CODON_START} == $seg1->{FIRST_CODON_START})
	  {
	    ##
	    ## Split to keep the right side of seg2 (but check for errors)
	    ##

	    #If seg2 IDEN_START < seg1 IDEN_START - this should not happen, so
	    #throw error
	    if($seg2->{IDEN_START} < $seg1->{IDEN_START})
	      {
		error("Identity starts not being processed in ascending ",
		      "order.  Unable to add segment.");
		$seg2->{FINAL_STOP} = 0;
	      }

	    #Else if seg2 IDEN_STOP <= seg1 IDEN_STOP - total overlap: set seg2
	    #FINAL_STOP to 0 to eliminate it
	    elsif($seg2->{IDEN_STOP} <= $seg1->{IDEN_STOP})
	      {$seg2->{FINAL_STOP} = 0}

	    #Else if seg2 IDEN_STOP > seg1 IDEN_STOP - (last codons don't
	    #overlap) keep left side of seg1 and right side of seg2
	    elsif($seg2->{IDEN_STOP} > $seg1->{IDEN_STOP})
	      {
		#If there is current overlap
		if($seg2->{FINAL_START} <= $seg1->{FINAL_STOP})
		  {
		    $seg1->{FINAL_STOP}  = $seg1->{LAST_CODON_STOP} - 3;
		    $seg2->{FINAL_START} = $seg1->{LAST_CODON_STOP} - 2;
		  }
	      }

	    #Else (unknown case) - error
	    else
	      {
		error("Unknown logic case.");
		$seg2->{FINAL_STOP} = 0;
	      }
	  }
	#else if seg2 last  codon stop  overlaps seg1 last codon stop
	#Example: seg1: ***---------***  (assumes seg2 start can't be < seg1's)
	#         seg2:       ***---***
	elsif($seg2->{LAST_CODON_STOP} == $seg1->{LAST_CODON_STOP})
	  {
	    ##
	    ## Split if the identity stops differ
	    ##

	    #It is assumed that seg2's first codon is inside the non-edge codon
	    #identity region of seg1

	    #If seg2 IDEN_STOP > seg1 IDEN_STOP - split to keep left side of
	    #seg1 and right side of seg2
	    if($seg2->{IDEN_STOP} > $seg1->{IDEN_STOP})
	      {
		#If there is current overlap
		if($seg2->{FINAL_START} <= $seg1->{FINAL_STOP})
		  {
		    $seg1->{FINAL_STOP}  = $seg1->{LAST_CODON_STOP} - 3;
		    $seg2->{FINAL_START} = $seg1->{LAST_CODON_STOP} - 2;
		  }
	      }

	    #Else set seg2 FINAL_STOP to 0 to eliminate it
	    else
	      {$seg2->{FINAL_STOP} = 0}
	  }
	#else if seg2 first codon start overlaps seg1 last codon stop
	#Example: seg1: ***------***
	#         seg2:          ***------***
	elsif($seg2->{FIRST_CODON_START} == ($seg1->{LAST_CODON_STOP} - 2))
	  {
	    #We will assume that in any case, the codon will be fixed by an
	    #alternate codon if there's a conflict, so we'll just clean up the
	    #ends to not overlap

	    #If there is current overlap
	    if($seg2->{FINAL_START} <= $seg1->{FINAL_STOP})
	      {
		$seg2->{FINAL_START} = $seg2->{FIRST_CODON_START} + 3;
	      }
	  }
	#else if seg2 last  codon stop  overlaps seg1 first codon start
	#Example: seg1:          ***------***  (SHOULD NOT BE POSSIBLE)
	#         seg2: ***------***
	elsif(($seg2->{LAST_CODON_STOP} - 2) == $seg1->{FIRST_CODON_START})
	  {
	    #This should not happen, so error
	    error("Identity starts not being processed in ascending order.  ",
		  "Unable to add segment.");
	    $seg2->{FINAL_STOP} = 0;
	  }
	#else if seg2 first codon start overlaps seg1 non-edge codon
	#Example: seg1: ***------------*** OR: ***------***
	#         seg2:       ***---***              ***------***
	elsif($seg2->{FIRST_CODON_START} < ($seg1->{LAST_CODON_STOP}   - 2) &&
	      $seg2->{FIRST_CODON_START} > ($seg1->{FIRST_CODON_START} + 2))
	  {
	    #It is assumed that seg2's last codon is to the right of seg1's
	    #last codon stop

	    #If seg2's LAST_CODON_STOP <= seg1's LAST_CODON_STOP - set
	    #FINAL_STOP to 0 to eliminate it for total overlap (bec. we assume
	    #seg1's start is <= seg2's start)
	    if($seg2->{LAST_CODON_STOP} <= $seg1->{LAST_CODON_STOP})
	      {$seg2->{FINAL_STOP} = 0}

	    #Else if seg2's LAST_CODON_STOP > seg1's LAST_CODON_STOP - keep
	    #seg1's left side up to the IDEN_STOP and start seg2 at seg1's
	    #IDEN_STOP + 1
	    elsif($seg2->{LAST_CODON_STOP} > $seg1->{LAST_CODON_STOP})
	      {
		#If there is current overlap
		if($seg2->{FINAL_START} <= $seg1->{FINAL_STOP})
		  {
		    $seg1->{FINAL_STOP}  = $seg1->{LAST_CODON_STOP} - 3;
		    $seg2->{FINAL_START} = $seg1->{LAST_CODON_STOP} - 2;
		  }
	      }

	    #Else (unknown case) - error
	    else
	      {
		error("Unknown logic case.");
		$seg2->{FINAL_STOP} = 0;
	      }
	  }
	#else - no overlap (ignore - keep both/no change)
	#Example: seg1: ***------***
	#         seg2:             ***---***
      }
  }

#This returns a hash containing starts and stops and which pair and sequence
#each segment should come from from beginning to end.  It finds a codon
#boundary closest to the midpoint between segments to determine where to get
#all the sequence from.  It also contains any recoded edge codons and the
#coordinates they go to
sub getFilledMap
  {
    my $soln      = $_[0];
    my $data      = $_[1];
    my $part_hash = $_[2];

    verbose("Selecting identity segment midpoints...");

    #Locate all the midpoints between which the sequence sources change
    #This will save the direct sources.  Indirect sources will be undefined
    my $filled_map = setMapDividers($soln,$data,$part_hash);

    #Fills in the stop coordinates and removes the last bogus divider (for the
    #sequence end)
    fillMapStops($filled_map);

    #Now we are left with a map that is only missing how to determine where
    #recoding comes from.  That will be determined during sequence construction
    #/"weaving".
    return($filled_map);
  }

#This basically finds all the start coordinates of sections of the map for each
#sequence where the source of the sequence changes (given each sequence's
#different encodings from the different alignments decided on by
#codonHomologizer).  The start coordinates are the second key in the hash (the
#first being the sequence ID).  The value is a hash with the keys PAIR_ID,
#SEQ_ID, TYPE, and STOP.
sub setMapDividers
  {
    my $soln      = $_[0];
    my $data      = $_[1];
    my $part_hash = $_[2];

    my $div_map = setStaticMapDividers($soln);

    createLeftDividers($soln,$data,$part_hash,$div_map);

    createRightDividers($soln,$data,$part_hash,$div_map);

    return($div_map);
  }

#This creates 'left' aka 'SOURCE' dividers in the divider map.  It looks to the
#left of each identity segment to find the closest nearby segment in all of the
#alignments each sequence is involved in and finds a midpoint in frame 1 to
#place a divider.  A divider tells the map which alignment a portion of a
#sequence should be obtained from or whether it should be re-encoded.  Left
#dividers always come directly from an alignment, which is why they are also
#known as 'source' dividers.  This is because downstream of the divider is an
#identity segment from that alignment.
sub createLeftDividers
  {
    my $soln      = $_[0]; #Identity segment solution
    my $data      = $_[1]; #Original data hash created from the input files
    my $part_hash = $_[2]; #Partner hash
    my $div_map   = $_[3]; #Partially finished map where dividers will be saved

    foreach my $seqid (sort {$a cmp $b} keys(%$soln))
      {
	my $lastind = -1;
	my $cnt     = 0;
	my @ordered_indexes = sort {$soln->{$seqid}->[$a]->{FINAL_START} <=>
				      $soln->{$seqid}->[$b]->{FINAL_START}}
	  0..$#{$soln->{$seqid}};

	#For each segment record index (in order of ascending FINAL_START)
	foreach my $pseudoind (0..$#ordered_indexes)
	  {
	    my $recind  = $ordered_indexes[$pseudoind];
	    my $rec     = $soln->{$seqid}->[$recind];

	    if(!defined($rec->{FINAL_STOP}))
	      {error("FINAL_STOP NOT DEFINED IN THE FOLLOWING HASH: [",
		     join(', ',map {"$_ => $rec->{$_}"} keys(%$rec)),"].")}

	    my $msg = "Determining left-side midpoint boundaries of " .
	      "[$seqid:$rec->{FINAL_START}-$rec->{FINAL_STOP}]";
	    isDebug() ? verbose($msg) : verboseOverMe($msg);

	    #Obtain all the pair IDs and Sequence IDs of all the sequences that
	    #$seqid has been aligned with
	    my $partners = $part_hash->{$seqid};

	    #Using the partner coordinates relative to this FINAL_START,
	    #look for any segment to the left among the partner sequences
	    #this sequence has been paired with.  If any exists, use the
	    #closest one's FINAL_STOP to determine this sequence's relative
	    #coord.  Then determine a mid-point on a codon boundary and
	    #record it in this sequence's map and copy/convert it to all
	    #partner sequences's div_maps.

	    my($divider_left,$closest_left) =
	      getClosestLeftDivider($seqid,
				    $rec->{FINAL_START},
				    $soln,
				    $data,
				    $partners,
				    $div_map);

	    #Set the default as the beginning of the sequence just in case
	    if(!defined($divider_left))
	      {$divider_left = 1}

	    #Now let's check the divider that we're going to add for conflicts
	    #it may have with a pre-existing divider.
	    #If the divider extends more than a codon in length...  We're going
	    #to re-use the dividerExists code to search for conflicts.
	    #Existing dividers at the divider location or 1 after the stop of
	    #this segment are OK.  We want to find anything in between that
	    #would be a problem.  So if the divider is only
	    my $overlapping_divs = [];
	    if(($divider_left + 2) != $rec->{FINAL_STOP})
	      {
		@$overlapping_divs = (dividerExists($seqid,
						    $divider_left + 3,
						    $rec->{FINAL_STOP} - 2,
						    undef,
						    $rec->{PAIR_ID},
						    $div_map,
						    $data,
						    1,
						    1));

		if(scalar(@$overlapping_divs))
		  {debug("New conflicting dividers ([",
			 scalar(@$overlapping_divs),
			 "] of them) detected at [$seqid:$divider_left].")}
	      }

	    #This divider might have been copied from another sequence, so
	    #check it for a conflict before recording it
	    if(exists($div_map->{$seqid}) &&
	       exists($div_map->{$seqid}->{$divider_left}))
	      {
		if(defined($div_map->{$seqid}->{$divider_left}->{PAIR_ID}) &&
		   $div_map->{$seqid}->{$divider_left}->{PAIR_ID} ne
		   $rec->{PAIR_ID})
		  {debug("Sequence [$seqid] has multiple source alignments: [",
			 $div_map->{$seqid}->{$divider_left}->{PAIR_ID},
			 "] and [$rec->{PAIR_ID}].  Both alignments are ",
			 "marked as being the source for a divider starting ",
			 "left of an identity segment at sequence-relative ",
			 "position [$seqid:segment:$rec->{FINAL_START} ",
			 "divider:$divider_left].  Ignoring alignment ",
			 "[$rec->{PAIR_ID}].  This should be checked ",
			 "manually to confirm that their identity segments ",
			 "are encoded the same in both alignments.  Even ",
			 "though completely overlapping redundant identical ",
			 "segments were aribtrarily removed from a sequence ",
			 "that was aligned independently with 2 other ",
			 "sequences, the segments still exist in the ",
			 "partners, which will trigger the copying of the ",
			 "same divider (from different sources) to the ",
			 "common sequence record - and is OK - as long as ",
			 "it's confirmed that this is what's happening.  One ",
			 "possible complication is if the identity in one ",
			 "pair is larger than in the other pair.  The other ",
			 "pair's segment would have been eliminated (since ",
			 "it was completely overlapping) and its partner ",
			 "would copy a unique divider at its start and 1 ",
			 "after the stop, changing the source unnecessarily, ",
			 "but since they MUST be identical anyway, it only ",
			 "results in making it more complicated than it ",
			 "needs to be.",{LEVEL => 3})}

		#If this was a recode or undefined divider
		if(!defined($div_map->{$seqid}->{$divider_left}->{TYPE}) ||
		   $div_map->{$seqid}->{$divider_left}->{TYPE} eq 'RECODE')
		  {
		    $div_map->{$seqid}->{$divider_left} =
		      {PAIR_ID    => $rec->{PAIR_ID},
		       TYPE       => 'SOURCE',
		       STOP       => undef,
		       HARD_START => $rec->{FINAL_START},#To swap bad REDOCEs
		       HARD_STOP  => $rec->{FINAL_STOP},#To prevent bad REDOCEs
		       IDEN_START => $rec->{IDEN_START}, #Might be outside the
		       IDEN_STOP  => $rec->{IDEN_STOP},  #coords, but that's OK
		       SEQ        => ''};
		  }
		#Else if this was a copied source divider, fill it in
		elsif($div_map->{$seqid}->{$divider_left}->{IDEN_START} == 0)
		  {
		    $div_map->{$seqid}->{$divider_left}->{PAIR_ID} =
		      $rec->{PAIR_ID};
		    $div_map->{$seqid}->{$divider_left}->{IDEN_START} =
		      $rec->{IDEN_START};
		    $div_map->{$seqid}->{$divider_left}->{IDEN_STOP} =
		      $rec->{IDEN_STOP};
		    $div_map->{$seqid}->{$divider_left}->{HARD_START} =
		      $rec->{FINAL_START};
		    $div_map->{$seqid}->{$divider_left}->{HARD_STOP} =
		      $rec->{FINAL_STOP};
		  }
	      }
	    else
	      {
		debug({LEVEL => 4},
		      (exists($div_map->{$seqid}) &&
		       exists($div_map->{$seqid}->{$divider_left}) ?
		       'Overwriting' : "Creating")," source divider ",
		      "$divider_left left of segment at $rec->{FINAL_START} ",
		      "for $seqid from $rec->{PAIR_ID} with identity coords ",
		      "$rec->{IDEN_START}-$rec->{IDEN_STOP}");

		$div_map->{$seqid}->{$divider_left} =
		  {PAIR_ID    => $rec->{PAIR_ID},
		   TYPE       => 'SOURCE',
		   STOP       => undef,
		   HARD_START => $rec->{FINAL_START},#To swap bad REDOCE divs
		   HARD_STOP  => $rec->{FINAL_STOP},#To prevent bad REDOCE divs
		   IDEN_START => $rec->{IDEN_START}, #Might be outside the
		   IDEN_STOP  => $rec->{IDEN_STOP},  #coords, but that's OK
		   SEQ        => ''};
	      }

	    my $partner_seqid = getPartnerSeqID($seqid,$rec->{PAIR_ID},$data);

	    #If we created conflicts by creating this new source divider
	    if(scalar(@$overlapping_divs))
	      {
		foreach my $conf_div (@$overlapping_divs)
		  {getFixDividerConflicts($seqid,
					  $conf_div,
					  (defined($div_map->{$seqid}
						   ->{$conf_div}
						   ->{HARD_STOP}) ?
					   $div_map->{$seqid}->{$conf_div}
					   ->{HARD_STOP} : $conf_div),
					  $partner_seqid,
					  $rec->{PAIR_ID},
					  $div_map,
					  $data,
					  defined($div_map->{$seqid}
						  ->{$conf_div}->{HARD_STOP}),
					  1)}
	      }

	    #Now let's check the identity segment's source partner
	    $overlapping_divs = [];
	    if(($divider_left + 2) != $rec->{FINAL_STOP})
	      {
		@$overlapping_divs =
		  #Convert the divider coordinates to the partner
		  map {convertDivider($seqid,
				      $_,
				      $rec->{PAIR_ID},
				      $partner_seqid,
				      $data)}
		    (dividerExists($seqid,
				   $divider_left + 3,
				   $rec->{FINAL_STOP} - 2,
				   $partner_seqid,
				   $rec->{PAIR_ID},
				   $div_map,
				   $data,
				   1,
				   0));

		if(scalar(@$overlapping_divs))
		  {debug("New conflicting dividers ([",
			 scalar(@$overlapping_divs),
			 "] of them) detected at [$partner_seqid:",
			 convertDivider($seqid,
					$divider_left,
					$rec->{PAIR_ID},
					$partner_seqid,
					$data),"].")}
	      }

	    #Copy this divider to the partner sequences
	    copyDivider($divider_left,
			$seqid,
			$div_map,
			$partners,
			$data,
			$rec->{PAIR_ID},
			'left',
			$soln,
			$closest_left,
			$rec->{FINAL_START},
		        $rec->{FINAL_STOP}); #This will be the HARD STOP

	    #If we created conflicts by copying the new source divider to its
	    #source partner
	    if(scalar(@$overlapping_divs))
	      {
		foreach my $conf_div (@$overlapping_divs)
		  {getFixDividerConflicts($partner_seqid,
					  $conf_div,
					  (defined($div_map->{$partner_seqid}
						   ->{$conf_div}
						   ->{HARD_STOP}) ?
					   $div_map->{$partner_seqid}
					   ->{$conf_div}->{HARD_STOP} :
					   $conf_div),
					  $partner_seqid,
					  $rec->{PAIR_ID},
					  $div_map,
					  $data,
					  defined($div_map->{$partner_seqid}
						  ->{$conf_div}->{HARD_STOP}),
					  1)}
	      }

	    #Check to make sure that the very beginning of the sequence is
	    #accounted for, because a divider could have been copied over and
	    #interrupted the sequence from ever having a divider at the
	    #beginning
	    if($pseudoind == 0 && $divider_left != 1 &&
	       (!exists($div_map->{$seqid}) ||
		!exists($div_map->{$seqid}->{1})))
	      {
		debug({LEVEL => 4},"Creating UNDEF left edge divider 1 (in ",
		      "addition to divider [$divider_left]) left of segment ",
		      "at $rec->{FINAL_START} for $seqid from ",
		      "$rec->{PAIR_ID} with identity coords ",
		      "$rec->{IDEN_START}-$rec->{IDEN_STOP}");
		$div_map->{$seqid}->{1} =
		  {PAIR_ID    => $rec->{PAIR_ID},
		   TYPE       => undef,
		   STOP       => undef,
		   HARD_START => undef,
		   HARD_STOP  => undef,
		   IDEN_START => 0,
		   IDEN_STOP  => 0,
		   SEQ        => ''}
		#No need to copy it, because it's not real.  It's just a place-
		#holder that could get copied over.
	      }

	    $lastind = $recind;
	    $cnt++;
	  }
      }

    return(0);
  }

#This creates 'right' aka 'RECODE' (or simply undefined) dividers in the
#divider map.  It looks to the right of each identity segment to find the
#closest nearby segment in all of the alignments each sequence is involved in
#(which doesn't already have a divider) and finds a midpoint in frame 1 to
#place a divider.  A divider tells the map which alignment a portion of a
#sequence should be obtained from or whether it should be re-encoded.  Right
#dividers (unless a left divider is already there), is always recoded.  This is
#because there may be another divider downstream and no identity segment in
#between because the closest identity segment comes from an unrelated
#alignment, e.g. the closest identity segment to A's segment with B is a
#segment B has with C.
sub createRightDividers
  {
    my $soln      = $_[0];
    my $data      = $_[1];
    my $part_hash = $_[2];
    my $div_map   = $_[3];

    foreach my $seqid (sort {$a cmp $b} keys(%$soln))
      {
	my $lastind = -1;
	my $cnt     = 0;
	my @ordered_indexes = sort {$soln->{$seqid}->[$a]->{FINAL_START} <=>
				      $soln->{$seqid}->[$b]->{FINAL_START}}
	  0..$#{$soln->{$seqid}};

	#For each segment record index (in order of ascending FINAL_START)
	foreach my $pseudoind (0..$#ordered_indexes)
	  {
	    my $recind  = $ordered_indexes[$pseudoind];
	    my $rec     = $soln->{$seqid}->[$recind];

	    if(!defined($rec->{FINAL_STOP}))
	      {error("FINAL_STOP NOT DEFINED IN THE FOLLOWING HASH: [",
		     join(', ',map {"$_ => $rec->{$_}"} keys(%$rec)),"].")}

	    my $msg = "Determining right-side midpoint boundaries of " .
	      "[$seqid:$rec->{FINAL_START}-$rec->{FINAL_STOP}]";
	    isDebug() ? verbose($msg) : verboseOverMe($msg);

	    #Obtain all the pair IDs and Sequence IDs of all the sequences that
	    #$seqid has been aligned with
	    my $partners = $part_hash->{$seqid};

	    #Using the partner coordinates relative to this FINAL_STOP,
	    #look for a any segment to the right among the partner
	    #sequences this sequence has been paired with.  If any exists,
	    #use the closest one's FINAL_START to determine this sequence's
	    #relative coord.  Then determine a mid-point on a codon
	    #boundary and record it in this sequence's map and copy/convert
	    #it to all partner sequences's div_maps

	    my($divider_right,$closest_right) =
	      getClosestRightDivider($seqid,
				     $rec->{FINAL_STOP},
				     $soln,
				     $data,
				     $partners,
				     $div_map);

	    #If the returned divider is indefined and we're at the end of the
	    #ordered segments, set the default as just after the end of the
	    #sequence just in case
	    if(!defined($divider_right) && $pseudoind == $#ordered_indexes)
	      {$divider_right =
		 length($data->{$rec->{PAIR_ID}}->{SEQS}->{$seqid}) + 1}
	    elsif(!defined($divider_right))
	      {
		$lastind = $recind;
		$cnt++;
		next;
	      }
	    #This divider might have been copied from another sequence, so
	    #check it for a conflict before recording it
	    if(!exists($div_map->{$seqid}) ||
	       !exists($div_map->{$seqid}->{$divider_right}))
	      {
		debug({LEVEL => 4},
		      (exists($div_map->{$seqid}) &&
		       exists($div_map->{$seqid}->{$divider_right}) ?
		       'Overwriting' : "Creating"),
		      " right divider $divider_right ",
		      "right of segment at $rec->{FINAL_STOP} for $seqid ",
		      "from $rec->{PAIR_ID} with identity coords ",
		      "$rec->{IDEN_START}-$rec->{IDEN_STOP}");

		$div_map->{$seqid}->{$divider_right} =
		  {PAIR_ID    => $rec->{PAIR_ID},
		   IDEN_START => 0,
		   IDEN_STOP  => 0,
		   HARD_START => undef,
		   HARD_STOP  => undef,
		   TYPE       => undef,  #Don't know seg this belongs to yet,
		   STOP       => undef,  #so allow it to be over-written with
		   SEQ        => ''};    #another pair's copy by setting =undef
	      }

            #Copy this divider to the partner sequences
	    copyDivider($divider_right,
			$seqid,
			$div_map,
			$partners,
			$data,
			$rec->{PAIR_ID},
			'right',
			$soln,
			$rec->{FINAL_STOP} + 1,
			$closest_right,
		        undef);

	    #Check to make sure that the very end of the sequence is accounted
	    #for, because a divider could have been copied over and interrupted
	    #the sequence from ever having a divider at the end (or rather, one
	    #past the end - the way it was coded to work)
	    my $dr = length($data->{$rec->{PAIR_ID}}->{SEQS}->{$seqid}) + 1;
	    if($pseudoind == $#ordered_indexes &&
	       $divider_right != length($data->{$rec->{PAIR_ID}}->{SEQS}
					->{$seqid}) + 1 &&
	       !exists($div_map->{$seqid}->{$dr}))
	      {
		debug({LEVEL => 4},"Creating UNDEF right edge divider $dr ",
		      "right of segment at $rec->{FINAL_STOP} for $seqid ",
		      "from $rec->{PAIR_ID} with identity coords ",
		      "$rec->{IDEN_START}-$rec->{IDEN_STOP} because the last ",
		      "divider [$divider_right] was not the sequence length ",
		      "+ 1 [$dr].");

		$div_map->{$seqid}->{$dr} =
		  {PAIR_ID    => $rec->{PAIR_ID},
		   TYPE       => undef,
		   STOP       => undef,
		   HARD_START => undef,
		   HARD_STOP  => undef,
		   IDEN_START => 0,
		   IDEN_STOP  => 0,
		   SEQ        => ''};
		#No need to copy it, because it's not real.  It's just a place-
		#holder that could get copied over.
	      }

	    $lastind = $recind;
	    $cnt++;
	  }
      }

    return(0);
  }

#This method pre-populates the map with those dividers that do not need a
#midpoint calculation because they are contiguous neighbors with a segment
#piece to the left
sub setStaticMapDividers
  {
    my $soln      = $_[0];
    my $div_map   = {};

    foreach my $seqid (sort {$a cmp $b} keys(%$soln))
      {
	my @ordered_indexes = sort {$soln->{$seqid}->[$a]->{FINAL_START} <=>
				      $soln->{$seqid}->[$b]->{FINAL_START}}
	  0..$#{$soln->{$seqid}};

	my $stop_hash = {};

	#For each segment record index (in order of ascending FINAL_START)
	foreach my $pseudoind (0..$#ordered_indexes)
	  {
	    my $recind  = $ordered_indexes[$pseudoind];
	    my $rec     = $soln->{$seqid}->[$recind];

	    if(exists($stop_hash->{$rec->{FINAL_START} - 1}))
	      {
		$div_map->{$seqid}->{$rec->{FINAL_START}} =
		  {PAIR_ID    => $rec->{PAIR_ID},
		   TYPE       => 'SOURCE',
		   STOP       => undef,
		   HARD_START => $rec->{FINAL_START},#To swap bad RECODE divs
		   HARD_STOP  => $rec->{FINAL_STOP},#To prevent bad RECODE divs
		   IDEN_START => $rec->{IDEN_START}, #Might be outside the
		   IDEN_STOP  => $rec->{IDEN_STOP},  #coords, but that's OK
		   SEQ        => ''};
	      }

	    $stop_hash->{$rec->{FINAL_STOP}} = 0;
	  }
      }

    return($div_map);
  }

#Given a pair of sequence coordinates, this sub returns the midpoint coordinate
#that is as close to the middle as is possible and where both sequences have a
#non-gap character.  Finding an alignment midpoint is necessary because if I
#found a sequence-based midpoint for each sequence in an alignment, depending
#on the number and distribution of gaps, the midpoints calculated for each
#sequence and converted to the other sequence would vary wildly.
#This sub assumes that the coordinates submitted are each a valid return value
#(say, if the are adjacent or the same coordinate).  The coordinate returned is
#a sequence coordinate for the $seqid sequence.  It is assumed that since the
#returned value is for a position where both sequences do not have a gap, that
#convertDivider will behave consistently regardless of which direction the
#coordinate conversion is going.  It is also assumed that the midpoint
#calculation will be the same regardless of which sequence is represented by
#$seqid versus $partner_seqid
sub getAlnMidpoint
  {
    my $pairid        = $_[0];
    my $seqid         = $_[1];
    my $left_coord    = $_[2];
    my $right_coord   = $_[3];
    my $partner_seqid = $_[4];
    my $data          = $_[5];

    my $query_pos       = 0;
    my $aln_pos         = 0;
    my $left_aln_coord  = 0;
    my $right_aln_coord = 0;
    my $query_aln       = '';
    my $partner_aln     = '';
    my $midpoint        = 0;

    while($query_pos < $right_coord)
      {
	my $query_char = substr($data->{$pairid}->{ALNS}->{$seqid},
				$aln_pos,1);
	if($query_char ne '-')
	  {$query_pos++}

	my $partner_char = substr($data->{$pairid}->{ALNS}->{$partner_seqid},
				  $aln_pos,1);

	if($query_pos >= $left_coord && $query_pos <= $right_coord)
	  {
	    $query_aln   .= $query_char;
	    $partner_aln .= $partner_char;
	  }

	$aln_pos++;
      }

    my $mid_aln = int(length($query_aln) / 2);
    my $last_both = 1;
    my $last_both_aln = 1;
    $aln_pos = 0;
    $query_pos = 0;

    while($aln_pos < length($query_aln))
      {
	my $query_char   = substr($query_aln,  $aln_pos,1);
	my $partner_char = substr($partner_aln,$aln_pos,1);
	if($query_char ne '-')
	  {$query_pos++}
	$aln_pos++;
	if($aln_pos < $mid_aln && $query_char ne '-' && $partner_char ne '-')
	  {
	    $last_both     = $query_pos;
	    $last_both_aln = $aln_pos;
	  }
	if($aln_pos >= $mid_aln && $query_char ne '-' && $partner_char ne '-')
	  {
	    if(abs($aln_pos - $mid_aln) < abs($last_both_aln - $mid_aln))
	      {$midpoint = $query_pos}
	    else
	      {$midpoint = $last_both}
	    last;
	  }
      }

    $midpoint = $left_coord + $midpoint - 1;

    return($midpoint);
  }

#Searches sequences that have been paired with a particular sequence to find
#the closest segment to its left and returns a midpoint on a codon boundary
#(frame 1)
sub getClosestLeftDivider
  {
    my $seqid    = $_[0];
    my $coord    = $_[1];
    my $soln     = $_[2];
    my $data     = $_[3];
    my $partners = $_[4];
    my $divs     = $_[5];

    #For each alignment $seqid has been aligned with
    my $closest = 0;
    my $closest_partner = $partners->[0];
    foreach my $partner (@$partners)
      {
	if($seqid eq $partner->[1])
	  {
	    error("Encountered sequence [$seqid] aligned with itself.",
		  {DETAIL => ("NOTE: The partners array may have been " .
			      "improperly constructed.")});
	    next;
	  }

	#Convert the coordinate to the partner sequence
	my $partner_coord = convertDivider($seqid,
					   $coord,
					   $partner->[0],  #pair ID
					   $partner->[1],  #partner seq ID
					   $data);

	#Find the closest coord to the left
	my $part_left_seg_coord = getClosestLeftSegCoord($partner->[1],
							 $partner_coord,
							 $soln);

	#Convert the coordinate back to this sequence
	my $orig_left_seg_coord = convertDivider($partner->[1],
						 $part_left_seg_coord,
						 $partner->[0],  #pair ID
						 $seqid,
						 $data);

	#This is intended to correct for end-gaps in the alignment and is
	#designed around the case where the query coord ($coord) is 1.  When
	#it's converted to the partner and aligns to something other than the
	#first base, the convert divider method is designed to return its first
	#base (accounting for the starting gap in the partner (i.e. all the
	#bases in the query sequence aligning with the starting gap return base
	#1 in the partner)).  But, in that scenario, when you convert back, it
	#returns the actual aligned base from the query sequence (e.g. base
	#16).  This catches that case and changes it back to 1.
	if($partner_coord == $part_left_seg_coord &&
	   $orig_left_seg_coord > $coord)
	  {$orig_left_seg_coord = $coord}

	#If this coordinate is closer and not overlapping, save it.
	if($orig_left_seg_coord > $closest)
	  {
	    debug({LEVEL => 4},"Closest left segment to [@$partner ",
		  "$partner_coord] is [@$partner $part_left_seg_coord] and ",
		  "converted back to [$partner->[0] $seqid ",
		  "$orig_left_seg_coord]");

	    #If the divider is inside the query coordinate, something went
	    #wrong.  (If it's the same coordinate, it's OK.  We'll just use the
	    #start of this identity segment as the location of the divider -
	    #though there may be an issue if this turns out to not be a valid
	    #codon boundary (i.e. frame 1).)
	    if($orig_left_seg_coord > $coord)
	      {
		error("Returned coord [$orig_left_seg_coord] overlaps query ",
		      "coord [$coord] in seq [$seqid] on its left side.  ",
		      "Obtained from a segment in [$partner->[1]] and ",
		      "alignment [$partner->[0]].  The query coord was ",
		      "converted from [$seqid:$coord] to [$partner->[1]:",
		      "$partner_coord] and the closest coord in that ",
		      "sequence was: [$partner->[1]:$part_left_seg_coord] - ",
		      "which should not overlap, but did overlap when ",
		      "converted back to [$seqid:$orig_left_seg_coord] ",
		      "coord.");
	      }
	    $closest = $orig_left_seg_coord;
	    $closest_partner = [@$partner,$part_left_seg_coord];
	    if($orig_left_seg_coord == $coord)
	      {return($coord,$coord)}
	  }
	else
	  {
	    debug({LEVEL => 4},"Further left segment to [@$partner ",
		  "$partner_coord] is [@$partner $part_left_seg_coord] and ",
		  "converted back to [$partner->[0] $seqid ",
		  "$orig_left_seg_coord]");
	  }
      }

    #This means that a segment abutts this one either directly or indirectly,
    #but in any case, the divider can't be anywhere else
    if($closest == $coord)
      {return($coord,$coord)}

    #See if there already exists a divider (or dividers) in this region in any
    #of the partners (or self).
    #Dividers with slightly different coordinates are common because their
    #position is calculated using different alignments and the alignments can
    #be very gappy.

    getFixDividerConflicts($seqid,
			   $closest,
			   $coord,
			   $closest_partner->[1], #SeqID
			   $closest_partner->[0], #PairID
			   $divs,
			   $data,
			   1,
			   1); #Looking in self/$seqid

    #Now check my partners
    #The divider we're going to create is going to be copied to all of $seqid's
    #partners, so we need to pre-emptively check for conflicts in the range in
    #which we're considering to create a divider.  Any source divider that
    #extends into that range is considered a conflict.
    foreach my $partner (@$partners)
      {
	getFixDividerConflicts($seqid,
			       $closest,
			       $coord,
			       $partner->[1],
			       $partner->[0],
			       $divs,
			       $data,
			       1, #Looking left to create a source div in seqid
			       0, #Looking in partner (seqid in $partner->[1])
			       $closest_partner->[1], #SeqID
			       $closest_partner->[0]);
      }

    #First check myself...
    my $self_divider =
      dividerExists($seqid,
		    $closest,
		    $coord,
		    $closest_partner->[1], #SeqID
		    $closest_partner->[0],
		    $divs,
		    $data,
		    1,
		    1);

    if($self_divider)
      {return($self_divider,$closest)}

    #Now check my partners
    my($right_most_divider);
    foreach my $partner (@$partners)
      {
	#This should be converted back
	my $converted_self_divider = dividerExists($seqid,
						   $closest,
						   $coord,
						   $partner->[1],
						   $partner->[0],
						   $divs,
						   $data,
						   1,
						   0);

	if(defined($converted_self_divider) && $converted_self_divider &&
	   (!defined($right_most_divider) ||
	    $converted_self_divider > $right_most_divider))
	  {$right_most_divider = $converted_self_divider}
      }

    if(defined($right_most_divider) && $right_most_divider)
	  {return($right_most_divider,$closest)}

    #If this is the left-most edge of the query sequence, then we can return
    #it, otherwise, we mist find a mid-way dividing point between this identity
    #segment and the nearest one on the left
    if($closest == 1)
      {return($closest,$closest)}

    #We now have the closest occupied coordinate to the left of $coord.  We
    #need to split the difference and then return the closest frame 1 codon
    #position.  We add 1 to the closest coord because we want to choose a
    #midpoint among available positions to put a divider.  That includes the
    #original $coord position, but not the last position of the nearest
    #neighbor to the left ($closest)
    debug({LEVEL => 4},"Calculating midpoint using alignment ",
	  "[$closest_partner->[0]]");
    my $midpoint = getAlnMidpoint($closest_partner->[0],  #pairID
				  $seqid,
				  $closest + 1,
				  $coord,
				  $closest_partner->[1],  #seqID
				  $data);

    #Determine the frame position (1, 2, or 3) of the coordinate
    my $frame_pos = ($midpoint - 1) % 3 + 1;

    if($frame_pos == 1)
      {return($midpoint,$closest)}
    elsif($frame_pos == 2)
      {
	if(($midpoint - 1) <= $closest)
	  {
	    #This should not happen, but just in case...
	    if(($midpoint + 2) > $coord)
	      {
		error("Bad frame boundary.  Frame 2.  Midpoint: ",
		      "[$midpoint].  Closest: [$closest].  Coord(/final ",
		      "start): [$coord].");
		return($coord,$closest);
	      }

	    return($midpoint + 2,$closest);
	  }

	return($midpoint - 1,$closest);
      }
    elsif($frame_pos == 3)
      {
	if(($midpoint + 1) > $coord)
	  {
	    error("Bad frame boundary for left divider.  Frame 3.  Midpoint: ",
		  "[$midpoint].  Closest: [$closest].  Coord(/final start): ",
		  "[$coord].");
	    return($coord,$closest);
	  }

	return($midpoint + 1,$closest);
      }
    else
      {error("Unexpected frame result.")}

    return($coord,$closest);
  }

sub getFixDividerConflicts
  {
    my $seqid         = $_[0];
    my $lesser_bound  = $_[1];
    my $greater_bound = $_[2];
    my $partner_seqid = $_[3];
    my $pair_id       = $_[4];
    my $div_map       = $_[5];
    my $data          = $_[6];
    my $looking_left  = $_[7]; #looking left = creating source divider in seqid
    my $look_in_self  = $_[8];
    my $from_seqid    = $_[9]; #This is the partner from which the closest
    my $from_pairid   = $_[10];#segment came from

    my $success = 1;

    #Check for divider conflicts and attempot to fix them, then re-check.
    #There is likely only to be a single conflict as long as we stay on top of
    #things, but I've coded this for the possibility that there could be
    #multiple conflicts
    my $tried = {};
    while(my $bad_div = getDividerConflict($seqid,
					   $lesser_bound,
					   $greater_bound,
					   $partner_seqid,
					   $pair_id,
					   $div_map,
					   $data,
					   $looking_left,#For ordering conf div
					   $look_in_self))
      {
	#The bad divider will probably be defined.  I don't think the loop
	#would be entered if it was not defined.
	if(!defined($bad_div) ||
	   exists($tried->{"$bad_div->{PAIR_ID}:$bad_div->{SEQID}:" .
			   "$bad_div->{DIVIDER}"}))
	  {
	    if(defined($bad_div))
	      {
		my $other_seqid =
		  getPartnerSeqID($bad_div->{SEQID},$bad_div->{PAIR_ID},$data);
	       error("Unable to fix ",($look_in_self ? 'self' : 'partner'),
		     " divider conflict: [$bad_div->{SEQID}:",
		     "$bad_div->{DIVIDER}] that ends at [",
		     $div_map->{$bad_div->{SEQID}}->{$bad_div->{DIVIDER}}
		     ->{HARD_START},"-",
		     $div_map->{$bad_div->{SEQID}}->{$bad_div->{DIVIDER}}
		     ->{HARD_STOP},"] and overlaps self-range [$seqid:",
		     "$lesser_bound-$greater_bound] when trying to create a ",
		     "[$seqid:",
		     ($look_in_self ? ($looking_left ? 'source' : 'recode') :
		      'unknown'),'] divider while looking [',
		     ($looking_left ? 'left' : 'right'),'], which interrupts ',
		     'a [',"$bad_div->{SEQID}:$bad_div->{DIVIDER}-",
		     $div_map->{$bad_div->{SEQID}}->{$bad_div->{DIVIDER}}
		     ->{HARD_STOP},":",
		     (defined($div_map->{$bad_div->{SEQID}}
			      ->{$bad_div->{DIVIDER}}->{TYPE}) ?
		      $div_map->{$bad_div->{SEQID}}->{$bad_div->{DIVIDER}}
		      ->{TYPE} : 'recode'),"] divider.  The search range is ",
		     "based on bounding segments in [$seqid] and [",
		     (defined($from_seqid) ? "$from_seqid in $from_pairid" :
		      'unknown'),"].",
		     {DETAIL => "Alignment that determines where the " .
		      "conflicting divider comes from:\n" .
		      getAlnSnippet($bad_div->{SEQID},
				    $bad_div->{DIVIDER},
				    $div_map->{$bad_div->{SEQID}}
				    ->{$bad_div->{DIVIDER}}
				    ->{HARD_STOP},
				    $bad_div->{PAIR_ID},
				    $other_seqid,
				    $data,
				    $div_map->{$bad_div->{SEQID}}
				    ->{$bad_div->{DIVIDER}}->{HARD_START},
				    $div_map->{$bad_div->{SEQID}}
				    ->{$bad_div->{DIVIDER}}->{HARD_STOP},
				    convertDivider($bad_div->{SEQID},
						   $div_map->{$bad_div
							      ->{SEQID}}
						   ->{$bad_div->{DIVIDER}}
						   ->{HARD_START},
						   $bad_div->{PAIR_ID},
						   $other_seqid,
						   $data),
				    convertAlnSeqCoord($bad_div->{SEQID},
						       $div_map->{$bad_div
								  ->{SEQID}}
						       ->{$bad_div->{DIVIDER}}
						       ->{HARD_STOP},
						       $bad_div->{PAIR_ID},
						       $other_seqid,
						       $data)) . "\n\n" .
		      "Alignment showing the partner where the conflicting " .
		      "divider will be copied:\n" .
		      getAlnSnippet($seqid,
				    $lesser_bound,
				    $greater_bound,
				    $pair_id,
				    $partner_seqid,
				    $data,
				    $lesser_bound,
				    $greater_bound,
				    convertDivider($seqid,
						   $lesser_bound,
						   $pair_id,
						   $partner_seqid,
						   $data),
				    convertDivider($seqid,
						   $greater_bound,
						   $pair_id,
						   $partner_seqid,
						   $data))
		     })
	     }
	    $success = 0;
	    last;
	  }

	debug("Bad dividers seem to be propagating [$bad_div->{PAIR_ID}:",
	      "$bad_div->{SEQID}:$bad_div->{DIVIDER}].")
	  if(scalar(keys(%$tried)));

	#We're going to send along the existing divider if one exists, to make
	#sure the fix addresses it specifically.
	fixDividerConflict($seqid,
			   $lesser_bound,
			   $greater_bound,
			   $partner_seqid,
			   $pair_id,
			   $div_map,
			   $data,
			   $look_in_self, #Whether we're fixing self or partner
			   $looking_left,
			   $bad_div);

	$tried->{"$bad_div->{PAIR_ID}:$bad_div->{SEQID}:$bad_div->{DIVIDER}"} =
	  $bad_div;
      }

    if(scalar(keys(%$tried)))
      {
	my $msg =
	  join('',(($success ? "Fixed" : "Failed to fix")," a ",
		   ($look_in_self ? 'self' : 'partner')," divider conflict: [",
		   join(',',keys(%$tried)),"] [",
		   join(',',map {defined($div_map->{$_->{SEQID}}
					 ->{$_->{DIVIDER}}->{HARD_STOP}) ?
					   "$_->{SEQID}:$_->{DIVIDER}-" .
					     $div_map->{$_->{SEQID}}
					       ->{$_->{DIVIDER}}->{HARD_START}
						 . "-" .
						   $div_map->{$_->{SEQID}}
						     ->{$_->{DIVIDER}}
						       ->{HARD_STOP} :
							 'changed to recode'}
			values(%$tried)),"]!"));
	$success ? debug($msg) : warning($msg);
      }
    else
      {debug("No conflicts were found in search range: [$seqid:$lesser_bound-",
	     "$greater_bound] with partner [$partner_seqid] from alignment ",
	     "[$pair_id].")}

    return($success);
  }

#This will fix conflicts between the range in which a divider will be created
#(or a pre-selected divider position) and an existing divider which overlaps
#the proposed divider area or specific divider location
sub fixDividerConflict
  {
    my $seqid         = $_[0];
    my $lesser_bound  = $_[1];
    my $greater_bound = $_[2];
    my $partner_seqid = $_[3];
    my $pair_id       = $_[4];
    my $div_map       = $_[5];
    my $data          = $_[6];
    my $look_in_self  = $_[7];
    my $looking_left  = $_[8];
    my $bad_div_obj   = $_[9]; #hash with keys: SEQID, DIVIDER, & PAIR_ID

    #We are assuming that only source dividers (i.e. with hard starts/stops)
    #will get in here, because getDividerConflict only finds conflicts for
    #these dividers.  But, we should make sure that the divider conflict object
    #has a value
    if(!defined($bad_div_obj))
      {return()}

    my $old_divider     = $bad_div_obj->{DIVIDER};
    my $new_divider     = $bad_div_obj->{DIVIDER};
    my $conflict_seqid  = $bad_div_obj->{SEQID};
    my $conflict_pairid = $bad_div_obj->{PAIR_ID};
    my $hard_start      = $div_map->{$conflict_seqid}->{$old_divider}
      ->{HARD_START};
    my $hard_stop       = $div_map->{$conflict_seqid}->{$old_divider}
      ->{HARD_STOP};
    my $conf_part_seqid =
      getPartnerSeqID($conflict_seqid,$conflict_pairid,$data);
    my $old_iden_start  = $div_map->{$conflict_seqid}->{$old_divider}
      ->{IDEN_START};
    my $old_iden_stop   = $div_map->{$conflict_seqid}->{$old_divider}
      ->{IDEN_STOP};

    my $partner_updated = 0;

    #We're fixing a divider in the partner, which may or may not be a source
    #divider
    if(!$look_in_self)
      {
	#This block is the only reason for passing the partner_seqid

	$lesser_bound  = convertDivider($seqid,
					$lesser_bound,
					$pair_id,
					$partner_seqid,
					$data);
	$greater_bound = convertDivider($seqid,
					$greater_bound,
					$pair_id,
					$partner_seqid,
					$data);
	$seqid = $partner_seqid;

	#The divider in the bad_div_obj should be the partner
	if($bad_div_obj->{SEQID} ne $seqid)
	  {error("Bad sequence ID in the bad divider object.  It doesn't ",
		 "match the partner sequence ID.")}

	if($div_map->{$conflict_seqid}->{$old_divider}->{PAIR_ID} ne
	   $conflict_pairid)
	  {
	    error("The pair ID that getDividerConflict put in the ",
		  "bad divider object does not match what was found ",
		  "in the divider map while processing the partner.");
	  }
      }
    elsif($div_map->{$conflict_seqid}->{$old_divider}->{PAIR_ID} ne
	  $conflict_pairid)
      {
	error("The pair ID that getDividerConflict put in the ",
	      "bad divider object does not match what was found ",
	      "in the divider map while processing self.");
      }

    my $conflict_with_self = 0;
    #If the conflicting source divider we found and the source divider we
    #are seeking to create in the search region are each being created
    #while looking left of the same segment from the same alignment, then
    #it's a "conflict with self".  When the first version of the divider
    #was created, it did not find a closer segment and created a divider
    #far off to the left.  Then, while looking left of the same segment in
    #the partner, we find a closer segment with a sequence with which the
    #first sequence was never aligned.  This must be handled differently
    #from a conflict from a different pair.
    #I shouldn't need to check that hard start == greater bound, but it doesn't
    #hurt.
    if($looking_left && $conflict_pairid eq $pair_id &&
       $hard_start == $greater_bound)
      {
	debug("Conflict with self detected in [$pair_id].  Conflicting ",
	      "divider: [$conflict_seqid:$old_divider-$hard_start-",
	      "$hard_stop].  Looking in self [$seqid] for a divider in the ",
	      "range: [$seqid:$lesser_bound-$greater_bound].  The greater ",
	      "bound should equal the hard start: [$greater_bound == ",
	      "$hard_start ? ",($greater_bound == $hard_start ? 'YES' : 'NO'),
	      "].");
	$conflict_with_self = 1;
      }

    #Error-check that the divider is in the same sequence as the one with the
    #conflicting search range
    if($seqid ne $bad_div_obj->{SEQID})
      {
	error("Conflicting divider is from a different sequence ",
	      "[$bad_div_obj->{SEQID}] than the one we are searching for ",
	      "conflicts in [$seqid].");
	return();
      }

    #1. If the HARD START of the conflicting divider occurs after the
    #   conflicting search range start
    if($hard_start > $lesser_bound)
      {
	#If the divider is not already at the hard start
	if($old_divider != $hard_start)
	  {
	    ##
	    ## Move the conflicting divider right-ward
	    ##

	    #1.1. If the hard start <= the end of the search range
	    if($hard_start <= $greater_bound)
	      {
                #If this is a conflict with self
                if($conflict_with_self)
                  {
		    #Mark the divider to be changed to undef/recode
		    $new_divider = undef;
		  }
		#Make sure we're not copying over something (which should not
		#happen)
		elsif(!exists($div_map->{$bad_div_obj->{SEQID}}
			      ->{$hard_start}))
		  {
		    #1.1.1. Copy the conflicting divider to its hard start
		    $new_divider = $hard_start;
		    $div_map->{$bad_div_obj->{SEQID}}->{$new_divider} =
		      {%{$div_map->{$bad_div_obj->{SEQID}}->{$old_divider}}};
		  }
		elsif($div_map->{$bad_div_obj->{SEQID}}->{$old_divider}
		      ->{PAIR_ID} ne
		      $div_map->{$bad_div_obj->{SEQID}}->{$hard_start}
		      ->{PAIR_ID})
		  {error("Unexpected pre-existing divider: ",
			 "[$bad_div_obj->{SEQID}:$hard_start].  It's a [",
			 (defined($div_map->{$bad_div_obj->{SEQID}}
				  ->{$hard_start}->{TYPE}) ?
			  $div_map->{$bad_div_obj->{SEQID}}->{$hard_start}
			  ->{TYPE} : 'recode'),"] divider from alignment [",
			 $div_map->{$bad_div_obj->{SEQID}}->{$hard_start}
			 ->{PAIR_ID},"].  Unable to copy [",
			 (defined($div_map->{$bad_div_obj->{SEQID}}
				  ->{$old_divider}->{TYPE}) ?
			  $div_map->{$bad_div_obj->{SEQID}}->{$old_divider}
			  ->{TYPE} : 'recode'),"] divider ",
			 "[$bad_div_obj->{SEQID}:$old_divider] from [",
			 $div_map->{$bad_div_obj->{SEQID}}->{$old_divider}
			 ->{PAIR_ID},"] over it.")}
                #It's still considered "changed" even though we didn't actually
                #need to move it.  This is so it will get edited to not overlap
                #anymore.
                else
		  {$new_divider = $hard_start}
	      }
	    #1.2. Else (the hard start > the end of the search range)
	    else
	      {
		#Make sure the end of the search range is in frame 1
		if(calcFrame($greater_bound) != 1)
		  {error("Larger search bound must be in frame 1.")}
		#Make sure we're not copying over something (which should not
		#happen)
		elsif(exists($div_map->{$bad_div_obj->{SEQID}}
			     ->{$greater_bound}) &&
		      $div_map->{$bad_div_obj->{SEQID}}->{$old_divider}
		      ->{PAIR_ID} ne
		      $div_map->{$bad_div_obj->{SEQID}}->{$greater_bound}
		      ->{PAIR_ID})
		  {error("Unexpected pre-existing divider: ",
			 "[$bad_div_obj->{SEQID}:$greater_bound].  It's a [",
			 (defined($div_map->{$bad_div_obj->{SEQID}}
				  ->{$greater_bound}->{TYPE}) ?
			  $div_map->{$bad_div_obj->{SEQID}}->{$greater_bound}
			  ->{TYPE} : 'recode'),"] divider from alignment [",
			 $div_map->{$bad_div_obj->{SEQID}}->{$greater_bound}
			 ->{PAIR_ID},"].  Unable to copy [",
			 (defined($div_map->{$bad_div_obj->{SEQID}}
				  ->{$old_divider}->{TYPE}) ?
			  $div_map->{$bad_div_obj->{SEQID}}->{$old_divider}
			  ->{TYPE} : 'recode'),"] divider ",
			 "[$bad_div_obj->{SEQID}:$old_divider] from [",
			 $div_map->{$bad_div_obj->{SEQID}}->{$old_divider}
			 ->{PAIR_ID},"] over it.")}
		elsif(!exists($div_map->{$bad_div_obj->{SEQID}}
			      ->{$greater_bound}))
		  {
		    #1.2.1. Copy the conflicting divider to the end of the
		    #       search range (making sure it's in frame 1 of
		    #       course)
		    $new_divider = $greater_bound;
		    $div_map->{$bad_div_obj->{SEQID}}->{$new_divider} =
		      {%{$div_map->{$bad_div_obj->{SEQID}}->{$old_divider}}};
		  }
                #It's still considered "changed" even though we didn't actually
                #need to move it.  This is so it will get edited to not overlap
                #anymore.
                else
		  {$new_divider = $hard_start}
	      }
	  }

	##
	## Re-use the location of the previously conflicting divider
	##

	#1.3. If the (previously) conflicting divider started at or after the
	#     start of the search range
	if($old_divider >= $lesser_bound)
	  {
	    #1.3.1. If the old divider was copied to a new location above
	    if(defined($new_divider) && $old_divider != $new_divider)
	      {
		#I check new_divider's defined state above, but given earlier
		#code, it should always be defined here.  Conflicts with self
		#(which cause the undefined new divider only occur here when
		#the old divider is less than or equal to the lesser bound

                #This is called to simply check that $partner_pair_id is valid
	        getPartnerSeqID($bad_div_obj->{SEQID},$pair_id,$data);

		#1.3.1.1. Change the original conflicting divider to undef/
		#         recode using the pair_id that is the source for the
		#         search range.  It will be selected as an existing
		#         divider anyway, if it's *inside* the search range.
		$div_map->{$bad_div_obj->{SEQID}}->{$old_divider}->{TYPE} =
		  undef;
		$div_map->{$bad_div_obj->{SEQID}}->{$old_divider}->{PAIR_ID} =
		  $pair_id;
		$div_map->{$bad_div_obj->{SEQID}}->{$old_divider}
		  ->{HARD_START} = undef;
		$div_map->{$bad_div_obj->{SEQID}}->{$old_divider}
		  ->{HARD_STOP} = undef;
		$div_map->{$bad_div_obj->{SEQID}}->{$old_divider}
		  ->{IDEN_START} = 0;
		$div_map->{$bad_div_obj->{SEQID}}->{$old_divider}
		  ->{IDEN_STOP} = 0;

		#Add/edit the conflicting source divider's partner divider
		updatePartnerSourceDivider($bad_div_obj->{SEQID},
					   $old_divider,
					   $conflict_pairid,
					   $conf_part_seqid,
					   $data,
					   $div_map,
					   $div_map->{$bad_div_obj->{SEQID}}
					   ->{$old_divider});

                $partner_updated = 1;
	      }
	    #1.3.2. Else nothing (leave the original divider alone)
	  }
	#1.4. Else (the previously conflicting divider started before the start
	#     of the search range)
	else
	  {
	    #If the divider changed
	    if(!defined($new_divider) || $old_divider != $new_divider)
	      {
		#1.4.1. Change it to a recode divider using the source
		#       alignment it already had
		$div_map->{$bad_div_obj->{SEQID}}->{$old_divider}->{TYPE} =
		  undef;
		$div_map->{$bad_div_obj->{SEQID}}->{$old_divider}
		  ->{HARD_START} = undef;
		$div_map->{$bad_div_obj->{SEQID}}->{$old_divider}
		  ->{HARD_STOP} = undef;
		$div_map->{$bad_div_obj->{SEQID}}->{$old_divider}
		  ->{IDEN_START} = 0;
		$div_map->{$bad_div_obj->{SEQID}}->{$old_divider}
		  ->{IDEN_STOP} = 0;

		#Add/edit the conflicting source divider's partner divider
		updatePartnerSourceDivider($bad_div_obj->{SEQID},
					   $old_divider,
					   $conflict_pairid,
					   $conf_part_seqid,
					   $data,
					   $div_map,
					   $div_map->{$bad_div_obj->{SEQID}}
					   ->{$old_divider});

                $partner_updated = 1;

                #Just let a new divider be selected in the search range and let
                #it copy back to the partner as it would have had it been found
                #first.  We get here when the search range, looking left of a
                #partner is smaller than looking left of the original due to a
                #segment the partner has found that's closer that the original
                #never aligned with
                if(!defined($new_divider))
		  {return()}
	      }

	    #1.4.2. Duplicate the source divider to every (hard start to hard
	    #       stop) overlapping frame 1 position within the search range
	    #       to protect those positions
	    for(#Start splitting at the point of overlap
		my $src_div = ($hard_start > $lesser_bound ?
			       $hard_start : $lesser_bound);

		#Only split at points where there is overlap
		$src_div >= $lesser_bound && $src_div <= $greater_bound &&
		$src_div < $hard_stop;

		#Split only in frame 1
		$src_div += 3)
	      {
		if($src_div == $greater_bound && $hard_stop > ($src_div + 2))
		  {debug("Fixing the end of conflict ",
			 "[$bad_div_obj->{SEQID}:$src_div].")}

		#If we're not splitting all the way to the hard stop, make sure
		#the last coordinates used go to the original end of the
		#segment
		my $new_hard_stop = ($src_div == $greater_bound &&
				     $hard_stop > ($src_div + 2) ?
				     $hard_stop : $src_div + 2);

		if(exists($div_map->{$bad_div_obj->{SEQID}}->{$src_div}))
		  {
		    if($new_divider == $src_div)
		      {
			if(defined($div_map->{$bad_div_obj->{SEQID}}
				   ->{$src_div}->{TYPE}) &&
			   $div_map->{$bad_div_obj->{SEQID}}->{$src_div}
			   ->{TYPE} eq 'SOURCE')
			  {
			    $div_map->{$bad_div_obj->{SEQID}}->{$src_div}
			      ->{HARD_STOP} = $new_hard_stop;

			    #If the iden stop needs to be updated as well
			    if($div_map->{$bad_div_obj->{SEQID}}->{$src_div}
			       ->{IDEN_STOP} > 0 &&
			       $div_map->{$bad_div_obj->{SEQID}}->{$src_div}
			       ->{IDEN_STOP} >
			       $div_map->{$bad_div_obj->{SEQID}}
			       ->{$src_div}->{HARD_STOP})
			      {$div_map->{$bad_div_obj->{SEQID}}->{$src_div}
				 ->{IDEN_STOP} =
				   $div_map->{$bad_div_obj->{SEQID}}
				     ->{$src_div}->{HARD_STOP}}
			  }
			else
			  {
			    debug("Encountered a recode divider ",
				  "overlapping a source divider.");
			    #Replace the recode divider with a source divider
			    $div_map->{$bad_div_obj->{SEQID}}->{$src_div} =
			      {%{$div_map->{$bad_div_obj->{SEQID}}
				   ->{$old_divider}}};
			    $div_map->{$bad_div_obj->{SEQID}}->{$src_div}
			      ->{HARD_STOP} = $new_hard_stop;

			    #If the iden stop needs to be updated as well
			    if($div_map->{$bad_div_obj->{SEQID}}->{$src_div}
			       ->{IDEN_STOP} > 0 &&
			       $div_map->{$bad_div_obj->{SEQID}}->{$src_div}
			       ->{IDEN_STOP} >
			       $div_map->{$bad_div_obj->{SEQID}}
			       ->{$src_div}->{HARD_STOP})
			      {$div_map->{$bad_div_obj->{SEQID}}->{$src_div}
				 ->{IDEN_STOP} =
				   $div_map->{$bad_div_obj->{SEQID}}
				     ->{$src_div}->{HARD_STOP}}
			  }
		      }
		    else
		      {
			my $msg =
			  join('',("Encountered a ",
				   (defined($div_map->{$bad_div_obj->{SEQID}}
					    ->{$src_div}->{TYPE}) &&
				    $div_map->{$bad_div_obj->{SEQID}}
				    ->{$src_div}->{TYPE} eq 'SOURCE' ?
				    'source' : 'recode')," divider ",
				   "[$bad_div_obj->{SEQID}:$src_div] ",
				   "overlapping a source divider from ",
				   ($div_map->{$bad_div_obj->{SEQID}}
				    ->{$src_div}->{PAIR_ID} eq
				    $div_map->{$bad_div_obj->{SEQID}}
				    ->{$old_divider}->{PAIR_ID} ?
				    'the same' : 'a different')," alignment [",
				   $div_map->{$bad_div_obj->{SEQID}}
				   ->{$src_div}->{PAIR_ID},"] vs [",
				   $div_map->{$bad_div_obj->{SEQID}}
				   ->{$old_divider}->{PAIR_ID},"]."));
			if($div_map->{$bad_div_obj->{SEQID}}->{$src_div}
			   ->{PAIR_ID} eq
			   $div_map->{$bad_div_obj->{SEQID}}->{$old_divider}
			   ->{PAIR_ID})
			  {debug($msg)}
			else
			  {warning($msg)}
			$div_map->{$bad_div_obj->{SEQID}}->{$src_div} =
			  {%{$div_map->{$bad_div_obj->{SEQID}}
			       ->{$old_divider}}};
			$div_map->{$bad_div_obj->{SEQID}}->{$src_div}
			  ->{HARD_START} = $src_div;
			$div_map->{$bad_div_obj->{SEQID}}->{$src_div}
			  ->{HARD_STOP} = $new_hard_stop;

			#If the iden start needs to be updated as well
			if($div_map->{$bad_div_obj->{SEQID}}->{$src_div}
			   ->{IDEN_START} > 0 &&
			   $div_map->{$bad_div_obj->{SEQID}}->{$src_div}
			   ->{IDEN_START} <
			   $div_map->{$bad_div_obj->{SEQID}}
			   ->{$src_div}->{HARD_START})
			  {$div_map->{$bad_div_obj->{SEQID}}->{$src_div}
			     ->{IDEN_START} =
			       $div_map->{$bad_div_obj->{SEQID}}
				 ->{$src_div}->{HARD_START}}
			#If the iden stop needs to be updated as well
			if($div_map->{$bad_div_obj->{SEQID}}->{$src_div}
			   ->{IDEN_STOP} > 0 &&
			   $div_map->{$bad_div_obj->{SEQID}}->{$src_div}
			   ->{IDEN_STOP} >
			   $div_map->{$bad_div_obj->{SEQID}}
			   ->{$src_div}->{HARD_STOP})
			  {$div_map->{$bad_div_obj->{SEQID}}->{$src_div}
			     ->{IDEN_STOP} =
			       $div_map->{$bad_div_obj->{SEQID}}
				 ->{$src_div}->{HARD_STOP}}
		      }
		  }
		else
		  {
		    $div_map->{$bad_div_obj->{SEQID}}->{$src_div} =
		      {%{$div_map->{$bad_div_obj->{SEQID}}->{$old_divider}}};
		    $div_map->{$bad_div_obj->{SEQID}}->{$src_div}
		      ->{HARD_START} = $src_div;
		    $div_map->{$bad_div_obj->{SEQID}}->{$src_div}
		      ->{HARD_STOP} = $new_hard_stop;

		    #If the iden start needs to be updated as well
		    if($div_map->{$bad_div_obj->{SEQID}}->{$src_div}
		       ->{IDEN_START} > 0 &&
		       $div_map->{$bad_div_obj->{SEQID}}->{$src_div}
		       ->{IDEN_START} <
		       $div_map->{$bad_div_obj->{SEQID}}
		       ->{$src_div}->{HARD_START})
		      {$div_map->{$bad_div_obj->{SEQID}}->{$src_div}
			 ->{IDEN_START} =
			   $div_map->{$bad_div_obj->{SEQID}}
			     ->{$src_div}->{HARD_START}}
		    #If the iden stop needs to be updated as well
		    if($div_map->{$bad_div_obj->{SEQID}}->{$src_div}
		       ->{IDEN_STOP} > 0 &&
		       $div_map->{$bad_div_obj->{SEQID}}->{$src_div}
		       ->{IDEN_STOP} >
		       $div_map->{$bad_div_obj->{SEQID}}
		       ->{$src_div}->{HARD_STOP})
		      {$div_map->{$bad_div_obj->{SEQID}}->{$src_div}
			 ->{IDEN_STOP} =
			   $div_map->{$bad_div_obj->{SEQID}}
			     ->{$src_div}->{HARD_STOP}}
		  }

		#Add/edit the conflicting source divider's partner divider
		updatePartnerSourceDivider($bad_div_obj->{SEQID},
					   $src_div,
					   $conflict_pairid,
					   $conf_part_seqid,
					   $data,
					   $div_map,
					   $div_map->{$bad_div_obj->{SEQID}}
					   ->{$src_div});
	      }
	  }
      }
    #2. Else (the HARD START of the conflicting divider occurs at or before the
    #   conflicting search range start)
    else
      {
	#If this is a conflict with self (in which case the lesser bound and
	#greater bound must each equal the hard start and the old divider must
	#be less than all of them), set the new divider to undef/recode
	if($hard_start == $lesser_bound && $hard_start == $greater_bound &&
	   $conflict_with_self && $old_divider < $hard_start)
	  {
	    $div_map->{$bad_div_obj->{SEQID}}->{$old_divider}->{TYPE} =
	      undef;
	    $div_map->{$bad_div_obj->{SEQID}}->{$old_divider}
	      ->{HARD_START} = undef;
	    $div_map->{$bad_div_obj->{SEQID}}->{$old_divider}
	      ->{HARD_STOP} = undef;
	    $div_map->{$bad_div_obj->{SEQID}}->{$old_divider}
	      ->{IDEN_START} = 0;
	    $div_map->{$bad_div_obj->{SEQID}}->{$old_divider}
	      ->{IDEN_STOP} = 0;

	    #Propagate the change back to the prior location in the partner
	    updatePartnerSourceDivider($bad_div_obj->{SEQID},
				       $old_divider,
				       $conflict_pairid,
				       $conf_part_seqid,
				       $data,
				       $div_map,
				       $div_map->{$bad_div_obj->{SEQID}}
				       ->{$old_divider});

	    $partner_updated = 1;

	    return();
	  }

	#2.1. Duplicate the source divider to every overlapping (hard start to
	#     hard stop) frame 1 position within the search range to protect
	#     those positions, (or we could somehow return the fact that the
	#     divider created should be a source divider matching the
	#     overlapping source divider).  Note, we are in this situation not
	#     because of indirect overlap of identity segments, but because a
	#     midpoint divider between segments was selected inside another
	#     segment.  Leave the conflicting divider's start unchanged.
	for(my $src_div = ($hard_start > $lesser_bound ?
			   $hard_start : $lesser_bound);
	    $src_div >= $lesser_bound && $src_div <= $greater_bound &&
	    $src_div < $hard_stop;
	    $src_div += 3)
	  {
	    if($src_div == $greater_bound && $hard_stop > ($src_div + 2))
	      {debug("Fixing the end of conflict ",
		     "[$bad_div_obj->{SEQID}:$src_div].")}

	    #If we're not splitting all the way to the hard stop, make sure the
	    #last coordinates used go to the original end of the segment
	    my $new_hard_stop = ($src_div == $greater_bound &&
				 $hard_stop > ($src_div + 2) ?
				 $hard_stop : $src_div + 2);

	    if(exists($div_map->{$bad_div_obj->{SEQID}}->{$src_div}))
	      {
		if($new_divider == $src_div)
		  {
		    if(defined($div_map->{$bad_div_obj->{SEQID}}
			       ->{$src_div}->{TYPE}) &&
		       $div_map->{$bad_div_obj->{SEQID}}->{$src_div}
		       ->{TYPE} eq 'SOURCE')
		      {
			$div_map->{$bad_div_obj->{SEQID}}->{$src_div}
			  ->{HARD_STOP} = $new_hard_stop;

			#If the iden stop needs to be updated as well
			if($div_map->{$bad_div_obj->{SEQID}}->{$src_div}
			   ->{IDEN_STOP} > 0 &&
			   $div_map->{$bad_div_obj->{SEQID}}->{$src_div}
			   ->{IDEN_STOP} >
			   $div_map->{$bad_div_obj->{SEQID}}
			   ->{$src_div}->{HARD_STOP})
			  {$div_map->{$bad_div_obj->{SEQID}}->{$src_div}
			     ->{IDEN_STOP} =
			       $div_map->{$bad_div_obj->{SEQID}}
				 ->{$src_div}->{HARD_STOP}}
		      }
		    else
		      {
			debug("Encountered a recode divider ",
			      "overlapping a source divider.");
			#Replace the recode divider with a source divider
			$div_map->{$bad_div_obj->{SEQID}}->{$src_div} =
			  {%{$div_map->{$bad_div_obj->{SEQID}}
			       ->{$old_divider}}};
			$div_map->{$bad_div_obj->{SEQID}}->{$src_div}
			  ->{HARD_STOP} = $new_hard_stop;

			#If the iden stop needs to be updated as well
			if($div_map->{$bad_div_obj->{SEQID}}->{$src_div}
			   ->{IDEN_STOP} > 0 &&
			   $div_map->{$bad_div_obj->{SEQID}}->{$src_div}
			   ->{IDEN_STOP} >
			   $div_map->{$bad_div_obj->{SEQID}}
			   ->{$src_div}->{HARD_STOP})
			  {$div_map->{$bad_div_obj->{SEQID}}->{$src_div}
			     ->{IDEN_STOP} =
			       $div_map->{$bad_div_obj->{SEQID}}
				 ->{$src_div}->{HARD_STOP}}
		      }
		  }
		else
		  {
		    my $msg =
		      join('',("Encountered a ",
			       (defined($div_map->{$bad_div_obj->{SEQID}}
					->{$src_div}->{TYPE}) &&
				$div_map->{$bad_div_obj->{SEQID}}->{$src_div}
				->{TYPE} eq 'SOURCE' ? 'source' : 'recode'),
			       " divider [$bad_div_obj->{SEQID}:$src_div",
			       (defined($div_map->{$bad_div_obj->{SEQID}}
					->{$src_div}->{TYPE}) &&
				$div_map->{$bad_div_obj->{SEQID}}->{$src_div}
				->{TYPE} eq 'SOURCE' ?
				"-" . $div_map->{$bad_div_obj->{SEQID}}
				->{$src_div}->{HARD_START} . "-" .
				$div_map->{$bad_div_obj->{SEQID}}->{$src_div}
				->{HARD_STOP} : ''),"] while splitting ",
			       "another source divider [$bad_div_obj->{SEQID}",
			       ":$old_divider-",
			       $div_map->{$bad_div_obj->{SEQID}}
			       ->{$old_divider}->{HARD_START},"-",
			       $div_map->{$bad_div_obj->{SEQID}}
			       ->{$old_divider}->{HARD_STOP},"] for overlap ",
			       "from ",
			       ($div_map->{$bad_div_obj->{SEQID}}->{$src_div}
				->{PAIR_ID} eq
				$div_map->{$bad_div_obj->{SEQID}}
				->{$old_divider}->{PAIR_ID} ?
				'the same' : 'a different')," alignment ",
			       "[encountered:",
			       $div_map->{$bad_div_obj->{SEQID}}->{$src_div}
			       ->{PAIR_ID},"] vs [splitting:",
			       $div_map->{$bad_div_obj->{SEQID}}
			       ->{$old_divider}->{PAIR_ID},"]."));
		    if(defined($div_map->{$bad_div_obj->{SEQID}}->{$src_div}
			       ->{TYPE}) &&
		       $div_map->{$bad_div_obj->{SEQID}}->{$src_div}->{TYPE} eq
		       'SOURCE' &&
		       $div_map->{$bad_div_obj->{SEQID}}->{$src_div}
		       ->{PAIR_ID} ne
		       $div_map->{$bad_div_obj->{SEQID}}->{$old_divider}
		       ->{PAIR_ID} &&

		       #This was added here when this warning was encountered
		       #in a dataset however, this likely needs to be added
		       #elsewhere as well
		       $src_div !=
		       $div_map->{$bad_div_obj->{SEQID}}->{$src_div}
		       ->{HARD_START})
		      {
			warning($msg);
			debug("The pre-existing alignment used to source ",
			      "[$bad_div_obj->{SEQID}:$src_div-",
			      $div_map->{$bad_div_obj->{SEQID}}->{$src_div}
			      ->{HARD_STOP},"] (",
			      $div_map->{$bad_div_obj->{SEQID}}->{$src_div}
			      ->{PAIR_ID},").  The identity coordinates [",
			      $div_map->{$bad_div_obj->{SEQID}}->{$src_div}
			      ->{IDEN_START},"-",
			      $div_map->{$bad_div_obj->{SEQID}}->{$src_div}
			      ->{IDEN_STOP},"] have been capitalized:\n",
			      getAlnSnippet($bad_div_obj->{SEQID},
					    $src_div,
					    $div_map->{$bad_div_obj->{SEQID}}
					    ->{$src_div}->{HARD_STOP},
					    $div_map->{$bad_div_obj->{SEQID}}
					    ->{$src_div}->{PAIR_ID},
					    getPartnerSeqID($bad_div_obj
							    ->{SEQID},
							    $div_map
							    ->{$bad_div_obj
							       ->{SEQID}}
							    ->{$src_div}
							    ->{PAIR_ID},
							    $data,
							    $div_map
							    ->{$bad_div_obj
							       ->{SEQID}}
							    ->{$src_div}
							    ->{IDEN_START},
							    $div_map
							    ->{$bad_div_obj
							       ->{SEQID}}
							    ->{$src_div}
							    ->{IDEN_STOP}),
					    $data),"\n\n",
			      "The alignment we want to use to source ",
			      "[$bad_div_obj->{SEQID}:$src_div-",
			      $div_map->{$bad_div_obj->{SEQID}}->{$src_div}
			      ->{HARD_STOP},"] (",
			      $div_map->{$bad_div_obj->{SEQID}}->{$old_divider}
			      ->{PAIR_ID},") - The portion of the pre-",
			      "existing identity coords [",
			      $div_map->{$bad_div_obj->{SEQID}}->{$old_divider}
			      ->{IDEN_START},"-",
			      ($new_hard_stop <
			       $div_map->{$bad_div_obj->{SEQID}}
			       ->{$old_divider}->{IDEN_STOP} ? $new_hard_stop :
			       $div_map->{$bad_div_obj->{SEQID}}
			       ->{$old_divider}->{IDEN_STOP}),"] we intend ",
			      "to replace have been capitalized:\n",
			      getAlnSnippet($bad_div_obj->{SEQID},
					    $src_div,
					    $new_hard_stop,
					    $div_map->{$bad_div_obj->{SEQID}}
					    ->{$old_divider}->{PAIR_ID},
					    getPartnerSeqID($bad_div_obj
							    ->{SEQID},
							    $div_map
							    ->{$bad_div_obj
							       ->{SEQID}}
							    ->{$old_divider}
							    ->{PAIR_ID},
							    $data),
					    $data,
					    $div_map->{$bad_div_obj->{SEQID}}
					    ->{$old_divider}->{IDEN_START},
					    ($new_hard_stop <
					     $div_map->{$bad_div_obj->{SEQID}}
					     ->{$old_divider}->{IDEN_STOP} ?
					     $new_hard_stop :
					     $div_map->{$bad_div_obj->{SEQID}}
					     ->{$old_divider}->{IDEN_STOP})));
		      }
		    else
		      {debug($msg)}

		    $div_map->{$bad_div_obj->{SEQID}}->{$src_div} =
		      {%{$div_map->{$bad_div_obj->{SEQID}}
			   ->{$old_divider}}};
		    $div_map->{$bad_div_obj->{SEQID}}->{$src_div}
		      ->{HARD_START} = $src_div;
		    $div_map->{$bad_div_obj->{SEQID}}->{$src_div}
		      ->{HARD_STOP} = $new_hard_stop;

		    #If the iden start needs to be updated as well
		    if($div_map->{$bad_div_obj->{SEQID}}->{$src_div}
		       ->{IDEN_START} > 0 &&
		       $div_map->{$bad_div_obj->{SEQID}}->{$src_div}
		       ->{IDEN_START} <
		       $div_map->{$bad_div_obj->{SEQID}}
		       ->{$src_div}->{HARD_START})
		      {$div_map->{$bad_div_obj->{SEQID}}->{$src_div}
			 ->{IDEN_START} =
			   $div_map->{$bad_div_obj->{SEQID}}
			     ->{$src_div}->{HARD_START}}
		    #If the iden stop needs to be updated as well
		    if($div_map->{$bad_div_obj->{SEQID}}->{$src_div}
		       ->{IDEN_STOP} > 0 &&
		       $div_map->{$bad_div_obj->{SEQID}}->{$src_div}
		       ->{IDEN_STOP} >
		       $div_map->{$bad_div_obj->{SEQID}}
		       ->{$src_div}->{HARD_STOP})
		      {$div_map->{$bad_div_obj->{SEQID}}->{$src_div}
			 ->{IDEN_STOP} =
			   $div_map->{$bad_div_obj->{SEQID}}
			     ->{$src_div}->{HARD_STOP}}
		  }
	      }
	    else
	      {
		$div_map->{$bad_div_obj->{SEQID}}->{$src_div} =
		  {%{$div_map->{$bad_div_obj->{SEQID}}->{$old_divider}}};
		$div_map->{$bad_div_obj->{SEQID}}->{$src_div}
		  ->{HARD_START} = $src_div;
		$div_map->{$bad_div_obj->{SEQID}}->{$src_div}
		  ->{HARD_STOP} = $new_hard_stop;

		#If the iden start needs to be updated as well
		if($div_map->{$bad_div_obj->{SEQID}}->{$src_div}
		   ->{IDEN_START} > 0 &&
		   $div_map->{$bad_div_obj->{SEQID}}->{$src_div}
		   ->{IDEN_START} <
		   $div_map->{$bad_div_obj->{SEQID}}
		   ->{$src_div}->{HARD_START})
		  {$div_map->{$bad_div_obj->{SEQID}}->{$src_div}
		     ->{IDEN_START} =
		       $div_map->{$bad_div_obj->{SEQID}}
			 ->{$src_div}->{HARD_START}}
		#If the iden stop needs to be updated as well
		if($div_map->{$bad_div_obj->{SEQID}}->{$src_div}
		   ->{IDEN_STOP} > 0 &&
		   $div_map->{$bad_div_obj->{SEQID}}->{$src_div}
		   ->{IDEN_STOP} >
		   $div_map->{$bad_div_obj->{SEQID}}
		   ->{$src_div}->{HARD_STOP})
		  {$div_map->{$bad_div_obj->{SEQID}}->{$src_div}
		     ->{IDEN_STOP} =
		       $div_map->{$bad_div_obj->{SEQID}}
			 ->{$src_div}->{HARD_STOP}}
	      }

	    #Add/edit the conflicting source divider's partner divider
	    updatePartnerSourceDivider($bad_div_obj->{SEQID},
				       $src_div,
				       $conflict_pairid,
				       $conf_part_seqid,
				       $data,
				       $div_map,
				       $div_map->{$bad_div_obj->{SEQID}}
				       ->{$src_div});
	  }

	#We started at the hard start or the lesser_bound, whichever is further
	#and that's where we started editing, so we still need to edit the
	#beginning of this conflicting divider
	#If beginning of the conflicting divider was left of the search range
	if($old_divider < $lesser_bound)
	  {
	    #If the hard start is right at the beginning of the search range
	    if($hard_start == $lesser_bound)
	      {
		#We can eliminate the hard and iden coords from the divider
		$div_map->{$bad_div_obj->{SEQID}}->{$old_divider}->{HARD_START}
		  = undef;
		$div_map->{$bad_div_obj->{SEQID}}->{$old_divider}->{HARD_STOP}
		  = undef;
		$div_map->{$bad_div_obj->{SEQID}}->{$old_divider}->{IDEN_START}
		  = 0;
		$div_map->{$bad_div_obj->{SEQID}}->{$old_divider}->{IDEN_STOP}
		  = 0;
	      }
	    else
	      {
		$div_map->{$bad_div_obj->{SEQID}}->{$old_divider}->{HARD_STOP}
		  = $lesser_bound - 1;
		$div_map->{$bad_div_obj->{SEQID}}->{$old_divider}->{IDEN_STOP}
		  = $lesser_bound - 1;
	      }
	  }
      }

    #If the old divider changed, copy the changes to the source partner
    if(!$partner_updated && $old_divider != $new_divider)
      {
	#Add/edit the conflicting source divider's partner divider
	updatePartnerSourceDivider($bad_div_obj->{SEQID},
				   $old_divider,
				   $conflict_pairid,
				   $conf_part_seqid,
				   $data,
				   $div_map,
				   $div_map->{$bad_div_obj->{SEQID}}
				   ->{$old_divider});
      }
  }

#This sub creates or overwrites source dividers in identity segment partners.
#It's intended to only be used to change an old divider that is not located at
#its hard start to an undef/recode divider and to add dividers at frame 1
#intervals inside an identity segment.  The changes are the result of resolving
#divider conflicts where a divider was moved or source dividers are created to
#"protect" a conflicting source divider from getting its identity segment
#interrupted by a recode divider.
sub updatePartnerSourceDivider
  {
    my $seqid      = $_[0];
    my $divider    = $_[1];
    my $pair_id    = $_[2];
    my $part_seqid = $_[3];
    my $data       = $_[4];
    my $div_map    = $_[5];
    my $div_hash   = $_[6];

    if(defined($div_hash->{STOP}))
      {
	error("Cannot call updatePartnerSourceDivider after updating the ",
	      "divider map stops.");
	return();
      }

    my $part_divider =
      convertDivider($seqid,
		     $divider,
		     $pair_id,
		     $part_seqid,
		     $data);

    #If the divider didn't convert properly due to a gap or weird overlap
    #issue, try to find it
    if(!exists($div_map->{$part_seqid}) ||
       !exists($div_map->{$part_seqid}->{$part_divider}))
      {
	my $tmp_part_divider =
	  convertAlnSeqCoord($seqid,
			     $divider,
			     $pair_id,
			     $part_seqid,
			     $data);
	if(calcFrame($tmp_part_divider) == 3)
	  {$tmp_part_divider++}
	if(calcFrame($tmp_part_divider) == 1)
	  {
	    #Now, if we can convert back and get the same coord, it's confirmed
	    my $check_div = convertDivider($part_seqid,
					   $tmp_part_divider,
					   $pair_id,
					   $seqid,
					   $data);
	    if($check_div != $divider)
	      {$check_div = convertAlnSeqCoord($part_seqid,
					       $tmp_part_divider,
					       $pair_id,
					       $seqid,
					       $data)}
	    if($check_div == $divider)
	      {$part_divider = $tmp_part_divider}
	  }
      }

    if(exists($div_map->{$part_seqid}) &&
       exists($div_map->{$part_seqid}->{$part_divider}) &&
       !defined($div_map->{$part_seqid}->{$part_divider}->{PAIR_ID}))
      {warning("Source partner [$part_seqid:$part_divider] to [$seqid:",
	       "$divider] doesn't have a defined pair ID in the map.  ",
	       "Overwriting source divider change created for ",
	       "[$div_hash->{PAIR_ID}] and converted using [$pair_id].")}
    elsif(exists($div_map->{$part_seqid}) &&
	  exists($div_map->{$part_seqid}->{$part_divider}) &&
	  $div_hash->{PAIR_ID} ne
	  $div_map->{$part_seqid}->{$part_divider}->{PAIR_ID} &&
	  defined($div_map->{$part_seqid}->{$part_divider}->{TYPE}) &&
	  $div_map->{$part_seqid}->{$part_divider}->{TYPE} eq 'SOURCE')
      {
	#This can legitimately happen when 1 sequence has overlapping identity
	#segments with 2 other sequences that do not overlap one another.  The
	#common sequence's segments get adjusted to cover all identity, but the
	#partner sequence's do not.  If they run into a divider conflict, they
	#will try and move their divider, but since the segment's start or stop
	#may have been adjusted to accommodate the overlap, the divider
	#encountered could be from a different pair.

	#If the divider we've adjusted (previously) extends over a larger
	#region, we'll allow the overwrite with a warning
	if($div_hash->{IDEN_STOP} >
	   $div_map->{$part_seqid}->{$part_divider}->{IDEN_STOP})
	  {
	    #We only need to warn about the adjustment being done if either
	    #the adjusting divider or the unexpected divider are not equal to
	    #their hard start
	    if($divider != $div_hash->{HARD_START} ||
	       $part_divider !=
	       $div_map->{$part_seqid}->{$part_divider}->{HARD_START})
	      {warning("The source partner has a different pair ID [",
		       $div_map->{$part_seqid}->{$part_divider}->{PAIR_ID},
		       "] than the one that was moved in the partner: ",
		       "[$div_hash->{PAIR_ID}].",
		       {DETAIL =>
			join('',("Divider found: [$part_seqid:$part_divider].",
				 "  Conflicting divider: [$seqid:$divider].  ",
				 "This can be caused by a pair of ",
				 "overlapping identity segments one sequence ",
				 "has with 2 other sequences, which causes ",
				 "potential collisions such as this one.  ",
				 "Keeping the adjusted larger/longer ",
				 "segment."))})}
	  }
	#If the divider we've adjusted (previously) extends over a smaller
	#region, we'll skip the overwrite with a warning
	elsif($div_hash->{IDEN_STOP} <
	      $div_map->{$part_seqid}->{$part_divider}->{IDEN_STOP})
	  {
	    #We only need to warn about the adjustment not being done if either
	    #the adjusting divider or the unexpected divider are not equal to
	    #their hard start
	    if($divider != $div_hash->{HARD_START} ||
	       $part_divider !=
	       $div_map->{$part_seqid}->{$part_divider}->{HARD_START})
	      {warning("The source partner has a different pair ID [",
		       $div_map->{$part_seqid}->{$part_divider}->{PAIR_ID},
		       "] than the one that was moved in the partner: ",
		       "[$div_hash->{PAIR_ID}].",
		       {DETAIL =>
			join('',("Divider found: [$part_seqid:$part_divider].",
				 "  Conflicting divider: [$seqid:$divider].  ",
				 "This can be caused by a pair of ",
				 "overlapping identity segments one sequence ",
				 "has with 2 other sequences, which causes ",
				 "potential collisions such as this one.  ",
				 "Keeping the original larger/longer ",
				 "segment."))})}
	    return();
	  }
	#Else the region is the same length, so we'll quietly skip it
      }
    #If we have not gotten to this sequence yet, let it quietly find the
    #adjusted partner divider
    elsif(!exists($div_map->{$part_seqid}))
      {return()}
    #If the expected partner source divider doesn't exist in the expected
    #place, make sure adjusting it won't cause a problem
    elsif(!exists($div_map->{$part_seqid}->{$part_divider}))
      {
	if((!defined($div_hash->{HARD_START}) ||
	    $divider != $div_hash->{HARD_START}) && isDebug())
	  {warning("Source divider [$seqid:$divider] doesn't exist in the ",
		   "partner's map in the expected location: [$part_seqid:",
		   "$part_divider].  Unable to adjust the divider position.  ",
		   "Associated region: [$seqid:$divider",
		   (defined($div_hash->{HARD_STOP}) ?
		    "-$div_hash->{HARD_STOP}" : ''),"] may not be encoded ",
		   "for optimal homology with [$part_seqid].",
		   {DETAIL =>
		    join('',("Divider [$seqid:$divider] for identity located ",
			     "at [",
			     (defined($div_hash->{HARD_STOP}) ?
			      "$div_hash->{HARD_START}-$div_hash->{HARD_STOP}"
			      : 'unknown'),"] based on alignment: ",
			     "[$div_hash->{PAIR_ID}] was moved due to a ",
			     "conflict with another divider (which happens ",
			     "when a sequence's partners have identity ",
			     "segments with sequences the first sequence has ",
			     "no observed identity with).  Each divider for ",
			     "an identity segment has a companion in its ",
			     "partner.  The expected divider in the partner ",
			     "was not found.  This can be caused by a pair ",
			     "of overlapping identity segments that one ",
			     "sequence has with 2 other sequences, when the ",
			     "partners have no corresponding segment overlap ",
			     "in their own map, which causes, as in this ",
			     "case, the dividers to become out of synch.  ",
			     "It's OK.  The sequence encoding won't be ",
			     "harmed, but it may just not have the highest ",
			     "homology in the region to the right of the ",
			     "divider to the next segment."))})}
	return();
      }

    #Copy the divider hash from the original source partner that changed
    $div_map->{$part_seqid}->{$part_divider} = {%$div_hash};

    #Error-check that the PAIR_ID is valid for this sequence
    getPartnerSeqID($part_seqid,
		    $div_map->{$part_seqid}->{$part_divider}->{PAIR_ID},
		    $data);

    my $new_type = $div_hash->{TYPE};

    #Update the coordinates to this partner if this will be a source divider
    if(defined($new_type) && $new_type eq 'SOURCE')
      {
	$div_map->{$part_seqid}->{$part_divider}->{HARD_START} =
	  convertAlnSeqCoord($seqid,
			     $div_hash->{HARD_START},
			     $pair_id,
			     $part_seqid,
			     $data);
	$div_map->{$part_seqid}->{$part_divider}->{HARD_STOP}  =
	  convertAlnSeqCoord($seqid,
			     $div_hash->{HARD_STOP},
			     $pair_id,
			     $part_seqid,
			     $data);

	#If the identity segment is stored in this divider
	if($div_map->{$part_seqid}->{$part_divider}->{IDEN_START} > 0)
	  {
	    $div_map->{$part_seqid}->{$part_divider}->{IDEN_START} =
	      convertAlnSeqCoord($seqid,
				 $div_hash->{IDEN_START},
				 $pair_id,
				 $part_seqid,
				 $data);
	    $div_map->{$part_seqid}->{$part_divider}->{IDEN_STOP}  =
	      convertAlnSeqCoord($seqid,
				 $div_hash->{IDEN_STOP},
				 $pair_id,
				 $part_seqid,
				 $data);
	  }
      }
    else
      {
	#Note: PAIR_ID will have already been changed in the passed div_hash
	#In fact, these are probably already changed as well, but just to be on
	#the safe side...
	$div_map->{$part_seqid}->{$part_divider}->{HARD_START} = undef;
	$div_map->{$part_seqid}->{$part_divider}->{HARD_STOP}  = undef;
	$div_map->{$part_seqid}->{$part_divider}->{IDEN_START} = 0;
	$div_map->{$part_seqid}->{$part_divider}->{IDEN_STOP}  = 0;
      }
  }

#Searches the segments of a solution for the closest right-side coordinate that
#is to the left of a given coordinate and returns the coordinate after it as a
#valid divider (assuming all end coordinates are in frame 3).  A linear search
#is performed.  Could improve this to be a binary search.  Note: Returns 1 if
#there are no segments to the left or if the coordinate passed in was 1.
sub getClosestLeftSegCoord
  {
    my $seqid = $_[0];
    my $coord = $_[1];
    my $soln  = $_[2];

    my $closest = 1;

    foreach my $rec (sort {$a->{FINAL_START} <=> $b->{FINAL_START}}
		     grep {$_->{TYPE} eq 'seg'}
		     @{$soln->{$seqid}})
      {
	#If we're back to the query coordinate, return the previous one
	if($rec->{FINAL_START} >= $coord)
	  {last}
	if($rec->{FINAL_STOP} < $coord)
	  {$closest = $rec->{FINAL_STOP} + 1}
	#If this is a previous segment (i.e. it starts before the query coord)
	#but it overlaps the query segment, set the divider at the query
	#coordinate
	elsif($rec->{FINAL_START} < $coord && $rec->{FINAL_STOP} >= $coord)
	  {
	    $closest = $coord;
	    last;
	  }
      }

    if((($closest - 1) % 3 + 1) != 1) #If the coord returned is not in frame 1
      {error("Records with stops not in frame 3 exist.")}

    return($closest);
  }

#Searches sequences that have been paired with a particular sequence to find
#the closest segment to its left and returns a midpoint on a codon boundary
#(frame 1)
sub getClosestRightDivider
  {
    my $seqid    = $_[0];
    my $coord    = $_[1];
    my $soln     = $_[2];
    my $data     = $_[3];
    my $partners = $_[4];
    my $divs     = $_[5];

    if($coord % 3 != 0) #If the coordinate is not in frame 3
      {error("Coordinate: [$coord] must be in frame 3.  Frame [",
	     ((($coord - 1) % 3) + 1),"] was sent in.")}

    #For each alignment $seqid has been aligned with
    my($closest,$seqlen);
    my $closest_partner = $partners->[0];
    foreach my $partner (@$partners)
      {
	if(!defined($closest))
	  {
	    $seqlen  = length($data->{$partner->[0]}->{SEQS}->{$seqid});
	    $closest = $seqlen + 1;
	  }

	#Convert the coordinate to the partner sequence
	#Note, we're adding 1 to the coordinate becasue it is assumed that
	#$coord is in frame 3 and the divider can only be in frame 1.
	#convertDivider therefor can only take a frame 1 coordinate.
	my $partner_coord = convertDivider($seqid,
					   $coord + 1,
					   $partner->[0],  #pair ID
					   $partner->[1],  #partner seq ID
					   $data);

	#Find the closest coord to the right
	my($part_rght_seg_coord,$part_right_seg_pair) =
	  getClosestRightSegCoord($partner->[1], #seqid
				  $partner_coord,
				  $soln,
				  length($data->{$partner->[0]}->{SEQS}
					 ->{$partner->[1]}));

	#Convert the coordinate back to this sequence
	#Note, I had originally used convertDivider here, but that causes the
	#following problem: If the query sequence has a segment from 1-6 and a
	#partner aligns with an end-gap at the beinning (aligning with the
	#query's coord 13) and the partner has a segment from 1-whatever, that
	#1 gets converted back to 1 here in the query sequence instead of 13.
	my $orig_rght_seg_coord =
	  convertAlnSeqCoord($partner->[1],  #partner seqID
			     $part_rght_seg_coord,
			     $partner->[0],  #pair ID
			     $seqid,
			     $data);
	#If $part_rght_seg_coord aligns inside a large gap, the returned value
	#will be in frame 3 because convertAlnSeqCoord returns the last
	#sequence coordinate counted.  We add 1 to assume that that is the
	#closest corresponding segment.
	if(calcFrame($orig_rght_seg_coord) == 3)
	  {$orig_rght_seg_coord++}

	#If this coordinate is closer, save it.
	if($orig_rght_seg_coord < $closest)
	  {
	    debug({LEVEL => 4},"Closest right segment to partner [@$partner ",
		  "$partner_coord] is [$part_right_seg_pair $partner->[1] ",
		  "$part_rght_seg_coord] (i.e. to self [$seqid:$coord+1]) is ",
		  "[$partner->[0] $seqid $orig_rght_seg_coord]");

	    #If the divider is inside the query coordinate, something went
	    #wrong.  If it's to the immediate right, it's OK to use as a
	    #divider - though there may be an issue if this turns out to not be
	    #a valid codon boundary (i.e. frame 1).)
	    if($orig_rght_seg_coord <= $coord)
	      {
		error("Returned coord [$orig_rght_seg_coord] overlaps query ",
		      "coord [$coord] in seq [$seqid] on its right side.  ",
		      "Obtained from a segment in [$partner->[1]] and ",
		      "alignment [$partner->[0]].  The query coord was ",
		      "converted from [$seqid:$coord+1] to [$partner->[1]:",
		      "$partner_coord] and the closest coord in that ",
		      "sequence was: [$partner->[1]:$part_rght_seg_coord] - ",
		      "which should not overlap, but did overlap when ",
		      "converted back to [$seqid:$orig_rght_seg_coord] ",
		      "coord.");
	      }
	    $closest = $orig_rght_seg_coord;
	    $closest_partner = [@$partner,$part_rght_seg_coord];
	    if(($orig_rght_seg_coord - 1) == $coord)
	      {return($orig_rght_seg_coord,$orig_rght_seg_coord)}
	  }
	else
	  {
	    debug({LEVEL => 4},"Further right segment to [@$partner ",
		  "$partner_coord] is [@$partner $part_rght_seg_coord] and ",
		  "converted back to [$partner->[0] $seqid ",
		  "$orig_rght_seg_coord]");
	  }
      }

    #This means that a segment abutts this one either directly or indirectly,
    #but in any case, the divider can't be anywhere else
    if(($closest - 1) == $coord)
      {return($closest,$closest)}

    #See if there already exists a divider (or dividers) in this region in any
    #of the partners (or self).
    #Dividers with slightly different coordinates are common because their
    #position is calculated using different alignments and the alignments can
    #be very gappy.

    getFixDividerConflicts($seqid,
			   $coord + 1,
			   $closest,
			   $closest_partner->[1],
			   $closest_partner->[0],
			   $divs,
			   $data,
			   0,
			   1);

    #Now check my partners
    foreach my $partner (@$partners)
      {
	getFixDividerConflicts($seqid,
			       $coord + 1,
			       $closest,
			       $partner->[1],
			       $partner->[0],
			       $divs,
			       $data,
			       0,
			       0,
			       $closest_partner->[1],
			       $closest_partner->[0]);
      }

    #First check myself...
    my $self_divider =
      dividerExists($seqid,
		    $coord,
		    $closest,
		    $closest_partner->[1],
		    $closest_partner->[0],
		    $divs,
		    $data,
		    0,
		    1);

    if($self_divider)
      {return($self_divider,$closest)}

    #Now check my partners
    my($left_most_divider);
    foreach my $partner (@$partners)
      {
	#This should be converted back
	my $converted_self_divider = dividerExists($seqid,
						   $coord,
						   $closest,
						   $partner->[1],
						   $partner->[0],
						   $divs,
						   $data,
						   0,
						   0);

	if(defined($converted_self_divider) && $converted_self_divider &&
	   (!defined($left_most_divider) ||
	    $converted_self_divider < $left_most_divider))
	  {$left_most_divider = $converted_self_divider}
      }

    if(defined($left_most_divider) && $left_most_divider)
      {return($left_most_divider,$closest)}

    #If the closest "segment" is the end of the sequence (plus 1), no need to
    #find the midpoint, because it's not a real segment - rather it's just the
    #end boundary
    if($closest == ($seqlen + 1))
      {return($closest,$closest)}

    #We now have the closest occupied coordinate to the right of $coord.  We
    #need to split the difference and then return the closest frame 1 codon
    #position.  We add 1 to $coord because we want to choose a midpoint among
    #available positions to put a divider.  That includes the $closest
    #position, but not the $coord position (the last position of the query
    #segment)
    debug({LEVEL => 4},"Calculating midpoint using alignment ",
	  "[$closest_partner->[0]]");
    my $midpoint = getAlnMidpoint($closest_partner->[0],  #pairID
				  $seqid,
				  $coord + 1,
				  $closest,
				  $closest_partner->[1],
				  $data);

    #Determine the frame position (1, 2, or 3) of the coordinate
    my $frame_pos = ($midpoint - 1) % 3 + 1;

    if($frame_pos == 1)
      {return($midpoint,$closest)}
    elsif($frame_pos == 2)
      {
	if(($midpoint + 2) > $closest)
	  {
	    #This should not happen, but just in case...
	    if(($midpoint - 1) <= $coord)
	      {
		error("Bad frame boundary.  Frame 2.  Midpoint: ",
		      "[$midpoint].  Closest: [$closest].  Coord(/final ",
		      "start): [$coord].");
		return($coord + 1,$closest);
	      }

	    return($midpoint - 1,$closest);
	  }

	return($midpoint + 2,$closest);
      }
    elsif($frame_pos == 3)
      {
	if(($midpoint + 1) > $closest)
	  {
	    error("Bad frame boundary for right divider.  Frame 3.  ",
		  "Midpoint: [$midpoint].  Closest: [$closest].  Coord(",
		  "/final start): [$coord].");
	    return($coord + 1,$closest);
	  }

	return($midpoint + 1,$closest);
      }
    else
      {error("Unexpected frame result.")}

    return($coord + 1,$closest);
  }

#Searches the segments of a solution for the closest left-most coordinate that
#is to the right of (or equal to) a given coordinate.  A linear search is
#performed.  Could improve this to be a binary search.  Note: Returns sequence
#length +1 if there are no segments to the right or of (or equal to) the
#coordinate passed in, or if it was the sequence length.
sub getClosestRightSegCoord
  {
    my $seqid  = $_[0];
    my $coord  = $_[1];
    my $soln   = $_[2];
    my $seqlen = $_[3];

    if($seqlen < ($coord - 1))
      {error("Coordinate: [$coord] submitted is larger than the sequence ",
	     "length submitted for [$seqid]: [$seqlen] (plus 1).")}

    my $closest = $seqlen + 1;
    my $from_pair = '';

    foreach my $rec (sort {$b->{FINAL_START} <=> $a->{FINAL_START}}
		     grep {$_->{TYPE} eq 'seg'}
		     @{$soln->{$seqid}})
      {
	if($rec->{FINAL_STOP} < $coord)
	  {last}
	if($rec->{FINAL_START} > $coord)
	  {
	    $closest   = $rec->{FINAL_START};
	    $from_pair = $rec->{PAIR_ID};
	  }
	elsif($rec->{FINAL_START} <= $coord && $rec->{FINAL_STOP} >= $coord)
	  {
	    $closest = $coord;
	    $from_pair = $rec->{PAIR_ID};
	    last;
	  }
      }

    return($closest,$from_pair);
  }

#This converts and copies the divider found when looking at one sequence's
#segments to the other sequences, though it only defines the alignment/pair the
#sequence should come from if it is a copy to the direct partner.  Indirect
#partners only create an undefined entry at the converted location.  An
#indirect partner is an alignment between the originating sequence and another
#sequence other than the one involved in the originating alignment.  What I
#mean by "originating", is where the segment came from the prompted the
#creation of the divider.  A divider is a starting coordinate and pair ID for
#a particular sequence indicating where sequence for that position forward
#should be obtained from.  Example: Sequence 1 is paired with sequence 2 in
#alignment A.  A segment in that alignment prompts the creation of a divider to
#the left and right of that segment.  Sequence 1 is also in alignment B with
#sequence 3.  The divider in sequence 1 has a corresponding location in
#sequence 3, so it is copied there, but is undefined because that segment of
#sequence needs to be recoded to best match sequence 1 as it exists in
#alignment A and not as it exists in alignment B anymore.
sub copyDivider
  {
    my $divider  = $_[0];
    my $seqid    = $_[1];
    my $div_map  = $_[2];
    my $partners = $_[3];
    my $data     = $_[4];
    my $pair_id  = $_[5];
    my $side     = $_[6]; #left side = we know the source or whether to recode
    my $soln     = $_[7]; #So see if a divider is copying into an identity
                          #segment, so the pairID & source type can be updated
    my $lesser_bound  = $_[8]; #Don't copy if a divider in this range exists
    my $greater_bound = $_[9]; #If $hard_stop defined, this will be hard_start
    my $hard_stop     = $_[10];

    my $uniq_seq_partners = {};
    foreach my $partner_seqid (grep {$_ ne $seqid} map {$_->[1]} @$partners)
      {$uniq_seq_partners->{$partner_seqid} = 0}

    foreach my $partner (@$partners)
      {
	my $partner_pair_id = $partner->[0];
	my $partner_seqid   = $partner->[1];
	my($partner_divider);

	#A divider can be in the middle of a gap and it doesn't convert
	#correctly, thus we need to check between the bounds to see if a
	#divider already exists here.  The coordinates of the bounds are
	#guaranteed/assumed to be the closest segments possible.  If a divider
	#exists in the bounds, it is assumed that when it's converted to the
	#query sequence, it's the same coordinate.  It's the copying here that
	#is assumed to be inaccurate.
	my $existing_divider = dividerExists($seqid,
					     $lesser_bound,
					     $greater_bound,
					     $partner_seqid,
					     $partner_pair_id,
					     $div_map,
					     $data,
					     $side,
					     0);
	if(defined($existing_divider) && $existing_divider &&
	   ($side ne 'left' || $pair_id ne $partner_pair_id ||
	    !exists($div_map->{$partner_seqid}) ||
	    !exists($div_map->{$partner_seqid}->{$existing_divider})))
	  {
	    debug("Skipping due to existing divider: ",
		  "[$seqid:$existing_divider].");
	    next;
	  }
	elsif(defined($existing_divider) && $existing_divider)
	  {$partner_divider = convertDivider($seqid,
					     $existing_divider,
					     $partner_pair_id,
					     $partner_seqid,
					     $data)}
	else
	  {$partner_divider = convertDivider($seqid,
					     $divider,
					     $partner_pair_id,
					     $partner_seqid,
					     $data)}

	#If we can't fix any divider conflicts that would result from this
	#copy, skip it.
	unless(getFixDividerConflicts($seqid,
				      $lesser_bound,
				      $greater_bound,
				      $partner_seqid,
				      $partner_pair_id,
				      $div_map,
				      $data,
				      0,
				      0))
	  {
	    warning("Conflicts found in divider copy: [$seqid:$divider to ",
		    "$partner_seqid:$partner_divider].");
	  }

	my $partner_bound = ($side ne 'left' ? undef :
			     convertDivider($seqid,
					    $greater_bound,
					    $partner_pair_id,
					    $partner_seqid,
					    $data) - 1);

	#Assuming that the hard-stop is always directly translatable when it's
	#a left/source divider.  It should be, because it's inside an identity
	#segment.
	my $partner_hard_stop = ($side ne 'left' || !defined($hard_stop) ?
				 undef : convertAlnSeqCoord($seqid,
							    $hard_stop,
							    $partner_pair_id,
							    $partner_seqid,
							    $data));

	#This divider might have been copied from another sequence, so
	#check it for a conflict before recording it

	#If this divider doesn't already exist in this sequence
	if(!exists($div_map->{$partner_seqid}) ||
	   !exists($div_map->{$partner_seqid}->{$partner_divider}))
	  {
	    #If this sequence is a part of the same pair as the originating
	    #sequence, set the pair ID in the div map
	    if($pair_id eq $partner_pair_id)
	      {
		debug({LEVEL => 4},"Copying divider $seqid:$divider in ",
		      "$pair_id to $partner_seqid:$partner_divider as ",
		      "$pair_id");

		#Error-check that the pair ID is correct:
		getPartnerSeqID($partner_seqid,$pair_id,$data);

		$div_map->{$partner_seqid}->{$partner_divider} =
		  {PAIR_ID    => $pair_id,
		   TYPE       => ($side eq 'left' ? 'SOURCE' : undef),
		   HARD_START => (defined($partner_bound) &&
				  defined($partner_hard_stop) ?
				  $partner_bound + 1 : undef),
		   HARD_STOP  => $partner_hard_stop,
		   STOP       => undef,
		   IDEN_START => 0, #This will hopefully get replaced when the
		   IDEN_STOP  => 0, #source is handled directly
		   SEQ        => ''};
	      }
	    else
	      {
		#We need to make sure that the divider will end up within the
		#expected bounds when copied back (like if it's found as an
		#existing divider, will it be in the expected region).  A
		#divider can end up being completely out of bounds if the
		#region we're looking in aligns with a lot of gaps with another
		#sequence.
		my $check_divider = convertDivider($partner_seqid,
						   $partner_divider,
						   $partner_pair_id,
						   $seqid,
						   $data);
		if(coordsOverlap($check_divider,$check_divider,
				 $lesser_bound,$greater_bound))
		  {
		    debug({LEVEL => 4},"Copying divider $seqid:$divider (in ",
			  "between $seqid:$lesser_bound and $seqid:",
			  "$greater_bound) in $pair_id to $partner_seqid:",
			  "$partner_divider as $partner_pair_id which is ",
			  "position: [$seqid:",
			  convertDivider($partner_seqid,
					 $partner_divider,
					 $partner_pair_id,
					 $seqid,
					 $data),"] when converted back.");

		    #Error-check that the pair ID is correct:
		    getPartnerSeqID($partner_seqid,$partner_pair_id,$data);

		    $div_map->{$partner_seqid}->{$partner_divider} =
		      {PAIR_ID    => $partner_pair_id,
		       TYPE       => undef,
		       STOP       => undef,
		       HARD_START => undef,
		       HARD_STOP  => undef,
		       IDEN_START => 0,
		       IDEN_STOP  => 0,
		       SEQ        => ''};
		  }
		else
		  {debug({LEVEL => 4},"Skipping divider copy of $seqid:",
			 "$divider (in between $seqid:$lesser_bound and ",
			 "$seqid:$greater_bound) in $pair_id to ",
			 "$partner_seqid:$partner_divider as ",
			 "$partner_pair_id because position: [$seqid:",
			  convertDivider($partner_seqid,
					 $partner_divider,
					 $partner_pair_id,
					 $seqid,
					 $data),"] is out of bounds when",
			 " converted back.")}
	      }
	  }
	elsif($side eq 'left' && $pair_id eq $partner_pair_id &&
	      (!defined($div_map->{$partner_seqid}->{$partner_divider}
			->{TYPE}) ||
	       $div_map->{$partner_seqid}->{$partner_divider}->{TYPE} eq
	       'RECODE'))
	  {
	    #This is called to simply check that $partner_pair_id is valid
	    getPartnerSeqID($partner_seqid,$partner_pair_id,$data);

	    debug("Over-writing recode divider with a source divider.");
	    $div_map->{$partner_seqid}->{$partner_divider}->{PAIR_ID} =
	      $partner_pair_id;
	    $div_map->{$partner_seqid}->{$partner_divider}->{TYPE} = 'SOURCE';
	    $div_map->{$partner_seqid}->{$partner_divider}->{HARD_START} =
	      (defined($partner_bound) && defined($partner_hard_stop) ?
	       $partner_bound + 1 : undef);
	    $div_map->{$partner_seqid}->{$partner_divider}->{HARD_STOP} =
	      $partner_hard_stop;
	  }
	else
	  {
	    debug("Source divider already exists in this location: ",
		  "[$partner_seqid:$partner_divider] and it came from [",
		  $div_map->{$partner_seqid}->{$partner_divider}->{PAIR_ID},
		  "].");
	  }
      }
  }

sub getPartnerSeqID
  {
    my $seqid   = $_[0];
    my $pair_id = $_[1];
    my $data    = $_[2];

    my $part_seqids = [grep {$_ ne $seqid} keys(%{$data->{$pair_id}->{ALNS}})];

    if(scalar(@$part_seqids) > 1)
      {error("Sequence ID [$seqid] is not a member of alignment: [$pair_id].")}

    return($part_seqids->[0]);
  }

sub dividersExist
  {
    my $seqid         = $_[0];
    my $lesser_bound  = $_[1];
    my $greater_bound = $_[2];
    my $div_map       = $_[3];
    my $data          = $_[4];
    my $looking_left  = $_[5]; #When choosing among multiple source dividers,
                               #choose the one closest to where you're looking
                               #from (if you're looking left, choose the right-
                               #most one.
    my $partners      = $_[6];
    my $pair_id       = $_[7]; #The pair for the segment starting at the
                               #greater bound that will be the basis for
                               #sourcing or recoding

    my $partner_seqid = getPartnerSeqID($seqid,$pair_id,$data);

    #First check myself...
    my $self_divider =
      dividerExists($seqid,
		    $lesser_bound,
		    $greater_bound,
		    $partner_seqid,
		    $pair_id,
		    $div_map,
		    $data,
		    $looking_left,
		    1);

    if($self_divider)
      {return($self_divider)}

    #Now check my partners
    foreach my $partner (@$partners)
      {
	#This should be converted back
	my $partner_divider = dividerExists($seqid,
					    $lesser_bound,
					    $greater_bound,
					    $partner->[1],
					    $partner->[0],
					    $div_map,
					    $data,
					    $looking_left,
					    0);

	if($partner_divider)
	  {return($partner_divider)}
      }

    return(0);
  }

#Checks the bounds of a given sequence against the partner sequence's converted
#coordinates to see if a divider already exists there.  Bounds are assumed to be the closest.  Only 1 or 0 dividers should exist, otherwise, warn.  It returns the divider coordinate OR (from the original sequence) or 0 (which is not a valid divider coordinate) so it can be used as a boolean.
sub dividerExists
  {
    my $seqid         = $_[0];
    my $lesser_bound  = $_[1];
    my $greater_bound = $_[2];
    my $partner_seqid = $_[3];
    my $pair_id       = $_[4];
    my $div_map       = $_[5];
    my $data          = $_[6];
    my $looking_left  = $_[7]; #When choosing among multiple source dividers,
                               #choose the one closest to where you're looking
                               #from (if you're looking left, choose the right-
                               #most one).
    my $check_self    = $_[8];

    if(calcFrame($lesser_bound) == 3)
      {$lesser_bound++}
    elsif(calcFrame($lesser_bound) == 2)
      {error("Bounds must be in either frame 1 or 3.")}

    if(calcFrame($greater_bound) == 3)
      {$greater_bound++}
    elsif(calcFrame($greater_bound) == 2)
      {error("Bounds must be in either frame 1 or 3.")}

    #If we're not checking a partner
    if($check_self)
      {
	my $existing_divs = [sort {$looking_left ? $b <=> $a : $a <=> $b}
			     grep {$_ >= $lesser_bound && $_ <= $greater_bound}
			     keys(%{$div_map->{$seqid}})];

	#If multiple dividers were found and at least 1 is nonsource and occurs
	#after a source divider (i.e. it interrupts a source divider), issue a
	#warning.
	if(scalar(@$existing_divs) > 1 &&
	   scalar(grep {(!defined($div_map->{$seqid}->{$existing_divs->[$_]}
				  ->{TYPE}) ||
			 $div_map->{$seqid}->{$existing_divs->[$_]}->{TYPE} ne
			 'SOURCE') && $_ > 0 &&
			   defined($div_map->{$seqid}
				   ->{$existing_divs->[$_ - 1]}->{TYPE}) &&
				     $div_map->{$seqid}
				       ->{$existing_divs->[$_ - 1]}->{TYPE} eq
					 'SOURCE'}
		  (0..$#{$existing_divs})))
	  {
	    my $divwarnkey = "$seqid:$lesser_bound-$greater_bound";
	    $divider_warnings->{$divwarnkey}->{MSG} =
	      join('',("Multiple dividers found between segment boundaries: ",
		       "[$seqid:$lesser_bound] and ",
		       "[$seqid:$greater_bound] that are deemed ",
		       "'closest': [",
		       join(' ',map {"$seqid:$_(" .
				       (defined($div_map->{$seqid}->{$_}
						->{TYPE}) ?
					$div_map->{$seqid}->{$_}->{TYPE} :
					'RECODE') . ")"}
			    @$existing_divs),"]."));
	    $divider_warnings->{$divwarnkey}->{SEQID}    = $seqid;
	    $divider_warnings->{$divwarnkey}->{DIVIDERS} = [@$existing_divs];
	    $divider_warnings->{$divwarnkey}->{ORDER}    =
	      scalar(keys(%$divider_warnings));
	    return(wantarray ? @$existing_divs : $existing_divs->[0]);
	  }
	elsif(scalar(@$existing_divs))
	  {return(wantarray ? @$existing_divs : $existing_divs->[0])}

	return(wantarray ? () : 0);
      }

    my $partner_lesser = convertDivider($seqid,
					$lesser_bound,
					$pair_id,
					$partner_seqid,
					$data);

    my $partner_lesser_gaps = $partner_lesser;

    if($lesser_bound > 1)
      {
	#If the query sequence has gaps to the immediate left in its
	#alignment with this partner sequence, include that partner
	#sequence's gap-aligned characters in the existing divider
	#search.  We can check this by converting the previous
	#coordinate and see if it is also 1 away
	$partner_lesser_gaps = convertAlnSeqCoord($seqid,
						  $lesser_bound - 1,
						  $pair_id,
						  $partner_seqid,
						  $data);
	if($partner_lesser_gaps < ($partner_lesser - 1))
	  {
	    #Add 1 to put it in frame 1
	    $partner_lesser_gaps++;
	  }
      }

    my $partner_greater = convertDivider($seqid,
					 $greater_bound,
					 $pair_id,
					 $partner_seqid,
					 $data);

    my $existing_divs = [sort {$looking_left ? $b <=> $a : $a <=> $b}
			 grep {my $c=0;if($_ >= $partner_lesser_gaps &&
					  $_ <= $partner_greater)
				 {$c=convertDivider($partner_seqid,  #Conv back
						    $_,
						    $pair_id,
						    $seqid,
						    $data)}
			       #(
				#A divider falls withing the search bounds
				#(
				 $_ >= $partner_lesser_gaps &&
				 $_ <= $partner_greater
				#) ||
				##A divider's extension to its hard-stop
				##overlaps the search bounds (if defined)
				#(defined($div_map->{$partner_seqid}->{$_}
				#	 ->{HARD_STOP}) &&
				# $_ < $partner_lesser_gaps &&
				# $div_map->{$partner_seqid}->{$_}
				# ->{HARD_STOP} >= $partner_lesser))
				   &&
				   #Make sure the divider returned is inside
				   #the original bounds
				   ($c == 0 || ($c >= $lesser_bound &&
						$c <= $greater_bound))}
			 keys(%{$div_map->{$partner_seqid}})];

    my $existing_orig_divs = [map {convertDivider($partner_seqid,  #Conv back
						  $_,
						  $pair_id,
						  $seqid,
						  $data)}
			      @$existing_divs];

    debug({LEVEL => 5},"Dividers already exist: [",
	  join(',',map {"$partner_seqid:$_"} @$existing_divs),"] ",
	  "between $seqid:$lesser_bound and $seqid:$greater_bound in ",
	  "$partner_seqid:$partner_lesser_gaps and $partner_seqid:",
	  "$partner_greater according to alignment $pair_id.  Converted back ",
	  "to: [",join(',',map {"$seqid:$_"} @$existing_orig_divs),"].");

    #If multiple dividers were found and at least 1 is nonsource and occurs
    #after a source divider (i.e. it interrupts a source divider), issue a
    #warning.
    if(scalar(@$existing_divs) > 1 &&
       scalar(grep {(!defined($div_map->{$partner_seqid}
			      ->{$existing_divs->[$_]}->{TYPE}) ||
		     $div_map->{$partner_seqid}->{$existing_divs->[$_]}
		     ->{TYPE} ne 'SOURCE') && $_ > 0 &&
		       defined($div_map->{$partner_seqid}
			       ->{$existing_divs->[$_ - 1]}->{TYPE}) &&
				 $div_map->{$partner_seqid}
				   ->{$existing_divs->[$_ - 1]}->{TYPE} eq
				     'SOURCE'}
		  (0..$#{$existing_divs})))
      {
	my $divwarnkey = "$partner_seqid:$partner_lesser-$partner_greater";
	$divider_warnings->{$divwarnkey}->{MSG} =
	  join('',("Multiple partner dividers found between segment ",
		   "boundaries: [$partner_seqid:$partner_lesser] and ",
		   "[$partner_seqid:$partner_greater] that are deemed ",
		   "'closest': [",
		   join(' ',map {"$partner_seqid:$_(" .
				   (defined($div_map->{$partner_seqid}
					    ->{$_}->{TYPE}) ?
				    $div_map->{$partner_seqid}->{$_}
				    ->{TYPE} : 'RECODE') . ")"}
			@$existing_divs),"]."));
	$divider_warnings->{$divwarnkey}->{DETAIL} =
	  join('',("Bounds were converted from [$seqid:$lesser_bound] and ",
		   "[$seqid:$greater_bound] and dividers were converted back ",
		   "to: [",
		   join(' ',map {"$seqid:$_"} @$existing_orig_divs),"]."));
	$divider_warnings->{$divwarnkey}->{SEQID} = $seqid;
	$divider_warnings->{$divwarnkey}->{DIVIDERS} =
	  [@$existing_orig_divs];
	$divider_warnings->{$divwarnkey}->{ORDER} =
	  scalar(keys(%$divider_warnings));
	return(wantarray ? @$existing_orig_divs : $existing_orig_divs->[0]);
      }
    elsif(scalar(@$existing_divs))
      {return(wantarray ? @$existing_orig_divs : $existing_orig_divs->[0])}

    return(wantarray ? () : 0);
  }

sub getDividerConflict
  {
    my $seqid         = $_[0];
    my $lesser_bound  = $_[1];
    my $greater_bound = $_[2];
    my $partner_seqid = $_[3];
    my $pair_id       = $_[4];
    my $div_map       = $_[5];
    my $data          = $_[6];
    my $looking_left  = $_[7]; #looking left = creating source divider in seqid
    my $look_in_self  = $_[8];

    if(calcFrame($lesser_bound) == 3)
      {$lesser_bound++}
    elsif(calcFrame($lesser_bound) == 2)
      {error("Bounds must be in either frame 1 or 3.")}

    if(calcFrame($greater_bound) == 3)
      {$greater_bound++}
    elsif(calcFrame($greater_bound) == 2)
      {error("Bounds must be in either frame 1 or 3.")}

    if(!defined($partner_seqid))
      {
	error("Partner sequence ID is required.");
	return(undef);
      }

    my($partner_lesser,$partner_lesser_gaps,$partner_greater);

    #If we're not checking a partner
    if($look_in_self)
      {
	#Note, the dividers' pair_ids may be the same.  It's still a conflict
	#and the divider to the left needs to be changed to a recode divider
	#based on the pair_id of the identity segment to the left
	my $conflicting_divs =  [sort {$looking_left ? $b <=> $a : $a <=> $b}
				 grep {#A divider's extension to its hard-stop
				   #overlaps the search bounds
				   defined($div_map->{$seqid}->{$_}
					   ->{HARD_STOP}) &&
					     $_ < $lesser_bound &&
					       $div_map->{$seqid}->{$_}
						 ->{HARD_STOP} >=
						   $lesser_bound &&
						     #The divider is for a
						     #region that is longer
						     #than 1 codon
						     ($div_map->{$seqid}->{$_}
						      ->{HARD_STOP} - $_) > 2}
				 keys(%{$div_map->{$seqid}})];

	if(scalar(@$conflicting_divs))
	  {
	    my $divwarnkey = "CONFLICT:$seqid:$lesser_bound-$greater_bound";

	    $partner_lesser = convertDivider($seqid,
					     $lesser_bound,
					     $pair_id,
					     $partner_seqid,
					     $data);

	    $partner_lesser_gaps = $partner_lesser;

	    if($lesser_bound > 1)
	      {
		#If the query sequence has gaps to the immediate left in its
		#alignment with this partner sequence, include that partner
		#sequence's gap-aligned characters in the existing divider
		#search.  We can check this by converting the previous
		#coordinate and see if it is also 1 away
		$partner_lesser_gaps = convertAlnSeqCoord($seqid,
							  $lesser_bound - 1,
							  $pair_id,
							  $partner_seqid,
							  $data);
		if($partner_lesser_gaps < ($partner_lesser - 1))
		  {
		    #Add 1 to put it in frame 1
		    $partner_lesser_gaps++;
		  }
	      }

	    $partner_greater = convertDivider($seqid,
					      $greater_bound,
					      $pair_id,
					      $partner_seqid,
					      $data);

	    my $other_pair = $div_map->{$seqid}->{$conflicting_divs->[0]}
	      ->{PAIR_ID};
	    my $other_seqid = getPartnerSeqID($seqid,$other_pair,$data);
	    $divider_warnings->{$divwarnkey}->{MSG} =
	      join('',
		   ("Conflicting self ",
		    (defined($div_map->{$seqid}->{$conflicting_divs->[0]}
			     ->{TYPE}) ?
		     $div_map->{$seqid}->{$conflicting_divs->[0]}->{TYPE} :
		     'RECODE')," divider(s: [",scalar(@$conflicting_divs),
		    "]) found overlapping segment ",
		    "boundaries: [$seqid:$lesser_bound] and ",
		    "[$seqid:$greater_bound] (while creating a ",
		    ($looking_left ? 'SOURCE' : 'RECODE')," divider in ",
		    "[$seqid]) that are deemed 'closest': [",
		    join(' ',map {"$seqid:$_-HardStop:" .
				    $div_map->{$seqid}->{$_}->{HARD_STOP} .
				      "(" .
					(defined($div_map->{$seqid}->{$_}
						 ->{TYPE}) ?
					 $div_map->{$seqid}->{$_}->{TYPE} :
					 'RECODE') . ")"}
			 @$conflicting_divs),"].  This can",
		    ($div_map->{$seqid}->{$conflicting_divs->[0]}->{HARD_START}
		     > $lesser_bound ? '' : 'not')," be fixed via a swap ",
		    "because the start of the next segment after the ",
		    "conflicting divider is located at [",
		    $div_map->{$seqid}->{$conflicting_divs->[0]}->{HARD_START},
		    "] needs to be greater than or equal to the beginning of ",
		    "the search range [$lesser_bound].  Alignment snippets ",
		    "from the query alignment and the first divider that had ",
		    "a conflict (capitalization represents the search bounds ",
		    "- all other sequence is added based on source dividers ",
		    "that overlap the search range):\n\nCapitalization in ",
		    "the first sequence is the conflicting divider's ",
		    "position and hard stop.  Capitalization in the second ",
		    "sequence is the converted search range.\n",
		    getAlnSnippet($seqid,
				  $conflicting_divs->[0],
				  $div_map->{$seqid}->{$conflicting_divs->[0]}
				  ->{HARD_STOP},
				  $pair_id,
				  $partner_seqid,
				  $data,
				  $lesser_bound,
				  $greater_bound,
				  $partner_lesser,
				  $partner_greater),"\n\nCapitalization in ",
		    "the first sequence is the search range.  Capitalization ",
		    "in the second sequence is the converted conflicting ",
		    "divider's position and hard stop.\n",
		    getAlnSnippet($seqid,
				  $conflicting_divs->[0],
				  $div_map->{$seqid}->{$conflicting_divs->[0]}
				  ->{HARD_STOP},
				  $other_pair,
				  $other_seqid,
				  $data,
				  $lesser_bound,
				  $greater_bound,
				  convertDivider($seqid,
						 $conflicting_divs->[0],
						 $other_pair,
						 $other_seqid,
						 $data),
				  convertAlnSeqCoord($seqid,
						     $div_map->{$seqid}
						     ->{$conflicting_divs->[0]}
						     ->{HARD_STOP},
						     $other_pair,
						     $other_seqid,
						     $data)),"\n\n"));
	    $divider_warnings->{$divwarnkey}->{DETAIL}   = '';
	    $divider_warnings->{$divwarnkey}->{SEQID}    = $seqid;
	    $divider_warnings->{$divwarnkey}->{DIVIDERS} =
	      [@$conflicting_divs];
	    $divider_warnings->{$divwarnkey}->{ORDER}    =
	      scalar(keys(%$divider_warnings));

	    #Error-check the conflicting divider
	    getPartnerSeqID($seqid,$other_pair,$data);

	    return({SEQID   => $seqid,
		    DIVIDER => $conflicting_divs->[0],
		    PAIR_ID => $other_pair});
	  }

	return(undef);
      }

    $partner_lesser = convertDivider($seqid,
				     $lesser_bound,
				     $pair_id,
				     $partner_seqid,
				     $data);

    $partner_lesser_gaps = $partner_lesser;

    if($lesser_bound > 1)
      {
	#If the query sequence has gaps to the immediate left in its
	#alignment with this partner sequence, include that partner
	#sequence's gap-aligned characters in the existing divider
	#search.  We can check this by converting the previous
	#coordinate and see if it is also 1 away
	$partner_lesser_gaps = convertAlnSeqCoord($seqid,
						  $lesser_bound - 1,
						  $pair_id,
						  $partner_seqid,
						  $data);
	if($partner_lesser_gaps < ($partner_lesser - 1))
	  {
	    #Add 1 to put it in frame 1
	    $partner_lesser_gaps++;
	  }
      }

    $partner_greater = convertDivider($seqid,
				      $greater_bound,
				      $pair_id,
				      $partner_seqid,
				      $data);

    #Note, the dividers' pair_ids may be the same.  It's still a conflict and
    #the divider to the left needs to be changed to a recode divider based on
    #the pair_id of the identity segment to the left
    my $conflicting_divs =  [sort {$looking_left ? $b <=> $a : $a <=> $b}
			     grep {#A divider's extension to its hard-stop
			       #overlaps the search bounds
			       defined($div_map->{$partner_seqid}->{$_}
				       ->{HARD_STOP}) &&
					 $_ < $partner_lesser &&
					   $div_map->{$partner_seqid}->{$_}
					     ->{HARD_STOP} >=
					       $partner_lesser &&
						 #The divider is for a region
						 #that is longer than 1 codon
						 ($div_map->{$partner_seqid}
						  ->{$_}->{HARD_STOP} - $_) >
						    2}
			     keys(%{$div_map->{$partner_seqid}})];

    my $conflicting_orig_divs =  [map {convertDivider($partner_seqid,#Conv back
						      $_,
						      $pair_id,
						      $seqid,
						      $data)}
				  @$conflicting_divs];

    if(scalar(@$conflicting_divs))
      {
	my $divwarnkey = "CONFLICT:$partner_seqid:$partner_lesser-" .
	  $partner_greater;
	my $other_pair = $div_map->{$partner_seqid}->{$conflicting_divs->[0]}
	  ->{PAIR_ID};
	my $other_seqid = getPartnerSeqID($partner_seqid,$other_pair,$data);
	$divider_warnings->{$divwarnkey}->{MSG} =
	  join('',
	       ("Conflicting partner ",
		(defined($div_map->{$partner_seqid}->{$conflicting_divs->[0]}
			 ->{TYPE}) ?
		 $div_map->{$partner_seqid}->{$conflicting_divs->[0]}->{TYPE} :
		 'RECODE')," divider(s: [",scalar(@$conflicting_divs),"]) ",
		"found overlapping segment boundaries: [$partner_seqid:",
		"$partner_lesser] and [$partner_seqid:$partner_greater] ",
		"(while creating an unknown divider while looking ",
		($looking_left ? 'left' : 'right')," of a segment in ",
		"[$seqid/$partner_seqid]) that are deemed 'closest': [",
		join(' ',map {$div_map->{$partner_seqid}->{$_}->{PAIR_ID} .
				":$partner_seqid:$_-HardStop:" .
				  $div_map->{$partner_seqid}->{$_}
				    ->{HARD_STOP} . "(" .
				      (defined($div_map->{$partner_seqid}
					       ->{$_}->{TYPE}) ?
				       $div_map->{$partner_seqid}->{$_}
				       ->{TYPE} : 'RECODE') . ")"}
		     @$conflicting_divs),"].  This can",
		($div_map->{$partner_seqid}->{$conflicting_divs->[0]}
		 ->{HARD_START} > $partner_lesser ? '' : 'not')," be fixed ",
		"via a swap because the start of the next segment after the ",
		"conflicting divider is located at [",
		$div_map->{$partner_seqid}->{$conflicting_divs->[0]}
		->{HARD_START},"] needs to be greater than or equal to the ",
		"beginning of the search range [$partner_lesser].  Alignment ",
		"snippets from the query alignment and the first partner ",
		"that had a conflict (capitalization represents the search ",
		"bounds - all other sequence is added based on source ",
		"dividers that overlap the search range):\n\n",
		getAlnSnippet($partner_seqid,
			      $conflicting_divs->[0],
			      $div_map->{$partner_seqid}
			      ->{$conflicting_divs->[0]}->{HARD_STOP},
			      $pair_id,
			      $seqid,
			      $data,
			      $partner_lesser,
			      $partner_greater,
			      $lesser_bound,
			      $greater_bound),"\n\n",
		getAlnSnippet($partner_seqid,
			      $conflicting_divs->[0],
			      $div_map->{$partner_seqid}
			      ->{$conflicting_divs->[0]}->{HARD_STOP},
			      $other_pair,
			      $other_seqid,
			      $data,
			      $partner_lesser,
			      $partner_greater),"\n\n"));
	$divider_warnings->{$divwarnkey}->{DETAIL} =
	  join('',("Bounds were converted from [$seqid:$lesser_bound] and ",
		   "[$seqid:$greater_bound] and the conflicting divider ",
		   "converts to: [",
		   join(' ',map {"$seqid:$_"} @$conflicting_orig_divs),"]."));
	$divider_warnings->{$divwarnkey}->{SEQID} = $seqid;
	$divider_warnings->{$divwarnkey}->{DIVIDERS} =
	  [@$conflicting_orig_divs];
	$divider_warnings->{$divwarnkey}->{ORDER} =
	  scalar(keys(%$divider_warnings));

	#Error-check the conflicting divider
	getPartnerSeqID($partner_seqid,$other_pair,$data);

	return({SEQID   => $partner_seqid,
		DIVIDER => $conflicting_divs->[0],
		PAIR_ID => $other_pair});
      }

    return(undef);
  }

#Builds a 2 line string from an alignment using the sequence coordinates of a
#query sequence.  Also capitalizes sequence in a given sequence coordinate
#range in both of the sequences.
sub getAlnSnippet
  {
    my $seqid          = $_[0];
    my $lesser_bound   = $_[1];
    my $greater_bound  = $_[2];
    my $pair_id        = $_[3];
    my $partner_seqid  = $_[4];
    my $data           = $_[5];
    my $cap_reg_start  = (defined($_[6]) ? $_[6] : 0);
    my $cap_reg_stop   = (defined($_[7]) ? $_[7] : 0);
    my $pcap_reg_start = (defined($_[8]) ? $_[8] : 0);
    my $pcap_reg_stop  = (defined($_[9]) ? $_[9] : 0);

    my $query_pos     = 0;
    my $partner_pos   = 0;

    my $query_aln     = '';
    my $partner_aln   = '';

    my $aln_pos       = 0;
    my $aln_start     = 0;
    my $aln_stop      = 0;

    my $partner_start = 0;
    my $partner_stop  = 0;

    #Expand to encompass all coords
    if($cap_reg_start > 0 && $cap_reg_start < $lesser_bound)
      {$lesser_bound = $cap_reg_start}
    if($cap_reg_stop > 0 && $cap_reg_stop > $greater_bound)
      {$greater_bound = $cap_reg_stop}

    while($query_pos < $greater_bound)
      {
	my $query_char =
	  substr($data->{$pair_id}->{ALNS}->{$seqid},$aln_pos,1);

	if($query_char ne '-')
	  {$query_pos++}

	$query_char = ($cap_reg_start > 0 &&
		       $query_pos >= $cap_reg_start &&
		       $query_pos <= $cap_reg_stop ?
		       uc($query_char) : lc($query_char));

	my $partner_char =
	  substr($data->{$pair_id}->{ALNS}->{$partner_seqid},$aln_pos,1);

	if($partner_char ne '-')
	  {$partner_pos++}

	$partner_char = ($pcap_reg_start > 0 &&
			 $partner_pos >= $pcap_reg_start &&
			 $partner_pos <= $pcap_reg_stop ?
			 uc($partner_char) : lc($partner_char));

	if($query_pos >= $lesser_bound && $query_pos <= $greater_bound)
	  {
	    $query_aln   .= $query_char;
	    $partner_aln .= $partner_char;

	    if($query_pos == $lesser_bound)
	      {
		$aln_start = $aln_pos + 1;
		$partner_start = $partner_pos + ($partner_char eq '-' ? 1 : 0);
	      }
	    $partner_stop = $partner_pos;

	    $aln_stop = $aln_pos + 1;
	  }

	$aln_pos++;
      }

    my $snippet = ("$query_aln $seqid:$lesser_bound-$greater_bound " .
		   "[$seqid:$cap_reg_start-$cap_reg_stop] " .
		   "$pair_id:$aln_start-$aln_stop\n$partner_aln " .
		   "$partner_seqid:$partner_start-$partner_stop" .
		   ($pcap_reg_start > 0 ?
		    " [$partner_seqid:$pcap_reg_start-$pcap_reg_stop]" : ''));

    return($snippet);
  }

#Takes a sequence ID, a coordinate, a pair ID, and the sequence ID of the
#partner and returns the corresponding coordinate in the partner according to
#the alignment.  This uses a linear search.  Could be better but would need a
#lot more work.  If the query coordinate is the full length of the query
#sequence, return the length of the target sequence.  If it's greater than the
#query sequence, return the length of the target sequence plus 1.  These rules
#are made so that if the alignment has a gap at the end, we will know how to
#encode the remaining sequence.
sub convertDivider
  {
    my $query_seqid  = $_[0];
    my $query_coord  = $_[1];
    my $pairid       = $_[2];
    my $target_seqid = $_[3];
    my $data         = $_[4];

    #The beginning of a sequence should always create a divider at the
    #beginning of the other sequence, so return 1 instead of counting through
    #bases that align with a beginning gap in the query.
    if($query_coord == 1)
      {return(1)}

    #A divider must be in frame 1 because if we divide the sequence source mid-
    #codon, it could change the encoded AA
    if($query_coord % 3 != 1) #If the coordinate is not in frame 1
      {error("Coordinate: [$query_coord] must be in frame 1.  Frame [",
	     ((($query_coord - 1) % 3) + 1),"] was sent in.")}

    debug("Checking if coord [$query_coord] is larger than [$query_seqid]'s ",
	  "length as it exists in [$pairid].  Data for [$pairid] ",
	  (exists($data->{$pairid}) ? 'exists' : 'does not exist'),".  Data ",
	  "for SEQS ",
	  (exists($data->{$pairid}) && exists($data->{$pairid}->{SEQS}) ?
	   'exists' : 'does not exist'),".  Data for [$query_seqid] ",
	  (exists($data->{$pairid}) && exists($data->{$pairid}->{SEQS}) &&
	   exists($data->{$pairid}->{SEQS}->{$query_seqid}) ?
	   'exists but is not defined' : 'does not exist'),'.')
      if(!defined($data->{$pairid}->{SEQS}->{$query_seqid}));

    if($query_coord >= length($data->{$pairid}->{SEQS}->{$query_seqid}))
      {
	if(length($data->{$pairid}->{SEQS}->{$query_seqid}) % 3 != 0)
	  {error("Length of sequence [$query_seqid]: [",
		 length($data->{$pairid}->{SEQS}->{$query_seqid}),
		 "] is not a multiple of 3!")}
	return(length($data->{$pairid}->{SEQS}->{$target_seqid}) + 1);
      }

    my $query_pos       = 0;
    my $target_pos      = 0;
    my $aln_pos         = 0;
    my $query_gaps      = 0;
    my $last_query_gaps = 0;

    while($query_pos < $query_coord)
      {
	if(substr($data->{$pairid}->{ALNS}->{$query_seqid},
		  $aln_pos,1) ne '-')
	  {
	    $query_pos++;
	    $query_gaps = 0;
	  }
	else
	  {$query_gaps++}

	my $target_char = substr($data->{$pairid}->{ALNS}->{$target_seqid},
				 $aln_pos,1);
	if($target_char ne '-')
	  {$target_pos++}

	#If we are at the query coordinate we were looking for and it
	#corresponds to a gap character in the target alignment sequence,
	#increment the target position, because the last target position was
	#aligned with a different query position.  Note, this only works for
	#the start of segments.  End coordinates should not be queried.
	if($query_pos == $query_coord && $target_char eq '-')
	  {$target_pos++}

	$last_query_gaps = $query_gaps;
	$aln_pos++;
      }

    if($target_pos % 3 != 1)
      {error("Coordinate computed: [$target_pos] is not in frame 1.  Frame [",
	     ((($target_pos - 1) % 3) + 1),"] was calculated.")}

    return($target_pos);
  }

#This subroutine takes a sequence coordinate (not including gap characters) of
#one sequence and finds the corresponding partner coordinate in an alignment.
#If a coordinate aligns with a gap, the coordinate of the last real sequence
#character of the partner is returned.  This differs from convertDivider, which
#attempts to include end-gaps, e.g. if convertDivider is given coordinate 1,
#which aligns with its partner at coordinate 15, it will return coordinate 1,
#to include the end-gap in front of coordinate 1 of the query sequence.  This
#method on the other hand, will return 15.  Also e.g., if convertDivider is
#given coordinate 1000, and that corresponds to a gap character in the partner,
#The returned target coordinate is incremented 1 more time so that the entire
#end gap is included.  Note, this assumes that convertDivider only takes frame
#1 positions and it's assumed that whole codons are aligned and the divider
#must always be in frame 1.
sub convertAlnSeqCoord
  {
    my $query_seqid  = $_[0];
    my $query_coord  = $_[1];
    my $pairid       = $_[2];
    my $target_seqid = $_[3];
    my $data         = $_[4];

    my $query_pos  = 0;
    my $target_pos = 0;
    my $aln_pos    = 0;

    while($query_pos < $query_coord)
      {
	if(substr($data->{$pairid}->{ALNS}->{$query_seqid},
		  $aln_pos,1) ne '-')
	  {$query_pos++}

	my $target_char = substr($data->{$pairid}->{ALNS}->{$target_seqid},
				 $aln_pos,1);
	if($target_char ne '-')
	  {$target_pos++}

	$aln_pos++;
      }

    return($target_pos);
  }

#Assumes the last divider is 1 past the sequence length
sub fillMapStops
  {
    my $div_map = $_[0];

    foreach my $seqid (sort {$a cmp $b} keys(%$div_map))
      {
	#Obtain an array of dividers(/start coordinates for each section)
	my $div_array = [sort {$a <=> $b} keys(%{$div_map->{$seqid}})];

	#For every index in the array except the last (which is assumed to be
	#a dummy marker for 1 past the end of the sequence)
	for(my $i = 0;$i < $#{$div_array};$i++)
	  {
	    my $start = $div_array->[$i];
	    my $stop  = $div_array->[$i + 1] - 1;
	    $div_map->{$seqid}->{$start}->{STOP} = $stop;
	  }

	#Remove the last divider, as we don't need that placeholder anymore
	delete($div_map->{$seqid}->{$div_array->[-1]});
      }
  }

#This returns a hash keyed on sequence ID whose values are arrays containing
#Each pair ID of the alignments the sequence is involved in along with the
#sequence ID it was aligned/optimized with (contained in 2-member sub-arrays)
sub getPartnersHash
  {
    my $soln  = $_[0];
    my $data  = $_[1];
    my $partners_hash = {};

    foreach my $seqid (sort {$a cmp $b} keys(%$soln))
      {$partners_hash->{$seqid} = getPartnerPairSeqIDs($seqid,$soln,$data)}

    return($partners_hash);
  }

#Extracts the alignments/pairIDs that were included in the solution and uses
#those pairIDs to obtain the sequence IDs from the data hash.  One could search
#all the segments to find this out, but that's more complicated
sub getPartnerPairSeqIDs
  {
    my $seqid = $_[0];
    my $soln  = $_[1];
    my $data  = $_[2];

    my $pair_ids = {};
    foreach my $pairid (map {$_->{PAIR_ID}} @{$soln->{$seqid}})
      {$pair_ids->{$pairid} = 0}

    my $partners = [sort {$a->[0] cmp $b->[0] || $a->[1] cmp $b->[1]}
		    grep {$_->[1] ne $seqid}
		    map {my $pid=$_;
			 map {[$pid,$_]} sort {$a cmp $b}
			   keys(%{$data->{$pid}->{SEQS}})}
		    keys(%$pair_ids)];

    return(wantarray ? @$partners : $partners);
  }

#This takes a map and fills in the sequence by grabbing it either directly from
#a sequence's version from an alignment (as generated by codonHomologizer) or
#it recodes the sequence to best match a partner sequence it was aligned with,
#but was encoded differently due to a segment match with another sequence.  The
#map structure is: $div_map->{$seqid}->{$divider} =
#{PAIR_ID => $id,
# TYPE    => 'SOURCE' or 'RECODE' (from what the partner sequence actually is),
# STOP    => $stop,
# HARD_START => $start,       #This is used to swap RECODE & SOURCE dividers
# HARD_STOP  => $stop,        #This is only used to make sure a RECODE divider
                              #after a SOURCE divider (interrupting it)
# IDEN_START => $iden_start,  #This is so we can represent the identity segs as
# IDEN_STOP  => $iden_stop,   #capital letters
# SEQ     => $seq} (This is what this subroutine fills in.)
#There should be a divider for the first coordinate and the coordinate after
#the last coordinate (i.e. sequence length + 1).  Forst coordinate is 1, not 0.
sub weaveSeqs
  {
    my $map           = $_[0];
    my $data          = $_[1];
    my $partners_hash = $_[2];
    my $matrix        = $_[3];
    my $codon_hash    = $_[4];
    my $soln          = $_[5];

    my $num_filled = 1;
    my $unfinished = 1;
    my $unf_str    = '';

    my $countdowns = {map {$_ => scalar(keys(%{$map->{$_}}))} keys(%$map)};

    #While all segments' alignment sources have not been determined and we're
    #still making progress
    while($unfinished && $num_filled)
      {
	#We will assume that this pass through the loop will complete the map
	#unless we hit something we cannot fill in
	$unfinished = 0;
	$unf_str    = '';
	#We will require that something has to have been added to the map on
	#each iteration.
	$num_filled = 0;

	#For each each sequence
	foreach my $seqid (sort {$a cmp $b} keys(%$map))
	  {
	    my $laststop = 0;
	    my $seqlen   = 0;

	    #Cycle through the subsequence records
	    foreach my $divider (sort {$a <=> $b} keys(%{$map->{$seqid}}))
	      {
		#The length of any version of the sequence will do, as long as
		#this sequence is a part of the pair the current divider was
		#copied from.  Note, this will fail if no identity segments
		#were found in a pair.
		if(exists($data->{$map->{$seqid}->{$divider}->{PAIR_ID}}
			  ->{SEQS}->{$seqid}) && $seqlen == 0)
		  {$seqlen = length($data->{$map->{$seqid}->{$divider}
					    ->{PAIR_ID}}->{SEQS}
				    ->{$seqid})}

		if($laststop + 1 != $divider)
		  {error("Map incomplete for sequence: [$seqid].  Missing ",
			 "section between [$laststop-$divider].")}
		$laststop = $map->{$seqid}->{$divider}->{STOP};

		if($map->{$seqid}->{$divider}->{SEQ} eq '')
		  {
		    #If the source type of the sequence is RECODE, we need to
		    #determine what sequence this sequence's recoding should be
		    #based on.  We will arbitrarily pick the forst completely
		    #overlapping and filled in partner
		    if(!defined($map->{$seqid}->{$divider}->{TYPE}) ||
		       $map->{$seqid}->{$divider}->{TYPE} eq 'RECODE')
		      {
			if(reCodeLeftoverSegment($seqid,
						 $divider,
						 $map,
						 $partners_hash->{$seqid},
						 $matrix,
						 $data,
						 $codon_hash))
			  {
			    $num_filled++;
			    $countdowns->{$seqid}--;
			    verboseOverMe("Sourcing.  Unfinished segments: [",
					  join(' ',map {"$_:$countdowns->{$_}"}
					       sort {$a cmp $b}
					       keys(%$countdowns)),
					  "] $seqid:$divider-",
					  $map->{$seqid}->{$divider}->{STOP},
					  " Recoding");
			  }
			#Could not recode this segment, so we'll need to start
			#the loop over after this round
			else
			  {
			    $unf_str .= ($unfinished ? ',' : '');
			    $unf_str .= "$seqid:$divider-" .
			      $map->{$seqid}->{$divider}->{STOP} .
				" (based on " .
				  $map->{$seqid}->{$divider}->{PAIR_ID} . ")";
			    $unfinished++;
			    verboseOverMe("Sourcing.  Unfinished segments: [",
					  join(' ',map {"$_:$countdowns->{$_}"}
					       sort {$a cmp $b}
					       keys(%$countdowns)),
					  "] ",
					  "$seqid:$divider-",
					  $map->{$seqid}->{$divider}->{STOP},
					  " TBD.  Should be based on [",
					  $map->{$seqid}->{$divider}
					  ->{PAIR_ID},"].");
			  }
		      }
		    else #TYPE is 'SOURCE' - so grab the sequence directly
		      {
			$map->{$seqid}->{$divider}->{SEQ} =
			  substr($data->{$map->{$seqid}->{$divider}->{PAIR_ID}}
				 ->{SEQS}->{$seqid},
				 ($divider - 1),
				 ($map->{$seqid}->{$divider}->{STOP} -
				  $divider + 1));
			$num_filled++;
			$countdowns->{$seqid}--;
			verboseOverMe("Sourcing.  Unfinished segments: [",
				      join(' ',map {"$_:$countdowns->{$_}"}
					   sort {$a cmp $b}
					   keys(%$countdowns)),
				      "] ",
				      "$seqid:$divider-",
				      $map->{$seqid}->{$divider}->{STOP},
				      " Sourcing");
		      }
		  }
	      }

	    if($laststop != $seqlen)
	      {error("Map incomplete for sequence: [$seqid].  Missing ",
		     "section at the end between [$laststop-$seqlen].")}
	  }
      }

    verbose("Sourcing.  Unfinished segments: [",
	    join(' ',map {"$_:$countdowns->{$_}"} sort {$a cmp $b}
		 keys(%$countdowns)),
	    "] Done");

    if($unfinished)
      {error("Unable to fill in sequence for [$unfinished] segments.  Unable ",
	     "to determine which alignment every segment of sequence should ",
	     "either be sourced from directly, or the basis for recoding (if ",
	     "it aligned to a stretch of identity with another sequence).",
	     {DETAIL => "Unfinished segments: [$unf_str]"})}

    #Assign letter codes to each sequence
    my $curcode     = 'a';
    my $sourcecodes = {};
    my $codecount   = 0;
    my $sourceseqs  = {};
    my $sourcelkup  = {};
    foreach my $seqid (sort {$a cmp $b} keys(%$map))
      {
	$codecount++;
	if($codecount > 26)
	  {
	    warning("Too many sequences to generate sourced sequences.");
	    last;
	  }
	$sourcecodes->{$seqid} = $curcode++;
	$sourceseqs->{LEGEND}->{$seqid} = $sourcecodes->{$seqid};
      }
    foreach my $pid (keys(%$data))
      {
	my($sid1,$sid2) = keys(%{$data->{$pid}->{SEQS}});
	if(!defined($sid1) || !defined($sid2))
	  {
	    warning("Unable to create source lookup.",
		    {DETAIL => ("2 sequence IDs are not defined for every " .
				"pairwise alignment.")});
	    $codecount = 27; #This turns off the sourced sequences feature
	    last;
	  }
	$sourcelkup->{$pid}->{$sid1} = $sid2;
	$sourcelkup->{$pid}->{$sid2} = $sid1;
      }

    #Now assemble the sequences
    my $seqhash = {};

    foreach my $seqid (sort {$a cmp $b} keys(%$map))
      {
	my $laststop = 0;
	foreach my $divider (sort {$a <=> $b} keys(%{$map->{$seqid}}))
	  {
	    my $start = $divider;
	    my $stop  = $map->{$seqid}->{$divider}->{STOP};
	    my $pair  = $map->{$seqid}->{$divider}->{PAIR_ID};
	    my $type  = (defined($map->{$seqid}->{$divider}->{TYPE}) &&
			 $map->{$seqid}->{$divider}->{TYPE} eq 'SOURCE' ?
			 'SOURCE' : 'RECODE');

	    verboseOverMe("Assembing.  [$seqid:$start-$stop]");

	    debug({LEVEL => 4},"$seqid:$start-$stop",
		  ($map->{$seqid}->{$divider}->{IDEN_START} > 0 ?
		   " $map->{$seqid}->{$divider}->{IDEN_START}-" .
		   "$map->{$seqid}->{$divider}->{IDEN_STOP}\t" : "\t\t"),
		  "$type $pair");

	    if(($laststop + 1) != $start)
	      {
		my $size = $start - $laststop - 1;
		error("Unable to determine sequence [$seqid] from [",
		      ($laststop + 1),"] ","to [",($start - 1),"].  This ",
		      "should not have happened.  Inserting 'Ns' so that you ",
		      "get something out.");
		$seqhash->{$seqid} .= ('N' x $size);
	      }

	    if(length($map->{$seqid}->{$divider}->{SEQ}) !=
	       ($stop - $start + 1))
	      {
		error("The length of the saved sequence [$seqid] (from ",
		      "[$start] to [$stop]): [",
		      length($map->{$seqid}->{$divider}->{SEQ}),"] is not ",
		      "the expected length [",($stop - $start + 1),"].  This ",
		      "should not have happened.  Inserting 'Ns'.",
		      {DETAIL =>
		       (join('',("PAIR_ID => $pair, TYPE => ",
				 (defined($map->{$seqid}->{$divider}->{TYPE}) ?
				  $map->{$seqid}->{$divider}->{TYPE} : 'undef')
				)))});
		$seqhash->{$seqid} .= ('N' x ($stop - $start + 1));
	      }
	    else
	      {$seqhash->{$seqid} .= $map->{$seqid}->{$divider}->{SEQ}}

	    #Update the sourced sequences
	    if($codecount < 27)
	      {
		if(!defined($sourceseqs->{SEQS}->{$seqid}))
		  {$sourceseqs->{SEQS}->{$seqid} = ''}
		if(!defined($pair))
		  {
		    warning("Alignment source not defined for all sequence ",
			    "portions.  Unable to generate sourced sequence ",
			    "string.");
		    $codecount = 27;
		    $sourceseqs = {}; #Clear out the hash, but keep it a hash
		  }
		elsif(!exists($sourcelkup->{$pair}->{$seqid}) ||
		      !defined($sourcelkup->{$pair}->{$seqid}))
		  {
		    error("Partner for sequence [$seqid:$start-$stop] not ",
			  "found for alignment [$pair].  Unable to generate ",
			  "sourced sequence string.");
		    $sourceseqs->{SEQS}->{$seqid} .=
		      '_' x ($stop - $start + 1);
		  }
		elsif($type eq 'SOURCE')
		  {$sourceseqs->{SEQS}->{$seqid} .=
		     uc($sourcecodes->{$sourcelkup->{$pair}->{$seqid}} x
			($stop - $start + 1))}
		else
		  {$sourceseqs->{SEQS}->{$seqid} .=
		     lc($sourcecodes->{$sourcelkup->{$pair}->{$seqid}} x
			($stop - $start + 1))}
	      }

	    $laststop = $stop;
	  }
      }

    if(isVerbose())
      {print STDERR "\n"}

    my $validation_output = "Result validation:\n";

    #Edit the sequences to add the alternate codons (and while we're at it,
    #we'll error-check that the sequences were created, the correct length,
    #encode the same AA sequence, and contain the included identity segments)
    foreach my $seqid (sort {$a cmp $b} keys(%$soln))
      {
	#The length of each version of a sequence should be the same in any
	#pair it's involved in, so let's arbitrarily grab one.
	my $any_pair_id  = $partners_hash->{$seqid}->[0]->[0];
	my $any_orig_seq = $data->{$any_pair_id}->{SEQS}->{$seqid};
	if(!exists($seqhash->{$seqid}) || !defined($seqhash->{$seqid}))
	  {
	    error("Sequence [$seqid] was not constructed!");
	    next;
	  }
	elsif(length($seqhash->{$seqid}) !=
	      length($any_orig_seq))
	  {error("The length of sequence [$seqid]: [",
		 length($seqhash->{$seqid}),"] is not the expected length: [",
		 length($any_orig_seq),"].")}

	#For each record that is of type 'alternate codon'
	foreach my $altrec (grep {$_->{TYPE} eq 'alt'} @{$soln->{$seqid}})
	  {
	    if(!exists($altrec->{ALT_CODON}) ||
	       !defined($altrec->{ALT_CODON}) ||
	       length($altrec->{ALT_CODON}) != 3)
	      {error("Mal-formed or missing alternate codon for sequence ",
		     "[$seqid].")}
	    elsif($altrec->{FINAL_STOP} > length($seqhash->{$seqid}))
	      {error("Woven sequence [$seqid] is too short to place ",
		     "alternate codon at position [",
		     "$altrec->{FINAL_START}-$altrec->{FINAL_STOP}].")}
	    else
	      {
		verbose("Inserting alternate codon: [$seqid:",
			"$altrec->{FINAL_START} $altrec->{ALT_CODON}]");
		substr($seqhash->{$seqid},($altrec->{FINAL_START} - 1),3,
		       $altrec->{ALT_CODON});
	      }
	  }

	#Check that all the identity segments in the solution were included
	my($bad_nts,$bad_segs) =
	  validateRepairSegments($seqid,$soln,$seqhash,$data);

	#Check that the sequences encode the same AA sequence
	my $any_orig_aa = uc(translate($any_orig_seq,     $codon_hash->{REV}));
	my $mixed_aa    = uc(translate($seqhash->{$seqid},$codon_hash->{REV}));
	if($any_orig_aa ne $mixed_aa)
	  {error("The new mixed hybrid version of sequence [$seqid]'s ",
		 "translation differs from the original alignment ",
		 "[$any_pair_id].",
		 {DETAIL => "ORIG: [$any_orig_aa]\nNEW:  [$mixed_aa]"})}
	if($bad_nts)
	  {warning("The new mixed hybrid version of sequence [$seqid] failed ",
		   "validation of the inclusion of all the identity segments ",
		   "and expected alternate codons, but was repaired.  There ",
		   "were [$bad_nts] nucleotides that were not as expected, ",
		   "and replaced to preserve the integrity of the identity ",
		   "segments.")}

	my($n,$real_nt_diffs) = numDiff($any_orig_seq,$seqhash->{$seqid},0);
	my $ambig_nts = $n - $real_nt_diffs;
	my($d,$real_aa_diffs) = numDiff($any_orig_aa,$mixed_aa,1);
	my $ambigs = $d - $real_aa_diffs;
	$validation_output .=
	  join('',("$seqid\n\t",($real_aa_diffs ? 'NOT ' : ''),
		   "Validated\tAA sequence correctness",
		   ($real_aa_diffs ? " ($real_aa_diffs aas differ)" : '^'),
		   ".\n\t",($ambigs ? 'NOT ' : ''),"Validated\tAA sequence ",
		   "completeness",($ambigs ? " ($ambigs ambiguous aas)" : ''),
		   ".\n\t",($bad_nts ? 'Repaired/' : ''),"Validated\t",
		   "Identity segment inclusion",
		   ($bad_nts ? " ([$bad_nts] nts among [$bad_segs] identity " .
		    "segments repaired)" : ''),
		   ".#\n\t",($ambig_nts ? 'NOT ' : ''),"Validated\t",
		   "NT sequence.",
		   ($ambig_nts ? "  There were [$ambig_nts] ambiguous nts."
		    : ''),"  (Note: [$real_nt_diffs/",length($any_orig_seq),
		   "] nts were changed*.)\n"));
      }

    $validation_output .=
      join('',("^ = Confirmed to encode the same amino acid sequence ",
	       "compared to 1 arbitrarily selected alignment.\n* = As ",
	       "compared to 1 arbitrarily selected alignment.\n# = All ",
	       "included identity segments in the solution are confirmed to ",
	       "be present and correctly encoded from their respective ",
	       "pairwise source alignments.\n"));

    #Now let's make the identity segments capitalized (using the solution
    #because there are known issues currently with the correctness of the map)
    foreach my $seqid (sort {$a cmp $b} keys(%$soln))
      {
	$seqhash->{$seqid} = lc($seqhash->{$seqid});
	foreach my $rec (@{$soln->{$seqid}})
	  {
	    #Next if this segment doesn't containg any identity
	    next if($rec->{IDEN_START} == 0 ||
		    $rec->{IDEN_STOP} < $rec->{IDEN_START});

	    my $iden_start = $rec->{IDEN_START};
	    my $iden_stop  = $rec->{IDEN_STOP};

	    my $tmp_seq = '';
	    if($iden_start > 1)
	      {$tmp_seq .= substr($seqhash->{$seqid},0,($iden_start - 1))}
	    $tmp_seq .= uc(substr($seqhash->{$seqid},
				  ($iden_start - 1),
				  ($iden_stop - $iden_start + 1)));
	    if($iden_stop < length($seqhash->{$seqid}))
	      {$tmp_seq .= substr($seqhash->{$seqid},$iden_stop)}
	    #Sanity check
	    if(uc($seqhash->{$seqid}) ne uc($tmp_seq))
	      {error("Identity-capitalized string construction problem.  ",
		     "Sequence changed.")}
	    $seqhash->{$seqid} = $tmp_seq;
	  }
      }

    return($seqhash,$sourceseqs,$validation_output);
  }

sub validateRepairSegments
  {
    my $seqid   = $_[0];
    my $soln    = $_[1];
    my $seqhash = $_[2];
    my $data    = $_[3];

    my $bad_nts  = 0;
    my $bad_segs = 0;
    my $rep_hash = {};

    #We need to pre-process the positions in the hash because some duplicates
    #can be left over because of the way ALT_CODONs were modified to work given
    #the discovery of indirect overlap.
    foreach my $solrec (@{$soln->{$seqid}})
      {
	my $type = $solrec->{TYPE};
	my($expected_nts,$idlen,$idstart,$idstop,$start,$stop,$repair_nts,
	   $len,$all_nts);

	if($type eq 'seg')
	  {
	    #Obtain the identity characters that we expect to find in the final
	    #re-encoded sequence from the source alignment
	    $idstart       = $solrec->{IDEN_START};
	    $idstop        = $solrec->{IDEN_STOP};
	    my $pair_id    = $solrec->{PAIR_ID};
	    my $source_seq = $data->{$pair_id}->{SEQS}->{$seqid};
	    $idlen         = $idstop - $idstart + 1;
	    $expected_nts  = uc(substr($source_seq,$idstart - 1,$idlen));
	    $start         = $solrec->{FIRST_CODON_START};
	    $stop          = $solrec->{LAST_CODON_STOP};
	    $len           = $stop - $start + 1;
	    $repair_nts    = substr($source_seq,$start - 1,$len);
	    $all_nts       = uc(substr($source_seq,$start - 1,$len));
	  }
	elsif($type eq 'alt')
	  {
	    #We're going to assume that the alternate codon was selected
	    #correctly.  There might be identity characters on each side and
	    #there's no way to know which characters are a part of the identity
	    #because it wasn't saved.  So just to make things easier yet still
	    #have a check, we well test to see if the alternate codon is what
	    #we expect it to be
	    $idstart      = $solrec->{FINAL_START};
	    $idstop       = $solrec->{FINAL_STOP};
	    $idlen        = 3;
	    $expected_nts = uc($solrec->{ALT_CODON});
	    $start        = $solrec->{FINAL_START};
	    $stop         = $solrec->{FINAL_STOP};
	    $len          = 3;
	    $repair_nts   = $solrec->{ALT_CODON};
	    $all_nts      = $solrec->{ALT_CODON};
	  }

	if(!defined($idstart) || !defined($idstop) || $idstop < $idstart)
	  {
	    debug("Identity start/stop invalid.");
	    next;
	  }

	if(!exists($rep_hash->{$start}) ||
	   ($rep_hash->{$start}->{DATA}->{TYPE} eq 'seg' && $type eq 'alt'))
	  {$rep_hash->{$start}->{DATA} = {IDSTART => $idstart,
					  IDSTOP  => $idstop,
					  IDLEN   => $idlen,
					  EXPNTS  => $expected_nts,
					  STOP    => $stop,
					  LEN     => $len,
					  REPNTS  => $repair_nts,
					  TYPE    => $type,
					  ALLNTS  => $all_nts}}
	else
	  {
	    if($type eq 'alt' &&
	       $rep_hash->{$start}->{DATA}->{TYPE} eq 'alt' &&
	       $expected_nts ne $rep_hash->{$start}->{DATA}->{EXPNTS})
	      {error("Expected identity sequence conflict: [$expected_nts] ",
		     "vs [$rep_hash->{$start}->{DATA}->{EXPNTS}] at ",
		     "position [$seqid:$idstart-$idstop].")}
	  }

	if($len == 3)
	  {$rep_hash->{$start}->{COUNT}->{$type}->{$all_nts}++}
	elsif(!exists($rep_hash->{$start}->{COUNT}) ||
	      !exists($rep_hash->{$start}->{COUNT}->{$type}))
	  {$rep_hash->{$start}->{COUNT}->{$type}->{$all_nts} = 0}
      }

    #If there are only segment records at a position and there are multiple
    #different codon values for that position
    my @segconflicts = map {[$_,keys(%{$rep_hash->{$_}->{COUNT}->{seg}})]}
      grep {my $k=$_;
	    #A type of 'seg' exists in the count hash at this position
	    exists($rep_hash->{$k}->{COUNT}->{seg}) &&
	      #Only one type exists in the count hash at this position
	      scalar(keys(%{$rep_hash->{$k}->{COUNT}})) == 1 &&
		#The one type that's present has multiple different codons
		scalar(grep {scalar(keys(%$_)) > 1}
		       values(%{$rep_hash->{$k}->{COUNT}}))}
	sort {$a <=> $b} keys(%$rep_hash);
    if(scalar(@segconflicts))
      {error("The following positions in the solution for [$seqid] have ",
	     "multiple 1-codon split segment records and no alt codon ",
	     "record: [",join(';',map {join(',',@$_)} @segconflicts),"].")}

    foreach my $start (sort {$a <=> $b} keys(%$rep_hash))
      {
	my $idstart      = $rep_hash->{$start}->{DATA}->{IDSTART};
	my $idstop       = $rep_hash->{$start}->{DATA}->{IDSTOP};
	my $idlen        = $rep_hash->{$start}->{DATA}->{IDLEN};
	my $expected_nts = $rep_hash->{$start}->{DATA}->{EXPNTS};
	my $stop         = $rep_hash->{$start}->{DATA}->{STOP};
	my $len          = $rep_hash->{$start}->{DATA}->{LEN};
	my $repair_nts   = $rep_hash->{$start}->{DATA}->{REPNTS};
	my $type         = $rep_hash->{$start}->{DATA}->{TYPE};

	my $actual_nts = uc(substr($seqhash->{$seqid},$idstart - 1,$idlen));
	my $new_bad    = 0;

	for(my $pos = 0;$pos < length($actual_nts);$pos++)
	  {
	    my $exp_nt = substr($expected_nts,$pos,1);
	    my $act_nt = substr($actual_nts,$pos,1);
	    if($exp_nt ne $act_nt)
	      {
		warning("Unexpected identity segment sequence at: [$seqid:",
			($idstart+$pos),"].");
		$bad_nts++;
		$new_bad = 1;
	      }
	  }

	#Repair problems
	if($new_bad)
	  {
	    $bad_segs++;
	    substr($seqhash->{$seqid},$start - 1,$len,$repair_nts);
	  }
      }

    return($bad_nts,$bad_segs);
  }

sub outputReport
  {
    my $score       = $_[0];
    my $data        = $_[1];
    my $max_size    = $_[2];
    my $soln_stats  = $_[3];
    my $sizes       = $_[4];
    my $valid_msg   = $_[5];
    my $report_file = $_[6];

    openOut(*RPT,$report_file);

    outputSolutionScore($score);

    outputSegmentTable($soln_stats,$sizes);

    outputValidation($valid_msg);

    closeOut(*RPT);
  }

sub outputSolutionScore
  {
    my $score_obj = $_[0];

    if(!defined($score_obj) || !exists($score_obj->{STRETCHES}))
      {print("No solutions found.\n")}
    else
      {
	my $multi = scalar(keys(%{$score_obj->{STRETCHES}})) > 1;
	print("\nSolution Score",($multi ? ":\n" : ''));
	foreach my $size (sort {$b <=> $a} keys(%{$score_obj->{STRETCHES}}))
	  {
	    if($multi)
	      {print("\tStretch Size $size:\n")}
	    else
	      {print(" (Stretch size $size):\n")}

	    my $score = $score_obj->{STRETCHES}->{$size};
	    print(($multi ? "\t" : ''),
		  "\t$score->[0]\tMax number of segments included from the ",
		  "alignment with the fewest contributed segments\n",
		  ($multi ? "\t" : ''),
		  "\t$score->[1]\tTotal number of unique segments included ",
		  "among all sequences\n",
		  ($multi ? "\t" : ''),
		  "\t$score->[2]\tSegment spacing score*.\n");
	  }
	print(($multi ? "\tOverall Coverage\n" : ''),
	      ($multi ? "\t" : ''),
	      "\t$score_obj->{BASESCORE}\tTotal number of unique bases ",
	      "included in an identity segment of any size.\n\n");
	print("* Int(sqrt()) of the sum of distances of the closest segment ",
	      "to its region's center, normalized by the number of regions.  ",
	      "(A region represents the most even spacing possible given the ",
	      "number of unique segments.)  Smaller spacing scores are ",
	      "better.\n");
      }
  }

sub scoreToString
  {
    my $score = $_[0];
    if(!defined($score) || !exists($score->{STRETCHES}))
      {return('')}
    return(join(',',
		(map {"$_\{" . join(',',@{$score->{STRETCHES}->{$_}}) . '}'}
		 sort {$b <=> $a}
		 keys(%{$score->{STRETCHES}})),$score->{BASESCORE}));
  }

sub outputSegmentTable
  {
    my $stats = $_[0];
    my $sizes = $_[1];

    my $seqids = [sort {$a cmp $b} keys(%$stats)];

    print("\n");

    foreach my $size (sort {$b <=> $a} @$sizes)
      {
	#Print the headers
	print("Stretch Inclusions for Size $size:\n",
	      "\t",join("\t",@$seqids),"\n");

	foreach my $rowid (@$seqids)
	  {
	    print($rowid);
	    foreach my $colid (@$seqids)
	      {
		if($rowid eq $colid)
		  {print("\t")}
		else
		  {print("\t",(exists($stats->{$rowid}->{$colid}->{$size}) ?
			       $stats->{$rowid}->{$colid}->{$size} : '0'))}
	      }
	    print("\n");
	  }
	print("\n");
      }
  }

sub outputValidation
  {print($_[0])}

sub numDiff
  {
    my $seq1 = $_[0];
    my $seq2 = $_[1];
    my $prot = $_[2];

    if(length($seq1) != length($seq2))
      {
	error("Sequence lengths differ.");
	return(-1);
      }

    my $d = 0;
    my $r = 0;
    foreach my $p (0..(length($seq1) - 1))
      {
	my $s1c = uc(substr($seq1,$p,1));
	my $s2c = uc(substr($seq2,$p,1));
	if($s1c eq '-' || $s2c eq '-')
	  {error("Alignment characters found.")}
	if($s1c ne $s2c)
	  {
	    $d++;
	    #If either is an ambiguous nucleotide
	    if(($prot && $s1c !~ /^x$/i && $s2c !~ /^x$/i) ||
	       (!$prot &&
		$s1c !~ /^[bdhvrykmswn]$/i && $s2c !~ /^[bdhvrykmswn]$/i))
	      {$r++}
	  }
      }

    return($d,$r);
  }

sub translate
  {
    my $nt         = $_[0];
    my $codon_hash = $_[1];
    my $len        = length($nt);
    my $aa         = '';

    if($len % 3)
      {
	error("Sequence length [$len] is not a multiple of 3.");
	$len -= ($len % 3);
      }

    for(my $i = 0;$i < $len;$i += 3)
      {
	my $codon = substr($nt,$i,3);
	if(!exists($codon_hash->{$codon}))
	  {
	    error("Codon [$codon] does not exist in the usage file.  ",
		  "Inserting an X.");
	    $aa .= 'X';
	  }
	elsif(!defined($codon_hash->{$codon}) || $codon_hash->{$codon} eq '')
	  {
	    error("Codon [$codon] exists in the usage file, but doesn't have ",
		  "an amino acid value.  Inserting an X.");
	    $aa .= 'X';
	  }
	elsif(length($codon_hash->{$codon}) != 1)
	  {
	    error("Bad amino acid value for codon [$codon] in the usage ",
		  "file.  Inserting an X.");
	    $aa .= 'X';
	  }
	else
	  {$aa .= $codon_hash->{$codon}}
      }

    return($aa);
  }

#This sub searches through a sequence's partners in the map to find one (or
#more) that has sequence for the supplied map section (i.e. the supplied
#divider).  If there does not exist a section (or sections) of sequence for the
#supplied divider, 0 is returned and it does not do any recoding.  If there is
#complete existing overlap, it obtains the sequence to be recoded and its
#partner sequence (for its gap characters).  It replaces the partner sequence
#with the sequence(s) saved in the map and then recodes each codon to best
#match (by homology) the new partner sequence (without changing the AA that was
#originally encoded).
sub reCodeLeftoverSegment
  {
    my $seqid      = $_[0];  #This is the sequence that needs to be recoded
    my $divider    = $_[1];  #This is the portion that needs recoded
    my $map        = $_[2];  #This is where the final sequences are saved
    my $partners   = $_[3];  #This is where to obtain the sequence to recode
    my $matrix     = $_[4];  #This has the encoding data
    my $data       = $_[5];
    my $codon_hash = $_[6];

    my $new_seq  = '';

    #Foreach partner sequence (that $seqid has been paired/aligned with)
    foreach my $partner (@$partners)
      {
	my $pairid        = $partner->[0];
	my $partner_seqid = $partner->[1];

	#We need to obtain the aligned sequence that we are going to recode and
	#its partner, which will have been replaced by the sequence that was
	#already saved in the map, either sourced from another alignment or
	#itself recoded
	my($fixed_partner_aln,$change_this_aln) =
	  getPartnerAlnSeqs($map,
			    $data,
			    $pairid,
			    $seqid,          #Obtain this seq from alignment
			    $divider,
			    $map->{$seqid}->{$divider}->{STOP},
			    $partner_seqid); #Replace this seq with seq in map

	if(defined($fixed_partner_aln) && $fixed_partner_aln ne '' &&
	   defined($change_this_aln)   && $change_this_aln ne '')
	  {
	    my $fixed_start = convertDivider($seqid,
					     $divider,
					     $pairid,
					     $partner_seqid,
					     $data);

	    $map->{$seqid}->{$divider}->{SEQ} =
	      reCodeOneSubseq($fixed_partner_aln,
			      $fixed_start,
			      $change_this_aln,
			      $divider,
			      $matrix,
			      $codon_hash);
	    return(1);
	  }
	elsif(defined($fixed_partner_aln) || defined($change_this_aln))
	  {debug("Will not change: $fixed_partner_aln\n",
		 "Should change:   $change_this_aln\n",
		 "But could not find a filled partner.")}
	else
	  {debug({LEVEL => 2},"No source sequence for [$partner_seqid] could ",
		 "be found to change sequence: [$seqid:",
		 "$divider-$map->{$seqid}->{$divider}->{STOP}] in [$pairid]");
	 }
      }

    #Could not find a completed partner sequence
    return(0);
  }

#This obtains a segment of an alignment (bounded by the coordinates of the
#sequence that will be recoded after it is returned).  The sequence that will
#be changed is obtained directly from the alignment.  The partner sequence is
#also obtained (for its alignment characters) but is replaced by the sequence
#saved in the map - if it's there.  If any sequence in the partner is missing
#from the map, undef is returned.
sub getPartnerAlnSeqs
  {
    my $map          = $_[0];
    my $data         = $_[1];
    my $pairid       = $_[2];
    my $change_seqid = $_[3]; #We will be recoding this sequence obtained
    my $change_start = $_[4]; #directly from the alignment
    my $change_stop  = $_[5];
    my $fixed_seqid  = $_[6]; #This sequence will be obtained from the
                              #alignment, but then replaced with the segment
                              #saved in the map

    #We're going to be looking through the partner sequences, because that's
    #the sequence we're going to be basing the recoding on.  We need to compile
    #all sections of the partner sequence in the map we need and check to see
    #if it's populated.  So we need the coordinates as a guide, thus we must
    #convert the coordinates we're looking for
    my $fixed_start = convertDivider($change_seqid,
				     $change_start,
				     $pairid,
				     $fixed_seqid,
				     $data);

    #This will get the coordinate of the last nt that aligned with the query
    #segment, even if it has an end gap (as opposed to how convertDivider
    #works)
    my $fixed_stop = convertAlnSeqCoord($change_seqid,
					$change_stop,
					$pairid,
					$fixed_seqid,
					$data);

    my $replace_seq = '';

    #These values will change
    my $cur_change_start = $change_start;
    my $cur_change_stop  = $change_stop;
    my($cur_fixed_stop);

    #For each section in order of start coordinate
    foreach my $cur_fixed_start (sort {$a <=> $b}
				 keys(%{$map->{$fixed_seqid}}))
      {
	$cur_fixed_stop = $map->{$fixed_seqid}->{$cur_fixed_start}->{STOP};

	#If there is overlap with the sequence slated to be recoded (we don't
	#care about the partner sequence except to know whether there is any
	#sequence saved in the map for it)
	if(coordsOverlap($cur_fixed_start,$cur_fixed_stop,
			 $fixed_start,$fixed_stop))
	  {
	    #If there is no sequence present yet for this section, we cannot
	    #proceed
	    if($map->{$fixed_seqid}->{$cur_fixed_start}->{SEQ} eq '')
	      {return(undef,undef)}
	    elsif(length($map->{$fixed_seqid}->{$cur_fixed_start}->{SEQ}) !=
		  ($map->{$fixed_seqid}->{$cur_fixed_start}->{STOP} -
		   $cur_fixed_start + 1))
	      {error("The sequence for [$fixed_seqid] stored in the map: ",
		     "[$map->{$fixed_seqid}->{$cur_fixed_start}->{SEQ}] of ",
		     "type [",
		     (defined($map->{$fixed_seqid}->{$cur_fixed_start}->{TYPE})
		      ? $map->{$fixed_seqid}->{$cur_fixed_start}->{TYPE} :
		      'undef'),"] and pair ",
		     "[$map->{$fixed_seqid}->{$cur_fixed_start}->{PAIR_ID}] ",
		     "is not the expected length: [",
		     ($map->{$fixed_seqid}->{$cur_fixed_start}->{STOP} -
		      $cur_fixed_start + 1),"].")}

	    my $substr_start = ($cur_fixed_start < $fixed_start ?
			        $fixed_start - $cur_fixed_start : 0);

	    #If this represents only a portion of the saved/recoded
	    #sequence, grab that subsequence from the map object
	    my $tmp_replace_seq = '';
	    if($cur_fixed_stop > $fixed_stop)
	      {
		my $remainder =
		  length($map->{$fixed_seqid}->{$cur_fixed_start}->{SEQ}) -
		    $substr_start;
		$remainder -= ($cur_fixed_stop - $fixed_stop);
		if(length($map->{$fixed_seqid}->{$cur_fixed_start}->{SEQ})
		   < ($substr_start + 1))
		  {error("Sequence: [",
			 $map->{$fixed_seqid}->{$cur_fixed_start}->{SEQ},
			 "] from [$fixed_seqid] is too short to be ",
			 "substringed starting at position [",
			 ($substr_start + 1),"].  cur_fixed_start:",
			 "$cur_fixed_start cur_fixed_stop:$cur_fixed_stop ",
			 "orig:$fixed_start-$fixed_stop")}
		elsif(length($map->{$fixed_seqid}->{$cur_fixed_start}
			     ->{SEQ}) < ($substr_start + $remainder))
		  {error("Sequence: [",
			 $map->{$fixed_seqid}->{$cur_fixed_start}->{SEQ},
			 "] is too short to obtain [$remainder] bases ",
			 "substringed starting from position [",
			 ($substr_start + 1),"].")}
		$tmp_replace_seq =
		  substr($map->{$fixed_seqid}->{$cur_fixed_start}->{SEQ},
			 $substr_start,$remainder);
	      }
	    else
	      {$tmp_replace_seq =
		 substr($map->{$fixed_seqid}->{$cur_fixed_start}->{SEQ},
			$substr_start)}

	    my $tmp_fixed_start = ($cur_fixed_start < $fixed_start ?
				   $fixed_start : $cur_fixed_start);
	    my $tmp_fixed_stop  = ($cur_fixed_stop > $fixed_stop ?
				   $fixed_stop  : $cur_fixed_stop);
	    #Sanity check 2 (after changes)
	    if(length($tmp_replace_seq) !=
	       ($tmp_fixed_stop - $tmp_fixed_start + 1))
	      {error("getAlnSegment failed sanity check 2 for the fixed aln ",
		     "section [$tmp_fixed_start-$tmp_fixed_stop]: ",
		     "[$tmp_replace_seq]")}

	    $replace_seq .= $tmp_replace_seq;

	    #If we've gotten all the sequence we were looking for, move on
	    last if($cur_fixed_stop >= $fixed_stop);
	  }
      }

    if(!defined($cur_fixed_stop) || $cur_fixed_stop < $fixed_stop)
      {
	error("Did not make it to the target coordinate: ",
	      "[$fixed_seqid:$fixed_stop] (converted from coordinate: ",
	      "[$change_seqid:$change_stop]).  The last coordinate processed ",
	      "was [",(defined($cur_fixed_stop) ? $cur_fixed_stop : 'undef'),
	      "].  Map incomplete.");
	return(undef,undef);
      }
    #Sanity check
    elsif(length($replace_seq) != ($fixed_stop - $fixed_start + 1))
      {error("Fixed aln [$replace_seq] failed sanity check [$fixed_start-",
	     "$fixed_stop]")}

    #Obtain the segment of the alignment containing the sequence of the segment
    #we're going to be changing upon return and basically just the sequence
    #positions/gaps from the partner sequence that we will be replacing from
    #the map
    my($fixed_aln,$change_aln) =
      getAlnSegment($pairid,
		    $fixed_seqid,
		    $change_seqid,
		    $change_start,
		    $change_stop,
		    $data);

    my $check_aln = $fixed_aln;
    $check_aln =~ s/-+//g;
    if(length($check_aln) < length($replace_seq))
      {error("The number of non-gap characters in the alignment sequence ",
	     "[$fixed_aln]: [",length($check_aln),"] was less than the ",
	     "unaligned replacement sequence: [$replace_seq]: [",
	     length($replace_seq),"].  Alignment: [$pairid]  Seq:Start-stop: ",
	     "[$fixed_seqid:$fixed_start-$fixed_stop]  Reference/change ",
	     "sequence: [$change_seqid:$change_start-$change_stop]")}

    #Replace the sequence in the fixed aligned string with what we got from the
    #map
    $fixed_aln = replaceAlnSegment($fixed_aln,
				   $replace_seq);

    #Now we can return the sequences
    return($fixed_aln,$change_aln);
  }

#Takes a target sequence's start & stop plus the partner sequence's start &
#stop (in case they align with gaps)
#Searches the alignment to find the first start (target or query) and records
#the alignment sequence from there to the last stop.
#Returns 2 things: the segment of the aligned target sequence and the segment
#of the aligned query sequence
#NOTE: This does a linear search, which is sub-optimal, but expect this to only
#be used with short sequences
sub getAlnSegment
  {
    my $pairid      = $_[0];
    my $getseqid    = $_[1];
    my $query_seqid = $_[2];
    my $query_start = $_[3];
    my $query_stop  = $_[4];
    my $data        = $_[5];

    my $target_aln = $data->{$pairid}->{ALNS}->{$getseqid};
    my $query_aln  = $data->{$pairid}->{ALNS}->{$query_seqid};
    my $query_len  = length($data->{$pairid}->{SEQS}->{$query_seqid});

    my $aln_pos        = 0;
    my $cur_target_pos = 0;
    my $cur_query_pos  = 0;
    my $max_pos = length($target_aln);

    my $target_aln_segment = '';
    my $query_aln_segment  = '';

    #While we are not yet at the end of the captured segment
    while(($cur_query_pos < $query_stop ||
	   ($cur_query_pos == $query_stop && $query_stop == $query_len)) &&
	  $aln_pos < $max_pos)
      {
	my $target_char = substr($target_aln,$aln_pos,1);
	my $query_char  = substr($query_aln, $aln_pos,1);

	$aln_pos++;

	if($target_char ne '-')
	  {$cur_target_pos++}
	if($query_char ne '-')
	  {$cur_query_pos++}

	if($cur_query_pos >= $query_start ||
	   ($cur_query_pos == 0 && $query_start == 1))
	  {
	    $target_aln_segment .= $target_char;
	    $query_aln_segment  .= $query_char;
	  }
      }

    return($target_aln_segment,$query_aln_segment);
  }

#Replace the sequence portion of an alignment with the characters from another
#unaligned sequence.  Note, we are basically inserting the characters from one
#alignment (optimized for one another) into another alignment.  Any gaps in the
#new alignment will cause the original sequence aligned with that gaps to not
#change, thus any extra sequence sent in is tacked on unchanged.
sub replaceAlnSegment
  {
    my $aln_seq     = $_[0];
    my $replace_seq = $_[1];

    my $seq_pos = 0;
    my $new_aln = '';

    for(my $aln_pos = 0;$aln_pos < length($aln_seq);$aln_pos++)
      {
	my $char = substr($aln_seq,$aln_pos,1);
	if($char ne '-')
	  {
	    #The replace sequence passed in is the sequence that best matched
	    #the sequence this sequence is being aligned with.  If that
	    #sequence that this sequence is aligned with had a gap at the end
	    #of this segment, then that portion of sequence can't change, so we
	    #just need to tac it on.
	    if(($seq_pos + 1) > length($replace_seq))
	      {$new_aln .= $char}
	    else
	      {$new_aln .= substr($replace_seq,$seq_pos,1)}
	    $seq_pos++;
	  }
	else
	  {$new_aln .= '-'}
      }

    if($seq_pos < length($replace_seq))
      {error("The number of non-gap characters in the alignment sequence ",
	     "[$aln_seq]: [$seq_pos] was less than the unaligned replacement ",
	     "sequence: [$replace_seq]: [",length($replace_seq),"].")}
    elsif($seq_pos > length($replace_seq))
      {debug("Replaced a portion of a larger sequence.  Assuming there was a ",
	     "gap.")}

    return($new_aln);
  }

sub coordsOverlap
  {
    my $start1 = $_[0];
    my $stop1  = $_[1];
    my $start2 = $_[2];
    my $stop2  = $_[3];

    if(($start1 >= $start2 && $start1 <= $stop2) ||
       ($stop1 >= $start2 && $stop1 <= $stop2) ||
       ($start1 < $start2 && $stop1 > $stop2))
      {return(1)}

    return(0);
  }

sub outputHybrids
  {
    my $seq_hash    = $_[0]->[0];
    my $source_hash = $_[0]->[1];
    my $soln_stats  = $_[0]->[2];
    my $valid_msg   = $_[0]->[3];
    my $outfile     = $_[1];
    my $mapfile     = $_[2];

    openOut(*OUT,$outfile) || return();

    foreach my $seqid (sort {$a cmp $b} keys(%$seq_hash))
      {print(">$seqid\n",$seq_hash->{$seqid},"\n")}

    closeOut(*OUT);

    if(exists($source_hash->{LEGEND}))
      {
	#If it doesn't open, we'll be printing to STDOUT, so print a newline to
	#visually separate it from the output above
	unless(openOut(*MAP,$mapfile))
	  {print("\n")}

	print("#Source sequence legend.  Each sequence is represented by a ",
	      "character code.  The sequences below are composed of ",
	      "character codes that indicate the source of the sequence, ",
	      "which can be directly extracted from a pairwise aligning ",
	      "(upper-case) or recoded to match a sequence that was ",
	      "extracted from an unrelated pairwise alignment (lower-case)\n");
	foreach my $seqid (sort {$a cmp $b} keys(%{$source_hash->{LEGEND}}))
	  {print("#$seqid\t",uc($source_hash->{LEGEND}->{$seqid}),"\n")}
	foreach my $seqid (sort {$a cmp $b} keys(%{$source_hash->{SEQS}}))
	  {print(">$seqid\n",$source_hash->{SEQS}->{$seqid},"\n")}

	closeOut(*MAP);
      }

    return($soln_stats,$valid_msg);
  }

#This creates a hash containing references to the hash with the codon keys
#pointing to the scores and data
sub makeReverseCodonHash
  {
    my $codon_hash = $_[0];
    my $rev_cdn_hash = {};

    foreach my $aa (keys(%$codon_hash))
      {foreach my $cdn (keys(%{$codon_hash->{$aa}}))
	 {$rev_cdn_hash->{$cdn} = $aa}}

    return($rev_cdn_hash);
  }

#Reads a given aligned fasta file and creates a data structure like the
#following:
#$hash->{$pairid}->{SEGS}     = HASH->{$sz}->{$sqid} = [[$start,$spc_scr],...]
#$hash->{$pairid}->{NUMUNIQSEGS} = HASH->{$sz} = # non-overlapping $sz ID segs
#$hash->{$pairid}->{NUMBASES} = HASH->{$sz} = # 100%ID bases in the segs
#$hash->{$pairid}->{SPACINGSCORE} = HASH->{$sz} = # where small score = better
#$hash->{$pairid}->{SEQS}     = HASH->{$sqid} = $seq (with no aln characters)
#$hash->{$pairid}->{ALNS}     = HASH->{$sqid} = $seq (with aln characters)
#Returns 0 upon success, a positive error code otherwise
sub loadAlignment
  {
    my $alnfile       = $_[0];
    my $stretch_sizes = $_[1];
    my $hash          = $_[2];
    my $recnum        = 0;
    my $codon_err     = 0;
    my $len           = -1;

    if(ref($stretch_sizes) ne 'ARRAY' || scalar(@$stretch_sizes) == 0)
      {
	error("Invalid stretch sizes array sent as second argument.");
	return(1);
      }

    my $pairid = $alnfile;
    $pairid =~ s%.*/%%;
    $pairid =~ s%\.[^\.]+$/%%;
    #Just in case 2 files have the same name but live in different directories
    if(exists($hash->{$pairid}))
      {$pairid = $alnfile}

    openIn(*ALN,$alnfile) || return(2);

    while(my $rec = getNextSeqRec(*ALN,0,$alnfile,'fasta'))
      {
	$recnum++;
	if($recnum > 2)
	  {
	    error("More than 2 sequence records were found in alignment file ",
		  "[$alnfile].  Only the first 2 will be used.  All other ",
		  "sequences will be ignored.",
		  {DETAIL => ("This script processes pairwise alignments " .
			      "and the file name is used as a pair ID.")});
	    last;
	  }

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

	#Set the length of the alignment
	if($len == -1)
	  {$len = length($rec->[1])}

	#Error-check the alignment length
	if(length($rec->[1]) == 0)
	  {
	    error("Length 0 sequence [$id] encountered in alignment file ",
		  "[$alnfile].  Unable to process.");
	    closeIn(*ALN);
	    return(3);
	  }
	elsif(length($rec->[1]) != $len)
	  {
	    error("Alignment strings not the same length [$len != ",
		  length($rec->[1]),"] in alignment file [$alnfile].  Unable ",
		  "to process.");
	    closeIn(*ALN);
	    return(4);
	  }

	#Error-check that the sequence is codon-aligned
	if(!codonAligned($rec->[1]))
	  {
	    error("Sequence [$id] from file [$alnfile] is not codon-aligned.",
		  {DETAIL => join('','This script operates on the nucleotide ',
				  '-o output of codonHomologizer.pl when -a ',
				  'is "on"/true (its default).  Any set of ',
				  'pair-wise alignments can be submitted, as ',
				  'long as they are codon-aligned, due to ',
				  'the assumption it makes that the codons ',
				  'are aligned.  Having aligned codons is ',
				  'ideal for crossover events resulting in ',
				  'frame-correct results.')});
	    $codon_err++;
	  }

	$hash->{$pairid}->{ALNS}->{$id} = $rec->[1];
	$hash->{$pairid}->{SEQS}->{$id} = $rec->[1];
	$hash->{$pairid}->{SEQS}->{$id} =~ s/-+//g;

	debug("DATA HASH: $pairid $id $hash->{$pairid}->{SEQS}->{$id}");
      }

    closeIn(*ALN);

    if($codon_err)
      {return(5)}

    if($recnum < 2)
      {
	error("Only [$recnum] aligned sequences were found in alignment file ",
	      "[$alnfile].  Two sequences are required in each alignment ",
	      "file.");
	return(6);
      }

    return(loadSegments($hash->{$pairid},
			$stretch_sizes,
			#Sequence IDs (sorted) - can only be 2 keys
			(sort(keys(%{$hash->{$pairid}->{ALNS}}))),
			#Sequences (sorted by IDs) - can only be 2 keys
			(map {$hash->{$pairid}->{ALNS}->{$_}}
			 sort(keys(%{$hash->{$pairid}->{ALNS}})))));
  }

#This checks to see if the supplied alignment sequence is codon-aligned
sub codonAligned
  {
    my $seq = $_[0];
    my $l   = length($seq);
    my $n   = 0;

    for(my $c = 0;$c < $l;$c++)
      {
	my $nt = substr($seq,$c,1);
	if($nt eq '-')
	  {if($n % 3)
	     {return(0)}}
	else
	  {$n++}
      }

    return(1);
  }

#Takes 2 alignment strings and finds identical segments to load into a hash for
#searching.  Sets the following values in the hash:
#$hash->{$pairid}->{SEGS}     = HASH->{$sz}->{$sqid} =[[$start,$spacing_score],
#                                                      ...]
#$hash->{$pairid}->{NUMUNIQSEGS}  = HASH->{$sz} = # non-overlapping $sz ID segs
#$hash->{$pairid}->{NUMBASES}     = HASH->{$sz} = # 100%ID bases in the segs
#$hash->{$pairid}->{SPACINGSCORE} = HASH->{$sz} = # where small score = better
#Returns 0 upon success, a positive error code otherwise
#NOTE: This sub assumes that only 2 sequences were aligned.  If more were
#aligned, it's possible that a segment of identity will be missed because a
#third sequence could cause gaps to be in the middle of otherwise identical
#segments.
#TODO: The above issue/assumption could be worked-around by preprocessing the
#      alignment strings to remove aligned gaps.
sub loadSegments
  {
    my $hash          = $_[0];
    my $stretch_sizes = $_[1];
    my $id1           = $_[2]; #SeqID1
    my $id2           = $_[3]; #SeqID2
    my $aln1          = $_[4];
    my $aln2          = $_[5];

    my $aln_len = length($aln1);
    if(length($aln2) != $aln_len)
      {error("Aligned sequences not the same length.",
	     {DETAIL => "$id1: $aln1\n$id2: $aln2"})}

    #For each segment size we're trying to optimize
    foreach my $sz (@$stretch_sizes)
      {
	my $spos1           = 0;
	my $spos2           = 0;
	my $last_apos       = -1; #Last aln pos saved as a unique segment
	my $num_unique_segs = 0;
	my $num_seg_bases   = 0;

	#For each position in the alignment
	foreach my $apos (0..($aln_len - $sz))
	  {
	    #Grab the subsequences of size $sz and track each sequence's coord
	    my $subaln1 = substr($aln1,$apos,$sz);
	    if($subaln1 !~ /^-/)
	      {$spos1++}
	    my $subaln2 = substr($aln2,$apos,$sz);
	    if($subaln2 !~ /^-/)
	      {$spos2++}

	    #If the subalns are the same and don't contain gaps
	    if($subaln1 eq $subaln2 && $subaln1 !~ /-/)
	      {
		debug({LEVEL => 4},
		      "Identity segment: $id1:$spos1/$id2:$spos2 $subaln1");

		#Save the segment and each sequence's respective position
		push(@{$hash->{SEGS}->{$sz}->{$id1}},[$spos1,-1]);
		push(@{$hash->{SEGS}->{$sz}->{$id2}},[$spos2,-1]);

		#If this is a non-overlapping segment
		if($last_apos == -1 || $apos >= ($last_apos + $sz))
		  {
		    $num_unique_segs++;
		    $num_seg_bases += $sz;
		    $last_apos = $apos;
		  }
		#Else if the segment is larger than the current segment size
		elsif($last_apos >= 0 && $apos < ($last_apos + $sz))
		  {$num_seg_bases++}
	      }
	  }

	#Determine a score for how spread out the segments are where a smaller
	#score is more evenly spaced/better
	my $period_len1  = int(length($aln1) / ($num_unique_segs + 1));
	$period_len1 = 1 if($period_len1 < 1);
	my $period_len2  = int(length($aln2) / ($num_unique_segs + 1));
	$period_len2 = 1 if($period_len2 < 1);
	my $dist_scores1 = {};
	my $dist_scores2 = {};
	my $regions      = {};
	foreach my $pos_index (0..$#{$hash->{SEGS}->{$sz}->{$id1}})
	  {
	    my $pos1    = $hash->{SEGS}->{$sz}->{$id1}->[$pos_index]->[0];
	    my $middle1 = int($pos1 + $sz / 2);
	    my $region1 = int($middle1 / $period_len1);
	    my $target1 = $period_len1 * ($region1 + 1);
	    #Find the closer evenly spaced coordinate/region/target
	    if(abs($middle1 - ($period_len1 * ($region1 + 2))) <
	       abs($middle1 - ($period_len1 * ($region1 + 1))))
	      {
		$region1++;
		$target1 = $period_len1 * ($region1 + 1);
	      }
	    my $offset1 = abs($target1 - $middle1);

	    #Record the spacing score for this segment in order to be able to
	    #compare the spacings of segments included in the final solution
	    #later
	    $hash->{SEGS}->{$sz}->{$id1}->[$pos_index]->[1] = $offset1;

	    #If there's no score for this region or this score's better/smaller
	    if(!exists($dist_scores1->{$region1}) ||
	       $offset1 < $dist_scores1->{$region1})
	      {$dist_scores1->{$region1} = $offset1}

	    my $pos2    = $hash->{SEGS}->{$sz}->{$id2}->[$pos_index]->[0];
	    my $middle2 = int($pos2 + $sz / 2);
	    my $region2 = int($middle2 / $period_len2);
	    my $target2 = $period_len2 * ($region2 + 1);
	    #Find the closer evenly spaced coordinate/region/target
	    if(abs($middle2 - ($period_len2 * ($region2 + 2))) <
	       abs($middle2 - ($period_len2 * ($region2 + 1))))
	      {
		$region2++;
		$target2 = $period_len2 * ($region2 + 1);
	      }
	    my $offset2 = abs($target2 - $middle2);

	    #Record the spacing score for this segment in order to be able to
	    #compare the spacings of segments included in the final solution
	    #later
	    $hash->{SEGS}->{$sz}->{$id2}->[$pos_index]->[1] = $offset2;

	    #If there's no score for this region or this score's better/smaller
	    if(!exists($dist_scores2->{$region2}) ||
	       $offset2 < $dist_scores2->{$region2})
	      {$dist_scores2->{$region2} = $offset2}

	    #Record the region the segment was found in
	    push(@{$regions->{$sz}->{$id1}->{$region1}},[$pos_index,$offset1]);
	    push(@{$regions->{$sz}->{$id2}->{$region2}},[$pos_index,$offset2]);
	  }

	my $spacing_score = 0;
	my $score1 = 0;
	my $score2 = 0;
	foreach my $reg (0..($num_unique_segs - 1))
	  {
	    my $score = exists($dist_scores1->{$reg}) ?
	      $dist_scores1->{$reg} : $period_len1;
	    $score1 += $score;
	    $spacing_score += $score;
	  }
	foreach my $reg (0..($num_unique_segs - 1))
	  {
	    my $score = exists($dist_scores2->{$reg}) ?
	      $dist_scores2->{$reg} : $period_len2;
	    $score2 += $score;
	    $spacing_score += $score;
	  }

	verbose({LEVEL => 2},"Populating regions.");

	#Record the regions as indexes of identity segments broken up by the
	#non-overlapping region they are closest to (given an ideal spacing
	#scheme).  These indexes will be used later to try and select segments
	#in a greedy fashion such that the first ones we attempt to add will
	#fit something from every alignment
	if(scalar(keys(%{$regions->{$sz}->{$id1}})) >
	   scalar(keys(%{$regions->{$sz}->{$id2}})) ||
	   (scalar(keys(%{$regions->{$sz}->{$id1}})) ==
	    scalar(keys(%{$regions->{$sz}->{$id2}})) && $score1 < $score2))
	  {
	    my $regidx = 0;
	    foreach my $regkey (sort {$a <=> $b}
				keys(%{$regions->{$sz}->{$id1}}))
	      {
		debug({LEVEL => 2},"Size [$sz], region [$regidx] is getting [",
		      scalar(@{$regions->{$sz}->{$id1}->{$regkey}}),
		      "] segment indexes sorted by score.");
		$hash->{REGS}->{$sz}->[$regidx] =
		  [map {$_->[0]} sort {$a->[1] <=> $b->[1]}
		   @{$regions->{$sz}->{$id1}->{$regkey}}];
		$regidx++;
	      }
	  }
	else
	  {
	    my $regidx = 0;
	    foreach my $regkey (sort {$a <=> $b}
				keys(%{$regions->{$sz}->{$id2}}))
	      {
		debug({LEVEL => 2},"Size [$sz], region [$regidx] is getting [",
		      scalar(@{$regions->{$sz}->{$id2}->{$regkey}}),
		      "] segment indexes sorted by score.");
		$hash->{REGS}->{$sz}->[$regidx] =
		  [map {$_->[0]} sort {$a->[1] <=> $b->[1]}
		   @{$regions->{$sz}->{$id2}->{$regkey}}];
		$regidx++;
	      }
	  }

	#If there are any regions with no segments, create empty arrays for
	#them
	foreach my $regidx (0..($num_unique_segs - 1))
	  {if(!defined($hash->{REGS}->{$sz}->[$regidx]))
	     {$hash->{REGS}->{$sz}->[$regidx] = []}}

	$hash->{NUMUNIQSEGS}->{$sz}  = $num_unique_segs;
	$hash->{NUMBASES}->{$sz}     = $num_seg_bases;
	$hash->{SPACINGSCORE}->{$sz} = $spacing_score;
      }

    return(0);
  }

#DEPRECATED: These warnings have been changed to debug messages.  There still
#exist conflict and multiple divider issues that are detected, but code was
#written to fix/deal-with conflicts, and the multiple divider issue is
#*expected* in gappy all versus all alignment pairs (but are innocuous as well,
#since the bad effects they used to have have been resolved as well).  There is
#a separate error that is issued if the attempt to fix a conflict fails.
#Globals used: $divider_warnings
sub issueDividerWarnings
  {
    my $div_map = $_[0];

    my $explanation =
      join('',('Here is what this warning means: Between each closest pair ',
	       'of included identity stretches, a midpoint is selected that ',
	       'determines which alignment the sequence will be extracted ',
	       'from or what regions need to be recoded.  Sometimes this ',
	       'region can be aligned with a gap in different pairwise ',
	       'alignments, which can confuse the determination of a ',
	       'midpoint.  This warning may or may not result in issues ',
	       'constructing the final sequences.  As long as the sequences ',
	       'are marked as validated, the issue is a minor one and will ',
	       'only affect homology in the regions between identity ',
	       'segments.  If the sequences do not validate, these warnings ',
	       'can be helpful in manually determining how to fix the ',
	       'encoding and mixing of identity segments.  Any divider ',
	       'labeled as "RECODE" indicates a position (from that to the ',
	       'next divider and/or identity segment that may not be ',
	       'optimally encoded.'));

    foreach my $divwarnkey (sort {$divider_warnings->{$a}->{ORDER} <=>
				    $divider_warnings->{$b}->{ORDER}}
			    keys(%$divider_warnings))
      {
	my $seqid    = $divider_warnings->{$divwarnkey}->{SEQID};
	my $dividers = $divider_warnings->{$divwarnkey}->{DIVIDERS};
	my $msg      = $divider_warnings->{$divwarnkey}->{MSG};
	my $detail   = (exists($divider_warnings->{$divwarnkey}->{DETAIL}) ?
			$divider_warnings->{$divwarnkey}->{DETAIL} : '');

	#First, determine whether the undefined dividers were updated to source
	#dividers.  No warning is necessary if all dividers are source dividers
	my $dowarn = 0;
	foreach my $divider (@$dividers)
	  {
	    #All these checks are expected to be true except the checks on TYPE
	    #which hopefully are all false
	    if(exists($div_map->{$seqid}) &&
	       exists($div_map->{$seqid}->{$divider}) &&
	       exists($div_map->{$seqid}->{$divider}->{TYPE}) &&
	       (!defined($div_map->{$seqid}->{$divider}->{TYPE}) ||
		$div_map->{$seqid}->{$divider}->{TYPE} ne 'SOURCE'))
	      {$dowarn = 1}
	  }

	if($dowarn)
	  {
	    #Put the 2 warning types on separate lines so that the run report
	    #shows them on separate lines.
	    if($msg =~ /^Conflicting/)
	      {debug($msg,{DETAIL => $detail . $explanation})}
	    else
	      {debug($msg,{DETAIL => $detail . $explanation})}
	  }
      }
  }




















