#!/usr/bin/perl

#USAGE: Run with no options to get usage or with --help for basic details

use CommandLineInterface;
use warnings;
use strict;
require('ch_lib.pl'); #TODO: I'll turn this into a module later

##
## Describe the script
##

our $VERSION = '1.000';

setScriptInfo(CREATED => '10/5/2017',
              VERSION => $VERSION,
              AUTHOR  => 'Robert William Leach',
              CONTACT => 'rleach@princeton.edu',
              COMPANY => 'Princeton University',
              LICENSE => 'Copyright 2017',
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

my $codon_file_type =
  addInfileOption(GETOPTKEY   => 'c|codon-file|codon-usage-file=s',
		  REQUIRED    => 0,
		  PRIMARY     => 0,
		  DEFAULT     => undef,
		  SMRY_DESC   => 'Codon usage file.',
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
		   DETAIL_DESC   => 'Name of the output DNA sequence file.',
		   FORMAT_DESC   => 'Nucleotide fasta file.',
		   COLLISIONMODE => 'merge',
		   PAIR_WITH     => $seq_file_type,
		   PAIR_RELAT    => 'ONETOMANY'
		  );

my $stretch_mins = [5,11];
addArrayOption(GETOPTKEY   => 'stretch-min=s',
	       GETOPTVAL   => $stretch_mins,
	       REQUIRED    => 0,
	       DEFAULT     => $stretch_mins,
	       HIDDEN      => 0,
	       INTERPOLATE => 1,
	       SMRY_DESC   => 'Length of contiguous identity to search & mix.',
	       DETAIL_DESC => << "END_DETAIL"

This is the number of contiguous identical nucleotides between the sequences in the supplied alignments (see -i) for searching and mixing into the output sequences.  Every identical stretch of this length will be a candidate for weaving into the final output sequence for each input sequence present among the alignments.  Multiple lengths can be supplied.  The algorithm will attempt to insert them in descending size order.  The default values of 11 and 5 are based on the characteristics of recombination.  See --help --extended for more information.

END_DETAIL
	 );

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

my $codon_data_hash  = {};
my $matrix_data_hash = {};
my $input_data_hash  = {}; #->{$ntOutFile} = {CDNHASH => $cuh, DATA => $hash}

#nextFileSet iterates sets of command-line-supplied files processed together
while(nextFileCombo())
  {
    my $codonUsageFile = getInfile($codon_file_type);
    my $alnFile        = getInfile($seq_file_type);
    my $ntOutFile      = getOutfile($seq_out_type);

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
      }

    loadAlignment($alnFile,
		  $stretch_mins,
		  $input_data_hash->{$ntOutFile}->{DATA});
  }

my $max_size = (sort {$b <=> $a} @$stretch_mins)[0];

foreach my $outfile (keys(%$input_data_hash))
  {
    my $solution = {};
    #Add segments to the solution starting with the largest segments and going
    #down in size (to squeeze less optimally sized segments in)
    foreach my $size (sort {$b <=> $a} @$stretch_mins)
      {
	my $sorted_pair_ids = sortPairIDs($input_data_hash->{$outfile}->{DATA},
					  $size);
	$solution = getMaxWeavedSegments($input_data_hash->{$outfile}->{DATA},
					 $size,
					 $solution,
					 $sorted_pair_ids,
					 $max_size,
					 $input_data_hash->{$outfile}
					 ->{CDNHASH},
					 $input_data_hash->{$outfile}
					 ->{MTXHASH});
      }

    outputHybrids(generateHybrids($solution,
				  $input_data_hash->{$outfile}->{DATA},
				  $input_data_hash->{$outfile}->{CDNHASH},
				  $input_data_hash->{$outfile}->{MTXHASH}),
		   $outfile);
  }










##
## Subroutines...
##

#This sub takes the input_data_hash (which I will refer to as the segments
#object) and tries to select the most non-overlapping segments for inclusion in
#the final sequence output.
#DEPRECATED & UNTESTED - After starting this sub, I realized that I needed an implementation that gives me a good answer in a short amount of time, so I appended NOTGREEDY to the name of this sub and rewrote it below.  I don't know if this sub works, but the one below was based on the same design.  It just added alternating regions too.
sub getMaxWeavedSegmentsNOTGREEDY
  {
    my $hash       = $_[0];
    my $size       = $_[1];
    my $soln       = $_[2];
    my $pairs      = $_[3];
    my $pair_start = defined($_[4]) ? $_[4] : 0; #DO NOT SUPPLY
    my $pos_start  = defined($_[5]) ? $_[5] : 0; #DO NOT SUPPLY
    my $best       = {};
    my $best_score;

    foreach my $pairindex ($pair_start..$#{$pairs})
      {
	my $pairid = $pairs->[$pairindex];

	#Arbitrarily (but consistently) select the first seqid from the pair
	my($seqid1,$seqid2) =
	  sort {$a cmp $b} keys(%{$hash->{$pairid}->{SEQS}});

	#For each identical segment position (they are already in order)
	foreach my $index ($pos_start..$#{$hash->{$pairid}->{SEGS}->{$size}
					  ->{$seqid1}})
	  {
	    my $pos1 =
	      $hash->{$pairid}->{SEGS}->{$size}->{$seqid1}->[$index]->[0];
	    my $pos2 =
	      $hash->{$pairid}->{SEGS}->{$size}->{$seqid2}->[$index]->[0];

	    #Determine the coordinates after expansion to include every touched
	    #codon
	    my $first_codon_start1 = $pos1 - ($pos1 % 3);
	    my $last_codon_end1    = ($pos1 + $size - 1) +
	      (($pos1 + $size - 1) % 3);
	    my $first_codon_start2 = $pos2 - ($pos2 % 3);
	    my $last_codon_end2    = ($pos2 + $size - 1) +
	      (($pos2 + $size - 1) % 3);

	    #If the solution will accept this segment
	    if(canAddSegment($soln,$seqid1,$seqid2,$pos1,$pos2,
			     $first_codon_start1,$last_codon_end1,
			     $first_codon_start2,$last_codon_end2))
	      {
		my $new_soln = copySolution($soln);

		#Add this new segment for the first sequence
		push(@{$new_soln->{$seqid1}},
		     [$pos1,
		      $pos1+$size-1,
		      $first_codon_start1,
		      $last_codon_end1,
		      $pairid]);

		#Add this new segment for the second sequence
		push(@{$new_soln->{$seqid2}},
		     [$pos2,
		      $pos2+$size-1,
		      $first_codon_start2,
		      $last_codon_end2,
		      $pairid]);

		#Determine where to start in the next recursion
		my $next_pair = $pairindex;
		my $next_pos  = $index + 1;
		if($next_pos >=
		   scalar(@{$hash->{$pairid}->{SEGS}->{$size}->{$seqid1}}))
		  {
		    $next_pair++;
		    $next_pos = 0;
		  }

		#Obtain a candidate solution
		my($cand);
		#If we're out of pairs, set the new solution as the one we were
		#going to send down the recursion (and don't do the recursion)
		if($next_pair >= scalar(@$pairs))
		  {$cand = $new_soln}
		else
		  {$cand = getMaxWeavedSegmentsNOTGREEDY($hash,$size,$new_soln,
							 $pairindex,$index)}

		#See if this new candidate solution is better
		my $cand_score = scoreSolution($cand);
		if(candidateBetter($cand_score,$best_score))
		  {
		    $best_score = $cand_score;
		    $best       = $cand;
		  }
	      }
	  }
      }

    return($best);
  }

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
    my $hash         = $_[0];
    my $size         = $_[1];
    my $soln         = $_[2];
    my $pairs        = $_[3];
    my $max_size     = $_[4];
    my $codon_hash   = $_[5];
    my $pair_least   = (defined($_[6]) ? $_[6] : 0);             #DO NOT SUPPLY
    my $pair_start   = (defined($_[7]) ? $_[7] : $pair_least);   #BELOW THIS
    my $region_first = (defined($_[8]) ? $_[8] : [(0) x scalar(@$pairs)]);
    my $region_start = (defined($_[9]) ? $_[9] : [@$region_first]);
    my $counts       = (defined($_[10]) ? $_[10] :
			[map {[(0) x $hash->{$pairs->[$_]}->{NUMUNIQSEGS}
			       ->{$size}]}
			 0..$#{$pairs}]);
    my $best         = {};
    my $best_score;

    #While the first pair index still containing segments to try and add to the
    #solution is less than the number of pairs
    while($pair_least < scalar(@$pairs))
      {
	foreach my $curpair ($pair_start..$#{$pairs})
	  {
	    my $pairid  = $pairs->[$curpair];

	    #Arbitrarily order the sequence IDs
	    my($seqid1,$seqid2) =
	      sort {$a cmp $b} keys(%{$hash->{$pairid}->{SEQS}});

	    my $num_regions =
	      $hash->{$pairs->[$curpair]}->{NUMUNIQSEGS}->{$size};

	    while($region_first->[$curpair] < $num_regions)
	      {
		foreach my $curregion ($region_start->[$curpair]..
				       $#{$region_start})
		  {
		    #If there are still more segments in this region
		    if($counts->[$curpair]->[$curregion] <
		       scalar(@{$hash->{REGS}->{$size}->[$curregion]}))
		      {
			#The counts keep track of indexes into the segment
			#indexes contained in each region in the REGS hash.
			#They are sorted by each segment's distance from the
			#region center.
			my $reg_seg_index = $counts->[$curpair]->[$curregion];
			#Each segmentUsing the region segment index, we can
			#obtain the actual segment index so we can get the
			#start position of the segment and its spacing score.
			my $index = $hash->{REGS}->{$size}->[$curregion]
			  ->[$reg_seg_index];

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
			my $first_codon_start1 = $pos1 - ($pos1 % 3);
			my $last_codon_end1    = ($pos1 + $size - 1) +
			  (($pos1 + $size - 1) % 3);
			my $first_codon_start2 = $pos2 - ($pos2 % 3);
			my $last_codon_end2    = ($pos2 + $size - 1) +
			  (($pos2 + $size - 1) % 3);

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
						     $codon_hash);

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
				 $region_first->[$next_pair]}

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
							     $pair_least,
							     [@$region_first],
							     $next_reg,
							     [map {[@$_]}
							      @$counts]);
			      }

			    #See if this new candidate solution is better
			    my $cand_score = scoreSolution($cand,
							   $hash,
							   $max_size);
			    if(candidateBetter($cand_score,$best_score))
			      {
				$best_score = $cand_score;
				$best       = $cand;
			      }
			  }

			$counts->[$curpair]->[$curregion]++;
		      }
		  }

		#While the segment number we are on in the first region is
		#greater than the segment array's last index (i.e. we are out
		#of segments to try and add from the first region), increment
		#the first region index
		while($counts->[$curpair]->[$region_first->[$curpair]] >
		      $#{$hash->{REGS}->{$size}->[$region_first->[$curpair]]})
		  {$region_first->[$curpair]++}

		#This exists to accommodate recursive calls that set
		#$start_region (where the search should continue from) but we
		#want to go back to the first region with segments in it once
		#we reach the end of this set of loop iterations
		$region_start->[$curpair] = $region_first->[$curpair];
	      }
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
      }

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
    my $alt_codon         = $_[8];

    my $final_start       = $_[9];
    my $final_stop        = $_[10];

    my $spacing_score     = $_[11];

    my $insert_index = findSegmentInsertPos($soln->{$seqid},$iden_start);
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
	    SPACING_SCORE     => $spacing_score});

    return(0);
  }

#This returns the index where the record should go.
sub findSegmentInsertPos
  {
    my $seg_array  = $_[0];
    my $iden_start = $_[1];
    my $index = 0;

    if(scalar(@$seg_array) == 0)
      {return($index)}

    my $mini = 0;
    my $maxi = $#{$seg_array};
    my $i    = int($maxi / 2);

    #There could be multiple records with the same start (from different
    #pairs), so let's close in on where it is
    while(($maxi - $mini) > 0)
      {
	#If $i's segment starts to the left of the stop coordinate, move the
	#left bound right
	if($iden_start < $seg_array->[$i]->{IDEN_START})
	  {$maxi = $i}
	elsif($iden_start > $seg_array->[$i]->{IDEN_START})
	  {$mini = $i}
	else
	  {return($i + 1)}

	$i = int(($maxi - $mini) / 2) + $mini;

	if($i == $mini || $i == $maxi)
	  {last}
      }

    if($iden_start > $seg_array->[$mini]->{IDEN_START} &&
       $iden_start < $seg_array->[$maxi]->{IDEN_START})
      {return($maxi)}
    elsif($iden_start < $seg_array->[$mini]->{IDEN_START})
      {return($mini)}
    elsif($iden_start > $seg_array->[$maxi]->{IDEN_START})
      {return($maxi + 1)}
    else
      {
	error("Binary search of a [",scalar(@$seg_array),"] member list ",
	      "failed for identity start coordinate [$iden_start].");
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
    my $max_size           = $_[12];
    my $hash               = $_[13];
    my $spacing_score1     = $_[14];
    my $spacing_score2     = $_[15];
    my $codon_hash         = $_[16];

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

	       undef,
	       undef,

	       $spacing_score1);

    addSegment($new_soln,
	       $pairid,

	       $seqid2,
	       $iden_start2,
	       $iden_stop2,
	       $first_codon_start2,
	       $last_codon_stop2,

	       'seg',
	       undef,

	       undef,
	       undef,

	       $spacing_score2);

    if(defined($fit_data1->{LEFT_CODON}))
      {
	addSegment($new_soln,
		   $pairid,

		   $seqid1,
		   $first_codon_start1,
		   $first_codon_start1 + 2,
		   $first_codon_start1,
		   $first_codon_start1 + 2,

		   'alt',
		   $fit_data1->{LEFT_CODON},

		   $first_codon_start1,
		   $first_codon_start1 + 2,

		   undef);
      }

    if(defined($fit_data1->{RIGHT_CODON}))
      {
	addSegment($new_soln,
		   $pairid,

		   $seqid1,
		   $last_codon_stop1 - 2,
		   $last_codon_stop1,
		   $last_codon_stop1 - 2,
		   $last_codon_stop1,

		   'alt',
		   $fit_data1->{RIGHT_CODON},

		   $last_codon_stop1 - 2,
		   $last_codon_stop1,

		   undef);
      }

    if(defined($fit_data2->{LEFT_CODON}))
      {
	addSegment($new_soln,
		   $pairid,

		   $seqid2,
		   $first_codon_start2,
		   $first_codon_start2 + 2,
		   $first_codon_start2,
		   $first_codon_start2 + 2,

		   'alt',
		   $fit_data2->{LEFT_CODON},

		   $first_codon_start2,
		   $first_codon_start2 + 2,

		   undef);
      }

    if(defined($fit_data2->{RIGHT_CODON}))
      {
	addSegment($new_soln,
		   $pairid,

		   $seqid2,
		   $last_codon_stop2 - 2,
		   $last_codon_stop2,
		   $last_codon_stop2 - 2,
		   $last_codon_stop2,

		   'alt',
		   $fit_data2->{RIGHT_CODON},

		   $last_codon_stop2 - 2,
		   $last_codon_stop2,

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
#will error out and return undef if they do not.
#Globals used: codon hash
sub reCodeMerge
  {
    my $codon_hash      = $_[0];
    my $new_codon       = uc($_[1]);
    my $new_fixed_poses = $_[2];
    my $cur_codon       = uc($_[3]);
    my $cur_fixed_poses = $_[4];

    #Do a quick error-check to make sure the codons each encode the same AA
    if(getAA($new_codon,$codon_hash) ne getAA($cur_codon,$codon_hash))
      {
	error("The codons [$new_codon,$cur_codon] from the same sequence ",
	      "optimized by different alignments do not encode the same ",
	      "amino acid.  Cannot recode.");
	return(undef);
      }

    my $new_codon_hash = codonToHashOfCodons($new_codon,$codon_hash);

    #Create an array of all possible alternate codons to choose from
    my $alternates = [grep {$_ ne $new_codon && $_ ne $cur_codon}
		      keys(%$new_codon_hash)];

    my $choices = [];
    foreach my $alt (@$alternates)
      {
	my $matches = 1;
	foreach my $new_pos (@$new_fixed_poses)
	  {
	    if(substr($new_codon,$new_pos - 1,1) ne
	       substr($alt,$new_pos - 1,1))
	      {
		$matches = 0;
		last;
	      }
	  }
	next if(!$matches);
	foreach my $cur_pos (@$cur_fixed_poses)
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
      {return(undef)}

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
	  [grep {length($_) == $flex_size}
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
      {error("Aligned sequence lengths are not the same.")}

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

    if(scalar(keys(%$soln)) == 0)
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
	#left bound right
	if($stop <= $soln->{$seqid}->[$i]->{FIRST_CODON_START})
	  {$maxi = $i}
	else
	  {$mini = $i}

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
	      $stop <=
	      ($soln->{$seqid}->[$i]->{FIRST_CODON_START} + $max_size - 1));

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

    foreach my $rec (@$seqlist)
      {
	if($rec->{PAIR_ID} ne $pairid)
	  {
	    if($rec->{TYPE} eq 'seg')
	      {
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
				      $left_codon_cur,$fixed_poses_cur);

			#If we couldn't select an alternative codon that
			#satisfies both sequences
			if(!defined($tmp_codon))
			  {return(undef)}
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
				      $right_codon_cur,$fixed_poses_cur);

			#If we couldn't select an alternative codon that
			#satisfies both sequences
			if(!defined($tmp_codon))
			  {return(undef)}
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

		   #the codons overlap (i.e. are in the same [frame] location)
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
				      $right_codon_cur,$fixed_poses_cur);

			#If we couldn't select an alternative codon that
			#satisfies both sequences
			if(!defined($tmp_codon))
			  {return(undef)}
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

		   #the codons overlap (i.e. are in the same [frame] location)
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
				      $left_codon_cur, $fixed_poses_cur);

			#If we couldn't select an alternative codon that
			#satisfies both sequences
			if(!defined($tmp_codon))
			  {return(undef)}
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
		if($rec->{IDEN_START} == $first_codon_start)
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
		      }
		  }

		#If the already added alternate codon overlaps the candidate's
		#right edge codon
		if($rec->{IDEN_START} == ($last_codon_stop - 2))
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
		      }
		  }

		#If the already added alternate codon overlaps the candidate's
		#identity region (not including the edge codons' identity).
		if($rec->{IDEN_START} >= ($first_codon_start + 3) &&
		   $rec->{IDEN_START} <= ($last_codon_stop - 3))
		  {
		    #If the candidate has any alternate codons for its edges,
		    #they do not matter in this case because they cannot
		    #overlap.

		    my $subseq =
		      substr($hash->{$pairid}->{SEQS}->{$seqid},
			     $rec->{IDEN_START} - 1,
			     3);

		    if($subseq ne $rec->{ALT_CODON})
		      {return(undef)}
		  }
	      }
	  }
      }

    return({LEFT_CODON  => ($no_left_cdn  ? undef : $alt_left_cdn),
	    RIGHT_CODON => ($no_right_cdn ? undef : $alt_right_cdn)});
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
#2 Sum of deviation distances from ideal spacing - minimize
#3 Number of bases (times number of sources) included in stretches of any size
#  - maximize
sub scoreSolution
  {
    my $soln      = $_[0];
    my $hash      = $_[1]; #Needed to get: ->{$pairid}->{NUMUNIQSEGS}->{$sz}
    my $max_size  = $_[2]; #The largest identity segment (i.e. "stretch") size
    my $score_obj = [];

    my $ubase_sums  = {}; #Hash used to sum of unique identity base positions
                          #included per pair/sequence
    my $dist_scores = {}; #A hash to keep track of the segment closest to a
                          #region center
    my $occupied    = {}; #A hash to count the number of segments included

    my($inclusion_score);
    my $spacing_score = 0;
    my $base_score    = 0;

    #For each sequence we are constructing
    foreach my $seqid (keys(%$soln))
      {
	foreach my $subseqrec (grep {$_->{TYPE} eq 'seg'} @{$soln->{$seqid}})
	  {
	    my $pairid = $subseqrec->{PAIR_ID};

	    ## Inclusion score calculations... Record the unique base locations

	    #For each segment of size $max_size included in the solution (not
	    #including alt. codons)
	    if(($subseqrec->{IDEN_STOP} - $subseqrec->{IDEN_START} + 1) ==
	       $max_size)
	      {
		#For each position of identity in this segment, record the
		#position as the key in the hash (then we'll count the number
		#of keys)
		foreach my $pos ($subseqrec->{IDEN_START}..
				 $subseqrec->{IDEN_STOP})
		  {$ubase_sums->{$pairid}->{$seqid}->{$pos} = 0}
	      }

	    ## Num score (number of non-identity-overlapped segments, any size)

	    foreach my $pos ($subseqrec->{IDEN_START}..$subseqrec->{IDEN_STOP})
	      {$occupied->{$pairid}->{$seqid}->{$pos} = 0}

	    ## Spacing score calculations for each region for each pair/seq...
	    ## Record the distance of the closest identity sequence from the
	    ## region center

	    my $middle = int($subseqrec->{IDEN_START} + $max_size / 2);
	    my $period_len =
	      int(length($hash->{$pairid}->{SEQS}->{$seqid}) /
		  ($hash->{$pairid}->{NUMUNIQSEGS}->{$max_size} + 1));
	    my $region = int($middle / $period_len);
	    if(!exists($dist_scores->{$pairid}) ||
	       !exists($dist_scores->{$pairid}->{$seqid}) ||
	       !exists($dist_scores->{$pairid}->{$seqid}->{$region}) ||
	       $subseqrec->{SPACING_SCORE} <
	       $dist_scores->{$pairid}->{$seqid}->{$region})
	      {$dist_scores->{$pairid}->{$seqid}->{$region} =
		 $subseqrec->{SPACING_SCORE}}
	  }
      }

    #For all possible pairs available
    foreach my $pairid (keys(%$hash))
      {
	## Spacing score calculations

	#For all possible sequences
	foreach my $seqid (keys(%{$hash->{$pairid}->{SEQS}}))
	  {
	    my $period_len =
	      int(length($hash->{$pairid}->{SEQS}->{$seqid}) /
		  ($hash->{$pairid}->{NUMUNIQSEGS}->{$max_size} + 1));
	    #For all possible regions
	    foreach my $reg (0..($hash->{$pairid}->{NUMUNIQSEGS}->{$max_size} -
				 1))
	      {
		if(exists($dist_scores->{$pairid}) &&
		   exists($dist_scores->{$pairid}->{$seqid}) &&
		   exists($dist_scores->{$pairid}->{$seqid}->{$reg}))
		  {$spacing_score += $dist_scores->{$pairid}->{$seqid}->{$reg}}
		else
		  {$spacing_score += $period_len}
	      }
	  }

	## Inclusion score calculations...

	if(!exists($ubase_sums->{$pairid}))
	  {
	    $inclusion_score = 0;
	    last;
	  }
	#For each sequence in the pair (Note: the score for each should
	#technically be the same, so doing it for both here is useless, though
	#it would be good to do a sanity check here in the future)
	foreach my $seqid (keys(%{$ubase_sums->{$pairid}}))
	  {
	    my $cur_size = 0;
	    my $tmp_inclusion_score = 0;
	    my($last_pos);
	    foreach my $pos (sort {$a <=> $b}
			     keys(%{$ubase_sums->{$pairid}->{$seqid}}))
	      {
		if(!defined($last_pos) || ($last_pos + 1) == $pos)
		  {$cur_size++}
		else
		  {
		    #Calculate the number of unique/idependent segments
		    my $usegs = int($cur_size / $max_size);
		    if($usegs > 0)
		      {$tmp_inclusion_score += $usegs}
		    $cur_size = 0;
		  }
		$last_pos = $pos;
	      }
	    if($tmp_inclusion_score <
	       $hash->{$pairid}->{NUMUNIQSEGS}->{$max_size} &&
	       (!defined($inclusion_score) ||
		$tmp_inclusion_score < $inclusion_score))
	      {$inclusion_score = $tmp_inclusion_score}
	  }
      }

    ## Bases score calculation

    foreach my $pairid (keys(%$occupied))
      {foreach my $seqid (keys(%{$occupied->{$pairid}}))
	 {$base_score += scalar(keys(%{$occupied->{$pairid}->{$seqid}}))}}

    push(@$score_obj,$inclusion_score,$spacing_score,$base_score);

    return(wantarray ? @$score_obj : $score_obj);
  }

sub candidateBetter
  {
    my $cand = $_[0];
    my $best = $_[1];

    if($cand->[0] > $best->[0] ||
       ($cand->[0] == $best->[0] && $cand->[1] < $best->[1]) ||
       ($cand->[0] == $best->[0] && $cand->[1] == $best->[1] &&
	$cand->[2] > $best->[2]))
      {return(1)}

    return(0);
  }

#Creates a hash of hybrid sequences
sub generateHybrids
  {
    my $soln   = $_[0];
    my $data   = $_[1];
    my $usage  = $_[2];
    my $matrix = $_[3];

    $soln = reduceSolution($soln);

    $soln = updateFinalCoords($soln,$data);

    my $partners_hash = getPartnersHash($soln,$data);

    my $map = getFilledMap($soln,$data,$partners_hash);

    my $woven_seqs = weaveSeqs($map,$data,$partners_hash,$matrix,$usage);

    return($woven_seqs);
  }

#Merges overlapping segments derived from the same pair
sub reduceSolution
  {
    my $soln = $_[0];
    my $new_soln = {};
    foreach my $seqid (keys(%$soln))
      {
	my @current = ();
	foreach my $currec (@{$soln->{$seqid}})
	  {
	    #Add any previous records that no longer overlap the current record
	    #to the new solution or else add them to the tmpcurrent array.  It
	    #is assumed that the records are in order of IDEN_START, and that
	    #that's always less than (or equal to) IDEN_STOP (and
	    #LAST_STOP_CODON.  We're also assuming that FIRST_START_CODON
	    #overlaps IDEN_START.  So once the IDEN_START is past
	    #LAST_CODON_STOP, then there is no more overlap with that record,
	    #thus we can add it to the reduced solution.
	    my @tmpcurrent = ();
	    foreach my $lastrec (@current)
	      {
		#If we have gone past this "last" record, it's good to add
		if($lastrec->{LAST_CODON_STOP} < $currec->{FIRST_CODON_START})
		  {push(@{$new_soln->{$seqid}},$lastrec)}
		else
		  {push(@tmpcurrent,$lastrec)}
	      }
	    @current = @tmpcurrent;

	    #If there are no remaining records in the current array or if this
	    #is an alternate codon (which does not get merged with other
	    #records), add this record to the current array
	    if(scalar(@current) == 0 || $currec->{TYPE} eq 'alt')
	      {push(@current,{%$currec})}
	    else
	      {
		my $add = 1;
		#The remaining records all overlap.  Search through them to
		#find the one from the same pair source and edit its record to
		#extend the stops (unless one from the same pair doesn't exist
		#- then "$add" it to the current array as its own thing to be
		#searched for overlap on the next iteration)
		foreach my $lastrec (grep {$_->{TYPE} eq 'seg' &&
					     $currec->{IDEN_STOP} >
					       $_->{IDEN_STOP} &&
						 $_->{PAIR_ID} eq
						   $currec->{PAIR_ID}}
				     @current)
		  {
		    $lastrec->{IDEN_STOP} = $currec->{IDEN_STOP};
		    $lastrec->{LAST_CODON_STOP} = $currec->{LAST_CODON_STOP};
		    $lastrec->{FINALN_STOP} = $currec->{FINAL_STOP};
		    #Keep the score of the better placed piece of identity
		    if($currec->{SPACING_SCORE} < $lastrec->{SPACING_SCORE})
		      {$lastrec->{SPACING_SCORE} = $currec->{SPACING_SCORE}}
		    $add = 0;
		    last;
		  }
		if($add)
		  {push(@current,{%$currec})}
	      }
	  }
	#All the records have been traversed.  Anything left in @current can be
	#added to the reduced solution
	push(@{$new_soln->{$seqid}},@current);
      }

    return($new_soln);
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

    foreach my $seqid (keys(%$soln))
      {
	my @current = ();
	#For each segment record (in order of ascending IDEN_START)
	foreach my $currec (@{$soln->{$seqid}})
	  {
	    #Error-check to make sure that the FINAL values are fresh
	    if($currec->{FIRST_CODON_START} != $currec->{FINAL_START})
	      {
		warning("Final start coordinate has an unexpected value.  ",
			"Resetting.",
			{DETAIL =>
			 join(',',("updateFinalCoords should only be called ",
				   "once on a solution.  It expects the ",
				   "final start and stop to be the same as ",
				   "the first codon start and the last codon ",
				   "stop."))});
		$currec->{FINAL_START} = $currec->{FIRST_CODON_START};
	      }
	    if($currec->{LAST_CODON_STOP} != $currec->{FINAL_STOP})
	      {
		warning("Final stop coordinate has an unexpected value.  ",
			"Resetting.",
			{DETAIL =>
			 join(',',("updateFinalCoords should only be called ",
				   "once on a solution.  It expects the ",
				   "final start and stop to be the same as ",
				   "the first codon start and the last codon ",
				   "stop."))});
		$currec->{FINAL_STOP} = $currec->{LAST_CODON_STOP};
	      }

	    #Add any previous records that no longer overlap the current record
	    #to the new solution or else add them to the tmpcurrent array.  It
	    #is assumed that the records are in order of IDEN_START, and that
	    #that's always less than (or equal to) IDEN_STOP (and
	    #LAST_STOP_CODON.  We're also assuming that FIRST_START_CODON
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
		    unless(($currec->{FINAL_START} >= $_->{FINAL_START} &&
			    $currec->{FINAL_START} <= $_->{FINAL_STOP}) ||
			   ($currec->{FINAL_STOP}  >= $_->{FINAL_START} &&
			    $currec->{FINAL_STOP}  <= $_->{FINAL_STOP}) ||
			   ($currec->{FINAL_START} <  $_->{FINAL_START} &&
			    $currec->{FINAL_STOP}  >  $_->{FINAL_STOP}))
		      {next}

		    #If lastrec is an alternate codon (and currec is not)
		    if($lastrec->{TYPE} eq 'alt' && $currec->{TYPE} eq 'seg')
		      {clipForAltCodon($lastrec,$currec)}
		    #Else if currec is an alternate codon (and lastrec is not)
		    elsif($lastrec->{TYPE} eq 'seg' &&
			  $currec->{TYPE} eq 'alt')
		      {clipForAltCodon($currec,$lastrec)}
		    #Else if neither are alternate codon
		    elsif($lastrec->{TYPE} eq 'seg' &&
			  $currec->{TYPE} eq 'seg')
		      {clipSegments($lastrec,$currec)}
		    #Else error
		    else
		      {error("Overlapping alternate codons encountered.  ",
			     "Skipping.")}
		  }

		#When segments are completely overlapped by other segments,
		#their FINAL_STOP is set to 0.  Other math in adjusting for
		#overlap can also result in the FINAL_START being larger than
		#the FINAL_STOP when total overlap exists.  These segments can
		#be skipped.
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

    #If the alternate codon overlaps the first codon and
    #does not overlap the second codon (assumes all are in
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

    #If the stop of segment 2 occurs before the stop of segment 1, set the stop
    #of segment 2 to 0 so that it will get eliminated
    if($seg2->{FINAL_STOP} <= $seg1->{FINAL_STOP} &&
       $seg2->{FINAL_STOP} >= $seg1->{FINAL_START})
      {$seg2->{FINAL_STOP} = 0}
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
		if($seg2->{FINAL_START} <= $seg1->{FINAL_STOP})
		  {
		    $seg1->{FINAL_STOP}  = $seg1->{IDEN_STOP};
		    $seg2->{FINAL_START} = $seg1->{IDEN_STOP} + 1;
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
		    $seg1->{FINAL_STOP}  = $seg1->{IDEN_STOP};
		    $seg2->{FINAL_START} = $seg1->{IDEN_STOP} + 1;
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
		    $seg1->{FINAL_STOP}  = $seg1->{IDEN_STOP};
		    $seg2->{FINAL_START} = $seg1->{IDEN_STOP} + 1;
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
	    #alternate codon if there's a conflict, so we'll just end the final
	    #stop of seg1 just before the identity start of seg2

	    #If there is current overlap
	    if($seg2->{FINAL_START} <= $seg1->{FINAL_STOP})
	      {
		$seg1->{FINAL_STOP}  = $seg2->{IDEN_START} - 1;
		$seg2->{FINAL_START} = $seg1->{IDEN_START};
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
		    $seg1->{FINAL_STOP}  = $seg1->{IDEN_STOP};
		    $seg2->{FINAL_START} = $seg1->{IDEN_STOP} + 1;
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
    my $div_map   = {};

    foreach my $seqid (keys(%$soln))
      {
	my $lastind = -1;
	my $cnt     = 0;
	my @ordered_indexes = sort {$soln->{$seqid}->[$a]->{FINAL_START} <=>
				      $soln->{$seqid}->[$b]->{FINAL_START}}
	  0..$#{$soln->{$seqid}};

	#For each segment record (in order of ascending FINAL_START)
	foreach my $pseudoind (@ordered_indexes)
	  {
	    my $recind  = $ordered_indexes[$pseudoind];
	    my $nextind = (($pseudoind + 1) > $#ordered_indexes ?
			   undef : $ordered_indexes[$pseudoind + 1]);
	    my $rec     = $soln->{$seqid}->[$recind];

	    #Obtain all the pair IDs and Sequence IDs of all the sequences that
	    #$seqid has been aligned with
	    my $partners = $part_hash->{$seqid};

	    #Using the partner coordinates relative to this FINAL_START,
	    #look for a any segment to the left among the partner sequences
	    #this sequence has been paired with.  If any exists, use the
	    #closest one's FINAL_STOP to determine this sequence's relative
	    #coord.  Then determine a mid-point on a codon boundary and
	    #record it in this sequence's map and copy/convert it to all
	    #partner sequences's div_maps.

	    my $divider_left = getClosestLeftDivider($seqid,
						     $rec->{FINAL_START},
						     $soln,
						     $data,
						     $partners);

	    #Set the default as the beginning of the sequence just in case
	    if(!defined($divider_left))
	      {$divider_left = 1}

	    #This divider might have been copied from another sequence, so
	    #check it for a conflict before recording it
	    if(exists($div_map->{$seqid}) &&
	       exists($div_map->{$seqid}->{$divider_left}) &&
	       $div_map->{$seqid}->{$divider_left}->{TYPE} eq 'SOURCE' &&
	       $div_map->{$seqid}->{$divider_left}->{PAIR_ID} ne
	       $rec->{PAIR_ID})
	      {error("Sequence [$seqid] starting at ",
		     "[$divider_left] has a conflict between alignments [",
		     $div_map->{$seqid}->{$divider_left}->{PAIR_ID},
		     "] and [$rec->{PAIR_ID}].  Unable to encode ",
		     "segment from alignment [$rec->{PAIR_ID}].")}
	    else
	      {$div_map->{$seqid}->{$divider_left} =
		 {PAIR_ID => $rec->{PAIR_ID},
		  TYPE    => 'SOURCE',
		  STOP    => undef,
		  SEQ     => ''}}

	    #Copy this divider to the partner sequences
	    copyDivider($divider_left,
			$seqid,
			$div_map,
			$partners,
			$data,
			$rec->{PAIR_ID});

	    #Using the partner coordinates relative to this FINAL_STOP,
	    #look for a any segment to the right among the partner
	    #sequences this sequence has been paired with.  If any exists,
	    #use the closest one's FINAL_START to determine this sequence's
	    #relative coord.  Then determine a mid-point on a codon
	    #boundary and record it in this sequence's map and copy/convert
	    #it to all partner sequences's div_maps

	    my $divider_right = getClosestRightDivider($seqid,
						       $rec->{FINAL_STOP},
						       $soln,
						       $data,
						       $partners);

	    #Set the default as just after the end of the sequence just in case
	    if(!defined($divider_right))
	      {$divider_right =
		 length($data->{$rec->{PAIR_ID}}->{SEQS}->{$seqid}) + 1}

	    #This divider might have been copied from another sequence, so
	    #check it for a conflict before recording it
	    if(exists($div_map->{$seqid}) &&
	       exists($div_map->{$seqid}->{$divider_right}) &&
	       $div_map->{$seqid}->{$divider_right}->{TYPE} eq 'SOURCE' &&
	       $div_map->{$seqid}->{$divider_right}->{PAIR_ID} ne
	       $rec->{PAIR_ID})
	      {error("Sequence [$seqid] starting at ",
		     "[$divider_right] has a conflict between alignments [",
		     $div_map->{$seqid}->{$divider_right}->{PAIR_ID},
		     "] and [$rec->{PAIR_ID}].  Unable to encode ",
		     "segment from alignment [$rec->{PAIR_ID}].")}
	    else
	      {$div_map->{$seqid}->{$divider_right} =
		 {PAIR_ID => $rec->{PAIR_ID},
		  TYPE    => 'SOURCE',
		  STOP    => undef,
		  SEQ     => ''}}

	    #Copy this divider to the partner sequences
	    copyDivider($divider_right,
			$seqid,
			$div_map,
			$partners,
			$data,
			$rec->{PAIR_ID});

	    $lastind = $recind;
	    $cnt++;
	  }
      }

    return($div_map);
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

    #For each alignment $seqid has been aligned with
    my $closest = 0;
    foreach my $partner (@$partners)
      {
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

	#If this coordinate is closer, save it.
	if($orig_left_seg_coord > $closest)
	  {
	    $closest = $orig_left_seg_coord;
	    if(($orig_left_seg_coord + 1) == $coord)
	      {return($coord)}
	  }
      }

    #We now have the closest occupied coordinate to the left of $coord.  We
    #need to split the difference and then return the closest frame 1 codon
    #position.  We add 1 to the closest coord because we want to choose a
    #midpoint among available positions to put a divider.  That includes the
    #original $coord position, but not the last position of the nearest
    #neighbor to the left ($closest)
    my $midpoint = int(($closest + 1 + $coord) / 2);

    #Determine the frame position (1, 2, or 3) of the coordinate
    my $frame_pos = ($midpoint - 1) % 3 + 1;

    if($frame_pos == 1)
      {return($midpoint)}
    elsif($frame_pos == 2)
      {
	if(($midpoint - 1) <= $closest)
	  {
	    #This should not happen, but just in case...
	    if(($midpoint + 2) > $coord)
	      {
		error("Bad frame boundary.");
		return($coord);
	      }

	    return($midpoint + 2);
	  }

	return($midpoint - 1);
      }
    elsif($frame_pos == 3)
      {
	if(($midpoint + 1) > $coord)
	  {
	    error("Bad frame boundary.");
	    return($coord);
	  }

	return($midpoint + 1);
      }
    else
      {error("Unexpected frame result.")}

    return($coord);
  }

#Searches the segments of a solution for the closest right-most coordinate that
#is to the left of a given coordinate.  A linear search is performed.  Could
#improve this to be a binary search.  Note: Returns 1 if there are no segments
#to the left or if the coordinate passed in was 1.
sub getClosestLeftSegCoord
  {
    my $seqid = $_[0];
    my $coord = $_[1];
    my $soln  = $_[2];

    my $closest = 1;

    foreach my $rec (sort {$a->{FINAL_START} <=> $b->{FINAL_START}}
		     @{$soln->{$seqid}})
      {
	if($rec->{FINAL_START} > $coord)
	  {last}
	if($rec->{FINAL_STOP} < $coord)
	  {$closest = $rec->{FINAL_STOP}}
	elsif($rec->{FINAL_START} < $coord && $rec->{FINAL_STOP} >= $coord)
	  {
	    $closest = $coord - 1;
	    last;
	  }
      }

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

    #For each alignment $seqid has been aligned with
    my($closest);
    foreach my $partner (@$partners)
      {
	if(!defined($closest))
	  {$closest = length($data->{$partner->[0]}->{SEQS}->{$seqid})}

	#Convert the coordinate to the partner sequence
	my $partner_coord = convertDivider($seqid,
					   $coord,
					   $partner->[0],  #pair ID
					   $partner->[1],  #partner seq ID
					   $data);

	#Find the closest coord to the left
	my $part_rght_seg_coord =
	  getClosestRightSegCoord($partner->[1], #seqid
				  $partner_coord,
				  $soln,
				  length($data->{$partner->[0]}->{SEQS}
					 ->{$partner->[1]}));

	#Convert the coordinate back to this sequence
	my $orig_rght_seg_coord = convertDivider($partner->[1],  #partner seqID
						 $part_rght_seg_coord,
						 $partner->[0],  #pair ID
						 $seqid,
						 $data);

	#If this coordinate is closer, save it.
	if($orig_rght_seg_coord < $closest)
	  {
	    $closest = $orig_rght_seg_coord;
	    if(($orig_rght_seg_coord - 1) == $coord)
	      {return($orig_rght_seg_coord)}
	  }
      }

    #We now have the closest occupied coordinate to the right of $coord.  We
    #need to split the difference and then return the closest frame 1 codon
    #position.  We add 1 to $coord because we want to choose a midpoint among
    #available positions to put a divider.  That includes the $closest
    #position, but not the $coordt position (the last position of the query
    #segment)
    my $midpoint = int(($closest + $coord + 1) / 2);

    #Determine the frame position (1, 2, or 3) of the coordinate
    my $frame_pos = ($midpoint - 1) % 3 + 1;

    if($frame_pos == 1)
      {return($midpoint)}
    elsif($frame_pos == 2)
      {
	if(($midpoint - 1) <= $coord)
	  {
	    #This should not happen, but just in case...
	    if(($midpoint + 2) > $closest)
	      {
		error("Bad frame boundary.");
		return($coord + 1);
	      }

	    return($midpoint + 2);
	  }

	return($midpoint - 1);
      }
    elsif($frame_pos == 3)
      {
	if(($midpoint + 1) > $closest)
	  {
	    error("Bad frame boundary.");
	    return($coord + 1);
	  }

	return($midpoint + 1);
      }
    else
      {error("Unexpected frame result.")}

    return($coord + 1);
  }

#Searches the segments of a solution for the closest left-most coordinate that
#is to the right of a given coordinate.  A linear search is performed.  Could
#improve this to be a binary search.  Note: Returns sequence length +1 if there
#are no segments to the right or if the coordinate passed in was the sequence
#length.
sub getClosestRightSegCoord
  {
    my $seqid  = $_[0];
    my $coord  = $_[1];
    my $soln   = $_[2];
    my $seqlen = $_[3];

    my $closest = $seqlen + 1;

    foreach my $rec (sort {$b->{FINAL_START} <=> $a->{FINAL_START}}
		     @{$soln->{$seqid}})
      {
	if($rec->{FINAL_STOP} < $coord)
	  {last}
	if($rec->{FINAL_START} > $coord)
	  {$closest = $rec->{FINAL_START}}
	elsif($rec->{FINAL_START} <= $coord && $rec->{FINAL_STOP} > $coord)
	  {
	    $closest = $coord + 1;
	    last;
	  }
      }

    return($closest);
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

    my $uniq_seq_partners = {};
    foreach my $partner_seqid (grep {$_ ne $seqid} map {$_->[1]} @$partners)
      {$uniq_seq_partners->{$partner_seqid} = 0}

    foreach my $partner (@$partners)
      {
	my $partner_pair_id = $partner->[0];
	my $partner_seqid   = $partner->[1];
	my $partner_divider = convertDivider($seqid,
					     $divider,
					     $partner_pair_id,
					     $partner_seqid,
					     $data);

	#This divider might have been copied from another sequence, so
	#check it for a conflict before recording it

	#If this divider doesn't already exist in this sequence
	if(!exists($div_map->{$partner_seqid}) ||
	   !exists($div_map->{$partner_seqid}->{$partner_divider}))
	  {
	    #If this sequence is a part of the same pair as the originating
	    #sequence, set the pair ID in the div map
	    if($pair_id eq $partner_pair_id)
	      {$div_map->{$partner_seqid}->{$partner_divider} =
		 {PAIR_ID => $pair_id,
		  TYPE    => 'SOURCE',
		  STOP    => undef,
		  SEQ     => ''}}
	    else
	      {$div_map->{$partner_seqid}->{$partner_divider} =
		 {PAIR_ID => undef,    #Not needed for type 'recode'. Obtained
		                       #from partners array
		  TYPE    => 'RECODE',
		  STOP    => undef,
		  SEQ     => ''}}
	  }
	elsif($div_map->{$partner_seqid}->{$partner_divider}->{TYPE} eq
	      'SOURCE' &&
	      $div_map->{$partner_seqid}->{$partner_divider}->{PAIR_ID} ne
	      $pair_id)
	  {error("Partner sequence [$partner_seqid] starting at ",
		 "[$partner_divider] has a conflict between alignments [",
		 $div_map->{$partner_seqid}->{$partner_divider}->{PAIR_ID},
		 "] and [$pair_id].  Unable to copy segment divider ",
		 "[$seqid:$divider] from alignment [$pair_id] to the ",
		 "partner sequence.")}
	#Else if this is a divider of type recode, but we're being copied from
	#the same pair, change it to source
	elsif($div_map->{$partner_seqid}->{$partner_divider}->{TYPE} ne
	      'SOURCE' && $pair_id eq $partner_pair_id)
	  {
	    $div_map->{$partner_seqid}->{$partner_divider}->{PAIR_ID} =
	      $partner_pair_id;
	    $div_map->{$partner_seqid}->{$partner_divider}->{TYPE} = 'SOURCE';
	  }
	#Else TYPE is 'SOURCE' and it's equal to $pair_id anyway, or it's not
	#'SOURCE' and this isn't the source alignment.  In the second case, it
	#means that there are multiple options to recode for (either for a
	#segment one of its partners from one of its alignments where it
	#matched a third sequence, or another of its partners with another
	#sequence from another alignment).  We will arbitrarily go with the
	#first one we encounter.  This is unlikely to happen if a user is
	#interested in one target sequence with multiple partners, but it could
	#happen if they submit all pairs of more than 3 sequences.  In either
	#case, we do not need to change anything.
      }
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

    if($query_coord > length($data->{$pairid}->{SEQS}->{$query_seqid}))
      {return(length($data->{$pairid}->{SEQS}->{$target_seqid}) + 1)}
    elsif($query_coord == length($data->{$pairid}->{SEQS}->{$query_seqid}))
      {return(length($data->{$pairid}->{SEQS}->{$target_seqid}))}

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

	#If we are at the query coordinate we were looking for and it
	#corresponds to a gap character in the target alignment sequence,
	#increment the target position, because the last target position was
	#aligned with a different query position.  Note, this only works for
	#the start of segments.  End coordinates should not be queried.
	if($query_pos == $query_coord && $target_char eq '-')
	  {$target_pos++}

	$aln_pos++;
      }

    return($target_pos);
  }

#Assumes the last divider is 1 past the sequence length
sub fillMapStops
  {
    my $div_map = $_[0];

    foreach my $seqid (keys(%$div_map))
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

sub getPartnersHash
  {
    my $soln  = $_[0];
    my $data  = $_[1];
    my $partners_hash = {};

    foreach my $seqid (keys(%$soln))
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

    my $partners = [grep {$_ ne $seqid}
		    map {my $pid=$_;
			 map {[$pid,$_]} keys(%{$data->{$pid}->{SEQS}})}
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

    my $num_filled = 1;
    my $unfinished = 1;

    #While all segments' alignment sources have not been determined and we're
    #still making progress
    while($unfinished && $num_filled)
      {
	#We will assume that this pass through the loop will complete the map
	#unless we hit something we cannot fill in
	$unfinished = 0;
	#We will require that something has to have been added to the map on
	#each iteration.
	$num_filled = 0;

	#For each each sequence
	foreach my $seqid (keys(%$map))
	  {
	    #Cycle through the subsequence records
	    foreach my $divider (keys(%{$map->{$seqid}}))
	      {
		if($map->{$seqid}->{$divider}->{SEQ} eq '')
		  {
		    #If the source type of the sequence is RECODE, we need to
		    #determine what sequence this sequence's recoding should be
		    #based on.  We will arbitrarily pick the forst completely
		    #overlapping and filled in partner
		    if($map->{$seqid}->{$divider}->{TYPE} eq 'RECODE')
		      {
			if(reCodeLeftoverSegment($seqid,
						 $divider,
						 $map,
						 $partners_hash->{$seqid},
						 $matrix,
						 $data,
						 $codon_hash))
			  {$num_filled++}
			else
			  {$unfinished = 0}
		      }
		    else #TYPE is 'SOURCE' - so grade the sequence directly
		      {
			$map->{$seqid}->{$divider}->{SEQ} =
			  substr($data->{$map->{$seqid}->{$divider}->{PAIR_ID}}
				 ->{SEQS}->{$seqid},
				 ($divider - 1),
				 ($map->{$seqid}->{$divider}->{STOP} -
				  $divider + 1));
			$num_filled++;
		      }
		  }
	      }
	  }
      }

    if($unfinished)
      {error("Unexpected case.  Unable to determine which alignment every ",
	     "segment of sequence (that will not be recoded) should come ",
	     "from.")}

    #Now assemble the sequences
    my $seqhash = {};

    foreach my $seqid (keys(%$map))
      {
	my $lastdiv = 0;
	foreach my $divider (sort {$a <=> $b} keys(%{$map->{$seqid}}))
	  {
	    if(($lastdiv + 1) != $divider)
	      {
		my $nstart = $lastdiv + 1;
		my $nstop  = $divider - 1;
		my $nsize = $nstop - $nstart + 1;
		error("Unable to determine sequence [$seqid] from [$nstart] ",
		      "to [$nstop].  This should not have happened.  ",
		      "Inserting 'Ns' so that you get something out.");
		$seqhash->{$seqid} .= ('N' x $nsize);
	      }
	    else
	      {$seqhash->{$seqid} .= $map->{$seqid}->{$divider}->{SEQ}}
	  }
      }

    return($seqhash);
  }

#This sub searches through a sequence's partners in the map to find one that
#has sequence that completely overlaps.  If there is not complete overlap, it
#returns 0 and does not do any recoding.  If there is complete existing
#overlap, it obtains the sequence to be recoded and then recodes each codon to
#best match (by homology) the partner sequence (without changing the AA that
#was originally encoded).
sub reCodeLeftoverSegment
  {
    my $seqid      = $_[0];
    my $divider    = $_[1];
    my $map        = $_[2];
    my $partners   = $_[3];
    my $matrix     = $_[4];
    my $data       = $_[5];
    my $codon_hash = $_[6];

    my $new_seq  = '';

    #Foreach partner sequence (that $seqid has been paired/aligned with)
    foreach my $partner (@$partners)
      {
	my $pairid        = $partner->[0];
	my $partner_seqid = $partner->[1];

	#We need to obtain the aligned sequence corresponding to the positions
	#in each sequence, if it's populated in the map
	my($fixed_partner_aln,$change_this_aln) =
	  getPartnerAlnSeqs($map,
			    $data,
			    $pairid,
			    $seqid,
			    $divider,
			    $map->{$seqid}->{$divider}->{STOP},
			    $partner_seqid);

	if(defined($fixed_partner_aln) && $fixed_partner_aln ne '' &&
	   defined($change_this_aln) && $change_this_aln ne '')
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
      }

    #Could not find a completed partner sequence
    return(0);
  }

#Compiles a portion of sequence (including alignment characters) from one or
#more sections from the map if those sections have all been populated with
#sequence.  This performs a linear search.  It could be improved to do a binary
#search.
sub getPartnerAlnSeqs
  {
    my $map          = $_[0];
    my $data         = $_[1];
    my $pairid       = $_[2];
    my $change_seqid = $_[3];
    my $change_start = $_[4];
    my $change_stop  = $_[5];
    my $fixed_seqid  = $_[6];

    my $fixed_start = convertDivider($change_seqid,
				     $change_start,
				     $pairid,
				     $fixed_seqid,
				     $data);
    my $fixed_stop  = convertDivider($change_seqid,
				     $map->{$change_seqid}->{$change_start}
				     ->{STOP} + 1,
				     $pairid,
				     $fixed_seqid,
				     $data) - 1;
    my $fixed_size  = $fixed_stop - $fixed_start + 1;

    #We need to first obtain the fixed sequence from the map.  It may actually
    #differ from its representation in the alignment/pair because it itself may
    #have been recoded.  We can determine whether it was recoded or not by its
    #type.  If its type is 'RECODE', then we need to obtain the alignment
    #sequence from the pair and alter its code to what the map contained
    #(preserving the gap characters/positions)
    my $fixed_aln = '';

    #Also obtain the sequence to change (based on the fixed sequence
    my $change_aln = '';

    #Keep track of whether any overlapping portion was sourced from an
    #alignment or whether it itself was recoded
    my($sourced);

    #These values will change
    my $cur_change_start = $change_start;
    my $cur_change_stop  = $change_stop;

    #For each section in order of start coordinate
    foreach my $cur_fixed_start (sort {$a <=> $b}
				 keys(%{$map->{$fixed_seqid}}))
      {
	my $cur_fixed_stop = $map->{$fixed_seqid}->{$cur_fixed_start}->{STOP};

	#If there is overlap with the target fixed start and stop
	if(coordsOverlap($cur_fixed_start,$cur_fixed_stop,
			 $fixed_start,$fixed_stop))
	  {
	    if($map->{$fixed_seqid}->{$cur_fixed_start}->{TYPE} eq 'RECODE')
	      {$sourced = 0}
	    else
	      {$sourced = 1}

	    #If there is no sequence present yet for this section, we cannot
	    #proceed
	    if($map->{$fixed_seqid}->{$cur_fixed_start}->{SEQ} eq '')
	      {return(0)}

	    my $substr_start = 0;

	    #If the sequence we're looking for starts somewhere in the middle
	    #of the current section, bump up the relative start
	    if($cur_fixed_start < $fixed_start)
	      {$substr_start = $fixed_start - $cur_fixed_start}

	    my $tmp_fixed_start = $cur_fixed_start + $substr_start;
	    my $tmp_fixed_stop  = ($cur_fixed_stop <= $fixed_stop ?
				   $cur_fixed_stop : $fixed_stop);

	    #Obtain the current change start/stop
	    my $cur_change_start = convertDivider($fixed_seqid,
						  $tmp_fixed_start,
						  $pairid,
						  $change_seqid,
						  $data);
	    if($cur_change_start < $change_start)
	      {$cur_change_start = $change_start}
	    my $cur_change_stop  = convertDivider($fixed_seqid,
						  $tmp_fixed_stop + 1,
						  $pairid,
						  $change_seqid,
						  $data) - 1;
	    if($cur_change_stop > $change_stop)
	      {$cur_change_stop = $change_stop}

	    my($next_fixed_seq,$next_change_seq) =
	      getAlnSegment($pairid,
			    $fixed_seqid,
			    $tmp_fixed_start,
			    $tmp_fixed_stop,
			    $change_seqid,
			    $cur_change_start,
			    $cur_change_stop,
			    $data);

	    #If the requence obtained has been previously recoded, replace the
	    #sequence characters in the alignment obtained with the recoded
	    #sequence that was saved in the map object (without alignment
	    #characters)
	    if(!$sourced)
	      {
		#Obtain the recoded sequence to use to replace the sequence in
		#the aligned version of this sequence
		my $replace_seq = '';

		#If this represents only a portion of the saved/recoded
		#sequence, grab that subsequence from the map object
		if($cur_fixed_stop > $fixed_stop)
		  {
		    my $remainder =
		      length($map->{$fixed_seqid}->{$cur_fixed_start}->{SEQ}) -
			$substr_start;
		    $remainder -= ($cur_fixed_stop - $fixed_stop);
		    $replace_seq =
		      substr($map->{$fixed_seqid}->{$cur_fixed_start}->{SEQ},
			     $substr_start,$remainder);
		  }
		else
		  {$replace_seq =
		     substr($map->{$fixed_seqid}->{$cur_fixed_start}->{SEQ},
			    $substr_start)}

		$next_fixed_seq = replaceAlnSegment($next_fixed_seq,
						    $replace_seq);
	      }

	    $fixed_aln  .= $next_fixed_seq;
	    $change_aln .= $next_change_seq;

	    #If we've gotten all the sequence we were looking for, move on
	    last if($cur_change_stop >= $change_stop);
	  }
      }

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
    my $get_start   = $_[2];
    my $get_stop    = $_[3];
    my $query_seqid = $_[4];
    my $query_start = $_[5];
    my $query_stop  = $_[6];
    my $data        = $_[7];

    my $target_aln = $data->{$pairid}->{ALNS}->{$getseqid};
    my $query_aln  = $data->{$pairid}->{ALNS}->{$query_seqid};

    my $aln_pos        = 0;
    my $cur_target_pos = 0;
    my $cur_query_pos  = 0;
    my $max_pos = length($target_aln);

    my $target_aln_segment = '';
    my $query_aln_segment  = '';

    #While we are not yet at the end of the captured segment
    while(($cur_target_pos < $get_stop || $cur_query_pos < $query_stop) &&
	  $aln_pos < $max_pos)
      {
	my $target_char = substr($target_aln,$aln_pos,1);
	my $query_char  = substr($query_aln, $aln_pos,1);

	$aln_pos++;

	if($target_char ne '-')
	  {$cur_target_pos++}
	if($query_char ne '-')
	  {$cur_query_pos++}

	if($cur_target_pos >= $get_start || $cur_query_pos >= $query_start)
	  {
	    if(($cur_target_pos + 1) < $get_start ||
	       ($cur_query_pos  + 1) < $query_start)
	      {error("The sequence positions (target: [$get_start] & query: ",
		     "[$query_start]) submitted for the 2 aligned sequences ",
		     "in alignment [$pairid] do not agree.",
		     {DETAIL =>
		      join('',("They should start at the same place in the ",
			       "alignment [$aln_pos] (or one less when there ",
			       "are gap characters) but the current ",
			       "positions are: (target: [$cur_target_pos] & ",
			       "query: [$cur_query_pos])."))})}
	    if(($cur_target_pos) > $get_stop ||
	       ($cur_query_pos)  > $query_stop)
	      {error("The sequence positions (target: [$get_stop] & query: ",
		     "[$query_stop]) submitted for the 2 aligned sequences ",
		     "in alignment [$pairid] do not agree.",
		     {DETAIL =>
		      join('',("They should stop at the same place in the ",
			       "alignment [$aln_pos] (or one less when there ",
			       "are gap characters) but the current ",
			       "positions are: (target: [$cur_target_pos] & ",
			       "query: [$cur_query_pos])."))})}
	    $target_aln_segment .= $target_char;
	    $query_aln_segment  .= $query_char;
	  }
      }

    return($target_aln_segment,$cur_query_pos);
  }

#Replace the sequence portion of an alignment with the characters from another
#unaligned sequence
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
	    $new_aln .= substr($replace_seq,$seq_pos,1);
	    $seq_pos++;
	  }
	else
	  {$new_aln .= '-'}
      }

    if($seq_pos != length($replace_seq))
      {error("The number of non-gap characters in the alignment sequence ",
	     "[$aln_seq] was not the same as the unaligned replacement ",
	     "sequence: [$replace_seq].")}

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
    my $seq_hash = $_[0];
    my $outfile  = $_[1];

    openOut(*OUT,$outfile) || return();

    foreach my $seqid (keys(%$seq_hash))
      {print(">$seqid\n",$seq_hash->{$seqid},"\n")}

    closeOut(*OUT);
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

    while(my $rec = getNextSeqRec(*ALN,0,$alnfile))
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

	$hash->{$pairid}->{ALNS}->{$id} = $rec->[1];
	$hash->{$pairid}->{SEQS}->{$id} = $rec->[1];
	$hash->{$pairid}->{SEQS}->{$id} =~ s/-+//g;
      }

    closeIn(*ALN);

    if($recnum < 2)
      {
	error("Only [$recnum] aligned sequences were found in alignment file ",
	      "[$alnfile].  Two sequences are required in each alignment ",
	      "file.");
	return(5);
      }

    return(loadSegments($hash->{$pairid},
			$stretch_sizes,
			#Sequence IDs (sorted) - can only be 2 keys
			(sort(keys(%{$hash->{$pairid}->{ALNS}}))),
			#Sequences (sorted by IDs) - can only be 2 keys
			(map {$hash->{$pairid}->{ALNS}->{$_}}
			 sort(keys(%{$hash->{$pairid}->{ALNS}})))));
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
    my $id1           = $_[1]; #SeqID1
    my $id2           = $_[2]; #SeqID2
    my $aln1          = $_[3];
    my $aln2          = $_[4];

    my $aln_len = length($aln1);
    if(length($aln2) != $aln_len)
      {error("Aligned sequences not the same length.")}

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

	#Record the regions as indexes of identity segments broken up by the
	#non-overlapping region they are closest to (given an ideal spacing
	#scheme).  These indexes will be used later to try and select segments
	#in a greedy fashion such that the first ones we attempt to add will
	#fit something from every alignment
	if(scalar(@{$regions->{$sz}->{$id1}}) >
	   scalar(@{$regions->{$sz}->{$id2}}) ||
	   (scalar(@{$regions->{$sz}->{$id1}}) ==
	    scalar(@{$regions->{$sz}->{$id2}}) && $score1 < $score2))
	  {
	    my $regidx = 0;
	    foreach my $regkey (sort {$a <=> $b}
				keys(%{$regions->{$sz}->{$id1}}))
	      {
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
		$hash->{REGS}->{$sz}->[$regidx] =
		  [map {$_->[0]} sort {$a->[1] <=> $b->[1]}
		   @{$regions->{$sz}->{$id2}->{$regkey}}];
		$regidx++;
	      }
	  }

	$hash->{NUMUNIQSEGS}->{$sz}  = $num_unique_segs;
	$hash->{NUMBASES}->{$sz}     = $num_seg_bases;
	$hash->{SPACINGSCORE}->{$sz} = $spacing_score;
      }

    return(0);
  }






















