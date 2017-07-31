#!/usr/bin/perl -w

#USAGE: Run with no options to get usage or with --help for basic details

use CommandLineInterface;
use warnings;
use strict;

##
## Describe the script
##

our $VERSION = '0.001';

setScriptInfo(CREATED => '7/31/2017',
              VERSION => $VERSION,
              AUTHOR  => 'Robert William Leach',
              CONTACT => 'rleach@princeton.edu',
              COMPANY => 'Princeton University',
              LICENSE => 'Copyright 2017',
              HELP    => << 'END_HELP'

This script, given a set of very dissimilar amino acid sequences, will locate all (sequentially ordered) stretches of identical residues between every possible pair and output a table of coordinates.  The coordinate table can then be supplied to the companion script "codonHomologizer.pl" to force the alignment to align these stretches of identity together.  It will preferentially select a combination of stretches that will create the fewest gaps in the final alignment produced by codonHomologizer (ignoring end gaps).

A word of caution: This is an exhaustive & recursive search.  There are combinations of parameters that can yield very long runtimes.  The search can be very fast by narrowing the relative coordinate window in which to search for identity.  A minimum of 3 idential AA residues is enforced, but you can use --force to over-ride the restriction.

END_HELP
	      ,
	      DETAILED_HELP => << 'END_DETAIL'

This script finds identical stretches of the requested length only, however if a stretch of identity is some multiple of the requested stretch size, the resulting coordinates in the output table will be merged.

A hash is used in the search.  The search is case-insensitive.

This script was created for very dissimilar sequences because existing alignment tools do not align to optimize stretches of identity, yet stretches of identity are required for crossover events to occur.*  If your sequences are very similar, this script is unnecessary.

* https://www.ncbi.nlm.nih.gov/pubmed/27671930

END_DETAIL
);

setDefaults(HEADER        => 1,
	    ERRLIMIT      => 3,
	    COLLISIONMODE => 'error', #,merge,rename (when outfiles conflict)
	    DEFRUNMODE    => 'usage');

##
## Configure the command line interface
##

my $mat_min  = 3;
my $mat_size = $mat_min;
addOption(GETOPTKEY   => 's|stretch-size=s',
	  GETOPTVAL   => \$mat_size,
	  REQUIRED    => 0,
	  DEFAULT     => $mat_size,
	  HIDDEN      => 0,
	  SMRY_DESC   => 'Stretch size.',
	  DETAIL_DESC => ('This is the length with which to look for ' .
			  'perfect matches between pairs of sequences.  ' .
			  'Only a series of matches of this length (or ' .
			  'multiples of this length) will be found.  The ' .
			  'script must be re-run to look for multiple ' .
			  'lengths.'));

my $global_play_min_def = 50;
my $global_play_def = 100;
my($global_play);
my($global_play_opt);      #How much the coordinates are allowed to slide apart
                          #Set to the difference in sequence length as a min
addOption(GETOPTKEY   => 'g|offset-max-global=i',
	  GETOPTVAL   => \$global_play_opt,
	  REQUIRED    => 0,
	  DEFAULT     => 'auto*',
	  HIDDEN      => 0,
	  SMRY_DESC   => 'Relative coordinate search range.',
	  DETAIL_DESC => ("Matches between sequences are only found if they " .
			  "occur within this relative coordinate range.  " .
			  "E.g. If this value is set to N, a match at " .
			  "coordinate 100 in sequence 1 will only be " .
			  "reported if it occurs (/starts) in sequence 2 " .
			  "between coordinates 100-N and 100+N." .
			  "However, matches are forced to be sequential.  " .
			  "In other words, once a match is located, the " .
			  "next match is restricted to everything to the " .
			  "right in each sequence.  Must be greater than or " .
			  "equal to --offset-max-local.  *If the difference " .
			  "between lengths of the sequences is greater than " .
			  "[$global_play_min_def], the default range is the " .
			  "lesser of either the sequence size difference or " .
			  "[$global_play_def]."));

my($local_play_opt);
my($local_play);          #How much the coordinates are allowed to slide apart
                          #from match to match.  The number of gaps must either
                          #be <= to this or the gap_ratio_max must be met.
addOption(GETOPTKEY   => 'l|offset-max-local=i',
	  GETOPTVAL   => \$local_play_opt,
	  REQUIRED    => 0,
	  DEFAULT     => 'auto*',
	  HIDDEN      => 0,
	  DETAIL_DESC => ("Each sequential match between sequences are only " .
			  "allowed to be found if they occur within this " .
			  "coordinate range relative to the positions of " .
			  "the previous match.  E.g. If this value is set " .
			  "to N and the previous match was at positions 100 " .
			  "and 200 for sequences 1 and 2 respectively, a " .
			  "match at coordinate 1000 in sequence 1 will only " .
			  "be reported if it occurs (/starts) in sequence 2 " .
			  "between coordinates 1100-N and 1100+N." .
			  "However, matches are forced to be sequential.  " .
			  "In other words, once a match is located, the " .
			  "next match is restricted to everything to the " .
			  "right in each sequence.  Must be less than or " .
			  "equal to --offset-max-global.  *The default " .
			  "value is the value of --offset-max-global."));

my $gap_ratio_max = 1;    #This is a local max (ratio of gaps allowed between 2
                          #matches)
addOption(GETOPTKEY   => 'r|gap-ratio-max-local=s',
	  GETOPTVAL   => \$gap_ratio_max,
	  REQUIRED    => 0,
	  DEFAULT     => $gap_ratio_max,
	  HIDDEN      => 0,
	  DETAIL_DESC => ('The maximum proportion of anticipated gaps to be ' .
			  'created between a pair of stretches (see -s).  ' .
			  'In other words, this is the difference in size ' .
			  'of the sequences between a pair of stretches ' .
			  'over the longer of the two.  A value between 0 ' .
			  'and 1 is recommended, but it can be greater than ' .
			  '1.  See -m for an alternative gap control.'));

my $global_gap_max = -1; #Max number of (internal) gaps allowed in the final
                          #alignment.  Note, this is a minimum number of gaps,
                          #assuming all the AAs in between the stretches align
                          #as tightly as possible.  -1 = off
addOption(GETOPTKEY   => 'm|gap-max=i',
	  GETOPTVAL   => \$global_gap_max,
	  REQUIRED    => 0,
	  DEFAULT     => 'infinity',
	  HIDDEN      => 0,
	  DETAIL_DESC => ('Stretch positions are restricted to a range (see ' .
			  '-g and -l), but many gaps can be created in the ' .
			  'resulting global alignment after submitting to ' .
			  'codonHomologizer.  While this maximum cannot ' .
			  'control preceisely how many gaps will be used in ' .
			  'the alignment, it can ensure that the sum total ' .
			  'of length differences between the stretches (see ' .
			  '-s) found is kept under this maximum.'));

my $time_limit = 60;  #Seconds
addOption(GETOPTKEY   => 't|time-limit=i',
	  GETOPTVAL   => \$time_limit,
	  REQUIRED    => 0,
	  DEFAULT     => $time_limit,
	  HIDDEN      => 0,
	  DETAIL_DESC => ('Report the best solution found for a pair of ' .
			  'sequences after searching for this number of ' .
			  'seconds.  Every possible pair will get this many ' .
			  'seconds to find the largest number of stretches ' .
			  'that meet all the requirements.  See -p to reset ' .
			  'this time restriction upon finding a better ' .
			  'solution.  Set to 0 to search all possible ' .
			  'solutions with no time limit.'));

my $improvement_reset = 1;
addOption(GETOPTKEY   => 'p|improvement-reset!',
	  GETOPTVAL   => \$improvement_reset,
	  REQUIRED    => 0,
	  DEFAULT     => 1,
	  HIDDEN      => 0,
	  DETAIL_DESC => ('As long as the solution keeps improving, reset ' .
			  'the --time-limit.  Only finding solutions with ' .
			  'more stretches (see -s) will reset the timer.  ' .
			  'Improved solutions with fewer alignment gaps, ' .
			  'but the same number of stretches will not reset ' .
			  'the timer.'));

my $seq_file_type =
  addInfileOption(GETOPTKEY   => 'i|aa-file=s',
		  REQUIRED    => 1,
		  PRIMARY     => 1,
		  DEFAULT     => undef,
		  SMRY_DESC   => 'Amino acid sequence file.',
		  FORMAT_DESC => << 'END_FORMAT'

The amino acid sequence file can be in fasta or fastq format.  There must be 2 or more sequences.  This is an unaligned sequence file.  Any alignment characters present such as gap characters ('-') or spaces will be removed.  It may be upper or lower case.

END_FORMAT
		 );

my $sequence_formats = ['auto','fasta','fastq'];
my $seq_format = $sequence_formats->[0];
addOption(GETOPTKEY   => 'f|aa-seq-format=s',
	  GETOPTVAL   => \$seq_format,
	  REQUIRED    => 0,
	  DEFAULT     => $seq_format,
	  HIDDEN      => 0,
	  DETAIL_DESC => ('Amino acid sequence file format.  Applies to ' .
			  'files submitted using -i.'),
	  ACCEPTS     => $sequence_formats);

my $extension = '.seg';
my $stretch_out_type =
  addOutfileSuffixOption(GETOPTKEY     => 'o|outfile-suffix=s',
			 FILETYPEID    => $seq_file_type,
			 GETOPTVAL     => \$extension,
			 REQUIRED      => 0,
			 PRIMARY       => 1,
			 DEFAULT       => $extension,
			 HIDDEN        => 0,
			 DETAIL_DESC   => ('The extension of the output tab-' .
					   'delimited coordinate file of ' .
					   'matching stretches (to be ' .
					   'supplied to codonHomologizer.pl ' .
					   'using the -d option).  An ' .
					   'outfile for every pair of ' .
					   'sequences in -i will be created ' .
					   'using a numeric ID per pair and ' .
					   'prepended to this suffix (e.g. ' .
					   ".1$extension, .2$extension,...)."),
			 FORMAT_DESC   => ('Tab-delimited text file.  The ' .
					   'first column is the sequence ' .
					   'ID, as parsed from the deflines ' .
					   'of -i.'),
			 COLLISIONMODE => 'error');

my $aa_suffix = '.faa';
my $aa_out_type =
  addOutfileSuffixOption(GETOPTKEY     => 'x|aa-pair-suffix=s',
			 FILETYPEID    => $seq_file_type,
			 GETOPTVAL     => \$aa_suffix,
			 REQUIRED      => 0,
			 PRIMARY       => 0,
			 DEFAULT       => $aa_suffix,
			 HIDDEN        => 0,
			 DETAIL_DESC   => << 'END_DETAIL'

Every pair of sequences used to find matching stretches (see -i and -s) generates a stretch coordinate file (see -o).  If there are more than 2 sequences in the supplied sequence file (see -i), a corresponding fasta file (per pair) will be generated with this extension in order to be supplied with the stretch coordinate file to codonHomologizer using codonHomologizer's -i option.  An outfile for every pair of sequences will be created using a numeric ID per pair and prepended to this suffix (e.g. .1$aa_suffix, .2$aa_suffix,...).  The pair ID will correspond to the pair ID of the --outfile-suffix.

END_DETAIL
			 ,
			 FORMAT_DESC   => ('Fasta format.'),
			 COLLISIONMODE => 'error');

my $expand_by = 0; #Coordinate pairs are expanded out by this number, e.g. 3-5
                   #is changed to 2-6.  Change to 0 to leave coords unchanged.
addOption(GETOPTKEY   => 'e|expand-coords-by=i',
	  GETOPTVAL   => \$expand_by,
	  REQUIRED    => 0,
	  DEFAULT     => $expand_by,
	  HIDDEN      => 0,
	  DETAIL_DESC => ('Pad the coordinates of the found stretches by ' .
			  'this amount.  This can help coerce the alignment ' .
			  'algorithm used in codonHomologizer to abutt non-' .
			  'matching amino acids with the stretch so that ' .
			  'selected codons could potentially lengthen the ' .
			  'stretch of homology by a nucleotide or 2.'));

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

if($mat_size < $mat_min)
  {
    error("Stretch size: [$mat_size] cannot be less than [$mat_min].",
	  {DETAIL => ('Setting a stretch size very small can lead to long ' .
		      'running times, so it is highly discouraged.')});
    quit(1);
  }

if(defined($global_play_opt) && defined($local_play_opt) &&
   $global_play_opt < $local_play_opt)
  {
    error("--offset-max-global [$global_play_opt] cannot be less than ",
	  "--offset-max-local [$local_play_opt].");
    quit(2);
  }
elsif(defined($global_play_opt))
  {$global_play = $global_play_opt}

if(defined($local_play_opt) && !defined($global_play_opt))
  {
    error("Setting --offset-max-local requires --offset-max-global to be ",
	  "set.");
    quit(6);
  }
elsif(defined($local_play_opt))
  {$local_play = $local_play_opt}

if($gap_ratio_max < 0 || $gap_ratio_max =~ /[^\d\.]/)
  {
    error("-r [$gap_ratio_max] must be a number and cannot be less than 0.");
    quit(3);
  }

if(defined($extension) && $extension eq '')
  {
    error('-o cannot be an empty string.');
    quit(4);
  }

if(defined($aa_suffix) && $aa_suffix eq '')
  {
    error('-x cannot be an empty string.');
    quit(5);
  }

my $global_max = 0; #Best solution thus far (per pair) - here so it's global
my $start_time = time(); #This is updated to current for each pair

#nextFileSet iterates sets of command-line-supplied files processed together
while(nextFileCombo())
  {
    my $aaFile = getInfile($seq_file_type);
    my $outStub = getOutfile($stretch_out_type);
    $outStub =~ s/\Q$extension\E$//;

    my $recs = [];
    openIn(*AA,$aaFile) || next;
    $recs = [getNextSeqRec(*AA,0,$aaFile)];
    closeIn(*AA);

    #Create a hash of all possible subsequences of size $mat_size
    my $hash = {};
    foreach my $rec (@$recs)
      {
	my $id = parseIDFromDef($rec->[0]);
	my $seq = $rec->[1];

	foreach my $pos (0..(length($seq) - $mat_size))
	  {
	    my $subseq = substr($seq,$pos,$mat_size);
	    push(@{$hash->{$id}->{$subseq}},$pos);
	  }
      }

    if(scalar(@$recs) < 2)
      {
	error("Too few sequences in file: [$aaFile]: [",scalar(@$recs),
	      "].  Skipping.");
	next;
      }

    my $pair_id = 0;
    foreach my $first (0..($#{$recs} - 1))
      {
	foreach my $second (($first + 1)..$#{$recs})
	  {
	    $start_time = time();

	    $pair_id++;
	    $global_max = 0; #Best solution thus far

	    my $len_diff = abs(length($recs->[$first]->[1]) -
			       length($recs->[$second]->[1]));
	    if(!defined($global_play_opt))
	      {
		if($len_diff < $global_play_min_def)
		  {$global_play = $global_play_def}
		else
		  {$global_play = $len_diff}
	      }
	    if($global_play < $len_diff)
	      {
		warning("--offset-max-global [$global_play] is less than the ",
			"difference in sequence lengths [$len_diff].");
	      }
	    if(!defined($local_play_opt))
	      {$local_play = $global_play}
	    if($global_play < $local_play)
	      {
		error("--offset-max-global [$global_play] is less than ",
		      "--offset-max-local [$local_play].");
		next;
	      }

	    my($max_num_stretches,$stretches,$junk1,$junk2,$abs_diffs_sum,
	       $searched) =
		 getMaxStretches(parseIDFromDef($recs->[$first]->[0]),
				 $recs->[$first]->[1],
				 parseIDFromDef($recs->[$second]->[0]),
				 $recs->[$second]->[1],
				 $hash,$mat_size,$global_play,$local_play,
				 $gap_ratio_max,$global_gap_max);

	    my $segfile = $outStub;
	    my $id1 = parseIDFromDef($recs->[$first]->[0]);
	    my $id2 = parseIDFromDef($recs->[$second]->[0]);
	    if($id1 =~ /^[A-Za-z0-9\-_\.]+$/ && $id2 =~ /^[A-Za-z0-9\-_\.]+$/)
	      {$segfile .= "_$id1\_vs_$id2"}
	    $segfile .= (scalar(@$recs) > 2 ? ".$pair_id" : '') . $extension;
	    openOut(*SEG,$segfile) || next;

	    print("#Stretch Size: $mat_size\n",
		  "#Offset Max Global: $global_play\n",
		  "#Offset Max Local: $local_play\n",
		  "#Gap Ratio Max Local: $gap_ratio_max\n",
		  "#Max Gaps Allowed: ",($global_gap_max < 0 ?
					 'ANY' : $global_gap_max),"\n",
		  "#Gaps Created: $abs_diffs_sum\n",
		  "#Max Stretches Found: $max_num_stretches\n",
		  "#Ratio of search space searched: $searched\n\n",
		  "#Stretches:\n#\t",
		  (parseIDFromDef($recs->[$first]->[0]),"\t",
		   parseIDFromDef($recs->[$second]->[0]),"\n#"),
		  join("\n#",map {join("\t",@$_)} @$stretches),"\n\n",
		  (scalar(@$stretches) ?
		   join("\t",
			(parseIDFromDef($recs->[$first]->[0]),
			 map {($_->[1] + 1 - $expand_by) . "\t" .
				($_->[1] + length($_->[0]) + $expand_by)}
			 @$stretches)) .
		   "\n" .
		   join("\t",
			(parseIDFromDef($recs->[$second]->[0]),
			 map {($_->[2] + 1 - $expand_by) . "\t" .
				($_->[2] + length($_->[0]) + $expand_by)}
			 @$stretches)) .
		   "\n" : ''));

	    closeOut(*SEG);

	    if(scalar(@$recs) != 2)
	      {
		my $aaoutfile = $outStub;
		if($id1 =~ /^[A-Za-z0-9\-_\.]+$/ &&
		   $id2 =~ /^[A-Za-z0-9\-_\.]+$/)
		  {$aaoutfile .= "_$id1\_vs_$id2"}
		$aaoutfile .= ".$pair_id$aa_suffix";
		openOut(*PAIR,$aaoutfile);
		my $def1 = $recs->[$first]->[0];
		my $def2 = $recs->[$second]->[0];
		$def1 =~ s/^\s*\@/>/;
		$def2 =~ s/^\s*\@/>/;
		print("$def1\n$recs->[$first]->[1]\n",
		      "$def2\n$recs->[$second]->[1]\n");
		closeOut(*PAIR);
	      }
	  }
      }
  }

sub getMaxStretches
  {
    my $id1            = $_[0];
    my $seq1           = $_[1];
    my $id2            = $_[2];
    my $seq2           = $_[3];
    my $hash           = $_[4];
    my $mat_size       = $_[5]  || 3;
    my $global_play    = $_[6]  || abs(length($seq2) - length($seq1));
    my $local_play     = defined($_[7]) ? $_[7] : $global_play;
    my $gap_ratio_max  = $_[8]  || 1;
    my $global_gap_max = $_[9]  || (length($seq2) > length($seq1) ?
				    length($seq2) : length($seq1));
    my $start1         = $_[10] || 0;
    my $start2         = $_[11] || 0;
    my $num_matches    = $_[12] || 0;
    my $abs_diffs_sum  = $_[13] || 0; #Sum of the local shifted distances (not
                                      #counting the first difference
    my $first          = $_[14] || 0;

    debug("getMaxStretches1(",
	  join(',',map {ref($_) eq '' ?
			  (length($_) <= 15 ? $_ : substr($_,0,12) . '...') :
			    ref($_)} @_),")\n");
    debug("getMaxStretches2(",
	  join(',',map {ref($_) eq '' ? (length($_) <= 15 ? $_ :
					 substr($_,0,12) . '...') : ref($_)}
	       ($id1,$seq1,$id2,$seq2,$hash,$mat_size,$global_play,
		$gap_ratio_max,$start1,$start2)),")\n");

    my $max_num_matches     = 0;
    my $matches             = [];
    my $contig              = 0;
    my $subseq2             = '';
    my $final_abs_diffs_sum = $abs_diffs_sum;
    my $pos1                = $start1;

    foreach($start1..(length($seq1) - $mat_size))
      {
	last if($time_limit && (time() - $start_time >= $time_limit));

	$pos1 = $_;
	my $subseq = substr($seq1,$pos1,$mat_size);
	if(exists($hash->{$id2}->{$subseq}))
	  {
	    foreach my $pos2 (grep {$_ >= $start2} @{$hash->{$id2}->{$subseq}})
	      {
		last if($time_limit && (time() - $start_time >= $time_limit));

		my $dist1       = $pos1 - $start1;
		my $dist2       = $pos2 - $start2;
		my $local_diff  = abs($dist2 - $dist1);
		my $longest     = ($dist1 > $dist2 ? $dist1 : $dist2);
		my $gap_ratio   = ($longest > 0 ? $local_diff / $longest : 0);
		my $global_diff = abs($pos1 - $pos2);
		my $gap_tally   = $abs_diffs_sum +
		  ($start1 == 0 || $start2 == 0 ? 0 : $local_diff);
		if($global_diff <= $global_play &&
		   ($local_diff <= $local_play ||
		    $gap_ratio <= $gap_ratio_max) &&
		   ($global_gap_max < 0 || $gap_tally <= $global_gap_max))
		  {
		    debug("getMaxStretches3(",
			  join(',',
			       map {ref($_) eq '' ?
				      (length($_) <= 15 ?
				       $_ : substr($_,0,12) . '...') :
					 ref($_)}
			       ($id1,$seq1,$id2,$seq2,$hash,$mat_size,
				$global_play,$gap_ratio_max,"$pos1+$mat_size",
				"$pos2+$mat_size")),")\n");
		    my($tmp_num_matches,$tmp_matches,$tmp_contig,$tmp_subseq,
		       $tmp_abs_diffs_sum) =
		      getMaxStretches($id1,$seq1,$id2,$seq2,$hash,$mat_size,
				      $global_play,$local_play,$gap_ratio_max,
				      $global_gap_max,$pos1 + $mat_size,
				      $pos2 + $mat_size,$num_matches + 1,
				      $gap_tally,
				      ($start1 == 0 ? $pos1 : $first));
		    if(#We got more matches...
		       $tmp_num_matches + 1 > $max_num_matches ||
		       #Or the number of matches is the same and the unevenness
		       #between all the matches is overall less
		       ($tmp_num_matches + 1 == $max_num_matches &&
		        $tmp_abs_diffs_sum < $final_abs_diffs_sum))
		      {
			$final_abs_diffs_sum = $tmp_abs_diffs_sum;
			$max_num_matches = $tmp_num_matches + 1;
			if($tmp_contig)
			  {
			    if(scalar(@$tmp_matches) == 1)
			      {@$matches =
				 (["$subseq$tmp_subseq",$pos1,$pos2])}
			    else
			      {
				my $mat = shift(@$tmp_matches);
				@$matches =
				  (["$subseq$tmp_subseq",$pos1,$pos2],
				   @$tmp_matches)
			      }
			  }
			else
			  {@$matches = ([$subseq,$pos1,$pos2],@$tmp_matches)}

			#If this match is a contiguous continuation of the last
			#match, mark it for concatenation upon return.
			if($start1 == $pos1 && $start2 == $pos2)
			  {
			    $contig = 1;
			    $subseq2 = $subseq;
			  }
			else
			  {
			    $contig = 0;
			    $subseq2 = '';
			  }
			if(($num_matches + $max_num_matches) > $global_max)
			  {
			    $global_max = $num_matches + $max_num_matches;
			    if($improvement_reset)
			      {
				verboseOverMe("IMPROVEMENT: Resetting timer")
				  if($time_limit);
				$start_time = time();
			      }
			  }
			verboseOverMe("Progress: ",int($first /
						       (length($seq1) -
							$mat_size) * 100),
				      "\% done.  Position: $first out of ",
				      (length($seq1) - $mat_size),
				      "  Max Stretches from positions ",
				      "$pos1/$pos2: $global_max",
				      " Recursion Depth: $num_matches ",
				      ($time_limit ? "Giving up in: " .
				       ($time_limit - (time() - $start_time)) :
				       ''),
				      "   \r");
		      }
		  }
	      }
	  }
      }
    return($max_num_matches,$matches,$contig,$subseq2,$final_abs_diffs_sum,
	   $pos1 / (length($seq1) - $mat_size));
  }

sub parseIDFromDef
  {
    my $def = $_[0];

    $def =~ s/^[>\@\+]\s*//;
    $def =~ s/^[>\@\+]\s*//;
    $def =~ s/ .*//;

    return($def);
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
		  {warning(scalar(@errors)," bad nucleotide characters were ",
			   "filtered out of your sequence: [",join('',@errors),
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
