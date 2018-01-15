#Library of methods used by codonHomologizer.pl and stretchMixer.pl

#Creates a hash from a tab delimited file like this:
#G	GGG	0.11	Gly
#G	GGA	0.21	Gly
#where:
#hash->{G}->{GGG} = {SCORE => 0.11,DATA => ['Gly']}
#hash->{G}->{GGG} = {SCORE => 0.21,DATA => ['Gly']}
#Globals used: $no_usage
sub readUsageFile
  {
    my $file     = $_[0];
    my $no_usage = $_[1];
    my $usg      = {};

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

	if(scalar(@data) < (!defined($no_usage) || $no_usage == 0 ? 3 : 2))
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
	elsif(!$no_usage && $data[2] !~ /\d/)
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

	$usg->{uc($data[0])}->{uc($data[1])} = {SCORE => ($no_usage ?
							  1.0 : $data[2]),
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
    my $usgobj           = $_[0];
    my $weighting_method = defined($_[1]) ? $_[1] : 'crossover1';

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

    verbose({LEVEL => 3},matrixToString($matrix));

    return($matrix);
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

#Returns a hash with keys SCORE and PAIRS.  The value for SCORE is a number.
#The value for PAIRS is a hash keyed on the category of flex_score keys that
#correspond to max_matches
sub getCrossover1ScoreHash
  {
    my $lesser           = $_[0];  #Lesser  (in alphabetical order) amino acid
    my $greater          = $_[1];  #Greater (in alphabetical order) amino acid
    my $matches          = $_[2];  #Maximum number of matches possible in cdns
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

    my $codons1  = [getCodons($lesser,$usgobj)];  #Assumed to be sorted by
    my $codons2  = [getCodons($greater,$usgobj)]; #descending usage value

    debug("Best homology among codons: [$matches].") if($testcase);
    debug("codons for $lesser: [",join(',',@$codons1),"].") if($testcase);
    debug("codons for $greater: [",join(',',@$codons2),"].") if($testcase);

    #This hash will track whether a codon pair with a matching pattern defined
    #by the flex scores exists or not, e.g. ATT/GTT is an MR match
    my $best_flex_min_usage  = {};
    my $best_pair            = [];
    my $best_min_usage_score = 0;

    foreach my $codon1 (@$codons1)
      {
	my $min_usage_score = getUsageScore($lesser,$codon1,$usgobj);

	#We're only going to look at codon pairs that have the maximum number
	#of matches possible given the 2 sets of codons (which was pre-
	#determined and supplied to this sub
	foreach my $codon2 (grep {getNumCodonMatches($codon1,$_) == $matches}
			    @$codons2)
	  {
	    #For each pair of codons possible between the sets defined by the 2
	    #amino acids, the score is influenced by the lesser of the 2 codon
	    #usages.  We want to select the pair with the largest lesser usage
	    #value.  This is selected from among the pairs that have the same
	    #flex code and the flex codes considered are those with the
	    #supplied number of matches.

	    #Determine the lesser of the 2 usage scores for this codon pair
	    if(getUsageScore($greater,$codon2,$usgobj) < $min_usage_score)
	      {$min_usage_score = getUsageScore($greater,$codon2,$usgobj)}

	    #Save the greater of the min usage scores among all pairs with the
	    #same flex code
	    my $flex_code = getCrossover1FlexCode($codon1,$codon2);
	    if(!exists($best_flex_min_usage->{$flex_code}) ||
	       $min_usage_score > $best_flex_min_usage->{$flex_code})
	      {$best_flex_min_usage->{$flex_code} = $min_usage_score}

	    #Save the greater of the min usage scores among all pairs
	    if($min_usage_score > $best_min_usage_score)
	      {
		$best_pair = [$codon1,$codon2];
		$best_min_usage_score = $min_usage_score;
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

sub matrixToString
  {
    my $matrix = $_[0];

    #Include the protein weight matrix
    my $string = "Alignment Weight Matrix (not used in alignment if -w was " .
      "supplied):\n\n" . pwMatrixToString($matrix) . "\n";

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

sub parseIDFromDef
  {
    my $def = $_[0];

    $def =~ s/^[>\@\+]\s*//;
    $def =~ s/^[>\@\+]\s*//;
    $def =~ s/ .*//;

    return($def);
  }

#Uses global variables: lastfiletype & seq_format
#Edited to take $seq_format instead of global 11/17/2017 -Rob
sub getNextSeqRec
  {
    #The first 2 parameters are passed along to the sub that reads the specific
    #file format detected (input handle and the no_format flag).  The third
    #parameter used by this sub is only for constructing errors and warnings
    #that mention the affected file name
    my $input_file = $_[2];
    my $filetype   = $_[3];

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
	    debug({LEVEL => 3},"System output from: [",
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

		debug({LEVEL => 3},"System output from: [",
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

1;
