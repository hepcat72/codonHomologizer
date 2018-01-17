#Makefile.PL version 1.1

use 5.010001;
use ExtUtils::MakeMaker;
my @scripts = qw(codonHomologizer.pl stretchAlign.pl stretchMixer.pl);

if(incompatibleMuscle(getExe('muscle')))
  {print STDERR ("WARNING: Muscle appears to be missing or potentially ",
                 "incompatible.  Installation will proceed, but in order to ",
                 "run filterIndels.pl without errors, you must install a ",
                 "compatible version of muscle and put it in your path.  ",
                 "Everything else will work fine.\n")}


WriteMakefile
  (
   NAME => 'codonHomologizer',
   ($[ >= 5.005 ?  #Add these new keywords supported since 5.005
    (AUTHOR        => 'Robert W. Leach (rleach@princeton.edu)',
     ABSTRACT_FROM => 'src/codonHomologizer.pl') : ()),
   VERSION_FROM    => 'src/codonHomologizer.pl', # finds $VERSION
   PREREQ_PM       => {'Getopt::Long' => 2.38,
		       'File::Glob'   => 1.17,
		       'strict'       => 1.07,
		       'warnings'     => 1.13},
   EXE_FILES       => [ map { "src/$_" } @scripts ]
  );




sub incompatibleMuscle
  {
    my $muscle = $_[0];

    if(!defined($muscle) || $muscle eq '' || !-e $muscle || !-x $muscle)
      {
        print STDERR ("ERROR: The muscle executable [$muscle] appears to ",
                      "either not be in your path, not exist, not have ",
                      "execute permissions, or you have not created a ",
                      "symbolic link named 'muscle' to the full name of the ",
                      "executable with version number.  If you have not ",
                      "installed muscle, you can find it here: ",
                      "http://www.drive5.com/muscle/downloads.htm\n");
        return(1);
      }

    my $version = `$muscle -version`;
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
        print STDERR ("WARNING: This script was tested with Muscle version ",
                      "3.8.31.  Your version appears to be [$vnum], thus it ",
                      "may not work properly.\n") unless($ok);
      }
    else
      {print STDERR ("WARNING: This script was tested with Muscle version ",
		     "3.8.31.  It may not work properly with the version you ",
                     "are using.\n")}

    return(0);
  }

sub getExe
  {
    my $command  = $_[0];
    my $sent_exe = $command;
    $sent_exe    =~ s/ .*//;
    my $exe      = '';

    if(eval("use File::Which;1;") ||
       eval("use local::lib;use File::Which;1;"))
      {
        $exe = which($sent_exe);
        if((!defined($exe) || $exe eq '') && -e $sent_exe && -x $sent_exe)
          {$exe = $sent_exe}
        elsif(!defined($exe))
          {$exe = ''}
      }
    else
      {
        print STDERR ("File::Which not found, switching to backup method.\n");
        $exe = `which $sent_exe`;
        chomp($exe);
        if($exe =~ /which: Command not found./ || $exe !~ /\S/)
          {
            print STDERR ("ERROR: System command 'which' does not appear to ",
			  "exist.  Please install the perl module ",
                          "File::Which.\n");
            $exe = '';
          }
        elsif($exe =~ /not found/i)
          {$exe = ''}
      }

    return($exe);
  }
