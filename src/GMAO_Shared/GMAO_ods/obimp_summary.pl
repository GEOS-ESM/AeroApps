#!/usr/bin/env perl

#-----------------------------------------------------------------------------
# !DESCRIPTION:
#
#    Summarize observation impact by running odsstats over files in given
#    location.
#
# !REVISION HISTORY:
#
#  16Apr2013  Todling   Initial code
#  24Feb2014  Todling   Expand capability to provide instrument (or other) summary
#
#-----------------------------------------------------------------------------

use Env;                 # make env vars readily available
use FindBin;             # so we can find where this script resides
use File::Basename;      # for basename(), dirname()
use File::Copy "cp";     # for cp()
use Shell qw(rm);        # rm commands
use Getopt::Long;        # command line options

# look for perl packages in the following locations
#--------------------------------------------------
use lib ( "$FindBin::Bin", "$FVROOT/bin", "$ESMADIR/$ARCH/bin" );

# FVROOT is where the binaries have been installed
# ------------------------------------------------
$fvroot  = dirname($FindBin::Bin);
$fvroot  =~ s|/u/.realmounts/share|/share|;   # for portability across
                                              # NAS machines
# Command line options
# --------------------
GetOptions( "h", "dir=s", "type=s", "rc=s", "ktsummary", "o=s" );

usage() if $opt_h;

$user = getlogin();

$init_status=init();

if (! $init_status ) {

   imp_summary();

} else {
   $rc = 1;
}

exit ($rc);

sub imp_summary {

   chdir("$filesdir");

   @files = glob("*${obimtyp}*${nymdhhz}.ods");
   if ( $kt_summary ) {

      @obscount = (0,0,0,0,0,0,0,0,0);
      @obsimp   = (0,0,0,0,0,0,0,0,0);
      @typs     = glob("spr temp uv hum spd pcp oz gps rad");
      foreach $fn ( @files ) {

        $cmd = "$fvroot/bin/odsstats -rc $rcfile -verbose $fn > /dev/null";
        print "$cmd \n";
        $rc = system($cmd);
        if (-e "odsstats_sum.txt") {
            $ic = 0;
            foreach $typ ( @typs ) {
               $vars = ( `grep $typ odsstats_sum.txt | cut -c40-` );
               $vars =~ s/^\s*//;
               @vars = split(/\s+/,$vars);
               $obscount[$ic] = $obscount[$ic] + $vars[0];
               $obsimp[$ic]   = $obsimp[$ic]   + $vars[1];
               $ic = $ic + 1;
            } # foreach typ
           } # sum exists
           rm("odsstats_sum.txt");
      } # imp files
   
      $ic = 0;
      print " obs-var     count      impact \n"; 
      foreach $typ ( @typs ) {
           printf " %8s %8d %11.4e \n", $typ, $obscount[$ic], $obsimp[$ic];
           $ic = $ic + 1;
      }

   } else { # not kt-summary
      $cmd = "$fvroot/bin/odsstats -rc $rcfile @files  > /dev/null";
      print "$cmd \n";
      $rc = system($cmd);
      if ( -e "odsstats_all.txt" ) {
         print "mv odsstats_all.txt $outfile \n";
         cp("odsstats_all.txt","$outfile");
         rm("odsstats_all.txt");
         rm("odsstats_sum.txt");
         rm("odsstats_numneg.txt");
      }
   } # kt-summary

} # end imp_summary

#=======================================================================
## name - init
##=======================================================================
sub init {


  if ( $#ARGV < 1 ) {
       print STDERR "missing nymd, nhms and/or expid; see usage";
       usage();
  } else {              # required command lile args
       $nymd    = $ARGV[0];
       $nhms    = sprintf("%6.6d",$ARGV[1]);
       $yyyy    = substr($nymd,0,4);
       $mm      = substr($nymd,4,2);
       $dd      = substr($nymd,6,2);
       $hh      = substr($nhms,0,2);
       $nymdhhz = "${nymd}_${hh}z";
  }

  $kt_summary = 0;
  if ( $opt_ktsummary ) {
      $kt_summary = 1;
  }

  if( $opt_rc ) { # rc file required by odsstats
      $rcfile = $opt_rc;
  } else {
      $kt_summary = 1;
      $rcfile = "$fvroot/etc/odsstats_ktonly.rc";
  }

  if( $opt_dir ) {  # location of input files
      $filesdir = $opt_dir;
  } else {
      $filesdir = "./";
  }

  $outfile = "obimp_summary.txt";
  if( $opt_o ) { $outfile = $opt_o };
  if( -e "$outfile" ) { rm("$outfile") };

  if( $opt_type ) {  # location of input files
      $obimtyp = $opt_type;
  } else {
      $obimtyp = "imp0hr";
  }

  $rc = 0;

} # end init

#=======================================================================
# name - usage
#=======================================================================
sub usage {

   print <<"EOF";

NAME
     obimp_summary - Summarizes observation impacts
          
SYNOPSIS

     obimp_summary [...options...]  nymd nhms
          
DESCRIPTION

     The following parameter are required 

     nymd     Year-month-day, e.g., 19990901  for 01 Sept 1999 
     nhms     Hour-minutes-seconds, e.g., 120000

OPTIONS
 
 -dir          full path of input files locations 
                (default: .)

 -h            prints this usage notice
 
 -rc RCFILE    full path of RC file for odsstats
                 (default: $fvroot/etc/odsstats_ktonly.rc)

 -type         allows choosing specific diag file, e.g., -type imp3_txe
                 (default: imp0hr)

ENVIRONMENT

AUTHOR
      R. Todling (ricardo.todling\@nasa.gov), NASA/GSFC/GMAO

EOF

  exit(1)

}
