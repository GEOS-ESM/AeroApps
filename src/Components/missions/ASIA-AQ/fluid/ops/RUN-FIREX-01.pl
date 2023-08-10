#!/usr/bin/perl -w
#
#########################################################################
#                                                                       # 
# RUN-FIREX-01.pl                                                       #
#                                                                       #
# Log:                                                                  #
# Mark Solomon SSAI 2019/08/02                                          #
#                                                                       #
#########################################################################

use integer;

# This module gives the routine access to the environment

use Env;

# This module contains the getopts() subroutine.
use Getopt::Long;

# Module contains "mkpath"
use File::Path;

# This module locates the full path name to the location of this file.  Variable
# $FindBin::Bin will contain that value.

use FindBin;

# This module contains the dirname() subroutine.

use File::Basename;
use File::Path;
use File::Copy;

my $env = "ops";
my $tab_argv;
my $fl_name;
my $sched_dir;
my $sched_sts_fl;
my $error;
my $die_away;
my $date = `date`;

# **********************************************
# section ()
# **********************************************
BEGIN {

    print `date`;
#$DB::single = 1;

# Keep track of errors within BEGIN block.

    $die_away = 0;

# Get options and arguments

   GetOptions ('env=s',\$env,
               'E=s',\$opt_E,
               't=s',\$opt_t,
               'O=s',\$opt_O,
               'R=s',\$opt_R,
               'I=s',\$opt_I,
               'L=s',\$opt_L,
               'F=s',\$opt_F,
               'sched_cnfg=s',\$sched_cnfg,
               'sched_id=s',\$sched_id,
               'sched_synp=s',\$sched_synp,
               'sched_c_dt=s',\$sched_c_dt,
               'sched_dir=s',\$sched_dir,
               'force',\$opt_force,
               'sched_sts_fl=s',\$sched_sts_fl,
               'sched_hs=s',\$sched_hs );


   if ( defined( $sched_hs ) ) {
       print "sched_hs: $sched_hs.\n";
   } 

   if ( defined( $sched_id ) )
   {
      $tab_argv = "$sched_cnfg, $sched_id, $sched_synp, $sched_c_dt";
      $fl_name = "RUN-FIREX-01.pl";
   }

# The pre-processing configuration file.

    if ( defined( $opt_E ) ) {
        $PREP_CONFIG_FILE = $opt_E;
        print "File: $PREP_CONFIG_FILE.\n";
    } else {
        print "File missing.\n";
        $error = "File missing.\n";
        $die_away = 1;
    }

# The Run_Config file to use.

    if ( defined( $opt_R ) ) {
        $RUN_CONFIG_FILE = $opt_R;
    } else {
        $RUN_CONFIG_FILE = "DEFAULT";
    }

# ID - default = "dao_ops".

    if ( defined( $opt_I ) ) {
        $proc_id = $opt_I;
    } else {
        $proc_id = "dao_ops";
    }

# The lag date command.

    if ( defined( $opt_L ) ) {
        $LAG_DATE = $opt_L;
    }
    else {
        $LAG_DATE = "DEFAULT";

    }

    if ( defined( $opt_t ) ) {
        $hour = $opt_t;
        print "Hour: $hour.\n";
    } else {
        $hour = 0;
        print "Hour set to: $hour.\n";
    }

    if ( defined( $opt_O ) ) {
        $output = $opt_O;
    } else {
        $output = "/discover/nobackup/dao_ops/";
        print "output directory missing set to: $output.\n";
    }

# Path to directory containing other GEOS DAS programs.
# Directory $GEOSDAS_PATH/bin will be searched for these
# programs.

    if ( defined( $opt_P ) ) {
        $GEOSDAS_PATH = $opt_P;
    } else {
        $GEOSDAS_PATH = "/home/dao_ops/GEOSadas-5_21/GEOSadas/Linux/";
    }

# The value for debug information

#


# The log file to use.

    if ( defined( $opt_F ) ) {
        $logfile = $opt_F;
    }
    else {
        $logfile = "DEFAULT";
    }


# Location for output listings
    $listing_file = 0;


    if ( $ARGV[0] eq 'help' || $ARGV[0] eq 'HELP' || $#ARGV < 0 ) {
        print STDERR <<'ENDOFHELP';
        
Usage:
        
        RUN-FIREX-01 [-E Prep_Config] [-P GEOSDAS_Path] [-O output_location] 
            [-D SYN_date] [-I Proc_ID] [-R Run_Config] Prep_ID
            
            Normal options and arguments:
            
            -E Prep_Config
            The full path to the preprocessing configuration file.  This file contains
            parameters needed by the preprocessing control programs. If not given, a
            file named $HOME/$prep_ID/Prep_Config is used.  get_map06_sst exits with an
            error if neither of these files exist. See below...

            -F Log_File
            The full path to the ERROR_LOG file ( Default is \~$user/GEOS_EVENT_LOG)

            -I Proc_ID
            The processing identifier. The default is "ops".
            
            -R Run_Config
            The full path and name of the Run_Config file to use. If not given,
            the default will be GEOSDAS_PATH/bin/Run_Config
            
            -D SYN_date
            The date for which the model is running for. FORMAT CCYYMMDD
            
            -P GEOSDAS_Path
            Path to directory containing other GEOS DAS programs.  The path is
            $GEOSDAS_PATH, where $GEOSDAS_PATH/bin is the directory containing these
            programs.  If -P GEOSDAS_Path is given, then other required programs not
            found in the directory where this program resides will be obtained from
            subdirectories in GEOSDAS_Path - the subdirectory structure is assumed
            to be the same as the operational subdirectory structure.  The default is
            to use the path to the subdirectory containing this program, which is what
            should be used in the operational environment.
            
            -O output_location
            If this option is specified, output listings (both STDERR and STDOUT) will be
            placed in the directory "output_location."

            -p    Use persisted anomaly instead of anomaly persistence.  This is done for sea-ice only.
            
            Prep_ID
            Identification tag for this run.
            
            
            Prep_Config setttings
            ---------------------
            
            The program will consult a Prep_Config file for configurable parameters.
            Many of these parameters can use the GrADS-style date and time templates,
            i.e., %y4, %y2, %m2, %d2.  The Prep_Config file is assumed to reside in
            ~/Prep_ID/Prep_Config, but any other path can be specified by the -E option.            
            Input files:
            
ENDOFHELP

            $die_away = 1;
        exit 0;
    }
# This module gives the routine access to the environment

    use Env;

# This module locates the full path name to the location of this file.  Variable
# $FindBin::Bin will contain that value.

    use FindBin;

# This module contains the dirname() subroutine.

    use File::Basename;

# If default GEOS DAS path, set path to parent directory of directory where this
# script resides.

    if ( $GEOSDAS_PATH eq "DEFAULT" ) {
        $GEOSDAS_PATH = dirname( $FindBin::Bin );
    }

# Set name of the bin directory to search for other programs needed by this one.

    $BIN_DIR = "$GEOSDAS_PATH/bin";
    $CORE_BIN_DIR = "$GEOSDAS_PATH/geosdas/Core/bin";
    $BIN_OPS = "$GEOSDAS_PATH/bin_ops";

# Get the name of the directory where this script resides.  If it is different
# than BIN_DIR, then this directory will also be included in the list of
# directories to search for modules and programs.

    $PROGRAM_PATH = $FindBin::Bin;

# Now allow use of any modules in the bin directory, and (if not the same) the
# directory where this program resides.  (The search order is set so that
# the program's directory is searched first, then the bin directory.)

    if ( $PROGRAM_PATH ne $BIN_DIR ) {
        @SEARCH_PATH = ( $PROGRAM_PATH, $BIN_DIR, $CORE_BIN_DIR, $BIN_OPS );
    } else {
        @SEARCH_PATH = ( $BIN_DIR, $CORE_BIN_DIR, $BIN_OPS );
    }
} # END 


# =================================================================
#
#
# =================================================================

# Any reason to exit found during the BEGIN block?

if ( $die_away == 1 ) {
    exit 1;
}

# Include the directories to be searched for required modules.

use lib ( @SEARCH_PATH );

# Set the path to be searched for required programs.

$ENV{'PATH'} = join( ':', @SEARCH_PATH, $ENV{'PATH'} );

# This module contains the extract_config() subroutine.

use Extract_config;

# Archive utilities: gen_archive

use Arch_utils;

# This module contains the z_time(), dec_time() and date8() subroutines.

use Manipulate_time;

# This module contains the mkpath() subroutine.

use File::Path;

use Conv_utils;

use Net::FTP;
use Remote_utils;

# Error logging utilities.
use Err_Log;

# Record FAILED to schedule status file.
use Recd_State;

# ID for the preprocessing run.
print "argv0 = $ARGV[0] , argv1 = $ARGV[1] , argv2 = $ARGV[2] \n";
$prep_ID = $ARGV[0];
print "prep_ID=$prep_ID\n";
$process_date = $ARGV[2];
print "process_date=$process_date\n";
$SYN_DATE = $process_date;



# Use Prep_Config file under the preprocessing run's directory in the user's home directory
# as the default.

if ( "$PREP_CONFIG_FILE" eq "DEFAULT" ) {
    $PREP_CONFIG_FILE = "$ENV{'HOME'}/$prep_ID/Prep_Config";
}

# Does the Prep_Config file exist?  If not, fatal_error.
  if ( ! -e "$PREP_CONFIG_FILE" ) {
    fatal_error("error $PREP_CONFIG_FILE not found.");
  }

if ( "$RUN_CONFIG_FILE" eq "DEFAULT" ) {
    $RUN_CONFIG_FILE = "$BIN_DIR/Run_Config";
}
#********#
# Set up #
#********#

($nextDate, $current_time) = inc_time ($process_date, 0, 1, 0);
print "nextDate=$nextDate\n";

# Defining an array which I will fill with found file names and fix some variables.


# Define job_id
$job_id = ${SYN_DATE};

# Extract working directory, staging, and working locations and the archive flag.

( $PROGRAM_NAME = extract_config( "PROGRAM_NAME", $PREP_CONFIG_FILE, "NONE" ) ) ne "NONE"
    or fatal_error ("RUN-FIREX-01 ERROR - can not set PROGRAM_NAME configuration value");

( $FIREX_BIN_DIR = extract_config( "FIREX_BIN_DIR", $PREP_CONFIG_FILE, "NONE" ) ) ne "NONE"
    or fatal_error ("RUN-FIREX-01 ERROR - can not set FIREX_BIN_DIR configuration value");

( $FIREX_YML = extract_config( "FIREX_YML", $PREP_CONFIG_FILE, "NONE" ) ) ne "NONE"
    or fatal_error ("RUN-FIREX-01 ERROR - can not set FIREX_YML configuration value");

    if ( defined( $opt_O ) ) {
        $listing_file = 1;
        if ( ! -d "$opt_O" ) {
            mkpath( "$opt_O" ) or fatal_error ("Can't make $opt_O: $!");
        }
        if ( -w "$opt_O" ) {
            $listing = "$opt_O/$PROGRAM_NAME.$$.listing";
            print "Standard output redirecting to $listing\n";
            open (STDOUT, ">$listing");
            open (STDERR, ">&" . STDOUT);

        } else {
            print "$0: WARNING: $listing is not writable for listing.\n";
        }
    }

# Write start message to Event Log

err_log (0, "RUN-FIREX-01", "$job_id", "$prep_ID", "-1",
         {'log_name' => "$logfile",
          'err_desc' => "$prep_ID has started - Standard output redirecting to $listing"});
$ENV{'FVROOT'} = "/home/dao_ops/jardizzo/FLUID/firex-aq/ops";
my $GITBIN = "/usr/local/other/git/1.7.3.4_GNU/bin";
$ENV{'GITBIN'} = "$GITBIN";
my $MYPYTHON = "/discover/nobackup/projects/gmao/share/dasilva/Python/epd-7.3-2-rh5-x86_64/bin";
$ENV{'MYPYTHON'} = "$MYPYTHON";
$ENV{'PATH'} = join( ':', "${FIREX_BIN_DIR}/", $ENV{'PATH'} );
$ENV{'PATH'} = join( ':', "/home/dao_ops/jardizzo/FLUID/firex-aq/utils", $ENV{'PATH'} );
$ENV{'PATH'} = join( ':', "/home/dao_ops/jardizzo/FLUID/firex-aq/bin", $ENV{'PATH'} );
$ENV{'PATH'} = join( ':', "$SHARE/dasilva/opengrads/Contents", $ENV{'PATH'} );
$ENV{'PATH'} = join( ':', "$GITBIN", $ENV{'PATH'} );
$ENV{'PATH'} = join( ':', "$MYPYTHON", $ENV{'PATH'} );
$ENV{'PATH'} = join( ':', "~adasilva/bin", $ENV{'PATH'} );
$ENV{'PATH'} = join( ':', "$SHARE/dasilva/bin", $ENV{'PATH'} );
$ENV{'PATH'} = join( ':', "/usr/local/bin", $ENV{'PATH'} );
$ENV{'PYTHONPATH'} = "/home/adasilva/workspace/QFED/Linux/lib/Python:/home/dao_ops/jardizzo/FLUID/firex-aq/lib";
$ENV{'GAVERSION'} = "2.1.0.oga.1";
do "/usr/share/modules/init/perl";
module ("purge");

print "$FIREX_BIN_DIR \n";
do "${FIREX_BIN_DIR}/pyg_modules_perl_wrapper";

foreach $var (sort(keys(%ENV))) {
    $val = $ENV{$var};
    $val =~ s|\n|\\n|g;
    $val =~ s|"|\\"|g;
    print "${var}=\"${val}\"\n";
}


$cmd="${FIREX_BIN_DIR}/firex_ops.py ${process_date} 0 ${FIREX_BIN_DIR}/models.rc ${FIREX_BIN_DIR}/models/${FIREX_YML}";
print "$cmd \n";
$result = system($cmd);
print "result = $result \n";
if( $result ) {
    fatal_error ("RUN-FIREX-01 ERROR - can not run /usr/bin/sh $PROGRAM_NAME".
                 " $process_date ...  ");
}



print `date`;

err_log (0, "RUN-FIREX-01", "$job_id","$prep_ID","-1",
         {'log_name' => "$logfile",
          'err_desc' => "RUN-FIREX-01: successfully finished "});
if ( $listing_file ) {
   close (STDERR);
   close (STDOUT);
   $out_listing = "$opt_O/${PROGRAM_NAME}.${process_date}.listing";
   system ("mv $listing $out_listing");
}
# Report success to D_BOSS Scheduler if run by scheduler
if ( defined( $sched_id ) ){ recd_state( $fl_name, "COMPLETE", $tab_argv, $sched_dir, $sched_sts_fl );}


$to = 'mark.s.solomon@nasa.gov,joseph.v.ardizzone@nasa.gov';
$subject = "${PROGRAM_NAME} for $process_date complete";
$message = $subject;

open(MAIL, "|/usr/lib/sendmail -t");

# Email Header
print MAIL "To: $to\n";
print MAIL "Subject: $subject\n\n";
# Email Body
print MAIL $message;

close(MAIL);
print "Email Sent Successfully\n";
#----------------------------------------------------------------------------
sub fatal_error{
    my ($message) = @_;
# Report failure to Error Log
    err_log (5, "RUN-FIREX-O1", "$prep_ID", "$env", "-1",
             {'err_desc' => "FATAL ERROR: $message"});

# Report failure to D_BOSS Scheduler if run by scheduler
    if ( defined( $sched_id ) ) {
        recd_state($fl_name, FAILED, $tab_argv, $sched_dir, $sched_sts_fl);
    }
    die "FATAL ERROR: $message \n";
}
