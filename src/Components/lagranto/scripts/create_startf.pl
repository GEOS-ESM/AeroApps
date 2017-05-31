#!/usr/bin/perl

# --------------------------------------------------------
# Separate different commands out (delimiter @)
# --------------------------------------------------------

# Get input command and remove all spaces
$cmd = $ARGV[0];
$_   = $cmd;s/\s+//g;$cmd=$_;
$len = length $cmd;

# Split the command string according to logical operators
$nline = 0;
while ( $len > 0 )
{
  # Get the length of the (remaining) command string
  $len = length $cmd; 

  # Get the position of the next command separator
  $at  = index($cmd,'@');
  if ( $at == -1 )
     { $next = $len+1; }
  elsif ( $at >= 0 )
     { $next = $at; }

  # A logical operator is not allowed to be at position 0
  if ( $next == 0) 
     {    
      die('Invalid expression... Check @ separator position ');
     }   

  # Extract the next substring
  $sub = substr($cmd,0,$next);
  $cmd = substr($cmd,$next+1,$len-$next-1);
  $len = length $cmd; 
 
  # Save the command in a new line
  @field[$nline] = $sub;
  $nline = $nline + 1;

}

# Number of lines ist too large by 1
$nline = $nline - 1;

# Set some defaults
if ( $nline == 2 ) 
   {
    $nline    = 3;
    @field[3] = 'nil';
   }
if ( $nline == 1 ) 
   {
    $nline    = 3;
    @field[2] = 'hPa';
    @field[3] = 'nil';
   }
if ( $nline == 0 ) 
   {
    $nline    = 3;
    @field[1] = 'nil()';
    @field[2] = 'hPa';
    @field[3] = 'nil';
   }

# --------------------------------------------------------
# Handle each command line separately
# --------------------------------------------------------

# ----- Horizontal grid ----------------------------------

# Now looking at the horizontal grid specifier
$entry = @field[0];

# Extract the command and the parameter list
$left  = index($entry,'(');
$right = index($entry,')');
if ( ($left != -1) && ($right != -1) )
   {
    $cmd  = substr($entry,0,$left);
    $list = substr($entry,$left+1,$right-$left-1);
    $len  = length $list; 
   }
else
   {    
    die('Invalid expression... Check horizontal [] separator position ');
   }   

# Now building the parametr list 
$len  = length $list; 
$npar = 0;
while ( $len > 0 ) 
   {
   $next = index($list,',');
   if ( $next != -1 )
      {
       @param[$npar] = substr($list,0,$next);
       $list         = substr($list,$next+1,$len-$next-1);
       $len          = length $list; 
       $npar         = $npar + 1;
      }
   else
      {
       @param[$npar] = substr($list,0,$len);
       $len          = 0;
       $npar         = $npar + 1;
      }
   }

# Check for syntax (needed number of parameters)
if ( ($cmd eq "file")      && ($npar != 1) ) 
    {  die('Invalid horizontal mode[file]... Check number of parameters ');     }
if ( ($cmd eq "line")      && ($npar != 5) ) 
    {  die('Invalid horizontal mode[line]... Check number of parameters ');     }  
if ( ($cmd eq "box.eqd")   && ($npar != 5) ) 
    {  die('Invalid horizontal mode[box.eqd]... Check number of parameters ');  }  
if ( ($cmd eq "box.grid")  && ($npar != 4) ) 
    {  die('Invalid horizontal mode[box.grid]... Check number of parameters '); }  
if ( ($cmd eq "point")     && ($npar != 2) ) 
    {  die('Invalid horizontal mode[point]... Check number of parameters ');    }  
if ( ($cmd eq "shift")     && ($npar != 4) ) 
    {  die('Invalid horizontal mode[shift]... Check number of parameters ');    }  
if ( ($cmd eq "poly.eqd")  && ($npar != 2) ) 
    {  die('Invalid horizontal mode[poly.eqd]... Check number of parameters '); }  
if ( ($cmd eq "poly.grid") && ($npar != 1) ) 
    {  die('Invalid horizontal mode[poly.grid]... Check number of parameters ');}  
if ( ($cmd eq "mask.grid") && ($npar != 2) ) 
    {  die('Invalid horizontal mode[mask.grid]... Check number of parameters ');}  
if ( ($cmd eq "mask.eqd") && ($npar != 3) ) 
    {  die('Invalid horizontal mode[mask.eqd]... Check number of parameters ');}  

# Write parameters
print "\"$cmd\"\n";
print "@param\n";

# ----- Vertical grid ----------------------------------------

# Now looking at the vertical grid specifier
$entry = @field[1];

# Extract the command and the parameter list
$left  = index($entry,'(');
$right = index($entry,')');
if ( ($left != -1) && ($right != -1) )
   {
    $cmd  = substr($entry,0,$left);
    $list = substr($entry,$left+1,$right-$left-1);
    $len  = length $list; 
   }
else
   {    
    die('Invalid expression... Check vertical [] separator position ');
   }   

# Now building the parametr list 
$len  = length $list; 
$npar = 0;
while ( $len > 0 ) 
   {
   $next = index($list,',');
   if ( $next != -1 )
      {
       @param[$npar] = substr($list,0,$next);
       $list         = substr($list,$next+1,$len-$next-1);
       $len          = length $list; 
       $npar         = $npar + 1;
      }
   else
      {
       @param[$npar] = substr($list,0,$len);
       $len          = 0;
       $npar         = $npar + 1;
      }
   }

# Check for syntax (needed number of parameters)
if ( ($cmd eq "file")      && ($npar != 1) ) 
    {  die('Invalid vertical mode[file]... Check number of parameters ');     }
if ( ($cmd eq "level")     && ($npar != 1) ) 
    {  die('Invalid vertical mode[level]... Check number of parameters ');    }
if ( ($cmd eq "list")      && ($npar == 0) ) 
    {  die('Invalid vertical mode[list]... Check number of parameters ');     }
if ( ($cmd eq "profile")   && ($npar != 3) ) 
    {  die('Invalid vertical mode[profile]... Check number of parameters ');  }
if ( ($cmd eq "grid")      && ($npar != 2) ) 
    {  die('Invalid vertical mode[grid]... Check number of parameters ');     }

# Write parameters
print "\"$cmd\"\n";
if ( $cmd eq "list") { print "$npar\n"; }
if ( $npar > 0 )
      { print "@param\n"; }

# ----- Vertical coordinate system  ----------------------------------

# Now looking at the vertical grid specifier
$cmd = @field[2];

# Check for allowed coordinate axes
if ( ($cmd ne "hPa") && ($cmd ne "hPa,agl") && ($cmd ne "K") && ($cmd ne "PVU") && ($cmd ne "INDEX") )
    {  die('Invalid vertical axis [allowed: hPa / hPa,agl / K / PVU / INDEX] ');     }

# Write command
print "\"$cmd\"\n";

# ----- Selection criteria --------------------------------------------

# Now looking at the selection specifier
$cmd = @field[3];

# Write command
print "$cmd\n";
