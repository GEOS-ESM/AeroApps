#!/usr/bin/perl

# --------------------------------------------------------
# Separate different commands out (delimiter & and |)
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
  $and  = index($cmd,'&');
  $or   = index($cmd,'|' );
  if ( ($and == -1) && ($or == -1) )
     { $next = $len+1; }
  elsif ( ($and == -1) && ($or >= 0) )
     { $next = $or;  @logop[$nline] = 'OR'; }
  elsif ( ($and >= 0) && ($or == -1) )
     { $next = $and; @logop[$nline] = 'AND'; }
  elsif ($and > $or)
     { $next = $or; @logop[$nline] = 'OR';  }
  else
      { $next = $and; @logop[$nline] = 'AND'; }
  
  # A logical operator is not allowed to be at position 0
  if ( $next == 0) 
   {    
      die('Invalid expression... Check logical operators & and |');
   }   

  # Extract the next substring
  $sub = substr($cmd,0,$next);
  $cmd = substr($cmd,$next+1,$len-$next-1);
  $len = length $cmd; 
 
  # Save the command in a new line
  @field[$nline] = $sub;
  $nline = $nline + 1;
}

# --------------------------------------------------------
# Handle each command line separately
# --------------------------------------------------------

# Write start marker
print "BEGIN \n";

foreach ( $i = 0; $i < $nline; $i++  )
{
    # Split the command into its four components
    @entry  = split /:/, @field[$i];
    $nentry = @entry;

    # Either four or three elements are needed
    if ( ($nentry < 3) | ($nentry > 4) )     
     {    
      die('Each expression needs either 3 or 4 fields...');
     }   

    # Write the command
    print "@entry[0] \n";

    # Get the variable and the 'mode' for this variable
    $left  = index(@entry[1],'(');
    $right = index(@entry[1],')');
    if ( ($left == -1) && ($right == -1) ) 
      { 
       $var  = @entry[1];
       $mode = 'VALUE';
      }
    elsif ( ($left > 0) && ($right > $left) ) 
      { 
       $var  = substr(@entry[1],0,$left);  
       $mode = substr(@entry[1],$left+1,$right-$left-1);
      }
    else
      {    
       die('Invalid variable specification...');
      }   
    print "$var  $mode \n";

    # Get the parameter list for this command
    @param  = split /,/, @entry[2];
    $nparam = @param;
    print "$nparam \n";
    print "@param \n";

    # If only three parameters are given, the time is assumed to be the first one
    if ( $nentry == 3 )
    {
	@entry[3]='FIRST';
    }    

    # Get the variable and the 'mode' for this variable
    $left  = index(@entry[3],'(');
    $right = index(@entry[3],')');
    if ( ($left == -1) && ($right == -1) ) 
      { 
       $time = @entry[3];
       $mode = 'ALL';
      }
    elsif ( ($left > 0) && ($right > $left) ) 
      { 
       $time = substr(@entry[3],0,$left);  
       $mode = substr(@entry[3],$left+1,$right-$left-1);
      }
    else
      {    
       die('Invalid time specification...');
      }   

    # Get the time list for this command
    if ( $time eq 'ALL' )
    {
	$time = -999;
        print "1 \n";
        print "$time \n";
    }
    elsif ( $time eq 'FIRST' )
    {
        $time = -996;
        print "1 \n";
	    print "$time \n";
    }
    elsif ( $time eq 'LAST' )
    {
        $time = -995;
        print "1 \n";
	    print "$time \n";
	}
	elsif ( $time eq 'TRIGGER' )
    {
        print "-993 \n";
    }
    else
    {
        @times  = split /,/, $time;
        $outstr = "";
	    $outlen = 0;
	    foreach $j ( @times )
        {
	    $to = index($j,'to');
            if ( $to == -1 )
               { 
		if ( $j eq 'FIRST' )
		  {
                      $outlen=$outlen+1;
		      $outstr=$outstr . " -996 ";
                  }
                elsif ( $j eq 'LAST' )
		  {
                      $outlen=$outlen+1;
		      $outstr=$outstr . " -995 ";
                  }
                else
		  {
                   $outlen=$outlen+1;
                   $outstr=$outstr . "$j "; 
	          }
               }
            else
	       { 
                $outlen = $outlen+3;
                $t1 = substr($j,0,$to);
		$t2 = substr($j,$to+2,length($j)-$to+1);
		if ( $t1 eq 'FIRST' ) 
		    { $t1='-996'; };
		if ( $t2 eq 'FIRST' ) 
		    { $t2='-996'; };
		if ( $t1 eq 'LAST'  ) 
		    { $t1='-995'; };
		if ( $t2 eq 'LAST'  ) 
		    { $t2='-995'; };
	        $outstr=$outstr . $t1 . " -994 " . $t2 . " "; 
	       }
        }
	print "$outlen \n";
        print "$outstr \n";
    }

    # Write the time mode
    if ( $mode eq 'ALL' )
    {
        print "$mode \n";
    }
    elsif ( $mode eq 'ANY' )
    {
        print "$mode \n";
    }
    elsif ( $mode eq 'NONE' )
    {
        print "$mode \n";
    }
    elsif ( $mode eq 'TRIGGER' )
    {
        print "$mode \n";
    }
    else
    {    
       die('Invalid time mode...');
    }   

	
    # Write the logical operator
    if ( $i < $nline-1 )
    {
      print "@logop[$i] \n";
    } 
}

# Write end marker
print "END \n";
