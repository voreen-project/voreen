#!/usr/bin/perl
use File::Copy;

sub usage() {
  print "\nUsage:  $0 <path> <win|unix>\n";
  print "\nSample: $0 voreen-src-2.6 win \n\n";
  exit 1;
}

usage() if (@ARGV < 2);
usage() if ((not $ARGV[1] eq "win") && (not $ARGV[1] eq "unix"));

die "Directory does not exist: '$ARGV[0]'" if (not -d $ARGV[0]);
        
# run convertEol function and print number of modified files
my @modifiedFiles;
convertEol($ARGV[0]);
print("convert-eol.pl: ".scalar(@modifiedFiles)." files modified\n");


# processes all .h/.cpp/.vert/.frag/.txt/.pro/.pri/.bat files
sub convertEol {
    my $dir = shift;
    
    for my $file (<$dir/*>) {
        if (-d $file) {
            convertEol($file);
        }   
        elsif (-f $file) {
            $_ = $file;
            if (m/\.(h|cpp|vert|frag|txt|pro|pri|bat)$/) {
                print "Processing '$file' ...\n";       
            
                # read lines from input file
                open F, $file or die "Can't open $file: $!";
                binmode F;
                my @f = <F>;
                close F;

                # open temporary output file
                my $tmpfile = $file . "-tmp";
                open F, ">", $tmpfile or die "Can't open $tmpfile for writing: $!";
                binmode F;
                
                # write convert input lines to output file
                my $changed = 0;
                my $l = '';
                foreach (@f) {
                    my $l = $_;
                    my $l_conv = $l;
                    
                    if ($ARGV[1] eq "win") {
                        $l_conv =~ s/\r\n|\n|\r/\r\n/g;  # Convert to win
                    }                    
                    elsif ($ARGV[1] eq "unix") {
                        $l_conv =~ s/\r\n|\n|\r/\n/g;    # Convert to unix
                    }  
                    else {
                        die "unknown mode: $ARGV[1]";
                    }
                    print F $l_conv;

                    $changed = 1 if (not $l_conv eq $l) 
                }
                close F;
                if ($changed) {
                    @modifiedFiles = (@modifiedFiles, $file);
                }

                # replace original file by tmp file
                move($tmpfile, $file);
            }
       }
    }
}

