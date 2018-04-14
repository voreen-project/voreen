#!/usr/bin/perl
#
# Recursively processes all *.h and *.cpp files in the passed directory $ARGV[0]
# and removes the preprocessor symbols listed in the files $ARGV[2:n].
#
use File::Copy;

sub usage() {
  print "\nUsage:  $0 <path> <symbolfile> [<symbolfile2> ... <symbolfileN>]\n";
  print "\nSample: $0 voreen-src-2.6 unifdef-symbols.txt\n\n";
  exit 1;
}

usage() if (@ARGV < 2);
die "Directory does not exist: '$ARGV[0]'" if (not -d $ARGV[0]);

# read symbol files
my $unifdef_params = "";
for (my $i = 1; $i < @ARGV; $i++) {
    print "Reading '$ARGV[$i]'...\n";
    open F, $ARGV[$i] or die "Can't open $ARGV[$i]: $!.";
    my @lines = <F>;
    foreach my $l (@lines) {
        chomp $l;
        
        # skip comments and lines with only whitespace
        next if $l =~ /^\s*(#.*)?$/;
        
        # break if invalid symbol is encountered
        if (not $l =~ m/\-(D|U)\w+/) {
            die "Invalid symbol '$l'. Expecting: -Dsym oder -Usym";
        }
        $unifdef_params .= " $l";
    }
    close F;    
}

print "Unifdef.pl: calling \'unifdef $unifdef_params <f>\'\n";

        
# run unifdef tool and list modified files
my @modifiedFiles;
unifdef($ARGV[0]);
print("Unifdef.pl: ".scalar(@modifiedFiles)." files modified\n");
foreach (@modifiedFiles) {
    print " * $_\n";
}


# Runs the unifdef tool with the passed parameters on all .h/.cpp files
sub unifdef {
    my $dir = shift;
    
    for my $file (<$dir/*>) {
        if (-d $file) {
            unifdef($file);
        }   
        elsif (-f $file) {
            $_ = $file;
            if (m/\.(h|cpp)$/) {
                print "Processing '$file' ...\n";       
            
                # run unifdef: write result to tmp file
                my $ret = system("unifdef $unifdef_params $file > $file-tmp");  
                
                # check return value
                die "Unifdef.pl: failed to execute unifdef" if ($ret == -1);
                $ret = $ret >> 8;
                die "Unifdef.pl: error on processing '$file' (return code: $ret)\n" if ($ret > 1);
                die "Unifdef.pl: failed to process '$file': '$file-tmp' not found" if (not -f "$file-tmp");                
                
                # replace original file by tmp file
                move("$file-tmp", $file);

                # note file if modified
                @modifiedFiles = (@modifiedFiles, $file) if ($ret == 1);
            }
       }
    }
}
