#!/usr/bin/perl
#
# Recursively removes files and directories from the passed directory ($ARGV[0])
# whose relative path matches any of the regular expressions listed in the 
# files $ARGV[2:n].
#
use strict;
use warnings;

use File::Path qw(rmtree);
use File::Spec;

# read work path
if (not $ARGV[0]) {
    die "error: no path specified\n";
}
elsif (not $ARGV[0] =~ m/\w/) {
    die "error: invalid path specified: $ARGV[0]\n";
}
elsif (not -d $ARGV[0]) {
    die "error: '$ARGV[0]' is not a directory\n";
}
my $workpath = $ARGV[0]; 

# read exclude files
my @excludelist;
for (my $i = 1; $i < @ARGV; $i++) {
    print "Reading '$ARGV[$i]'...\n";
    open F, $ARGV[$i] or die "Can't open $ARGV[$i]: $!.";
    my @lines = <F>;
    foreach my $l (@lines) {
        chomp $l;
        
        # skip comments and lines with only whitespace
        next if $l =~ /^\s*(#.*)?$/;
        
        @excludelist = (@excludelist,$l);
    }
    close F;    
}

cleanpath($workpath);

# Recursively removes all files and directories whose path (relative to $workpath) 
# match any of the expressions in the global @excludelist. 
sub cleanpath {
    my $dir = shift;
    
    local *DIR;
    opendir DIR, $dir or return;
    
    for (readdir DIR) {
        next if /^\.{1,2}$/;
        my $path = "$dir/$_";
        my $relpath = File::Spec->abs2rel($path, $workpath);
        if (-d $path) {  # directory
            if (listmatch($relpath) == 1) {
                print "Delete dir:   $path \n";
                rmtree($path);
            }
            else {
                cleanpath($path);
            }
        }
        elsif (-f $path) { # file
            if (listmatch($relpath) == 1) {
                print "Delete file:  $path \n";
                unlink($path);
            }
        }
        else {
            die "error: not a file or directory: $path\n";
        }
    }
    
    closedir DIR;
}

# Returns 1, if the passed string matches any 
# of the expressions in the global @excludelist.
sub listmatch {
    my $path = $_[0];
    foreach my $exp (@excludelist) {
        return 1 if ($path =~ m/($exp)/);
    }    
    return 0;
}
