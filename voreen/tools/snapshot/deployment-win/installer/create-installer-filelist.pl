#!/usr/bin/perl
#
use strict;
use warnings;

use File::Path qw(make_path remove_tree rmtree);
use File::Spec::Win32;

# read work path
if (not $ARGV[0]) {
    die "error: no path specified\n";
}
elsif (not -d $ARGV[0]) {
    die "error: '$ARGV[0]' is not a directory\n";
}
my $workpath = $ARGV[0]; 

#open (MYFILE, '>voreen-installer-files.nsi');
listfiles($workpath);
#close (MYFILE); 

# Recursively removes all files and directories whose path (relative to $workpath) 
# match any of the expressions in the global @excludelist. 
sub listfiles {
    my $dir = shift;
 
    my $reldir = File::Spec::Win32->abs2rel($dir, $workpath);
    $reldir =~ s/\//\\/g;

    #print MYFILE "\n\${AddItem} \"\$INSTDIR\\$reldir\"\n" if not ($reldir eq ".");
    #print MYFILE "\${SetOutPath} \"\$INSTDIR\\$reldir\"\n";
    print "\n\${AddItem} \"\$INSTDIR\\$reldir\"\n" if not ($reldir eq ".");
    print "\${SetOutPath} \"\$INSTDIR\\$reldir\"\n";
    
    local *DIR;

    # first run: files only
    opendir DIR, $dir or return;
    for (readdir DIR) {
        next if /^\.{1,2}$/;
        my $path = "$dir/$_";
        next if (-d $path);
        die "error: not a file or directory: $path\n" if (not -f $path);
        
        my $relpath = File::Spec::Win32->abs2rel($path, $workpath);
        my ($volume, $filepath, $filename) = File::Spec->splitpath($path);
        $filepath =~ s/\//\\/g;
        #print MYFILE "\${File} \"$filepath\" \"$filename\"\n";
        print "\${File} \"$filepath\" \"$filename\"\n";
    }
    closedir DIR;
    
    # second run: directories only
    opendir DIR, $dir or return;
    for (readdir DIR) {
        next if /^\.{1,2}$/;
        my $path = "$dir/$_";
        next if (-f $path);
        die "error: not a file or directory: $path\n" if (not -d $path);
        listfiles($path);
    }
    closedir DIR;
}
