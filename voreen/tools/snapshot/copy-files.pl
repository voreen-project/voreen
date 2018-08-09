#!/usr/bin/perl
#
# Creates a snapshot archive from the files and directories specified in 'snapshot-include-common.txt'
# and the platform-specific 'snapshot-include-unix.txt' and 'snapshot-include-win32.txt' include files.
#
use strict;
use warnings;

use File::Copy;
use File::Copy::Recursive qw(fcopy rcopy fmove rmove dirmove);
use File::Path qw(rmtree mkpath);

# read output archive name
my $output = $ARGV[0];
if (not $output) {
    die "error: no output path specified\n";
}
elsif (not $ARGV[0] =~ m/\w/) {
    die "error: invalid output path specified: $ARGV[0]\n";
}
elsif (-f $output) {
    die "error: output path '$output' is a file\n";
}
elsif (@ARGV < 2) {
    die "error: no include files specified.\n";
}
my $outPath = $output;

my $force = 0;
$force = 1 if ($ARGV[@ARGV-1] eq "--force");

# read -r parameter
my $really = 1;

# read include files
my @filelist;
for (my $i = 1; $i < @ARGV; $i++) {
    
    # skip modifier parameters (--force)
    next if $ARGV[$i] =~ /^--/;
    
    print "Reading '$ARGV[$i]'...\n";
    open F, $ARGV[$i] or die "Can't open $ARGV[$i]: $!.";
    my @filelist_temp = <F>;
    close F;
    @filelist = (@filelist,@filelist_temp);
}
if (@filelist == 0) {
    print "\nNo files to copy.\n";
    exit;
} 
    
# create output directory
if (-d $outPath and $really) {
    if (not $force) {
        print "Directory '$outPath' already exists. Overwrite? (y=yes, c=cancel) ";
        my $userinput =  <STDIN>;
        chomp ($userinput);
        if (not ($userinput eq "y" or $userinput eq "yes")) {
            die "Cancel.\n";
        }
    }
    print "Delete directory '$outPath'.\n";
    rmtree($outPath);
}
print "Create directory '$outPath'.\n";
mkpath($outPath) if $really;
   
#
# Copy specified files/dirs to output directory
#   
print "\nCopying specified files and directories to directory '$outPath'...\n";
foreach my $line (@filelist) {
    chomp $line;

    # skip comments and lines with only whitespace
    next if $line =~ /^\s*(#.*)?$/;

    my $file;
    my $destpath;
    
    if ($line =~ m/:/) {
        # destination path is specified (substring behind ':')
        my @splitline = split(/:/, $line);
        $file = $splitline[0];
        $destpath = $splitline[1];
    }
    else {
        # lines contains only file path
        $file = $line;
    }
        
    if (-f $file) { # file
        my $volume, my $filepath, my $filename;
        ($volume, $filepath, $filename) = File::Spec->splitpath($file);
        $filepath = $destpath if $destpath;
        if ($filepath and not (-d '$outPath/$filepath')) {
            print "create dir:  $outPath/$filepath\n";
            mkpath("$outPath/$filepath") if $really;
        }
        print "copy file:   $file\n             -> $outPath/$filepath\n"; 
        copy($file, "$outPath/$filepath") if $really;
    }
    elsif (-d $file) { # directory
        my $volume, my $dirpath, my $dir;
        ($volume, $dirpath, $dir) = File::Spec->splitpath($file);
        $dirpath = $destpath if $destpath;
        if ($dirpath and not (-d '$outPath/$dirpath')) {
            print "create dir:  $outPath/$dirpath\n";
            mkpath("$outPath/$dirpath") if $really;
        }
        # copy directory recursively, unless it has been specified as "dir/"
        if ($dir) {
            print "copy dir:    $file\n             -> $outPath/$dirpath\n";
            rcopy($file, "$outPath/$dirpath/$dir") if $really;
        }
    }
    else { # directory files (without subdirs)
        my $volume, my $filepath, my $filename;
        ($volume, $filepath, $filename) = File::Spec->splitpath($file);
        my $srcpath = $filepath;
        $filepath = $destpath if $destpath;
        if ($filename and ($filename eq "*")) {
            if ($filepath and not (-d '$outPath/$filepath')) {
                print "create dir:  $outPath/$filepath\n";
                mkpath("$outPath/$filepath") if $really;
            }
            # copy each file in directory
            print "copy files:  $srcpath*\n             -> $outPath/$filepath\n";  
            for my $tfile (<$srcpath*>) {
                if (-f $tfile) {
                    my $tdirs, my $tfilename;
                    ($volume, $tdirs, $tfilename) = File::Spec->splitpath($tfile);
                    print "                $tfilename\n";  
                    copy($tfile, "$outPath/$filepath") if $really;
                }
            }
        }
        else {
            die "error: not a file or directory: $file\n";
        }
    }
}

print "\nFinished copy.\n";

exit;
