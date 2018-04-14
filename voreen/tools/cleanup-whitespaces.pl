#!/usr/bin/perl
# Replaces tabs by 4 spaces, removes trailing whitespace on line end, adds newline to EOF

use strict;
use warnings;
use File::Find;
use File::Copy;

my @dirs = qw(src include ext/tgt);
@dirs = (@dirs, qw(apps/voreenve apps/voltool apps/voreentool apps/simple apps/itk_wrapper apps/tests));
@dirs = (@dirs, qw(modules/core modules/base modules/experimental modules/deprecated modules/staging));
@dirs = (@dirs, qw(modules/advancedlighting modules/camino modules/connexe modules/devil modules/dti));
@dirs = (@dirs, qw(modules/dynamicglsl modules/ffmpeg modules/flowreen modules/gdcm modules/hdf5));
@dirs = (@dirs, qw(modules/itk modules/itk_generated modules/onscreengui modules/opencl modules/openmp modules/pfskel));
@dirs = (@dirs, qw(modules/plotting modules/pvm modules/python modules/randomwalker modules/roi modules/sample modules/segy modules/tiff modules/volumelabeling));
@dirs = (@dirs, qw(modules/shadowmapping modules/stereoscopy modules/tiff modules/touchtable modules/volumelabeling modules/voreenbiology));
@dirs = (@dirs, qw(modules/workflow modules/zip));
my $interesting = 'cpp|h|frag|vert|geom|cu|cl';
my $skip_dirs = '\.svn|\.moc|\.obj|\.ui|ext';
$skip_dirs .= '|modules\/\w+\/ext';

find(\&wanted, @dirs);

sub wanted {
    if ($_ =~ /^($skip_dirs)$/) {
        $File::Find::prune = 1;
        return;
    };
    if ((-f $_) and ($_ =~ /\.($interesting)$/)) {
        update_untabify($_);
    }
}

sub update_untabify {
    my $file = shift;
    open F, $file or die "Can't open $file: $!";
    my @f = <F>;
    close F;

    my $tmpfile = $file . ".tmp";
    open F, ">", $tmpfile or die "Can't open $tmpfile for writing: $!";
    my $line = 1;
    my $changed = 0;
    my $tabs = 0;
    my $trailing = 0;

    my $l = '';
    foreach (@f) {
        $l = $_;
        if ($l =~ s/\t/    /g) {
            $changed = 1;
            $tabs++;
        }
        if ($l =~ s/\s+\n$/\n/) {
            $changed = 1;
            $trailing++;
        }

        print F $l;
        $line++;
    }
    my $changeinfo = '';
    $changeinfo .= "$tabs tabs" if ($tabs > 0);
    if ($trailing > 0) {
        $changeinfo .= ', ' if $changeinfo ne '';
        $changeinfo .= "$trailing trailing" ;
    }
    if ($l !~ /\n$/) {
        print F "\n";
        $changed = 1;
        $changeinfo .= ', ' if $changeinfo ne '';
        $changeinfo .= 'added newline to end';
    }

    close F;
    if ($changed) {
        move($tmpfile, $file);
        print $File::Find::name . ": $changeinfo\n";
    } else {
        unlink($tmpfile);
    }
}
