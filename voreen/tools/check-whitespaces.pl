#!/usr/bin/perl
# Checks source files for use of tabs and outputs results in JUnit format for use with Hudson.

use strict;
use warnings;
use File::Find;

my @dirs = qw(src include ext/tgt);
@dirs = (@dirs, qw(apps/voreenve apps/voltool apps/voreentool apps/simple apps/itk_wrapper apps/tests));
@dirs = (@dirs, qw(modules/core modules/base modules/experimental modules/deprecated modules/staging));
@dirs = (@dirs, qw(modules/advancedlighting modules/camino modules/connexe modules/devil modules/dti));
@dirs = (@dirs, qw(modules/dynamicglsl modules/ffmpeg modules/flowanalysis modules/gdcm modules/hdf5));
@dirs = (@dirs, qw(modules/itk modules/itk_generated modules/onscreengui modules/opencl modules/openmp modules/pfskel));
@dirs = (@dirs, qw(modules/plotting modules/pvm modules/python modules/randomwalker modules/roi modules/sample modules/segy modules/tiff modules/volumelabeling));
@dirs = (@dirs, qw(modules/shadowmapping modules/stereoscopy modules/tiff modules/touchtable modules/volumelabeling modules/voreenbiology));
@dirs = (@dirs, qw(modules/workflow modules/zip));
my $interesting = 'cpp|h|frag|vert|geom|cu|cl';
my $skip_dirs = '\.svn|\.moc|\.obj|\.ui|ext';
$skip_dirs .= '|modules\/\w+\/ext';
my $total_tabs = 0;
my @messages = ();

chdir('..');
find({ wanted => \&wanted, no_chdir => 1 }, @dirs);

my $result = "";

if ($total_tabs > 0) {
  $result = qq(<error message="found $total_tabs lines with tabs in source code, use spaces instead!">);
  foreach (@messages) {
      $result .= $_ . "\n";
  }
  $result .= '</error>';
}

    print << "(EOF)";
<testsuite name="CodeTests">
   <testcase classname="Source.Style" name="tabs">
      $result
   </testcase>
</testsuite>
(EOF)

sub wanted {
    if ($_ =~ /^($skip_dirs)$/) {
        $File::Find::prune = 1;
        return;
    };
    if ((-f $_) and ($_ =~ /\.($interesting)$/)) {
        check_tabs($_);
    }
}

sub check_tabs {
    my $file = shift;
    open F, $file or die "Can't open $file: $!";
    my @f = <F>;
    close F;

    my $line = 1;
    my $changed = 0;
    my $tabs = 0;
    my $trailing = 0;

    my $l = '';
    foreach (@f) {
        $l = $_;
        if ($l =~ /\t/) {
            $changed = 1;
            $tabs++;
        }
        $line++;
    }
    close F;

    if ($tabs > 0) {
        $total_tabs += $tabs;
        push @messages, "$file: found $tabs lines with tabs";
    }
}
