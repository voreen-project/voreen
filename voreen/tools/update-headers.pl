#!/usr/bin/perl

use strict;
use warnings;
use File::Find;
use File::Copy;

my $header_file_voreen = "header.txt";
my $header_file_tgt = "ext/tgt/header.txt";

my @dirs_voreen = qw(src include modules custommodules apps tools);
my @dirs_tgt =    qw(ext/tgt);
my $interesting = 'cpp|h|frag|vert|geom|cu|cl';
my $skip_dirs = '\.moc|\.obj|\.ui|ext';

open HEADER_VOREEN, $header_file_voreen or die "Can't open $header_file_voreen: $!";
my @header_voreen = <HEADER_VOREEN>;
close HEADER_VOREEN;

open HEADER_TGT, $header_file_tgt or die "Can't open $header_file_tgt: $!";
my @header_tgt = <HEADER_TGT>;
close HEADER_TGT;

print "1. Updating license headers in voreen files ...\n";
find(\&wanted_voreen, @dirs_voreen);
print "\n";

print "2. Updating license headers in tgt files ...\n";
find(\&wanted_tgt, @dirs_tgt);

sub wanted_voreen {
    if ($_ =~ /^($skip_dirs)$/) {
        $File::Find::prune = 1;
        return;
    };
    if ((-f $_) and ($_ =~ /\.($interesting)$/)) {
        update_header($_, "voreen");
    }
}

sub wanted_tgt {
    if ($_ =~ /^($skip_dirs)$/) {
        $File::Find::prune = 1;
        return;
    };
    if ((-f $_) and ($_ =~ /\.($interesting)$/)) {
        update_header($_, "tgt");
    }
}

sub update_header {
    my $file = shift;
    my $mode = shift;
    print $File::Find::name . "\n";
    open F, $file or die "Can't open $file: $!";
    my @f = <F>;
    close F;

    my $tmpfile = $file . ".tmp";
    open F, ">", $tmpfile or die "Can't open $tmpfile for writing: $!";
    if ($mode =~ "voreen") {
        print F @header_voreen;
    }
    elsif ($mode =~ "tgt") {
        print F @header_tgt;
    }
    else {
        die "ERROR: unknown mode: $mode";
    }
    print F "\n";
    my $inblock = 0;
    my $wait_for_start = 1;
    my $line = 1;
    my $do_replace = 1;
    foreach (@f) {
        if ($line < 15 and (not $inblock) and ($_ =~ m/auto-generated|copyright|(\(c\))/i)) {
            print "  * Copyright or auto-generated marker found in $file:$line, skipping this file\n";
            $do_replace = 0;
            last;
        }
        if ($line < 5 and not $inblock) {
            if ($_ =~ m!/\*{70,85}$!) {
              $inblock = 1;
              $wait_for_start = 0;
            }
        }
        if ($wait_for_start && $_ !~ m/^\s*$/) {
            $wait_for_start = 0;
        }
        print F $_ if not ($inblock or $wait_for_start);
        $line++;

        if ($inblock) {
          if ($_ =~ m!\*{30,}/!) {
                $inblock = 0;
                $wait_for_start = 1;
          }
        }
    }
    close F;
    if ($do_replace) {
        move($tmpfile, $file);
    } else {
        unlink($tmpfile);
    }
}
