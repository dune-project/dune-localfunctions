#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;

#######################################################################
#
# Help text
#
my $help = <<EOF;
make-localfunctions2vtu -- compile a VTK output executable for a specific
                           local basis

SYNOPSIS:
    make-localfunction2vtu --help|-h|-?
    make-localfunction2vtu --basis-help|-H LOCAL-BASIS
    make-localfunction2vtu --list
    make-localfunction2vtu --source LOCAL-BASIS
    make-localfunction2vtu LOCAL-BASIS

PARAMETERS:
    --basis-help LOCAL-BASIS
        Give help for the given LOCAL-BASIS.
    --help|-h|-?
        Print this help.
    --list
        List available local basis templates.
    --source LOCAL-BASIS
        Generate source for the given LOCAL-BASIS.  The header is written to a
        file named after the program name with ".cc" appended.
    LOCAL-BASIS
        Call the makefile with the apropriate parameters to build the program
        for this basis.
EOF
    ;
#######################################################################
#
# List of supported bases
#

#hash with key shortname -- full name of the base excluding template params
#each value is a hash with keys
# tparams  -- descriptive string of template params, if any
# help     -- help text
# headers  -- string holding required headers
# progname -- subroutine taking a list of template params and returning a
#             suffix for the resulting program
my %bases = (

#======================================================================
'Dune::P0LocalBasis' => {
tparams  => '<typename DomainFieldType, typename RangeFieldType, int dimension>',

help     => <<EOH,
Defines the constant scalar shape function in d dimensions. Is valid on any
type of reference element.
 * DomainFieldType: Type to represent the field in the domain.
 * RangeFieldType: Type to represent the field in the range.
 * dimension: Domain dimension.
EOH

headers  => <<EOH,
#include <dune/finiteelements/p0/p0localbasis.hh>
EOH

progname => sub {
    my $D = shift;
    my $R = shift;
    my $d = shift;
    return sprintf "p0-%s-%s-%dd", lc $D, lc $R, $d;
}},

#======================================================================
);

sub decode_basis($ ) {
    my $basis = $_[0];
    $basis =~ s/[[:space:]]+/ /g;
    $basis =~ s/^ //;
    $basis =~ s/ $//;

    $basis =~ s/^([[:alnum:]:_]+)[[:space:]]*//;
    my $name = $1;
    $name =~ s/ //g;
    my @result = ($name,);
    if($basis =~ s/^<(.*)>$/$1/) {
        push @result, map { s/^ //; s/ $//; $_ } split m/,/, $basis;
    }
    elsif($basis ne "") {
        die "Invalid basis specification \"$_[0]\"\n";
    }
    if(wantarray) { return @result }
    else          { return $result[0] }
}

sub find_basis($ ) {
    my $name = shift;
    if(defined $bases{$name}) {
        return $name;
    }
    if(defined $bases{"Dune::$name"}) {
        return "Dune::$name";
    }
    die "Can't find basis $name\n";
}

sub basis_help($ ) {
    my $name = shift;
    my $help = $name;
    $help .= $bases{$name}{tparams}
        if defined $bases{$name}{tparams};
    $help .= "\n\n";
    $help .= $bases{$name}{help}
        if defined $bases{$name}{help};
    return $help;
}

sub list_bases() {
    my @bases;
    for(sort keys %bases) {
        if(defined $bases{$_}{tparams}) {
            push @bases,  $_.$bases{$_}{tparams};
        } else {
            push @bases, $_;
        }
    }
    if(wantarray) { return @bases }
    else {
        return join "", map { "$_\n" } @bases;
    } 
}

sub gen_source($ ) {
    my($basis, @tparams) = decode_basis shift;
    $basis = find_basis $basis;
    my $fullname = $basis;
    if(@tparams) {
        $fullname .= "<".join(", ", @tparams).">";
    }

    my $progname = "functions2vtu-".$bases{$basis}{progname}(@tparams);
    my $sourcename = "$progname.cc";

    open(SOURCE, '>', $sourcename) or die "Can't open $sourcename for writing\n";

    print SOURCE <<EOF;
#ifdef HAVE_CONFIG_H
# include "config.h"     
#endif

EOF
    print SOURCE $bases{$basis}{headers}
        if defined $bases{$basis}{headers};
    print SOURCE <<EOF;

#define LOCAL_BASIS_TYPE $fullname
#define LOCAL_BASIS_TYPE_S "$fullname"
#define PROG_NAME "$progname"

#include "localfunctions2vtu.hh"
EOF

    close SOURCE;
}

#######################################################################
#
#  Parse Options
#

my $result = GetOptions(
    "help|l|?" => sub {
        print $help;
        exit 0;
    },
    "basis-help|H=s" => sub {
        print basis_help find_basis $_[1];
        exit 0;
    },
    "list" => sub {
        print scalar list_bases;
        exit 0;
    },
    "source=s" => sub {
        gen_source $_[1];
        exit 0;
    },
    );
if(not $result) { exit 1 }

if(@ARGV != 1) {
    die "Need a basis to generate the Program for.  Try --help.\n";
}

my($basis, @tparams) = decode_basis shift;
$basis = find_basis $basis;
my $fullname = $basis;
if(@tparams) {
    $fullname .= "<".join(", ", @tparams).">";
}

my $progname = "functions2vtu-".$bases{$basis}{progname}(@tparams);
print "Generating $progname\n";

my $command = "";
if(defined $ENV{MAKE} and $ENV{MAKE} ne "") {
    $command .= $ENV{MAKE};
} else {
    $command .= "make";
}

$command .= " -f localfunctions2vtu.make BASIS='$fullname' PROG_NAME=$progname $progname";

print "$command\n";
exec $command;
die "could not execute $command\n"
