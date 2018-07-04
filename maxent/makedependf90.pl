#!/usr/bin/perl

# Generic parser module, which allows you to tokenize fortran source files
package FParser;
use strict; use warnings;
my $tokenizer = qr/(?>\G\s*(?> !.*                                   # comments
      | ((?> "(?>""|&\s*[^\s]|[^"&])*"? | '(?>''|&\s*[^\s]|[^'&])*'? # strings
        | \d[+\-\.\w]* | \w+                                         # words 
        | ::? | =>? | [\<>]=? | \.\w+\. | \(\/? | \/[=\)]? | \*\*?   # operators
        | .))\s*                                                     # any,pp 
        \s*) )/x;
my $preprocessor = qr/(?>^\s*#)/;   # preprocessor statements
my $fixedcomment = qr/(?>^[*cC!])/; # fixed-format comment
my $contsym = qr/(?>^\s*&)?/;       # continuation symbol on line 
my $fortfixed = qr/\.f(?>77|ix|or|pp|tn|)$/i;  # fixed-format extension
my $fortfree = qr/\.f(?>90|95|03|08)$/i;       # free-format extension
 
our $debug;        
our $srcform;
our $ignorepp;
our $linelen = 72;

# Tokenize free or fixed-format Fortran source file into statements consisting 
# of a list of lexical tokens. No syntax check is done. Expects two arguments:
#  - 0th: the filename. The tokenizer distinguishes between fixed and free format
#         based on common extension conventions if $srcform is not specified.
#  - 1st: a callback subroutine, which is called once for every statement 
#         (semicolons and line continuations are handled correctly). The 0-th
#         argument is the line number (or -1 if none is present), the others
#         correspond to the lexical tokens of the statement.
sub tokenizef {
    my ($fname, $parser) = @_;
    my $fixed;

    open(my $SRC, '<', $fname) or die "Error opening $fname: ".$!."\n";
    if($srcform ? $srcform eq 'free' : $fname =~ /$fortfree/) { 
        $fixed = undef; 
    } elsif($srcform ? $srcform eq 'fixed' : $fname =~ /$fortfixed/) {
        $fixed = 1;
    } else{ die "$fname: Illegal file format and no valid fallback -f <form>"; }
    
    print STDERR "==== FILE: $fname (", $fixed?'fixed':'free', ")\n" if $debug;
    my @toks = (); my $contd = undef;
    LINE: while(my $line = <$SRC> || '') {
        if($line =~ /$preprocessor/) {    # preprocessor statements (ignore)
            while($line =~ /\\$/) { $line .= <$SRC>; }   # continuation lines
            print STDERR "$fname($.): Warning: Preprocessor directive (ignored)"
                            ."\n$line" unless $ignorepp;
            next LINE;
        }
        if($fixed) {   # fixed-format fortran source readout
        	next LINE if $line =~ /$fixedcomment/;
            chomp($line);
            if(length($line) > 6) {
                $contd = substr($line, 5, 1) ne ' ';
                $line = ($contd ? '' : substr($line,0,5))    # cut off line
                           . substr($line, 6, $linelen-6)."\n";
            }
        }
        PARSELINE: if($contd) {   # remove & sign, handle continued strings
            die "$fname($.): Illegal line continuation\n$line" unless @toks;
            $line =~ s/$contsym//;
            $line = pop(@toks).$line if $toks[-1] =~ /^['"]/;
            $contd = undef;
        } elsif(@toks) {   # at non-continued lines (fixed-form)
            print STDERR 'Tokens: ', join("|",@toks), "\n" if $debug;
            unshift(@toks, -1) unless $toks[0] =~ /^\d/;  # add line number 
            $parser->(@toks);  # invoke callback parsing routinge
            @toks = ();
        }
        last LINE if eof($SRC);
        TOKEN: while(($line =~ /$tokenizer/g) and defined($1)) {
            if($1 eq ';') {            # multiple statements
                goto PARSELINE;
            } elsif($1 eq '&') {       # continued statements
                $contd = 1;
            } else {                   # push rest on token list'
                push @toks, $1;        
            }
        } 
    } continue {
        print STDERR "Line $.: $line" if $debug;   # print afterwards
    }
    close $SRC;
}

# Actual perl script mimicing the behaviour of gcc -M.
package main;
use strict; use warnings;
use Getopt::Long qw(:config bundling);
use Pod::Usage;

my $outfile;
my $out_ext = 'o';
my $moddir = '';
my @exclude_mod_array = ();
my @infiles = ();

GetOptions(
    "o|output-file=s" => \$outfile,
    "f|file-format=s" => \$FParser::srcform,
    "v|verbose!"      => \$FParser::debug,
    "l|line-length=i" => \$FParser::linelen,
    "i|ignore-pp!"    => \$FParser::ignorepp,
    "x|exclude-mod=s" => \@exclude_mod_array,
    "e|output-ext=s"  => \$out_ext,
    "M|module-dir=s"  => sub{ $_[0] =~ s/[^\/]$/$&\//; $moddir = $_[0]; },
    "help"            => sub{ pod2usage(1); },
    "version"         => sub{ print "makedependf90  1.0\n"; exit 1; },
    man               => sub{ pod2usage(-exitval=>0, -verbose=>2); }
) or pod2usage(2);
pod2usage("No input files (--man for details)") unless @ARGV; 

my %exclude_mods = map { $_ => 1 } @exclude_mod_array;

open(OUT, '>'.($outfile||'-')) or die "Error opening '$outfile': $!\n";
for my $fname (@ARGV) {  # for each input file
    my (%deps, %mods) = ();   # array of required/provided dependencies
    &FParser::tokenizef($fname, sub {
    	    # gets line number as zeroth, tokens as other arguments  
            my $disc = lc($_[1]);  
            if($disc eq 'use') { # required modules
                $deps{ lc($_[2]) } = 1;
            } elsif($disc eq 'module') {  # provided modules
                my $name = lc($_[2]);
                if($name ne 'procedure') {  # ignore "module procedure"
                    $mods{ $name } = 1;
                }
            }
        });
    $fname =~ s/\..+?$/.$out_ext/;   # name mangling
    print OUT "\n# Dependency rules generated for $fname:\n";
    if(%mods) {
        print OUT "#   modules provided (empty production rules):\n";
        for my $mod (keys(%mods)) {
            print OUT "$mod.mod: $fname\n\t\@true\n";
        }}
    if(%deps) {
        print OUT "#   modules required (dependencies):\n";
        for my $mod (keys(%deps)) {
            # avoid dependency loop and do not include explicit excludes
            if(not exists($mods{$mod}) and not exists($exclude_mods{$mod})) { 
                print OUT "$fname: $moddir$mod.mod\n";
            }
        }}
}
close OUT;

# Man page
__END__

=head1 NAME

makedependf90 - Generate Fortran 90 module dependencies for Makefiles

=head1 SYNOPSIS

makedependf90 [B<-o> I<outputfile>] [B<-f> free|fixed] I<sourcefile>...

=head1 OPTIONS

 -o,--output-file <file>   Write to <outfile> instead of STDOUT
 -f,--file-format <form>   Enforces free or fixed form source files
 -l,--line-length <len>    Specifies maximum fixed-form line length
 -M,--module-dir  <dir>    Modules are located in <moddir>
 -i,--ignore-pp            Explicity ignore preprocessor directives
 -x,--exclude-mod <mod>    Modules that should be ignored
 -e,--output-ext <ext>     Assumed extension for object files
 -v,--verbose              Prints debug output (tokens, lines read) 
 --help                    Prints help and exits
 --man                     Displays the man page, then exits

=head1 DESCRIPTION

Parses free-form and fixed-form Fortran 90/95 files and extracts C<use> and 
C<module> statements, assembling a list of dependencies (used modules) and 
provided dependencies (declared modules).  These dependencies are then printed
as makefile rules to standard output (or the file specified in B<-o>).  This is 
essentially the behaviour of B<gcc -M> for Fortran files.

This allows you to use the C<include> feature of GNU Make to automatically
generate the correct dependency tree for compiling Fortran source files by
adding something like that to your B<makefile>: 

    SOURCES = $(wildcard *.f90)
    
    $(SOURCES): %.dep: %.f90
        perl makedependf90.pl -o $@ $<

    # Make reinvokes itself on updated dependency files
    ifneq ($(MAKECMDGOALS),clean)  
       -include $(SOURCES:.f90=.dep)
    endif
 
=head1 EXTENDED OPTIONS DESCRIPTION

=over 

=item B<-f, --file-format> free|fixed

Enforce either free-format or fixed-format Fortran source files. If not
supplied, the program tries to guess the format by the file extension by the 
following conventions:

   Fixed-format:  .f .f77 .for .fix .ftn .fpp  (case-insensitive)
   Free-format:   .f90 .f95 .f03 .f08          (case-insensitive)

=item B<-i, --ignore-pp>

Explicitly ignore pre-processor directives. Pre-processor directives are always
ignored by the program, but usually a warning is issued if such a directive is 
encountered.  Precede with B<fpp> in your toolchain to support preprocessing.  

=back

=head1 KNOWN BUGS AND LIMITATIONS

All preprocessor directives are ignored.  There is currently no support for
Fortran 2008 submodules.  The moronic makedepend Makefile meddling is also not
supported (DO NOT DELETE THIS LINE)

=head1 AUTHOR AND COPYRIGHT

Copyright (c) 2012 Markus Wallerberger (L<wallerberger@ifp.tuwien.ac.at>)

This program is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software 
Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with 
this program.  If not, see L<http://www.gnu.org/licenses/>.

=cut
