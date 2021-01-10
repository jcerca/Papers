#!/usr/bin/perl

# USAGE: perl Split_STACKS_fasta_file.pl NAME-OF-STACKS-SAMPLE-FILE
# Torsten H. Struck, Jose Cerca
use strict;

my $file = $ARGV[0];
my $locusname = ();
my $locusname_Allele0 = ();
my $locusname_Allele1 = ();
my $newheader = ();

open FILE, $file;
while (my $line = <FILE>) {
 chomp $line;
 if ($line =~ m/^>C(Locus_\d+)/) {
  $locusname = "$1.fas";
  $locusname_Allele0 = "$1_Allele_0.fas";
  $locusname_Allele1 = "$1_Allele_1.fas";
  ($newheader = $line) =~ s/>CLocus_\d+_Sample_\d+_Locus_\d+_(Allele_\d+)\s+\[(\w+)\]/>$2_$1/;
  (my $newheader_sep = $line) =~ s/>CLocus_\d+_Sample_\d+_Locus_\d+_Allele_\d+\s+\[(\w+)\]/>$1/;
  open (OUTFILE, '>>', $locusname);
   print OUTFILE "$newheader\n";
  close OUTFILE;
  if ($newheader =~ m/_Allele_0$/) {
   open (OUTFILE1, '>>', $locusname_Allele0);
    print OUTFILE1 "$newheader_sep\n";
   close OUTFILE1;
  }
  if ($newheader =~ m/_Allele_1$/) {
   open (OUTFILE2, '>>', $locusname_Allele1);
    print OUTFILE2 "$newheader_sep\n";
   close OUTFILE2;
  }
 } else {
  open (OUTFILE, '>>', $locusname);
   print OUTFILE "$line\n";
  close OUTFILE;
  if ($newheader =~ m/_Allele_0$/) {
   open (OUTFILE1, '>>', $locusname_Allele0);
    print OUTFILE1 "$line\n";
   close OUTFILE1;
  }
  if ($newheader =~ m/_Allele_1$/) {
   open (OUTFILE2, '>>', $locusname_Allele1);
    print OUTFILE2 "$line\n";
   close OUTFILE2;
  }
 }
}
close FILE;
