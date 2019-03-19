#!/software/bin/perl -w

use strict;
use Bio::SeqIO;
use DateTime;
$|=1;

my $file = shift;
my $outfile = shift;
my $seqin  = Bio::SeqIO->new('-format' => 'fasta',
                             '-file'   => $file);
$outfile = "./trimreport.txt" unless ($outfile);

# parameters
my $minleftover = 100;   # after trimmning start/end, at least this many bp should be left
my $maxperc     = 0.5;   # the sequence should contain less than this fraction of Ns in total
my $winsize     = 5000;  # for sliding window analysis
my $maxslidingN = 0.6;   # maximum fraction of Ns in sliding window before alarm sets off

# START

open(OUT,">${outfile}");
while (my $seqobj = $seqin->next_seq()) {    

    # total N count
    my $count = () = $seqobj->seq() =~ /N|n/g;
    my $Nperc = $count/$seqobj->length;
    print OUT "# WARNING: ",$seqobj->id(),"\t", int($Nperc*10000)/100," % Ns of total ",$seqobj->length,"\n" if ($Nperc > 0.8);
    
    # Ns at start/end
    my ($startseq) = (uc($seqobj->seq()) =~ /^(N+)/);
    my ($endseq) = (uc($seqobj->seq()) =~ /(N+)$/);
    if(!defined ($startseq)) {$startseq = ''}
    if(!defined ($endseq)) {$endseq = ''}
    my $realseq = substr(uc($seqobj->seq()), length($startseq), length($seqobj->seq()) - length($startseq) - length($endseq) );

    if ((length($startseq) > 0) || (length($endseq) > 0)) {
        if (length($startseq) > 0) {
            print OUT join "\t",("TRIM:",$seqobj->id(),1,length($startseq)), "\n";
        }
        if (length($endseq) > 0) {
            print OUT join "\t",("TRIM:",$seqobj->id(),length($startseq)+length($realseq)+1,length($seqobj->seq)),"\n";
        }
    }
    print OUT "REMOVE: ",$seqobj->id(),"\t",length($realseq)," bp leftover after trimming\n" if (length($realseq) <= $minleftover);

    my @seq = split(//, $seqobj->seq);

    # sliding window N analysis fwd
    my $Ncount;
    for(my $i = 1; $i <= $seqobj->length() - ($winsize-1); $i++) {
        if($i == 1) {
            my $window = $seqobj->subseq($i,$i+($winsize-1));
            $Ncount = () = $window =~ /N|n/g;
            last if ($Ncount/$winsize < $maxslidingN);
        }
        else {
            if($seq[$i-2] =~ /N/i) {$Ncount--}
            if(!defined($seq[$i+$winsize-3]) or $seq[$i+$winsize-3] =~ /N/i) {$Ncount++}
        }

        
        if (($i>1) && ($Ncount/$winsize < $maxslidingN)) {
            my $window = $seqobj->subseq($i,$i+($winsize-1));
            my ($clipseq,$realseq) = (uc($window) =~ /(N*)(\S*)$/);
            
            print OUT "FWDCLIP:\t",$seqobj->id(),"\t1\t",$i - 1 + length($clipseq) ,"\n";
            last;
        }
    }
    
    # sliding window N analysis rev
    my $revobj = $seqobj->revcom;
    my @revseq = split(//, $revobj->seq);
    for(my $i = 1; $i <= $revobj->length() - ($winsize-1); $i++) {
        if($i == 1) {
            my $window = $revobj->subseq($i,$i+($winsize-1));
            $Ncount = () = $window =~ /N|n/g;
            last if ($Ncount/$winsize < $maxslidingN);
        }
        else {
            if($revseq[$i-2] =~ /N/i) {$Ncount--}
            if(!defined($revseq[$i+$winsize-3]) or $revseq[$i+$winsize-3] =~ /N/i) {$Ncount++}
        }

        if (($i>1) && ($Ncount/$winsize < $maxslidingN)) {
            my $window = $revobj->subseq($i,$i+($winsize-1));
            my ($clipseq,$realseq) = (uc($window) =~ /(N*)(\S*)$/);
            print OUT "REVCLIP:\t",$seqobj->id(),"\t",$seqobj->length - $i + 2 - length($clipseq),"\t",$seqobj->length,"\n";
            last;
        }
    }

    undef $seqobj;
}