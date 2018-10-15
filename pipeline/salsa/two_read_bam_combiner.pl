#!/usr/bin/perl

use strict;

MAIN : {

    my ($read1_bam, $read2_bam) = @ARGV;
    
	if ((not defined $read1_bam) ||
	(not defined $read2_bam)) {
	die ("Usage: ./two_read_bam_combiner.pl <read 1 bam> <read 2 bam>\n");
    }
    
    open(FILE1,"samtools view -h $read1_bam |");
    open(FILE2,"samtools view -h $read2_bam |");
	
    my $line1 = <FILE1>;
    my $line2 = <FILE2>;

    my $counter = 0;
    my $new_counter = 0;

    while (defined $line1) {
	## the line directly below this needed to be modified slightly from the genotyping pipeline ##
	if ($line1 =~ /^(\@)SQ/){
		if ($line1 ne $line2){print $line1;print $line2; die ("inconsistent bam header.");}
		else{
			print $line1;
		}
		$line1 = <FILE1>;
		$line2 = <FILE2>;
		next;
	}

	$counter++;
	if ($counter == ($new_counter + 1000000)) {
	    print STDERR $counter . "\n";
	    $new_counter = $counter;
	}

	chomp $line1;
	chomp $line2;
	
	my ($id1, $flag1, $chr_from1, $loc_from1, $mapq1, $cigar1, $d1_1, $d2_1, $d3_1, $read1, $read_qual1, @rest1) = split(/\t/,$line1);
	my ($id2, $flag2, $chr_from2, $loc_from2, $mapq2, $cigar2, $d1_2, $d2_2, $d3_2, $read2, $read_qual2, @rest2) = split(/\t/,$line2);

	if ($id1 ne $id2) {
	    die ("The id's of the two files do not match up at line number $counter\n");
	}
	
	my $bin1 = reverse(dec2bin($flag1));
	my $bin2 = reverse(dec2bin($flag2));

	my @binary1 = split(//,$bin1);
	my @binary2 = split(//,$bin2);

	my $trouble = 0;
	if (($binary1[2] == 1) || ($mapq1 < 10)) {
	    $trouble = 1;
	}
	if (($binary2[2]== 1) || ($mapq2 < 10)) {
            $trouble = 1;
    }

	my $proper_pair1;
	my $proper_pair2;
	my $dist1;
	my $dist2;

	if (($binary1[2] == 0) && ($binary2[2] == 0)) {
	    $proper_pair1 = 1;
	    $proper_pair2 = 1;
	    if ($chr_from1 eq $chr_from2) {
		my $dist = abs($loc_from1 - $loc_from2);
		if ($loc_from1 >= $loc_from2) {
		    $dist1 = -1*$dist;
		    $dist2 = $dist;
		} else {
		    $dist1 = $dist;
		    $dist2 = -1*$dist;
		}
	    } else {
		$dist1 = 0;
		$dist2 = 0;
	    }
	} else {
	    $proper_pair1 = 0;
            $proper_pair2 = 0;
	    $dist1 = 0;
	    $dist2 = 0;
	}

	my $new_bin1 = join("","000000000000000000000",$binary1[10],$binary1[9],"0","0","1",$binary2[4],$binary1[4],$binary2[2],$binary1[2],$proper_pair1,"1");
	my $new_bin2 = join("","000000000000000000000",$binary2[10],$binary2[9],"0","1","0",$binary1[4],$binary2[4],$binary1[2],$binary2[2],$proper_pair2,"1");

	my $new_flag1 = bin2dec($new_bin1);
	my $new_flag2 = bin2dec($new_bin2);

	unless ($trouble == 1) {

	    print(join("\t",$id1,$new_flag1,$chr_from1,$loc_from1,$mapq1,$cigar1,$chr_from2,$loc_from2,$dist1,$read1,$read_qual1,@rest1) . "\n");
	    print(join("\t",$id2,$new_flag2,$chr_from2,$loc_from2,$mapq2,$cigar2,$chr_from1,$loc_from1,$dist2,$read2,$read_qual2,@rest2) . "\n");

	}

	$line1 = <FILE1>;
	$line2 = <FILE2>;

    }

}



sub dec2bin {

    my $str = unpack("B32", pack("N", shift));
    return $str;

}

sub bin2dec {
    return unpack("N", pack("B32", substr("0" x 32 . shift, -32)));
}
    
