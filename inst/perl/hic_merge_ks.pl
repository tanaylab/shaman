
use strict;

open(DATA1, $ARGV[0]) || die "Cannot open data1\n";
open(DATA2, $ARGV[1]) || die "Cannot open data2\n";

my($l1, $l2);
my($skip) = 0;
while(<DATA1>) {
	$l1 = $_;
	$l2 = <DATA2>;


	my(@x1) = split("\t", $l1);
	my(@x2) = split("\t", $l2);

	if($#x2 == -1) {
		if($skip != 0) {
			die "empty lines, terminating\n";
		}
		$skip++;
		next;
	}

	my($i1) = 1;
	my($i2) = 1;
	my($n1) = $#x1;
	my($n2) = $#x2;
	my($c1) = 1/$n1;
	my($c2) = 1/$n2;
	$x1[$n1] = $x2[$n2]+1;
	$x2[$n2] = $x1[$n2]+1;
	my($min_d) = 0;
	my($max_d) = 0;
	my($d) = 0;
	if($x2[2] < $x2[1]) {
		die "non monotonic distance sequence at x2!\n";
	}
	if($x1[2] < $x1[1]) {
		die "non monotonic distance sequence at x1!\n";
	}
	while($i1 < $n1 && $i2 < $n2) {
		if($x1[$i1] < $x2[$i2]) {
			$d += $c1;
			$i1++;
		} else {
			$d -= $c2;
			$i2++;
		}
		if($d > $max_d) {
			$max_d = $d;
		} elsif($d < $min_d) {
			$min_d = $d;
		}
	}
	$min_d = int($min_d*1000)/1000;
	$max_d = int($max_d*1000)/1000;
	print "$min_d\t$max_d\n";
}
