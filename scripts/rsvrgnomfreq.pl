#! /usr/bin/perl 

if ($ARGV[0] eq "v37") {
	$par = (2699520, 154931044);
} elsif ($ARGV[0] eq "v38") {
	$par = (2781479, 155701383);
} else {
	die("Unrecognised genome version")
}

while (<STDIN>) {
	chomp;
	s/^chr//;
	@F=split /\t/;
	@anno=split /\|/, $F[8];
	$inpar = ($F[2] < $par[0]) or ($F[2] > $par[1]);
	$suf= ($F[1] eq "X" and not $inpar) ? "_[mM]ale" : "";
	$id = $F[0];
	@pops=(/AC_([a-zA-Z]{3})$suf=\d+/g);
	@acs=(/AC_([a-zA-Z]{3})$suf=(\d+)/g);
	@ans=(/AN_([a-zA-Z]{3})$suf=(\d+)/g);
	while (@acs) { $acm{pop @acs} = pop @acs }
	while (@ans) { $anm{pop @ans} = pop @ans }
	for (@pops) {
		if ($anm{$_} > 0) { print "$id $acm{$_} $anm{$_} $_\n" if ($acm{$_} > 0) }
	}
}


