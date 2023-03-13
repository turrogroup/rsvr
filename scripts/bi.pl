#!/usr/bin/perl
my %h = (Q=>8,L=>4,f=>2,d=>4,S=>2,C=>1);
my @codes=split //, $ARGV[1];
my @lens=map { $h{$_} } @codes;
if ($ARGV[0] eq 'u') {
	my $v;
	binmode(STDIN); 
	while (read(STDIN, $v, $lens[0])) {
		print (unpack($codes[0], $v));
		if ((scalar @codes) > 1) {
			for (1..$#lens) {
				read(STDIN, $v, $lens[$_]); 
				print "\t" . unpack($codes[$_], $v);
			}
		}
		print "\n";
	}
} else {
	while (<STDIN>) {
		chomp;
		@F=split /[^\d]/;
		for (0..$#lens) {
			print pack($codes[$_], $F[$_]);
		}
	}
}
