#! /usr/bin/perl
# arguments should be supplied in order: 0. if there are four arguments, the first one will be the input file - otherwise stdin is used. 1. <START|BEGIN>, for MySQL/SQLite resp., 2. table name, 3. rows per insert
$a2= pop @ARGV;
$a1= pop @ARGV;
$a0= pop @ARGV;
if (scalar @ARGV > 0) {
	open($f, "<", $ARGV[0]);
} else {
	$f = STDIN;
}
print "$a0 TRANSACTION;\n";until(scalar (my @details = grep { defined } map { scalar <$f> } 1 .. $a2) eq 0) { print "INSERT INTO `$a1` VALUES \n", join(",\n", map { chomp; "\t(".$_.")" } @details), ";\n"; } print "COMMIT;\n"
