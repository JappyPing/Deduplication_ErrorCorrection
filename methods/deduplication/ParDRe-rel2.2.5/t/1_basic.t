# perl test script; generally run using "prove t/1_basic.t" after make
# or just run "make test"

use strict;
use warnings;
use Test::More tests => 25;

require "t/test_subs.pl";
my $EMPTY = q{};

my $prog = "./ParDRe";
my ($out, $err, $rc) = run($prog);
like $err, qr/Usage: ParDRe/, 'program runs';
reset_output_dir();

my $in = read_file("t/in/100.1.gz");
write_file("t/out/100.1.fastq", $in);
ok -s "t/out/100.1.fastq", 'unzipped input';

($out, $err, $rc) = run("$prog -i t/out/100.1.fastq");
is $rc, 0, 'single read - no error code';
is $out, $EMPTY, 'single read - no stdout';
my $expect = <<'END';
Process 0/1: Initialized
Command: ./ParDRe -i t/out/100.1.fastq 
FASTQ format identified
Process 0/1: To print in t/out/100.1.fastq.NonDup_1
Process 0/1: finished the clustering of 25 reads in [\d.]+ seconds
Process 0/1: Time to analyze 25 reads with 1 threads: [\d.]+ seconds
SUMMARY: Non duplicated paired reads 24/25 \(96.00 %\)
SUMMARY: Min/Max analyzed reads per process: 25/25 \(0.00 % of imbalance\)
SUMMARY: Min/Max runtime per process: [\d.]+/[\d.]+ \([\d.]+ % of imbalance\)
Process 0/1: Overall time: [\d.]+ seconds
END
like $err, qr/$expect/, 'single read - stderr';
ok -s "t/out/100.1.fastq.NonDup_1", 'single read - output exists';

($out, $err, $rc) = run("$prog -i t/in/100.1.gz -o t/out/100.1.gz -z");
is $rc, 0, 'single gz read - no error code';
is $out, $EMPTY, 'single gz read - no stdout';
$expect = <<'END';
Process 0/1: Initialized
Command: ./ParDRe -i t/in/100.1.gz -o t/out/100.1.gz -z 
FASTQ format identified
Process 0/1: To print in t/out/100.1.gz
Process 0/1: finished the clustering of 25 reads in [\d.]+ seconds
Process 0/1: Time to analyze 25 reads with 1 threads: [\d.]+ seconds
SUMMARY: Non duplicated paired reads 24/25 \(96.00 %\)
SUMMARY: Min/Max analyzed reads per process: 25/25 \(0.00 % of imbalance\)
SUMMARY: Min/Max runtime per process: [\d.]+/[\d.]+ \([\d.]+ % of imbalance\)
Process 0/1: Overall time: [\d.]+ seconds
END
like $err, qr/$expect/, 'single gz read - stderr';
my $size_in = -s "t/in/100.1.gz";
my $size = -s "t/out/100.1.gz";
ok $size > 0 && $size < $size_in, 'single read - output size is approx right';

($out, $err, $rc) = run("$prog -i t/in/100.1.gz -p t/in/100.2.gz -o t/out/100.1.gz  -r t/out/100.2.gz -z");
is $rc, 0, 'paired read - no error code';
is $out, $EMPTY, 'paired read - no stdout';
$expect =~ s{Process 0/1: To print in t/out/100.1.gz}
            {FASTQ format identified\nTo print in t/out/100.1.gz and t/out/100.2.gz};
$expect =~ s{-o}{-p t/in/100.2.gz -o};
$expect =~ s{-z}{-r t/out/100.2.gz -z};
$expect =~ s{25 reads}{25 paired reads}g;
$expect =~ s{24}{25};
$expect =~ s{96}{100};
like $err, qr/$expect/, 'paired read - stderr';
$size_in = -s "t/in/100.2.gz";
$size = -s "t/out/100.2.gz";
ok $size > 0 && $size < $size_in * 1.01, 'paired read - output size';

$in = sorted_deflines("t/in/100.2.gz");
$out = sorted_deflines("t/out/100.2.gz");
is $out, $in, 'paired read - same read deflines';

# test output name generation
($out, $err, $rc) = run("$prog -i t/out/100.1.gz -p t/out/100.2.gz -z");
ok -s "t/out/100.1.NonDup_1.gz" && -s "t/out/100.2.NonDup_1.gz", "default zip name";
is $rc, 0, 'paired gz default files - no error code';
is $out, $EMPTY, 'paired gz default files - no stdout';
$expect =~ s{/in/}{/out/}g;
$expect =~ s{-o.+-z}{-z};
$expect =~ s{print in.+}{print in t/out/100.1.NonDup_1.gz and t/out/100.2.NonDup_1.gz};
like $err, qr/$expect/, "paired gz default files - stderr";

# even if start weird
unlink "t/out/100.1.NonDup_1.gz", "t/out/100.2.NonDup_1.gz";
system "mv t/out/100.1.gz t/out/100.1";
system "mv t/out/100.2.gz t/out/100.2";
($out, $err, $rc) = run("$prog -i t/out/100.1 -p t/out/100.2 -z");
is $rc, 0, 'weird gz default files - no error code';
$expect =~ s{100.([12]).gz}{100.$1}g;
like $err, qr/$expect/, "weird gz default files - stderr";

# test no-print option
unlink "t/out/100.1.NonDup_1.gz", "t/out/100.2.NonDup_1.gz";
($out, $err, $rc) = run("$prog -i t/out/100.1 -p t/out/100.2 -z -np");
is $rc, 0, 'no-print - no error code';
is $out, $EMPTY, 'no-print - no stdout';
ok !-f "t/out/100.1.NonDup_1.gz" && !-f "t/out/100.2.NonDup_1.gz",
   'no-print - no files';
$expect =~ s{-z}{-z -np};
like $err, qr/$expect/, 'no-print - stderr reports duplicates';

sub sorted_deflines {
    my ($file) = @_;
    my $fh = Open($file);
    my @dl;
    while (<$fh>) {
        push @dl, $_;
        <$fh>;
        <$fh>;
        <$fh>;
    }
    return join($EMPTY, sort @dl);
}
