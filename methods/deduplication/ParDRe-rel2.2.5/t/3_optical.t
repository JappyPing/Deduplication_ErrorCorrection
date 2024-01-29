# perl test script; generally run using "prove t/1_basic.t" after make
# or just run "make test"

use strict;
use Test::More tests => 13;

require "t/test_subs.pl";
my $EMPTY = q{};

my $prog = "./ParDRe";
reset_output_dir();

my ($out, $err, $rc) = run("$prog -d 2500 -i t/in/100.1.gz -o t/out/100.1.gz -z 1");
is $rc, 0, 'single read - no error code';
is $out, $EMPTY, 'single read - no stdout';
my $expect = <<'END';
Process 0/1: Initialized
Command: ./ParDRe -d 2500 -i t/in/100.1.gz -o t/out/100.1.gz -z 1 
FASTQ format identified
Process 0/1: To print in t/out/100.1.gz
Process 0/1: finished the clustering of 25 reads in [\d.]+ seconds
Process 0/1: Time to analyze 25 reads with 1 threads: [\d.]+ seconds
SUMMARY: Non duplicated paired reads 24/25 \(96.00 %\)
SUMMARY: Min/Max analyzed reads per process: 25/25 \(0.00 % of imbalance\)
SUMMARY: Min/Max runtime per process: [\d.]+/[\d.]+ \([\d.]+ % of imbalance\)
Process 0/1: Overall time: [\d.]+ seconds
END
like $err, qr/$expect/, 'single read - stderr';
my $size_in = -s "t/in/100.1.gz";
my $size = -s "t/out/100.1.gz";
ok $size > 0 && $size < $size_in, 'single read - output size is approx right';

# same thing, fast compression level
($out, $err, $rc) = run("$prog -d 2500 -i t/in/100.1.gz -o t/out/100.1.gz -z 2");
is $rc, 0, 'single read fast compress - no error code';
is $out, $EMPTY, 'single read fast compress - no stdout';
$expect =~ s{-z 1}{-z 2};
like $err, qr/$expect/, 'single read fast compress - stderr';
$expect =~ s{-z 2}{-z 1};
ok -s "t/out/100.1.gz" > $size, 'single read fast compress - larger';

($out, $err, $rc) = run("$prog -d 2500 -i t/in/100.1.gz -p t/in/100.2.gz -o t/out/100.1.gz  -r t/out/100.2.gz -z 1");
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

my $in = sorted_deflines("t/in/100.2.gz");
$out = sorted_deflines("t/out/100.2.gz");
is $out, $in, 'paired read - same read deflines';

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
