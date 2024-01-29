# subroutines useful for running tests

use strict;
use File::Temp qw(tempfile);

my $EMPTY = q{};
our ($Gzip, $GzipIn, $GzipOut);
$Gzip = $GzipIn = $GzipOut = 'gzip';

sub run {
    my ($cmd, $input) = @_;
    my $infile = $EMPTY;
    if (defined $input) {
        my $fh;
        ($fh, $infile) = tempfile();
        print $fh $input;
    }
    my ($fh, $stderr) = tempfile();
    close $fh;
    print "Running: $cmd $infile\n";
    my $output = `$cmd $infile 2>$stderr`;
    my $rc = $?;
    my $error = $EMPTY;
    local *ERROR1;
    if (open ERROR1, $stderr) {
        local $/;
        $error = <ERROR1>;
        close ERROR1;
    }
    unlink $stderr;
    if ($infile) {
        unlink $infile;
    }
    return wantarray ? ($output, $error, $rc) : $output;
}

sub today {
    my @lt = localtime;
    return sprintf '%02d-%02d-%d', $lt[3], $lt[4]+1, $lt[5]+1900;
}

sub write_file {
    my ($file, $data) = @_;
    my $fh = Open($file, "w");
    print $fh $data;
    close $fh or die "Failed writing to '$file', $!\n";
}

sub read_file {
    my ($file) = @_;
    my $fh = Open($file);
    local $/;
    my $data = <$fh>;
    return $data;
}

sub reset_output_dir {
    my $dir = "t/out";
    if (-d $dir) {
        system "rm -rf t/out";
    }
    mkdir $dir;
    return $dir;
}

# cut and pasted from GTB::File, to avoid having non-core dependencies
# then: s/croak/die/ (all of the "or die" were originally die, others croak)
sub Open {
    my ($file, $mode) = @_;
    if (!$file && $file ne '0') {
        die "Open: no filename provided";
    }
    if ($mode) {
        $mode = lc $mode;
    }
    elsif ($file =~ /^\s*\|/) {
        $mode = 'w';
    }
    else {
        $mode = 'r';
    }
    my $fh;
    if ($file =~ /\|/) {
        if ($mode eq 'r') {
            if ($file =~ /\|\s*$/) {
                open $fh, $file or die "Can't open pipe '$file', $!\n";
            }
            else {
                die "To open pipe for reading, pipe character must "
                    . "appear at end of command";
            }
        }
        elsif ($mode eq 'w') {
            if ($file =~ /^\s*\|/) {
                open $fh, $file or die "Can't open pipe '$file', $!\n";
            }
            else {
                die "To open pipe for writing, pipe character must "
                    . "appear at beginning of command";
            }
        }
        else { # pipe, but not first or last in sequence
            die << "END_MSG";
If a pipe character is present in the open string, there must be a pipe at
the beginning or end of the string, depending upon whether you plan to
write or read to the filehandle; '$file' is not valid.  If you need to read
and write to a program, try IPC::Open2 or IPC::Open3.
END_MSG
        }
    }
    elsif ($file =~ /\.(b?gz|bz2|zip|Z)$/) {
        if ($mode eq 'r') {
            my $prog = $1 eq 'bz2' ? 'bzip2' : ($GzipIn || $Gzip);
            die "File ($file) not found" unless (-e $file);
            die "File ($file) was not readable" unless (-r $file);
            open $fh, "$prog -dc $file |"
                or die "Can't read $file with $prog, $!\n";
        }
        elsif ($mode eq 'w') {
            my $prog = $1 eq 'bz2' ? 'bzip2' : ($GzipOut || $Gzip);
            open $fh, "| $prog > $file"
                or die "Can't create $prog file $file, $!\n";
        }
        elsif ($mode eq 'a') {
            if ($1 eq 'bz2') {
                die "Open: mode 'a' not supported for bzip2 file $file";
            }
            my $prog = $GzipOut || $Gzip;
            open $fh, "| $prog >> $file"
                or die "Can't append $prog output to $file, $!\n";
        }
        else {
            die "Open: mode '$mode' not supported; use 'r', 'w' or 'a'";
        }
    }
    elsif ($file eq '-') {
        if ($mode eq 'r') {
            open $fh, '-' or die "Can't read from STDIN, $!\n";
        }
        elsif ($mode eq 'w' || $mode eq 'a') {
            open $fh, '>-' or die "Can't write to STDOUT, $!\n";
        }
        else {
            die "Open: mode '$mode' not supported; use 'r', 'w' or 'a'";
        }
    }
    elsif ($mode eq 'r') {
        open $fh, '<', $file or die "Can't open $file, $!\n";
    }
    elsif ($mode eq 'w') {
        open $fh, '>', $file or die "Can't create $file, $!\n";
    }
    elsif ($mode eq 'a') {
        open $fh, '>>', $file or die "Can't append to $file, $!\n";
    }
    else {
        die "Open: mode '$mode' not supported; use 'r', 'w' or 'a'";
    }
    return $fh;
}

1;
