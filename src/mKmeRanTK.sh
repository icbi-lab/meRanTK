#!/bin/sh


/usr/local/bioinf/perl/perl-5.24.3/bin/pp -B -M Data::Dumper -M diagnostics -M GD -M GD::Text -M GD::Text::Align -M GD::Graph::lines -l /usr/lib64/libgd.so.2 -c -vv -o meRanGs meRanGs.pl
/usr/local/bioinf/perl/perl-5.24.3/bin/pp -B -M Data::Dumper -M diagnostics -M GD -M GD::Text -M GD::Text::Align -M GD::Graph::lines -l /usr/lib64/libgd.so.2 -c -vv -o meRanGh meRanGh.pl
/usr/local/bioinf/perl/perl-5.24.3/bin/pp -B -M Data::Dumper -M diagnostics -M GD -M GD::Text -M GD::Text::Align -M GD::Graph::lines -l /usr/lib64/libgd.so.2 -c -vv -o meRanT meRanT.pl
/usr/local/bioinf/perl/perl-5.24.3/bin/pp -B -M Data::Dumper -l /lib64/libm.so.6 -c -vv -o meRanCall meRanCall.pl
/usr/local/bioinf/perl/perl-5.24.3/bin/pp -B -M Data::Dumper -l /lib64/libm.so.6 -c -vv -o meRanCompare meRanCompare.pl
/usr/local/bioinf/perl/perl-5.24.3/bin/pp -B -M Data::Dumper -M MCE::Loop -c -vv -o meRanAnnotate meRanAnnotate.pl
