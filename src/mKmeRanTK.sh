#!/bin/sh


pp -B -M Data::Dumper -M diagnostics -M GD -M GD::Text -M GD::Text::Align -M GD::Graph::lines -l /usr/lib64/libgd.so.2 -c -vv -o meRanGs meRanGs.pl
pp -B -M Data::Dumper -M diagnostics -M GD -M GD::Text -M GD::Text::Align -M GD::Graph::lines -l /usr/lib64/libgd.so.2 -c -vv -o meRanGt meRanGt.pl
pp -B -M Data::Dumper -M diagnostics -M GD -M GD::Text -M GD::Text::Align -M GD::Graph::lines -l /usr/lib64/libgd.so.2 -c -vv -o meRanGh meRanGh.pl
pp -B -M Data::Dumper -M diagnostics -M GD -M GD::Text -M GD::Text::Align -M GD::Graph::lines -l /usr/lib64/libgd.so.2 -c -vv -o meRanT meRanT.pl
pp -B -M Data::Dumper -l /lib64/libm.so.6 -c -vv -o meRanCall meRanCall.pl
pp -B -M Data::Dumper -l /lib64/libm.so.6 -c -vv -o meRanCompare meRanCompare.pl
pp -B -M Data::Dumper -M MCE::Loop -c -vv -o meRanAnnotate meRanAnnotate.pl
