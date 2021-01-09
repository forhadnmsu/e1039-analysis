#!/bin/sh

#awk '($2=="H1XT")|| ($2=="H1XB")||($2=="H2XT")|| ($2=="H2XB")||($2=="H3XT")|| ($2=="H3XB") ||($2=="H4XT-u")||($2=="H4XT-d")||($2=="H4XB-u")||($2=="H4XB-d"){print $2$3" "$4}' /data2/e1039/daq/slowcontrols/lecroy/settings/current/H1 /data2/e1039/daq/slowcontrols/lecroy/settings/current/H2 /data2/e1039/daq/slowcontrols/lecroy/settings/current/H3 > outputx.txt

awk '($2=="H1YL")|| ($2=="H1YR")||($2=="H2YL")|| ($2=="H2YR")||($2=="H4Y1L-l")|| ($2=="H4Y1L-r") ||($2=="H4Y1R-l")||($2=="H4Y1R-r")||($2=="H4Y2L-l")||($2=="H4Y2L-r")||($2=="H4Y2R-l")||($2=="H4Y2R-r"){print $2$3" "$4}' /data2/e1039/daq/slowcontrols/lecroy/settings/current/H1 /data2/e1039/daq/slowcontrols/lecroy/settings/current/H2 /data2/e1039/daq/slowcontrols/lecroy/settings/current/H3 > outputy.txt
