#!/bin/sh

awk '($2=="H1XT")|| ($2=="H1XB")||($2=="H2XT")|| ($2=="H2XB")||($2=="H3XT")|| ($2=="H3XB") ||($2=="H4XT-u")||($2=="H4XT-d")||($2=="H4XB-u")||($2=="H4XB-d"){print $2$3" "$4}' /data2/e1039/daq/slowcontrols/lecroy/settings/current/H* > output.txt
#awk '($2=="H1XT")|| ($2=="H1XB")||($2=="H2XT")|| ($2=="H2XB")||($2=="H3XT")|| ($2=="H3XB") ||($2=="H4XT-u")||($2=="H4XT-d")||($2=="H4XB-u")||($2=="H4XB-d"){print $2" "$4}' /data2/e1039/daq/slowcontrols/lecroy/settings/current/H* > output.txt

