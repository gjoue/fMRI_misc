# takes 1st folder = sub*, strips "sub" and prints gaps in numerical sequence
 find sub*/func/reinfL/run1/ -name PhysIO_output.ps | sort | awk -F "/" '{print $1}' | sed -e "s/sub//" | awk 'p && p != $1 { for( i = p; i < $1; i++ ) print "missing " i; } {p = $1 + 1 }'


