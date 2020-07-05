## search which Ss are missing whatever, collapsing across series

find sub*/GLM2/aligned*/r*_L/ -name con_0007.nii | sort | awk -F "/" '{print $1}' | sed -e "s/sub//" | awk 'p && p != $1 { for( i = p; i < $1; i++ ) print "missing " i; } {p = $1 + 1 }'
find sub*/GLM2/aligned*/r*_R/ -name con_0007.nii | sort | awk -F "/" '{print $1}' | sed -e "s/sub//" | awk 'p && p != $1 { for( i = p; i < $1; i++ ) print "missing " i; } {p = $1 + 1 }'

