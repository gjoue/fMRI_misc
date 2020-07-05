find . -type f -name "*str2match*" -exec bash -c 'mv $0 ${0/match/newstring}' {} \;

