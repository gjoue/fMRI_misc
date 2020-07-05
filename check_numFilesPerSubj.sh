find . -maxdepth 4 -mindepth 3 -type d -exec sh -c "echo {}; ls -1 {} | sort | wc -l" \; | xargs -n 2 | awk '{print $2"\t"$1}'   > filecheck.out
