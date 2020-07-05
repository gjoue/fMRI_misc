find . -type f -name "*.txt" -exec sh -c 'mv "$0" "${0%.*}.JUL.ErC_R.txt"' {} \;

