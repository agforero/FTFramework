if [ $# -ne 1 ]; then
    echo "Usage: ./addTest.sh <directory>"
    exit 1
fi

cd source
if [ ! -d $1 ]; then
    echo "/$1/ does not exist in /source/"
    exit 2
fi

echo "Enter preemptive environmental commands to be included in .bats test, e.g. \`module load gcc\`."
echo -e "Enter \`ctrl + D\` to conclude your input.\n"
echo "---------------------------------------- BEGIN INPUT ----------------------------------------"
echo "# Environmental commands:" >> test.bats	
while read comm
do
    echo "$comm" >> test.bats
done

echo -e "----------------------------------------- END INPUT -----------------------------------------\n"
echo >> test.bats
echo "# Generic Makefile test" >> test.bats
echo "@test \"$d make\" {" >> test.bats
echo "  make" >> test.bats
echo "}" >> test.bats
echo "test.bats generated in $(pwd)."