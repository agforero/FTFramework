export OMP_PLACES=cores
export OMP_PROC_BIND=close

make
cd runs
mkdir $$
cd $$

# if no argument is provided
if [ $# -eq 0 ]; then
	threads=$(nproc)
	echo "running benchmark..."
	for i in $(seq 100000 100000 12500000)
	do
		../../nstream.exe $i $threads >> output.txt
	done
	echo "benchmark complete."

	echo "drawing graph..."
	cp ../../ograph.py .
	./ograph.py 
	rm ograph.py

	echo "done; results saved in /$$/"

elif [ $# -eq 1 ]; then
		echo "running benchmark..."
	for i in $(seq 100000 100000 12500000)
	do
		../../nstream.exe $i $2 >> output.txt
	done
	echo "benchmark complete."

	echo "drawing graph..."
	cp ../../ograph.py .
	./ograph.py 
	rm ograph.py

	echo "done; results saved in /$$/"

else
	echo "Usage: ./generate.sh <optional-threads>"
	exit 1
fi