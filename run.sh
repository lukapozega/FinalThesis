if [[ "$OSTYPE" == "linux-gnu" ]]; then
        type=linux_64
elif [[ "$OSTYPE" == "darwin"* ]]; then
	type=MacOS
fi
if [[ $1 != *.fa ]]
then
	echo "Reference file should be uncompressed and have .fa extension"
	exit
fi
file=${1%.*} 
./external_tools/minimap2/$type/minimap2 $1 $2 > build/bin/overlaps.paf
mkdir external_tools/RED/$type/src
mkdir external_tools/RED/$type/des
mv $1 external_tools/RED/$type/src
./external_tools/RED/$type/Red -gnm external_tools/RED/$type/src -rpt external_tools/RED/$type/des
mv external_tools/RED/$type/src/$1 .
mv external_tools/RED/$type/des/$file.rpt build/bin
rm -r external_tools/RED/$type/src
rm -r external_tools/RED/$type/des
printf "\n\nStarting assembly\n"
./build/bin/assembly build/bin/overlaps.paf $1 build/bin/$file.rpt $2