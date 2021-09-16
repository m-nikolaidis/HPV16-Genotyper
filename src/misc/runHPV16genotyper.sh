args=("$@")
if [ "${args[1]}" != "" ]; then
echo "The script takes at most 1 argument"
exit 1
fi

libflt=$(ls /usr/lib/x86_64-linux-gnu/|grep "libfltk_images.so.1.3")

if [ "${args[0]}" == "-i" ]; then
	sudo apt-get install libfltk1.3*;
	exit 0
fi

if [ "${args[0]}" == "-h" ]; then
	echo "HPV16genotyper basic usage"
	echo ""
	echo "Syntax: ./runHPV16genotyper [-h|-d|-i]"
	echo ""
	echo "options:"
	echo " -d   DEBUG   When used the default terminal output is not silenced"
	echo " -h   HELP    Display this message and exit"
	echo " -i   INSTALL Installs the required library (libfltk1.3*) (Needs sudo priviledges)"
	exit 0
fi

if ["$libflt" == ""]; then
	echo "Error libfltk_images.so.1.3 is not found in /usr/lib/x86_64-linux-gnu/";
	echo "";
	echo "You need to install it manually or with the help of this script"
	echo "Use the -h flag for more information"
	exit 1
fi

if [ "${args[0]}" == "-d" ]; then
	cd application
	./HPV16genotyper
fi
if [ "${args[0]}" == "" ]; then
	cd application
	./HPV16genotyper &>/dev/null
fi