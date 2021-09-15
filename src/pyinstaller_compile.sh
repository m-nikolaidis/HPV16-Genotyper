pyinstaller --distpath HPV16genotyper \
	--workpath tmpBuild \
	--name  HPV16genotyper \
	--add-data /media/marios/Data/Lab/HPV16/HPV16_genotyping_tool/src/resources:resources/ \
	--add-data /media/marios/Data/Lab/HPV16/HPV16_genotyping_tool/src/annot.py:. \
	--add-data /media/marios/Data/Lab/HPV16/HPV16_genotyping_tool/src/blast.py:. \
	--add-data /media/marios/Data/Lab/HPV16/HPV16_genotyping_tool/src/env_setup.py:. \
	--add-data /media/marios/Data/Lab/HPV16/HPV16_genotyping_tool/src/files_rc.py:. \
	--add-data /media/marios/Data/Lab/HPV16/HPV16_genotyping_tool/src/phylogeny.py:. \
	--add-data /media/marios/Data/Lab/HPV16/HPV16_genotyping_tool/src/simplot.py:. \
	--add-data /media/marios/Data/Lab/HPV16/HPV16_genotyping_tool/src/wrapper.py:. \
	-d all \
	--windowed \
	app.py;

rm -r tmpBuild;
mv *.pyo HPV16genotyper/;
mv *.spec HPV16genotyper/;