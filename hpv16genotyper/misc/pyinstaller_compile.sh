pyinstaller --distpath HPV16-GenotyperUb20 \
	--workpath tmpBuild \
	--name  HPV16-Genotyper \
	--add-data PathTo/resources:resources/ \
	--add-data PathTo/files_rc.py:. \
	--add-data PathTo/simplot.py:. \
	--add-data PathTo/appFunctions.py:. \
	-d all \
	--windowed \
	../app.py;

rm -r tmpBuild;
mv ../*.pyo HPV16-GenotyperUb20/HPV16-Genotyper/;
mv *.spec HPV16-GenotyperUb20/HPV16-Genotyper/;
