pyi-makespec \
	--name RetSynth \
	--icon='rs.icns' \
	--windowed \
	--add-binary='/usr/local/lib/python3.7/site-packages/libsbml/_libsbml.cpython-37m-darwin.so':'.' \
	--add-binary='/usr/local/Cellar/glpk/4.65/bin/glpsol':'.' \
	--add-binary='indigopython130_mac/lib/Mac/10.7/lib*':'indigopython130_mac/lib/Mac/10.7' \
	--add-data='ConstructedDatabases/DBINCHIECOLIDH1*.tar.gz':'ConstructedDatabases' \
	--add-data='Database/data/*':'Database/data' \
	--add-data='GeneCompatibility/BLAST/ncbi-blast-2.9.0+_mac/bin/*':'GeneCompatibility/BLAST/ncbi-blast-2.9.0+_mac/bin' \
	--add-data='GeneCompatibility/retrievegeneseqs/data/*_defaultdb*':'GeneCompatibility/retrievegeneseqs/data' \
	--add-data='GeneCompatibility/NCBI_SSU/data/*_defaultdb*':'GeneCompatibility/NCBI_SSU/data' \
	--add-data='GeneCompatibility/NCBI_SSU/data/*.zip':'GeneCompatibility/NCBI_SSU/data' \
	--add-data='GeneCompatibility/BLAST/blastdbs/*_defaultdb*':'GeneCompatibility/BLAST/blastdbs' \
	--add-data='GeneCompatibility/D_Tailor/CAI_Tables/*':'GeneCompatibility/D_Tailor/CAI_Tables' \
	rs_gui.py \

gsed -i 's|block_cipher = None|import sys'$"\n"'sys.setrecursionlimit(5000)'$"\n\n"'block_cipher = None|' RetSynth.spec \

pyinstaller \
	--clean \
	--noconfirm \
	--distpath dist \
	--workpath build \
	RetSynth.spec
