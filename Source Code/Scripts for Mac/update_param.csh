rsync -v -u '/Volumes/Software_on_K20/PARAM/Source Code/MOPAC and PARAM source code/MOPAC_Source_code/'*.F90 ~/PARAM/PARAM_Source_code  
rsync -v -u '/Volumes/Software_on_K20/PARAM/Source Code/MOPAC and PARAM source code/PARAM_Source_code/'*.F90 ~/PARAM/PARAM_Source_code  
rsync -v -u '/Volumes/Software_on_K20/PARAM/Source Code/MOPAC and PARAM source code/Makefile_for_PARAM_Mac.txt' ~/PARAM
rsync -v -u '/Volumes/Software_on_K20/PARAM/Source Code/MOPAC and PARAM source code/PARAM_Makefile_files.txt'  ~/PARAM
cd ~/PARAM
make -f Makefile_for_PARAM_Mac.txt 
