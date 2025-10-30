rsync -v -u  "/media/psf/Software_on_K20/PARAM/Source code/MOPAC and PARAM source code/MOPAC_source_code/"*.F90 ~/PARAM/PARAM_Source_code  
rsync -v -u  "/media/psf/Software_on_K20/PARAM/Source code/MOPAC and PARAM source code/PARAM_source_code/"*.F90 ~/PARAM/PARAM_Source_code  
rsync -v -u  "/media/psf/Software_on_K20/PARAM/Source code/MOPAC and PARAM source code/Makefile_for_PARAM_Linux.txt"  ~/PARAM 
rsync -v -u  "/media/psf/Software_on_K20/PARAM/Source code/MOPAC and PARAM source code/PARAM_Makefile_files.txt"  ~/PARAM
cd ~/PARAM
source  /opt/intel/bin/compilervars.sh intel64
make -f Makefile_for_PARAM_Linux.txt 
