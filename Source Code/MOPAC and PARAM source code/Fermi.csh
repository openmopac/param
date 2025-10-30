cd /home/jstewart/MOPAC2016GPU
source  /opt/intel/bin/compilervars.sh intel64
if [ -f /opt/mopac/MOPAC2016.exe ];
then
 rm /opt/mopac/MOPAC2016.exe
fi
rsync -v -u -r "/media/Software/PARAM/Source Code/MOPAC and PARAM source code/MOPAC_Source_code/" /home/jstewart/MOPAC_Source_code
rsync -v -u -r "/media/Software/PARAM/Source Code/MOPAC and PARAM source code/MOPAC_CUDA_files/" /home/jstewart/MOPAC_CUDA_files
make -f "/media/Software/PARAM/Source Code/MOPAC and PARAM source code/MOPAC_Source_code/Makefile_for_GPU_MOPAC"
cp /opt/mopac/MOPAC2016.exe  /media/Software/Utility/Linux_64_bit_GPU_MOPAC2016.exe
cd ~/work
/opt/mopac/MOPAC2016.exe Crambin_1SCF
cp Crambin_1SCF.out /media/Working/Crambin_1SCF_GPU.out
todos    /media/Working/Crambin_1SCF_GPU.out
objdump -f /opt/mopac/MOPAC2016.exe | grep ^archit
file  /opt/mopac/MOPAC2016.exe
