extract_pdbs.mpi.linuxgccrelease -in::file::silent /home/ictsun/dyf/rosetta/result/1e74.A.silent -in:auto_setup_metals 
for file in result_*.pdb; do mv -v "$file" "${file/result/1e74.A}"; done;
extract_pdbs.mpi.linuxgccrelease -in::file::silent /home/ictsun/dyf/rosetta/result/2ifi.A.silent -in:auto_setup_metals 
for file in result_*.pdb; do mv -v "$file" "${file/result/2ifi.A}"; done;
extract_pdbs.mpi.linuxgccrelease -in::file::silent /home/ictsun/dyf/rosetta/result/2m3i.A.silent -in:auto_setup_metals 
for file in result_*.pdb; do mv -v "$file" "${file/result/2m3i.A}"; done;
extract_pdbs.mpi.linuxgccrelease -in::file::silent /home/ictsun/dyf/rosetta/result/6ota.A.silent -in:auto_setup_metals 
for file in result_*.pdb; do mv -v "$file" "${file/result/6ota.A}"; done;
extract_pdbs.mpi.linuxgccrelease -in::file::silent /home/ictsun/dyf/rosetta/result/6kmy.B.silent -in:auto_setup_metals 
for file in result_*.pdb; do mv -v "$file" "${file/result/6kmy.B}"; done;
extract_pdbs.mpi.linuxgccrelease -in::file::silent /home/ictsun/dyf/rosetta/result/6kn3.A.silent -in:auto_setup_metals 
for file in result_*.pdb; do mv -v "$file" "${file/result/6kn3.A}"; done;