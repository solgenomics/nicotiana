import os
import subprocess

cwd = os.getcwd()
processes = []
init = True
for (dirpath, dirnames, filenames) in os.walk(cwd):
    if init:
        init = False
        continue
    name = "SRR636" + dirpath[dirpath.rfind('/')+1:len(dirpath)]
    p = subprocess.Popen(f'htseq-count -t gene -i Name -f bam {dirpath}/Aligned.out.bam ~/work/2019_Summer/Nitab-v5_annotation/Nitab-v5_gene_models.gff --stranded=no -r name > {dirpath}/{name}.htseq', shell=True)
    processes.append(p)

for p  in processes:
    p.wait()


print('Done!')