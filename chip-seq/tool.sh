#！ /bin/bash
mamba create -n chipseq python=3.9
mamba activate chipseq
mamba install -y sra-tools trim-galore samtools deeptools homer meme macs2 bowtie bowtie2 