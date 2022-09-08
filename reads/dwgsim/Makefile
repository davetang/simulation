all: test/hg38.eval.txt.gz

.PHONY: all clean

data/hg38.fa:
	data/download.sh

bin/samtools:
	bin/install_samtools.sh

bin/minimap2:
	bin/install_minimap2.sh

bin/dwgsim:
	bin/install_dwgsim.sh

test/hg38.bwa.read1.fastq.gz: bin/dwgsim data/hg38.fa
	bin/sim.sh

test/hg38.bam: test/hg38.bwa.read1.fastq.gz bin/minimap2 bin/samtools
	bin/map.sh

test/hg38.eval.txt.gz: test/hg38.bam
	bin/eval.sh

clean:
	rm -rf test/hg38.* bin/minimap2* bin/samtools* bin/dwgsim* bin/DWGSIM
