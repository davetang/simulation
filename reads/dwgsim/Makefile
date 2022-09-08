all: test/hg38.eval.txt.gz

.PHONY: all clean

data/hg38.fa: data/download.sh
	$<

bin/samtools: bin/install_samtools.sh
	$<

bin/minimap2: bin/install_minimap2.sh
	$<

bin/dwgsim: bin/install_dwgsim.sh
	$<

test/hg38.bwa.read1.fastq.gz: bin/sim.sh bin/dwgsim data/hg38.fa
	$<

test/hg38.bam: bin/map.sh test/hg38.bwa.read1.fastq.gz bin/minimap2 bin/samtools
	$<

test/hg38.eval.txt.gz: bin/eval.sh test/hg38.bam
	$<

clean:
	rm -rf test/hg38.* bin/minimap2* bin/samtools* bin/dwgsim* bin/DWGSIM
