all: bam eval
eval: test/hg38.mm.eval.txt.gz test/hg38.bwa.eval.txt.gz
bam: test/hg38.mm.bam test/hg38.bwa.bam

.PHONY: all clean

data/hg38.fa: data/download.sh
	$<

bin/samtools: bin/install_samtools.sh
	$<

bin/minimap2: bin/install_minimap2.sh
	$<

bin/dwgsim: bin/install_dwgsim.sh
	$<

bin/bwa: bin/install_bwa.sh
	$<

data/hg38.fa.sa: bin/bwa data/hg38.fa
	$< index $(word 2,$^)

test/hg38.bwa.read1.fastq.gz: bin/sim.sh bin/dwgsim data/hg38.fa
	$<

test/hg38.mm.bam: bin/map_mm.sh test/hg38.bwa.read1.fastq.gz bin/minimap2 bin/samtools
	$<

test/hg38.bwa.bam: bin/map_bwa.sh test/hg38.bwa.read1.fastq.gz bin/bwa bin/samtools
	$<

test/hg38.mm.eval.txt.gz: bin/eval.sh test/hg38.mm.bam
	$^

test/hg38.bwa.eval.txt.gz: bin/eval.sh test/hg38.bwa.bam
	$^

clean:
	rm -rf test/hg38.* bin/minimap2* bin/samtools* bin/dwgsim* bin/DWGSIM bin/bwa*
