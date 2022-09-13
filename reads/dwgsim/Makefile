all: tool data bam eval
eval: test/hg38.mm.eval.txt.gz test/hg38.bwa.eval.txt.gz test/hg38.dr.eval.txt.gz
bam: test/hg38.mm.bam test/hg38.bwa.bam test/hg38.dr.bam
tool: bin/samtools bin/minimap2 bin/dwgsim bin/bwa bin/dragen-os
data: data/hg38.fa data/hg38_dragmap data/hg38.fa.sa

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

bin/dragen-os: bin/install_dragmap.sh
	$<

data/hg38_dragmap: bin/dragen-os data/hg38.fa
	mkdir -p $@
	$< --num-threads 8 --build-hash-table true --ht-reference $(word 2,$^) --output-directory $@

data/hg38.fa.sa: bin/bwa data/hg38.fa
	$< index $(word 2,$^)

test/hg38.bwa.read1.fastq.gz: bin/sim.sh bin/dwgsim data/hg38.fa
	$<

test/hg38.mm.bam: bin/map_mm.sh test/hg38.bwa.read1.fastq.gz bin/minimap2 bin/samtools
	$<

test/hg38.bwa.bam: bin/map_bwa.sh test/hg38.bwa.read1.fastq.gz bin/bwa bin/samtools
	$<

test/hg38.dr.bam: bin/map_dragen.sh test/hg38.bwa.read1.fastq.gz bin/dragen-os bin/samtools
	$<

test/hg38.mm.eval.txt.gz: bin/eval.sh test/hg38.mm.bam
	$^

test/hg38.bwa.eval.txt.gz: bin/eval.sh test/hg38.bwa.bam
	$^

test/hg38.dr.eval.txt.gz: bin/eval.sh test/hg38.dr.bam
	$^

clean:
	rm -rf test/hg38.* bin/minimap2* bin/samtools* bin/dwgsim* bin/DWGSIM bin/bwa*
