THREADS := 8
PREFIX := hg38

# settings for simulating reads
# -N INT number of read pairs (-1 to disable) [-1]
N := 10000000

# -1 INT length of the first read [70]
# -2 INT length of the second read [70]
READ_LEN := 100

# -e FLOAT per base/color/flow error rate of the first read [from 0.020 to 0.020 by 0.000]
# -E FLOAT per base/color/flow error rate of the second read [from 0.020 to 0.020 by 0.000]
ERROR_RATE := 0

# note that the mutation rate is with respect to the reference FASTA. The
# # default rate of 0.0010 will generate around 3 million mutations in a human
# # reference genome
# # -r FLOAT rate of mutations [0.0010]
MUT_RATE := 0.002

# -R FLOAT fraction of mutations that are indels [0.10]
INDEL_RATE := 1

# -X FLOAT probability an indel is extended [0.30]
INDEL_EXT := 0.75

# -I INT the minimum length indel [1]
# -y FLOAT probability of a random DNA read [0.05]
RAND_PROB = 0

# -z INT random seed (-1 uses the current time) [-1]
SEED := 1984

OUTDIR := test/dwgsim_$(PREFIX)_$(N)_$(READ_LEN)_$(ERROR_RATE)_$(MUT_RATE)_$(INDEL_RATE)_$(INDEL_EXT)_$(RAND_PROB)_$(SEED)

all: $(OUTDIR) tool data bam eval
eval: $(OUTDIR)/$(PREFIX).mm.eval.txt.gz $(OUTDIR)/$(PREFIX).bwa.eval.txt.gz $(OUTDIR)/$(PREFIX).dr.eval.txt.gz
bam: $(OUTDIR)/$(PREFIX).mm.bam $(OUTDIR)/$(PREFIX).bwa.bam $(OUTDIR)/$(PREFIX).dr.bam
tool: bin/samtools bin/minimap2 bin/dwgsim bin/bwa bin/dragen-os
data: data/$(PREFIX).fa data/$(PREFIX)_dragmap data/$(PREFIX).fa.sa

READ1 := $(OUTDIR)/$(PREFIX).bwa.read1.fastq.gz
READ2 := $(OUTDIR)/$(PREFIX).bwa.read2.fastq.gz

.PHONY: all clean

$(OUTDIR):
	mkdir $@

data/$(PREFIX).fa: data/download.sh
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

data/$(PREFIX)_dragmap: bin/dragen-os data/$(PREFIX).fa
	mkdir -p $@
	$< --num-threads 8 --build-hash-table true --ht-reference $(word 2,$^) --output-directory $@

data/$(PREFIX).fa.sa: bin/bwa data/$(PREFIX).fa
	$< index $(word 2,$^)

$(READ1): bin/dwgsim data/$(PREFIX).fa
	$< -N $(N) -1 $(READ_LEN) -2 $(READ_LEN) -e $(ERROR_RATE) -E $(ERROR_RATE) -r $(MUT_RATE) -R $(INDEL_RATE) -X $(INDEL_EXT) -y $(RAND_PROB) -z $(SEED) $(word 2,$^) $(OUTDIR)/$(PREFIX)

$(OUTDIR)/$(PREFIX).mm.bam: bin/minimap2 $(READ1) bin/samtools
	$< -t $(THREADS) -ax sr data/$(PREFIX).fa $(READ1) $(READ2) | bin/samtools sort -@ $(THREADS) -O BAM | tee $(OUTDIR)/$(PREFIX).mm.bam | bin/samtools index - $(OUTDIR)/$(PREFIX).mm.bam.bai

$(OUTDIR)/$(PREFIX).bwa.bam: bin/bwa $(READ1) bin/samtools
	$< mem -t $(THREADS) data/$(PREFIX).fa $(READ1) $(READ2) | bin/samtools sort -@ $(THREADS) -O BAM | tee $(OUTDIR)/$(PREFIX).bwa.bam | bin/samtools index - $(OUTDIR)/$(PREFIX).bwa.bam.bai

$(OUTDIR)/$(PREFIX).dr.bam: bin/dragen-os $(READ1) bin/samtools
	$< --num-threads $(THREADS) -r data/$(PREFIX)_dragmap -1 $(READ1) -2 $(READ2) | bin/samtools sort -@ $(THREADS) -O BAM | tee $(OUTDIR)/$(PREFIX).dr.bam | bin/samtools index - $(OUTDIR)/$(PREFIX).dr.bam.bai

$(OUTDIR)/$(PREFIX).mm.eval.txt.gz: bin/eval.sh $(OUTDIR)/$(PREFIX).mm.bam
	$^

$(OUTDIR)/$(PREFIX).bwa.eval.txt.gz: bin/eval.sh $(OUTDIR)/$(PREFIX).bwa.bam
	$^

$(OUTDIR)/$(PREFIX).dr.eval.txt.gz: bin/eval.sh $(OUTDIR)/$(PREFIX).dr.bam
	$^

clean:
	rm -rf $(OUTDIR)/$(PREFIX).* bin/minimap2* bin/samtools* bin/dwgsim* bin/DWGSIM bin/bwa*
