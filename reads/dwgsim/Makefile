include tool_versions
.PHONY: all clean data tool index vcf bam eval

THREADS := 8
PREFIX := hg38

# settings for simulating reads
# -N INT number of read pairs (-1 to disable) [-1]
N := 200000000

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
MUT_RATE := 0.005

# -R FLOAT fraction of mutations that are indels [0.10]
INDEL_RATE := 0.50

# -X FLOAT probability an indel is extended [0.30]
INDEL_EXT := 0.75

# -I INT the minimum length indel [1]
# -y FLOAT probability of a random DNA read [0.05]
RAND_PROB = 0

# -z INT random seed (-1 uses the current time) [-1]
SEED := 1984

OUTDIR := test/dwgsim_$(PREFIX)_$(N)_$(READ_LEN)_$(ERROR_RATE)_$(MUT_RATE)_$(INDEL_RATE)_$(INDEL_EXT)_$(RAND_PROB)_$(SEED)

all: $(OUTDIR) data tool index vcf bam eval
eval: $(OUTDIR)/$(PREFIX).mm.eval.txt.gz $(OUTDIR)/$(PREFIX).bwa.eval.txt.gz $(OUTDIR)/$(PREFIX).dr.eval.txt.gz $(OUTDIR)/$(PREFIX).bwa.bamreadcount.tsv
bam: $(OUTDIR)/$(PREFIX).mm.bam $(OUTDIR)/$(PREFIX).bwa.bam $(OUTDIR)/$(PREFIX).dr.bam
vcf: $(OUTDIR)/$(PREFIX).mutations.vcf.gz.csi
index: data/$(PREFIX)_dragmap/hash_table.cfg data/$(PREFIX).fa.sa
tool: bin/samtools-$(samtools_version)/samtools bin/minimap2-$(minimap2_version)/minimap2 bin/dwgsim-$(dwgsim_version)/dwgsim bin/bwa-$(bwa_version)/bwa bin/DRAGMAP-$(dragmap_version)/build/release/dragen-os bin/bcftools-$(bcftools_version)/bcftools bin/bam-readcount-$(bamreadcount_version)/build/bin/bam-readcount
data: data/$(PREFIX).fa

READ1 := $(OUTDIR)/$(PREFIX).bwa.read1.fastq.gz
READ2 := $(OUTDIR)/$(PREFIX).bwa.read2.fastq.gz

$(OUTDIR):
	mkdir $@

data/$(PREFIX).fa: data/download.sh
	$<

bin/samtools-$(samtools_version)/samtools: bin/install_samtools.sh
	$< $(samtools_version)

bin/bcftools-$(bcftools_version)/bcftools: bin/install_bcftools.sh
	$< $(bcftools_version)

bin/minimap2-$(minimap2_version)/minimap2: bin/install_minimap2.sh
	$< $(minimap2_version)

bin/dwgsim-$(dwgsim_version)/dwgsim: bin/install_dwgsim.sh
	$< $(dwgsim_version)

bin/bwa-$(bwa_version)/bwa: bin/install_bwa.sh
	$< $(bwa_version)

bin/DRAGMAP-$(dragmap_version)/build/release/dragen-os: bin/install_dragmap.sh
	$< $(dragmap_version)

bin/bam-readcount-$(bamreadcount_version)/build/bin/bam-readcount: bin/install_bamreadcount.sh
	$< $(bamreadcount_version)

data/$(PREFIX)_dragmap/hash_table.cfg: bin/DRAGMAP-$(dragmap_version)/build/release/dragen-os data/$(PREFIX).fa
	mkdir -p $(dir $@)
	$< --num-threads $(THREADS) --build-hash-table true --ht-reference $(word 2,$^) --output-directory $(dir $@)

data/$(PREFIX).fa.sa: bin/bwa-$(bwa_version)/bwa data/$(PREFIX).fa
	$< index $(word 2,$^)

$(READ1): bin/dwgsim-$(dwgsim_version)/dwgsim data/$(PREFIX).fa
	$< -N $(N) -1 $(READ_LEN) -2 $(READ_LEN) -e $(ERROR_RATE) -E $(ERROR_RATE) -r $(MUT_RATE) -R $(INDEL_RATE) -X $(INDEL_EXT) -y $(RAND_PROB) -z $(SEED) $(word 2,$^) $(OUTDIR)/$(PREFIX)

$(OUTDIR)/$(PREFIX).mm.bam: bin/minimap2-$(minimap2_version)/minimap2 $(READ1) bin/samtools-$(samtools_version)/samtools
	$< -t $(THREADS) -ax sr data/$(PREFIX).fa $(READ1) $(READ2) | bin/samtools sort -@ $(THREADS) -O BAM | tee $(OUTDIR)/$(PREFIX).mm.bam | bin/samtools index - $(OUTDIR)/$(PREFIX).mm.bam.bai

$(OUTDIR)/$(PREFIX).bwa.bam: bin/bwa-$(bwa_version)/bwa data/$(PREFIX).fa.sa $(READ1) bin/samtools-$(samtools_version)/samtools
	$< mem -t $(THREADS) data/$(PREFIX).fa $(READ1) $(READ2) | bin/samtools sort -@ $(THREADS) -O BAM | tee $(OUTDIR)/$(PREFIX).bwa.bam | bin/samtools index - $(OUTDIR)/$(PREFIX).bwa.bam.bai

$(OUTDIR)/$(PREFIX).dr.bam: bin/DRAGMAP-$(dragmap_version)/build/release/dragen-os data/$(PREFIX)_dragmap/hash_table.cfg $(READ1) bin/samtools-$(samtools_version)/samtools
	$< --num-threads $(THREADS) -r $(dir $(word 2,$^)) -1 $(READ1) -2 $(READ2) | bin/samtools sort -@ $(THREADS) -O BAM | tee $(OUTDIR)/$(PREFIX).dr.bam | bin/samtools index - $(OUTDIR)/$(PREFIX).dr.bam.bai

$(OUTDIR)/$(PREFIX).mm.eval.txt.gz: bin/eval.sh $(OUTDIR)/$(PREFIX).mm.bam
	$^

$(OUTDIR)/$(PREFIX).bwa.eval.txt.gz: bin/eval.sh $(OUTDIR)/$(PREFIX).bwa.bam
	$^

$(OUTDIR)/$(PREFIX).dr.eval.txt.gz: bin/eval.sh $(OUTDIR)/$(PREFIX).dr.bam
	$^

$(OUTDIR)/$(PREFIX).bwa.bamreadcount.tsv: bin/bam-readcount-$(bamreadcount_version)/build/bin/bam-readcount $(OUTDIR)/$(PREFIX).bwa.bam
	$< -w 0 -f data/$(PREFIX).fa -l data/site_list $(word 2,$^) > $@

$(OUTDIR)/$(PREFIX).mutations.vcf.gz: bin/bcftools-$(bcftools_version)/bcftools $(READ1)
	$< view -O z $(OUTDIR)/$(PREFIX).mutations.vcf > $(OUTDIR)/$(PREFIX).mutations.vcf.gz

$(OUTDIR)/$(PREFIX).mutations.vcf.gz.csi: bin/bcftools-$(bcftools_version)/bcftools $(OUTDIR)/$(PREFIX).mutations.vcf.gz
	$< index $(word 2,$^)

clean:
	rm -rf $(OUTDIR)/$(PREFIX).* bin/minimap2* bin/samtools* bin/bcftools* bin/dwgsim* bin/bwa* bin/dragen-os bin/DRAGMAP-* data/$(PREFIX)_dragmap bin/bam-readcount*
