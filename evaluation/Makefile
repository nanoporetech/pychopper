READS="data/SIRV_E0_pcs109_25k.fq"
REF="data/sirv_transcriptome.fas"
CORES=8
TUN_SAMPLE=10000
IN=$(READS)
IN_EDLIB=edlib_classified_reads.fq
RES_EDLIB=edlib_rescued_reads.fq
IN_PHMM=phmm_classified_reads.fq
RES_PHMM=phmm_rescued_reads.fq

all:
	@echo Reference is at: $(READS)
	@echo Reads are at: $(REF)
	# seqkit sample -n $(NR_READS) $(READS) > $(IN)
	./scripts/process_reads.sh $(REF) $(IN) $(CORES) RAW
	pychopper -S phmm_report.tsv -A phmm_hits.bed -w $(RES_PHMM) -t $(CORES) -m phmm -Y $(TUN_SAMPLE) -r phmm_pychopper_report.pdf $(IN) $(IN_PHMM)
	./scripts/process_reads.sh $(REF) $(IN_PHMM) $(CORES) TRIM_PHMM
	./scripts/process_reads.sh $(REF) $(RES_PHMM) $(CORES) RES_PHMM
	pychopper -S edlib_report.tsv -A edlib_hits.bed -w $(RES_EDLIB) -t $(CORES) -m edlib -Y $(TUN_SAMPLE) -r edlib_pychopper_report.pdf $(IN) $(IN_EDLIB)
	./scripts/process_reads.sh $(REF) $(IN_EDLIB) $(CORES) TRIM_EDLIB
	./scripts/process_reads.sh $(REF) $(RES_EDLIB) $(CORES) RES_EDLIB

	./scripts/compare.sh wspace_RAW RAW wspace_TRIM_PHMM PHMM
	./scripts/compare.sh wspace_RAW RAW wspace_TRIM_EDLIB EDLIB

clean:
	rm -fr wspace_* *.pdf *.fq cmp_*
