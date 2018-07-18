#!/usr/bin/env nextflow

//+++++++++++++++++ SETUP++++++++++++++++++++++++
params.workdir = "/home/johannes/rdrive/PPG_SEQ_DATA-LICHTJ-SE00182/johannes/notebook/bcinerea_test"
params.refs = "${params.workdir}/reference/*.fasta"
params.input = "${params.workdir}/input/*.fasta"
params.outdir = "${params.workdir}/output"
//+++++++++++++++++++++++++++++++++++++++++++++++


log.info "====================================================================="
log.info "genome coverage compared to a reference"
log.info "Output    : ${params.outdir}"
log.info "====================================================================="

refs = Channel
.fromPath(params.refs)
.map{file -> tuple(file.simpleName, file)}

qrys = Channel
.fromPath(params.input)
.map{file -> tuple(file.simpleName, file)}


refs
.tap { genomes1 }
.combine(qrys)
.filter{ it[0] != it[2] }
.tap { genomePairs }

process makeGenomeFiles {
  tag { id }

  input:
  set id, "genome.fasta" from genomes1

  output:
  set id, "genome.fasta.fai" into genomeFiles

  """
samtools faidx genome.fasta
  """
}

process nucmer {
  tag { "${idR} vs ${idQ}" }

  input:
  set idR, seqR, idQ, seqQ from genomePairs

  output:
  set idR, idQ, val("raw"), "out.delta" into rawDeltas
  set idR, idQ, val("filtered"), "filtered.delta" into fltDeltas

  """
cat $seqR > ref.fasta
cat $seqQ > qry.fasta
nucmer --maxmatch ref.fasta qry.fasta
delta-filter -1 out.delta > filtered.delta
  """
}

process showcoords {
  tag { "${idR} vs ${idQ}" }

  input:
  set idR, idQ, type, "in.delta" from rawDeltas.mix(fltDeltas)

  output:
  set idR, idQ, "out.bed" into rawBed

"""
show-coords -T in.delta \
 | tail -n +5 \
 | awk 'BEGIN{OFS="\\t"} {print(\$8, \$1-1, \$2, ".", "+", \$7)}' \
 | sort -k1,1 -k2,2n \
 | bedtools merge -i - \
 | sed 's/\$/\\t$type/g' \
 > out.bed
"""
}

process genomeCoverage {
  tag { "${idR} vs ${idQ}" }

  input:
  set idR, idQ, "in.bed", "in.genome" from rawBed.combine(genomeFiles, by: 0)

  output:
  set idR, idQ, "out.txt" into genomeCoverages, genomeCoverages2

  """
sed 's/\$/\\t$idQ/g' in.bed > out.txt
awk 'BEGIN{OFS="\\t"} {print(\$1, 0, \$2, "base", "$idQ")}' in.genome >> out.txt
  """
}

genomeCoverages
.groupTuple()
.set { plottingInputs }

genomeCoverages2
.groupTuple()
.set { plottingInputs2 }

process plotCoverages1 {
  publishDir "${params.outdir}/nucmer", mode: 'copy'
  tag { "${idR} vs ${idQ}" }

  input:
  set idR, queryIDs, "matches.*.txt" from plottingInputs

  output:
  file("genomeCoverage.svg") into debug

  """
#!/usr/bin/env Rscript

library(readr)
library(dplyr)
library(ggplot2)
library(magrittr)

list.files(".", "matches.*.txt") %>%
  lapply(read_tsv, col_names = c("seqid", "start", "stop", "type", "query")) %>%
  do.call(what=rbind) -> data

  data %>%
    filter(type == "base") %>%
    ggplot(aes(x = query, y=start + (stop-start)/2, height=(stop-start))) +
    geom_tile(fill="white") +
    geom_tile(data=filter(data, type == "raw"), fill="#E6A0C4") +
    geom_tile(data=filter(data, type == "filtered"), fill="#7294D4") +
    theme_minimal() +
    facet_grid(. ~ seqid) +
    theme(axis.text.x = element_text(angle = 90)) -> plot

    ggsave(filename="genomeCoverage.svg", plot=plot, width=40)
    """
}
process plotCoverages2 {
  publishDir "${params.outdir}/nucmer", mode: 'copy'
  tag { "${idR} vs ${idQ}" }

  input:
  set idR, queryIDs, "matches.*.txt" from plottingInputs2

  output:
  file("genomeCoverage2.svg") into debug2

  """
#!/usr/bin/env Rscript

library(readr)
library(dplyr)
library(ggplot2)
library(magrittr)

list.files(".", "matches.*.txt") %>%
  lapply(read_tsv, col_names = c("seqid", "start", "stop", "type", "query")) %>%
  do.call(what=rbind) -> data

data %>%
  filter(type == "base") %>%
  ggplot(aes(xmin=start, xmax=stop, ymin=0, ymax=1)) +
  geom_rect(fill="white") +
  geom_rect(data=filter(data, type == "raw"), fill="#E6A0C4") +
  geom_rect(data=filter(data, type == "filtered"), fill="#7294D4") +
  facet_grid(seqid ~ query) +
  theme(strip.text.y = element_text(angle = 0),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_blank()) -> plot

  ggsave(filename="genomeCoverage2.svg", plot=plot, width=20)
  """
}

debug.println()
debug2.println()
