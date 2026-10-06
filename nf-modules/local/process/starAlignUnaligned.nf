/*
 * STAR reads alignment
 * Keeping unaligned reads in output FASTQ
 */

process starAlignUnaligned {
  tag "$meta.id"
  label 'star'
  label 'highCpu'
  label 'extraMem'
  conda "bioconda::star=2.7.11b"

  input:
  tuple val(meta), path(reads)
  path index
  path gtf

  output:
  tuple val(meta), path('*Aligned.out.bam'), emit: bam
  tuple val(meta), path("*_unaligned*.fastq.gz"), emit: unaligned_fastq
  path ("*out"), emit: logs
  path ("versions.txt"), emit: versions
  tuple val(meta), path("*ReadsPerGene.out.tab"), optional: true, emit: counts
  path("*out.tab"), optional: true, emit: countsLogs
  tuple val(meta), path("*Aligned.toTranscriptome.out.bam"), optional: true, emit: transcriptsBam

  when:
  task.ext.when == null || task.ext.when

  script:
  def args = task.ext.args ?: ''
  def prefix = task.ext.prefix ?: "${meta.id}"
  def gtfOpts = gtf.size() > 0 ? "--sjdbGTFfile ${gtf}" : ""
  """
  echo "STAR "\$(STAR --version 2>&1) > versions.txt
  STAR --genomeDir $index \\
       --readFilesIn $reads  \\
       --runThreadN ${task.cpus} \\
       --runMode alignReads \\
       --outSAMtype BAM Unsorted  \\
       --readFilesCommand zcat \\
       --runDirPerm All_RWX \\
       --outTmpDir "star_\$(date +%d%s%S%N)"\\
       --outFileNamePrefix ${prefix}  \\
       --outSAMattrRGline ID:$meta.id SM:$meta.id LB:Illumina PL:Illumina  \\
       --outSAMunmapped Within \\
       --outReadsUnmapped Fastx \\
       ${gtfOpts} \\
       ${args}

  # Post-processing pour harmoniser les noms de fichiers unmapped avec les autres aligneurs
  if [ -f ${prefix}Unmapped.out.mate2 ]; then
      mv ${prefix}Unmapped.out.mate1 ${prefix}_unaligned_1.fastq
      mv ${prefix}Unmapped.out.mate2 ${prefix}_unaligned_2.fastq
      gzip ${prefix}_unaligned_1.fastq ${prefix}_unaligned_2.fastq
  elif [ -f ${prefix}Unmapped.out.mate1 ]; then
      mv ${prefix}Unmapped.out.mate1 ${prefix}_unaligned.fastq
      gzip ${prefix}_unaligned.fastq
  fi
  """
}