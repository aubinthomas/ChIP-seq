/*
 * Alignement on reference genome with Bwa-mem
 * Keeping unaligned reads in output FASTQ
 */

process bwaMemUnaligned {
  tag "${meta.id}"
  label 'bwa'
  label 'highCpu'
  label 'highMem'
  conda "bioconda::bwa=0.7.19 bioconda::samtools=1.21"

  input:
  tuple val(meta), path(reads)
  path(index)

  output:
  tuple val(meta), path("*.bam")                  , emit: bam
  tuple val(meta), path("*_unaligned*.fastq.gz")  , emit: unaligned_fastq
  path("*.log")                                   , emit: logs
  path("versions.txt")                            , emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def args = task.ext.args ?: ''
  def prefix = task.ext.prefix ?: "${meta.id}"
  def readList = reads instanceof List ? reads : [reads]
  """
  localIndex=`find -L ./ -name "*.amb" | sed 's/.amb//'`

  bwa \\
    mem \\
    $args \\
    -t $task.cpus \\
    \${localIndex} \\
    $reads | samtools view -bS -@ $task.cpus -o ${prefix}.bam -

  if [ ${readList.size()} -eq 2 ]; then
      samtools view -f 4 -@ $task.cpus ${prefix}.bam | samtools fastq -1 ${prefix}_unaligned_1.fastq -2 ${prefix}_unaligned_2.fastq -0 /dev/null -s /dev/null -n
      gzip ${prefix}_unaligned_1.fastq ${prefix}_unaligned_2.fastq
  else
      samtools view -f 4 -@ $task.cpus ${prefix}.bam | samtools fastq - > ${prefix}_unaligned.fastq
      gzip ${prefix}_unaligned.fastq
  fi

  getBWAstats.sh -i ${prefix}.bam -p ${task.cpus} > ${prefix}_bwa.log
  echo "Bwa-mem "\$(bwa 2>&1 | grep Version | cut -d" " -f2) &> versions.txt
  """
}