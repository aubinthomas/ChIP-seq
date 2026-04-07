/*
 * Alignment on reference genome with Bowtie2
 * Keeping unaligned reads in output FASTQ
 */

process bowtie2Unaligned {
  tag "${meta.id}"
  label 'bowtie2'
  label 'highCpu'
  label 'highMem'
  conda "bioconda::bowtie2=2.5.4 bioconda::samtools=1.21"

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
  script:
  def args = task.ext.args ?: ''
  def prefix = task.ext.prefix ?: "${meta.id}"
  def inputOpts = meta.singleEnd ? "-U ${reads[0]}" : "-1 ${reads[0]} -2 ${reads[1]}"
  def readList = reads instanceof List ? reads : [reads]
  """
  localIndex=`find -L ./ -name "*.rev.1.bt2" | sed 's/.rev.1.bt2//'`
  echo \$(bowtie2 --version | awk 'NR==1{print "bowtie2 "\$3}') > versions.txt

  # Note : L'option -S de samtools view est obsolète, -b suffit
  bowtie2 -p ${task.cpus} \\
          ${args} \\
          -x \${localIndex} \\
          $inputOpts 2> ${prefix}_bowtie2.log | samtools view -b -@ $task.cpus -o ${prefix}.bam -

  if [ ${readList.size()} -eq 2 ]; then
      # Utilisation de samtools fastq directement sur le BAM avec le filtre -f 4
      samtools fastq -f 4 -@ $task.cpus -1 ${prefix}_unaligned_1.fastq -2 ${prefix}_unaligned_2.fastq -0 /dev/null -s /dev/null -n ${prefix}.bam
      gzip ${prefix}_unaligned_1.fastq ${prefix}_unaligned_2.fastq
  else
      # Pareil pour le single-end
      samtools fastq -f 4 -@ $task.cpus ${prefix}.bam > ${prefix}_unaligned.fastq
      gzip ${prefix}_unaligned.fastq
  fi
  """
}