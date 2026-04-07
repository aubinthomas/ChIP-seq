/* 
 * Mapping Worflow 2
 * Strategy: Spike-first subtraction.
 */

if (params.aligner == "bwa-mem"){
  include { bwaMemUnaligned as mappingSpike } from '../process/bwaMemUnaligned'
  include { bwaMem as mapping } from '../../common/process/bwa/bwaMem'
}else if (params.aligner == "bowtie2"){
  include { bowtie2Unaligned as mappingSpike } from '../process/bowtie2Unaligned'
  include { bowtie2 as mapping } from '../../common/process/bowtie2/bowtie2'
}else if (params.aligner == "star"){
  include { starAlignUnaligned as mappingSpike } from '../process/starAlignUnaligned'
  include { starAlign as mapping } from '../../common/process/star/starAlign'
}

include { compareBams } from '../../local/process/compareBams'
include { samtoolsSort as samtoolsSort } from '../../common/process/samtools/samtoolsSort'
include { samtoolsSort as samtoolsSortSpike } from '../../common/process/samtools/samtoolsSort'
include { samtoolsIndex as samtoolsIndex } from '../../common/process/samtools/samtoolsIndex'
include { samtoolsIndex as samtoolsIndexSpike } from '../../common/process/samtools/samtoolsIndex'
include { samtoolsFlagstat } from '../../common/process/samtools/samtoolsFlagstat'


process FIX_BAM_PREFIX {
    tag "$meta.id"
    executor 'local'

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${meta.id}_${params.genome}.bam"), emit: bam

    script:
    """
    mv $bam ${meta.id}_${params.genome}.bam
    """
}

workflow mappingFlow2 {

  take:
  reads
  indexRef
  indexSpike

  main:
  
  chVersions = Channel.empty()
  
  // Align on spike genome
  if (params.aligner == "star"){
    mappingSpike(
      reads,
      indexSpike.collect(),
      Channel.empty().collect().ifEmpty([])
    )
  }else{
    mappingSpike(
      reads,
      indexSpike.collect()
    )
  }
  chVersions = chVersions.mix(mappingSpike.out.versions)

  // Now align on reference genome the unaligned reads from the spike alignment
  if (params.aligner == "star"){
    mapping(
      mappingSpike.out.unaligned_fastq,
      indexRef.collect(),
      Channel.empty().collect().ifEmpty([])
    )
  }else{
    mapping(
      mappingSpike.out.unaligned_fastq,
      indexRef.collect()
    )
  }
  chVersions = chVersions.mix(mapping.out.versions)

  // Dans la stratégie de soustraction, la comparaison est implicite (filtre au niveau FASTQ).
  // compareBams n'est pas compatible ici car les sets de lectures sont asymétriques.
  // chBam = mapping.out.bam
  FIX_BAM_PREFIX(mapping.out.bam)
  chBam = FIX_BAM_PREFIX.out.bam
  chSpikeBam = mappingSpike.out.bam

  // Post-processing Reference BAMs
  samtoolsSort(chBam)
  samtoolsIndex(samtoolsSort.out.bam)
  samtoolsFlagstat(samtoolsSort.out.bam)
  chVersions = chVersions.mix(samtoolsSort.out.versions, samtoolsIndex.out.versions, samtoolsFlagstat.out.versions)

  // Post-processing Spike BAMs
  samtoolsSortSpike(chSpikeBam)
  samtoolsIndexSpike(samtoolsSortSpike.out.bam)

  emit:
  bam      = samtoolsSort.out.bam.join(samtoolsIndex.out.bai)
  logs     = mapping.out.logs
  flagstat = samtoolsFlagstat.out.stats
  spikeBam = samtoolsSortSpike.out.bam.join(samtoolsIndexSpike.out.bai)
  spikeLogs = mappingSpike.out.logs
  compareBamsMqc = Channel.empty()
  versions = chVersions
}