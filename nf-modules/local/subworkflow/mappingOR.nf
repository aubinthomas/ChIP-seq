/* 
 * Mapping Worflow 2
 * Strategy: Strict Double Subtraction (Cross-mapping prevention)
 */

if (params.aligner == "bwa-mem"){
  // Voie 1 : Isoler la Référence
  include { bwaMemUnaligned as mappingSpikeUnaligned } from '../process/bwaMemUnaligned'
  include { bwaMem as mappingRefFinal } from '../../common/process/bwa/bwaMem'
  
  // Voie 2 : Isoler le Spike
  include { bwaMemUnaligned as mappingRefUnaligned } from '../process/bwaMemUnaligned'
  include { bwaMem as mappingSpikeFinal } from '../../common/process/bwa/bwaMem'

} else if (params.aligner == "bowtie2"){
  // Voie 1 : Isoler la Référence
  include { bowtie2Unaligned as mappingSpikeUnaligned } from '../process/bowtie2Unaligned'
  include { bowtie2 as mappingRefFinal } from '../../common/process/bowtie2/bowtie2'
  
  // Voie 2 : Isoler le Spike
  include { bowtie2Unaligned as mappingRefUnaligned } from '../process/bowtie2Unaligned'
  include { bowtie2 as mappingSpikeFinal } from '../../common/process/bowtie2/bowtie2'

} else if (params.aligner == "star"){
  // Voie 1 : Isoler la Référence
  include { starAlignUnaligned as mappingSpikeUnaligned } from '../process/starAlignUnaligned'
  include { starAlign as mappingRefFinal } from '../../common/process/star/starAlign'
  
  // Voie 2 : Isoler le Spike
  include { starAlignUnaligned as mappingRefUnaligned } from '../process/starAlignUnaligned'
  include { starAlign as mappingSpikeFinal } from '../../common/process/star/starAlign'
}

include { compareBams } from '../../local/process/compareBams'
include { samtoolsSort as samtoolsSort } from '../../common/process/samtools/samtoolsSort'
include { samtoolsSort as samtoolsSortSpike } from '../../common/process/samtools/samtoolsSort'
include { samtoolsIndex as samtoolsIndex } from '../../common/process/samtools/samtoolsIndex'
include { samtoolsIndex as samtoolsIndexSpike } from '../../common/process/samtools/samtoolsIndex'
include { samtoolsFlagstat } from '../../common/process/samtools/samtoolsFlagstat'

// NETTOIE TOUS LES NOMS AVANT L'EXPORT
process CLEAN_REF_OUTPUTS {
    tag "$meta.id"
    label 'unix'
    label 'minCpu'
    label 'minMem'

    input:
    tuple val(meta), path(bam), path(bai), path(flagstat)

    output:
    tuple val(meta), path("${meta.id}_${params.genome}.bam"), path("*.bai"), emit: bam
    tuple val(meta), path("${meta.id}_${params.genome}.flagstats"), emit: flagstat

    script:
    """
    mv $bam ${meta.id}_${params.genome}.bam
    mv $bai ${meta.id}_${params.genome}.bam.bai
    mv $flagstat ${meta.id}_${params.genome}.flagstats
    """
}

workflow mappingFlow2 {

  take:
  reads
  indexRef
  indexSpike

  main:
  chVersions = Channel.empty()
  
  // =========================================================================
  // VOIE 1 : OBTENIR LE BAM RÉFÉRENCE PUR
  // =========================================================================
  
  // 1a. Aligner sur le Spike et récupérer les non-alignés (unmapped)
  if (params.aligner == "star"){
    mappingSpikeUnaligned( reads, indexSpike.collect(), Channel.empty().collect().ifEmpty([]) )
  }else{
    mappingSpikeUnaligned( reads, indexSpike.collect() )
  }
  chVersions = chVersions.mix(mappingSpikeUnaligned.out.versions)

  // 1b. Aligner ces restes sur la Référence
  if (params.aligner == "star"){
    mappingRefFinal( mappingSpikeUnaligned.out.unaligned_fastq, indexRef.collect(), Channel.empty().collect().ifEmpty([]) )
  }else{
    mappingRefFinal( mappingSpikeUnaligned.out.unaligned_fastq, indexRef.collect() )
  }
  chVersions = chVersions.mix(mappingRefFinal.out.versions)

  // Add genome information (needed for consistent file naming downstream)
  chBam = mappingRefFinal.out.bam.map{ meta, bam ->
    def newMeta = [ id: meta.id, name: meta.name, singleEnd: meta.singleEnd, genome: params.genome ]
    [newMeta, bam]
  }

  // =========================================================================
  // VOIE 2 : OBTENIR LE BAM SPIKE PUR
  // =========================================================================

  // 2a. Aligner sur la Référence et récupérer les non-alignés (unmapped)
  if (params.aligner == "star"){
    mappingRefUnaligned( reads, indexRef.collect(), Channel.empty().collect().ifEmpty([]) )
  }else{
    mappingRefUnaligned( reads, indexRef.collect() )
  }
  chVersions = chVersions.mix(mappingRefUnaligned.out.versions)

  // 2b. Aligner ces restes sur le Spike
  if (params.aligner == "star"){
    mappingSpikeFinal( mappingRefUnaligned.out.unaligned_fastq, indexSpike.collect(), Channel.empty().collect().ifEmpty([]) )
  }else{
    mappingSpikeFinal( mappingRefUnaligned.out.unaligned_fastq, indexSpike.collect() )
  }
  chVersions = chVersions.mix(mappingSpikeFinal.out.versions)

  // Add genome information (needed for consistent file naming downstream)
  chSpikeBam = mappingSpikeFinal.out.bam.map{ meta, bam ->
    def newMeta = [ id: meta.id, name: meta.name, singleEnd: meta.singleEnd, genome: params.spike ]
    [newMeta, bam]
  }

  // =========================================================================
  // POST-PROCESSING
  // =========================================================================

  // Post-processing Reference BAMs
  samtoolsSort(chBam)
  samtoolsIndex(samtoolsSort.out.bam)
  samtoolsFlagstat(samtoolsSort.out.bam)

  ch_to_clean = samtoolsSort.out.bam
    .join(samtoolsIndex.out.bai)
    .join(samtoolsFlagstat.out.stats)

  CLEAN_REF_OUTPUTS(ch_to_clean)

  chVersions = chVersions.mix(samtoolsSort.out.versions, samtoolsIndex.out.versions, samtoolsFlagstat.out.versions)

  // Post-processing Spike BAMs
  samtoolsSortSpike(chSpikeBam)
  samtoolsIndexSpike(samtoolsSortSpike.out.bam)

  emit:
  bam      = CLEAN_REF_OUTPUTS.out.bam
  logs     = mappingRefFinal.out.logs
  flagstat = CLEAN_REF_OUTPUTS.out.flagstat
  
  spikeBam = samtoolsSortSpike.out.bam.join(samtoolsIndexSpike.out.bai)
  spikeLogs = mappingSpikeFinal.out.logs
  compareBamsMqc = Channel.empty()
  versions = chVersions
}

