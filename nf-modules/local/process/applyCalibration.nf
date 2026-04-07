process APPLY_CALIBRATION {
    tag "$meta.id"
    label 'process_high'

    input:
    tuple val(meta), path(bam), path(bai), val(scaleFactor)
    path blacklist
    val effGenomeSize

    output:
    tuple val(meta), path("*.bigwig"), emit: bigwig
    path "versions.txt"              , emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def blacklist_option = blacklist ? "--blackListFileName $blacklist" : ""
    """
    bamCoverage \\
        --bam $bam \\
        --outFileName ${prefix}_calibrated.bigwig \\
        --outFileFormat bigwig \\
        --binSize 10 \\
        --numberOfProcessors $task.cpus \\
        --scaleFactor $scaleFactor \\
        $blacklist_option \\
        --effectiveGenomeSize $effGenomeSize

    echo \$(bamCoverage --version | awk '{print \$2}') > versions.txt
    """
}
