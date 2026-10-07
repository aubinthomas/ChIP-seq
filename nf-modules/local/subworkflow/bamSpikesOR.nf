/* * Occupancy Ratio (OR) Analysis 
 * Strategy: Cross-normalization between IP and Input using Spike-in counts.
 */

// Ce process génère les statistiques et exporte correctement COUNT
process SAMTOOLS_IDXSTATS_REF {
    tag "$meta.id"
    label 'samtools'
    label 'minCpu'
    label 'lowMem'

    input:
    tuple val(meta), path(bam), path(bai)
    
    output:
    tuple val(meta), env(COUNT), emit: count
    tuple val(meta), path("*.idxstats"), emit: stats
    path "versions.txt"                , emit: versions
    
    script:
    """
    samtools idxstats $bam > ${meta.id}.idxstats
    export COUNT=\$(awk '{s+=\$3} END {printf "%d", s+0}' ${meta.id}.idxstats)
    echo "samtools \$(samtools --version | head -n 1 | awk '{print \$2}')" > versions.txt
    """
}

// Process pour les BAMs Spike
process SAMTOOLS_IDXSTATS_SPIKE {
    tag "$meta.id"
    label 'samtools'
    label 'minCpu'
    label 'lowMem'

    input:
    tuple val(meta), path(bam), path(bai)
    
    output:
    tuple val(meta), env(COUNT), emit: count
    tuple val(meta), path("*.idxstats"), emit: stats
    path "versions.txt"                , emit: versions
    
    script:
    """
    samtools idxstats $bam > ${meta.id}_spike.idxstats
    export COUNT=\$(awk '{s+=\$3} END {printf "%d", s+0}' ${meta.id}_spike.idxstats)
    echo "samtools \$(samtools --version | head -n 1 | awk '{print \$2}')" > versions.txt
    """
}

// Ce process applique le facteur calculé
process APPLY_CALIBRATION {
    tag "$meta.id"
    label 'deeptools'
    label 'medCpu'
    label 'lowMem'

    input:
    tuple val(meta), path(bam), path(bai), val(scaleFactor)
    path blacklist
    val effGenomeSize

    output:
    tuple val(meta), path("*.bigwig"), emit: bigwig
    path "versions.txt"              , emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    // Sécurisation de l'option blacklist au cas où le fichier serait absent
    def blacklist_option = blacklist && blacklist.name != 'NO_FILE' ? "--blackListFileName $blacklist" : ""
    
    """
    bamCoverage \\
        --bam $bam \\
        --outFileName ${prefix}_calibrated.bigwig \\
        --outFileFormat bigwig \\
        --binSize 10 \\
        --numberOfProcessors $task.cpus \\
        --scaleFactor $scaleFactor \\
        $blacklist_option \\
        --effectiveGenomeSize ${effGenomeSize[0]}

    echo "bamCoverage \$(bamCoverage --version | awk '{print \$2}')" > versions.txt
    """
}

// Ce process calcule le facteur de normalisation Occupancy Ratio (OR).
// Il sauvegarde le facteur dans un fichier texte pour permettre la traçabilité des calculs.
process CALCULATE_OR {
    tag "${metaIP.id}"
    label 'unix'
    label 'minCpu'
    label 'minMem'
    

    input:
    // Ajout de refCtrl dans le tuple d'entrée
    tuple val(metaIP), val(spikeIP), val(spikeCtrl), val(refCtrl)

    output:
    tuple val(metaIP.id), env(OR_FACTOR), emit: factor
    path "*_or_factor.txt"              , emit: txt

    script:
    """
    # Calcul du facteur avec la formule de l'Occupancy Ratio
    export OR_FACTOR=\$(awk -v Wc="$spikeCtrl" -v IPc="$spikeIP" -v Wx="$refCtrl" 'BEGIN { 
        if (IPc == 0 || Wx == 0) { 
            printf "0" 
        } else { 
            # Facteur = SpikeCtrl / (SpikeIP * RefCtrl)
            printf "%.6e", Wc / (IPc * Wx) 
        } 
    }')
    
    # Écriture dans le fichier pour la traçabilité
    echo "\$OR_FACTOR" > ${metaIP.id}_or_factor.txt
    """
}



workflow bamSpikesORflow {
    take:
    bamsRef       // [meta, bam, bai] - Génome Référence
    bamsSpikes    // [meta, bam, bai] - Génome Spike
    chDesign      // [metaIP, metaControl] - Paires définies dans le design.csv
    blacklist
    effGenomeSize

    main:
    chVersions = Channel.empty()

    // 1. Extraction des statistiques de comptage
    ch_ref_out   = SAMTOOLS_IDXSTATS_REF(bamsRef)
    ch_spike_out = SAMTOOLS_IDXSTATS_SPIKE(bamsSpikes)

    chVersions = chVersions.mix(ch_ref_out.versions)

    // CORRECTION: Préparation des canaux avec l'ID en première position (index 0) 
    // pour garantir une jointure stricte et sans erreur sur les métadonnées.
    ch_ref_counts = ch_ref_out.count.map { meta, count -> [ meta.id, meta, count ] }
    ch_spike_counts = ch_spike_out.count.map { meta, count -> [ meta.id, meta, count ] }

    // [id, meta, ref_count, spike_count]
    ch_sample_counts = ch_ref_counts
        .join(ch_spike_counts)
        .map { id, meta, ref, meta_spike, spike -> [ id, meta, ref, spike ] }

    // 2 & 3. Recoupement IP / Control et calcul Groovy natif
    ch_inputs_for_or = chDesign
        // CORRECTION 1 : On lit la liste brute [ip, ctrl, rep, peak]
        .map { row -> 
            def ip_id = row[0]
            def ctrl_id = row[1]
            [ ip_id, ctrl_id ] 
        }
        // CORRECTION 2 : Jointure avec l'IP (la clé est l'id de l'IP)
        // ch_sample_counts contient : [id, meta, ref_count, spike_count]
        .join( ch_sample_counts ) 
        // Résultat du join : [ip_id, ctrl_id, metaIP, refIP, spikeIP]
        .map { ip_id, ctrl_id, metaIP, refIP, spikeIP -> 
            // On met l'id du Control en 1ère position pour la prochaine jointure
            [ ctrl_id, metaIP, spikeIP ] 
        }
        // CORRECTION 3 : Jointure avec le Control (la clé est l'id du Control)
        .join( ch_sample_counts ) 
        // Résultat du join : [ctrl_id, metaIP, spikeIP, metaCtrl, refCtrl, spikeCtrl]
        .map { ctrl_id, metaIP, spikeIP, metaCtrl, refCtrl, spikeCtrl -> 
            // On garde l'Input cible (refCtrl) pour la formule complète
            [ metaIP, spikeIP, spikeCtrl, refCtrl ] 
        }
    
    CALCULATE_OR(ch_inputs_for_or)
    ch_paired_counts = CALCULATE_OR.out.factor

    // 4. Application de la calibration sur le BAM de l'IP
    // On mappe bamsRef pour mettre l'ID en clé de jointure
    ch_bams_for_join = bamsRef.map { meta, bam, bai -> [ meta.id, meta, bam, bai ] }
    
    ch_to_calibrate = ch_bams_for_join
        .join(ch_paired_counts)
        // Reformatage propre pour le process APPLY_CALIBRATION : [meta, bam, bai, factor]
        .map { id, meta, bam, bai, factor -> [ meta, bam, bai, factor ] }

    APPLY_CALIBRATION(ch_to_calibrate, blacklist, effGenomeSize)
    chVersions = chVersions.mix(APPLY_CALIBRATION.out.versions)

    emit:
    bigwig = APPLY_CALIBRATION.out.bigwig
    versions = chVersions
}