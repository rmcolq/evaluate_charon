#!/usr/bin/env nextflow

/*
 * Limits how many samples are in flight at once, so that large inputs do not fill the disk. A slot is a directory
 * <slot_dir>/<unique_id>, removed by release_slot, or by the errorStrategy in conf/base.config when a sample fails.
 * Slots older than --slot_timeout minutes are reclaimed (e.g. after a killed run).
 */

process acquire_slot {
    label "process_single"
    tag "${unique_id}"
    maxForks 1
    cache false
    errorStrategy 'terminate'

    input:
    tuple val(unique_id), val(source)

    output:
    tuple val(unique_id), val(source)

    script:
    """
    slots="${params.slot_dir}"
    mkdir -p "\$slots"
    while true; do
      find "\$slots" -mindepth 1 -maxdepth 1 -type d -mmin +${params.slot_timeout} -exec rm -rf {} +
      n=\$(find "\$slots" -mindepth 1 -maxdepth 1 -type d | wc -l)
      if [ "\$n" -lt ${params.max_concurrent_samples} ]; then
        mkdir -p "\$slots/${unique_id}"
        break
      fi
      sleep 30
    done
    """
}

process release_slot {
    label "process_single"
    tag "${unique_id}"
    cache false

    input:
    val(unique_id)

    script:
    """
    rm -rf "${params.slot_dir}/${unique_id}"
    """
}

// Remote paths pass as strings, so they are not staged while a sample waits
workflow throttle_samples {
    take:
    sample_ch

    main:
    sample_ch.map{ id, f -> [id, f instanceof List ? f.collect{ it.toUriString() } : f.toUriString()] }
             .set{ source_ch }
    acquire_slot(source_ch)
    acquire_slot.out.map{ id, uri -> [id, uri instanceof List ? uri.collect{ file(it, checkIfExists: true) } : file(uri, checkIfExists: true)] }
                    .set{ gated_ch }

    emit:
    gated_ch
}

workflow release_finished_samples {
    take:
    done_ch
    n_methods

    main:
    done_ch.groupTuple(size: n_methods)
           .map{ id, classifiers -> id }
           .set{ finished_ch }
    release_slot(finished_ch)
}
