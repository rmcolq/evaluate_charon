#!/usr/bin/env nextflow

/*
 * Limits how many samples are in flight at once so that large inputs do not fill the disk.
 *
 * A "slot" is a directory <slot_dir>/<unique_id>. acquire_slot only lets a sample through when
 * fewer than params.max_concurrent_samples slots exist. The slot is removed by release_slot once
 * every enabled method has finished the sample, or by fail_sample() in conf/base.config when a
 * task for that sample fails for good. Slots older than params.slot_timeout minutes are reclaimed
 * as a safety net (e.g. after a killed run). Runs on the host (no container) so that slot_dir is
 * shared with the Nextflow head process.
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

// [id, path] -> [id, path], emitted only as slots become free. Remote paths are not
// staged here: the source is passed as a string and re-wrapped as a path afterwards.
workflow throttle_samples {
    take:
    sample_ch

    main:
    sample_ch.map{ id, f -> [id, f.toUriString()] }
             .set{ source_ch }
    acquire_slot(source_ch)
    acquire_slot.out.map{ id, uri -> [id, file(uri, checkIfExists: true)] }
                    .set{ gated_ch }

    emit:
    gated_ch
}

// Takes [id, classifier] "done" events and frees the sample's slot once all enabled
// methods have reported.
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
