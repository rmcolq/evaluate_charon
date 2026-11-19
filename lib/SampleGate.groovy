// Called by the errorStrategy of per-sample processes (conf/base.config) when a task has failed for good:
// frees the sample's slot (see modules/gate.nf) so a failed sample cannot block the queue, then ignores
// the error. Task names look like "run_charon (sample_id)", so every per-sample process must set
// tag "${unique_id}".
//
// This lives in lib/ rather than in the config because newer Nextflow versions (strict config syntax)
// do not allow functions to be defined in config files.
class SampleGate {
    static String failSample(task, String slotDir) {
        def m = (task.name =~ /\((.+)\)\s*$/)
        if (m) {
            def slot = new File("${slotDir}/${m[0][1]}")
            if (slot.exists()) slot.deleteDir()
        }
        return 'ignore'
    }
}
