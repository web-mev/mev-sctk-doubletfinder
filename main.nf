process run_doublet_finder {

    tag "Run SCTK DoubletFinder tool"
    publishDir "${params.output_dir}/SctkDoubletFinder.doublet_removed_counts", mode:"copy", pattern:"${output_name_prefix}*"
    publishDir "${params.output_dir}/SctkDoubletFinder.doublet_ids", mode:"copy", pattern:"${doublet_file_prefix}*"
    container "ghcr.io/web-mev/mev-sctk-doubletfinder"
    cpus 2
    memory '16 GB'

    input:
        path raw_counts

    output:
        path "${output_name_prefix}*"
        path "${doublet_file_prefix}*"

    script:
        output_name_prefix = "sctk_doublet_finder_reduced_counts"
        doublet_file_prefix = "sctk_doublet_finder_ids"
        """
        Rscript /opt/software/doubletfinder_qc.R \
            -f ${raw_counts} \
            -o ${output_name_prefix} \
            -d ${doublet_file_prefix}
        """
}

workflow {
    run_doublet_finder(params.raw_counts)
}