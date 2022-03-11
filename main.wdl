workflow SctkDoubletFinder {
    
    # An integer matrix of counts
    File raw_counts

    call runDoubletFinder {
        input:
            raw_counts = raw_counts
    }

    output {
        File doublet_removed_counts = runDoubletFinder.output_counts
        File doublet_ids = runDoubletFinder.output_ids
    }
}

task runDoubletFinder {
    File raw_counts

    String output_name_prefix = "sctk_doublet_finder_reduced_counts"
    String doublet_file_prefix = "sctk_doublet_finder_ids"

    Int disk_size = 20

    command {
        Rscript /opt/software/doubletfinder_qc.R \
        -f ${raw_counts} \
        -o ${output_name_prefix} \
        -d ${doublet_file_prefix}
    }

    output {
        File output_counts = glob("${output_name_prefix}*")[0]
        File output_ids = glob("${doublet_file_prefix}*")[0]
    }

    runtime {
        docker: "ghcr.io/web-mev/mev-sctk-doubletfinder"
        cpu: 2
        memory: "16 G"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 0
    }
}
