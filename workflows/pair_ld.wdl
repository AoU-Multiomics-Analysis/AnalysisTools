version 1.0

# Calculate unphased dosage r2 in one batch for the full cohort and each population.
workflow PairLDByPopulation {
    input {
        File pgen
        File pvar
        File psam
        File variant_pairs
        File sample_populations
        File? pgi
        String docker_image
        Int cpu = 2
        Int memory_gb = 8
        Int extra_disk_gb = 10
        Boolean allow_small_groups = false
    }

    call CalculatePairLD {
        input:
            pgen = pgen,
            pvar = pvar,
            psam = psam,
            variant_pairs = variant_pairs,
            sample_populations = sample_populations,
            pgi = pgi,
            docker_image = docker_image,
            cpu = cpu,
            memory_gb = memory_gb,
            extra_disk_gb = extra_disk_gb,
            allow_small_groups = allow_small_groups
    }

    output {
        File ld_results = CalculatePairLD.ld_results
        File sample_counts = CalculatePairLD.sample_counts
        File plink_logs = CalculatePairLD.plink_logs
        File plink_version = CalculatePairLD.plink_version
        File command_stdout = CalculatePairLD.command_stdout
        File command_stderr = CalculatePairLD.command_stderr
    }
}

task CalculatePairLD {
    input {
        File pgen
        File pvar
        File psam
        File variant_pairs
        File sample_populations
        File? pgi
        String docker_image
        Int cpu
        Int memory_gb
        Int extra_disk_gb
        Boolean allow_small_groups
    }

    # Allow room for localized inputs and one subset of the PGEN data.
    Int disk_gb = ceil(2.0 * (size([pgen, pvar, psam, variant_pairs, sample_populations], "GB") + size(select_all([pgi]), "GB"))) + extra_disk_gb
    Int plink_memory_mb = memory_gb * 750

    command <<<
        set -euo pipefail
        echo "[PairLD] Start input checks and batched population LD calculation"
        python /opt/pair_ld/pair_ld.py \
            --pgen '~{sub(pgen, "'", "'\"'\"'")}' \
            --pvar '~{sub(pvar, "'", "'\"'\"'")}' \
            --psam '~{sub(psam, "'", "'\"'\"'")}' \
            --pairs '~{sub(variant_pairs, "'", "'\"'\"'")}' \
            --populations '~{sub(sample_populations, "'", "'\"'\"'")}' \
            ~{if defined(pgi) then "--pgi '" + sub(select_first([pgi]), "'", "'\"'\"'") + "'" else ""} \
            ~{if allow_small_groups then "--allow-small-groups" else ""} \
            --threads ~{cpu} \
            --memory-mb ~{plink_memory_mb} \
            --out-dir results
        echo "[PairLD] Complete. Results: results/ld.tsv"
    >>>

    runtime {
        docker: docker_image
        cpu: cpu
        memory: memory_gb + " GB"
        disks: "local-disk " + disk_gb + " SSD"
        preemptible: 1
        maxRetries: 1
    }

    output {
        File ld_results = "results/ld.tsv"
        File sample_counts = "results/sample_counts.tsv"
        File plink_logs = "results/plink_logs.tar.gz"
        File plink_version = "results/plink_version.txt"
        File command_stdout = stdout()
        File command_stderr = stderr()
    }
}
