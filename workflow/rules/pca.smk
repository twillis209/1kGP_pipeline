rule pca:
    input:
        rules.prune_variants.output
    output:
        "results/1kG/{assembly}/{relatedness}/{ancestry}/{variant_type}/{maf}/qc/{variant_set}/pca/{window_size}_1_{r2}/merged.pca.eigenvec",
        "results/1kG/{assembly}/{relatedness}/{ancestry}/{variant_type}/{maf}/qc/{variant_set}/pca/{window_size}_1_{r2}/merged.pca.eigenval"
