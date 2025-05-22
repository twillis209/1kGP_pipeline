variant_set_flags = "|".join(["sans_mhc", "sans_pars", "sans_at_gc"])

wildcard_constraints:
    snp_set = 'with_mhc|sans_mhc',
    chr = "chr[0-9XY]{1,2}",
    chr_no = "[0-9XY]{1,2}",
    assembly = "hg19|hg38",
    relatedness = "all|unrelated",
    ancestry = "eur|afr|amr|eas|sas|all",
    #variant_set = rf"^(?:all|(?:{variant_set_flags})(?:_and_(?:{variant_set_flags}))*)$",
    variant_set = "all|sans_pars_and_sans_long_range_ld",
    variant_type = "all|snps_only",
    window_size = "\\d+",
    r2 = "0_\\d+",
    maf = "\\d+",
    vmiss = "\\d+",
    post_vmiss = "\\d+",
    seed = "\\d+"
