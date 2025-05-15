library(data.table)

if (grepl("hg38", snakemake@input[["ped"]])) {
  dat <- fread(snakemake@input[["ped"]], header = TRUE)

  if (snakemake@wildcards$relatedness == "unrelated") {
    dat <- dat[FatherID == "0" & MotherID == "0"]
  }

  for (x in c("eur", "afr", "amr", "eas", "sas")) {
    fwrite(dat[Superpopulation == toupper(x), .(SampleID)], file = file.path(snakemake@params$out_dir, sprintf("%s.samples", x)), sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
  }

  fwrite(dat[, .(SampleID)], file = file.path(snakemake@params$out_dir, "all.samples"), sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
} else if (grepl("hg19", snakemake@input[["ped"]])) {
  ped_dat <- fread(snakemake@input[["ped"]], header = T, sep = "\t")

  panel_dat <- fread(snakemake@input[["panel"]], header = F)

  names(panel_dat) <- c("ID", "Population", "Superpopulation", "Sex")

  merged_dat <- merge(panel_dat, ped_dat[, c("Individual ID", "Paternal ID", "Maternal ID")], by.x = "ID", by.y = "Individual ID", all.x = T)

  if (snakemake@wildcards$relatedness == "unrelated") {
    merged_dat <- merged_dat[`Paternal ID` == "0" & `Maternal ID` == "0"]
  }

  for (x in c("eur", "afr", "amr", "eas", "sas")) {
    fwrite(merged_dat[Superpopulation == toupper(x), .(ID)], file = file.path(snakemake@params$out_dir, sprintf("%s.samples", x)), sep = " ", col.names = FALSE, row.names = FALSE, quote = FALSE)
  }

  fwrite(merged_dat[, .(ID)], file = file.path(snakemake@params$out_dir, "all.samples"), sep = " ", col.names = FALSE, row.names = FALSE, quote = FALSE)
} else {
  stop("Input file contains neither hg38 nor hg19.")
}
