### For etc
file <- "samtools2SIC"
if(file.exists(file)){
	etcarch <- if (nzchar(R_ARCH)) paste("etc", R_ARCH, sep='') else "etc"
	dest <- file.path(R_PACKAGE_DIR, etcarch)
	dir.create(dest, recursive = TRUE, showWarnings = FALSE)
	file.copy(file, dest, overwrite = TRUE)
}

