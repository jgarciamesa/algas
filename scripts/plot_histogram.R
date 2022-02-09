
basic_plot = function(counts_file) {
	counts = read.csv(file = counts_file, header = FALSE)
	counts$V2 = counts$V2/sum(counts$V2)
	pdf(file = paste0("freq_plot_", basename(counts_file), ".pdf"))
	plot(counts, type = 'p', main = 'Gap Size Frequency', xlab = 'Gap lengths', ylab = 'Frequency')
	dev.off()
}

if(!interactive()) {
	ARGS = commandArgs(trailingOnly = TRUE)
	basic_plot(ARGS[1])
}

