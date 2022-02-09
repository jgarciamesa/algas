
basic_plot = function(counts_file) {
	counts = read.csv(file = counts_file, header = FALSE)
	pdf(file = "gap_histogram.pdf")
	plot(counts, type = 'h', main = 'Gap histogram', xlab = 'Gap lengths', ylab = 'Frequency')
	dev.off()
}

if(!interactive()) {
	ARGS = commandArgs(trailingOnly = TRUE)
	basic_plot(ARGS[1])
}

