
basic_plot = function(counts_file) {
	counts = read.csv(file = counts_file, header = FALSE)

	# count total length of gaps (total number of "-") and number of gaps
	total_len_gaps = sum(counts$V2)
	num_gaps = nrow(counts)

	# convert gap counts to percentages
	counts$V2 = counts$V2/sum(counts$V2)

	 # plotting
	pdf(file = paste0("freq_plot_", basename(counts_file), ".pdf"))
	plot(counts, type = 'p', main = 'Gap Size Frequency', xlab = 'Gap lengths', ylab = 'Frequency')
	mtext(text = paste("Number of gaps", num_gaps, ", total gap length", total_len_gaps), side = 3)
	dev.off()
}

len3 = function(counts_file) {
    counts = read.csv(file = counts_file, header = FALSE)

    # count number of gaps that are of length multiple of 3 and which aren't
    no_len3 = 0
    len3 = 0
    for(i in 1:nrow(counts)) {
        if(counts$V1[i] %% 3 == 0) {
            no_len3 = no_len3 + counts$V2[i]
        } else {
            len3 = len3 + counts$V2[i]
        }
    }

    # convert counts to percentages
    len3 = round((len3/sum(counts$V2))*100)
    no_len3 = round((no_len3/sum(counts$V2))*100)

    # display results
    print(paste0(len3, "% length multiple of 3 gaps and ",
                no_len3, "% length NOT multiple of 3 gaps"),
          quote = FALSE)
}

if(!interactive()) {
	ARGS = commandArgs(trailingOnly = TRUE)
	if(ARGS[1] == "histogram") {
	    basic_plot(ARGS[2])
	} else if(ARGS[1] == "len3") {
	    len3(ARGS[2])
	} else {
	    print("Command not supported.")
	}
}

