bams <- dir("../../results/RNAseq/sortAlignment/v1/", pattern = "bam$", full.names = TRUE)
samples <- gsub(".bam", "", basename(bams), fixed = TRUE)
df <- data.frame(id = samples, path = bams, stringsAsFactors = FALSE)
df$Status <- ifelse(grepl("G3", df$id), "Control", "Case")
df[df$id == "G2-CS10-CO", "Status"] <- "G2CS10CO"
write.table(df, file = "ggsashimi_tab_G2CS10CO.tsv",
            quote = FALSE, row.names = FALSE, 
            col.names = FALSE, sep = "\t")
## Docker command
##docker run -v /home/SHARED/PROJECTS/CHD_MARATO/:/home/SHARED/PROJECTS/CHD_MARATO/ -w $PWD -v $PWD:$PWD guigolab/ggsashimi -b ggsashimi_tab_G2CS10CO.tsv \
-c 4:84515747-84519350 -M 10 -C 3 -O 3 -A median -F png -o G2CS10CO.png