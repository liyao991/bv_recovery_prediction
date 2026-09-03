# test_shortread_parse.R - how ShortRead parses different header formats
suppressWarnings(Sys.setlocale("LC_CTYPE", "Chinese (Simplified)_China.936"))
suppressPackageStartupMessages({ library(ShortRead) })

# write two test fasta files: (a) accession + space + taxonomy (DADA2 standard)
#                              (b) taxonomy only (no accession)
fa1 <- tempfile(fileext = ".fa")
writeLines(c(">AB000000.1.1234 Bacteria;Proteobacteria;Gammaproteobacteria;Enterobacterales;Enterobacteriaceae;Escherichia;Escherichia coli",
             "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT",
             ">CD000000.1.5678 Bacteria;Firmicutes;Bacilli;Lactobacillales;Lactobacillaceae;Lactobacillus;Lactobacillus crispatus",
             "TGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC"), fa1)
ref1 <- readFasta(fa1)
cat("=== format A: accession + taxonomy ===\n")
print(as.character(id(ref1)))
print(as.character(description(ref1)))

fa2 <- tempfile(fileext = ".fa")
writeLines(c(">Bacteria;Proteobacteria;Gammaproteobacteria;Enterobacterales;Enterobacteriaceae;Escherichia;Escherichia coli",
             "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT",
             ">Bacteria;Firmicutes;Bacilli;Lactobacillales;Lactobacillaceae;Lactobacillus;Lactobacillus crispatus",
             "TGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC"), fa2)
ref2 <- readFasta(fa2)
cat("=== format B: taxonomy only ===\n")
print(as.character(id(ref2)))
print(as.character(description(ref2)))
cat("== PARSE_TEST_DONE ==\n")
