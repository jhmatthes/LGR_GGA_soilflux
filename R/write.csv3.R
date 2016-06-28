#write.csv function with the ability to write a second header line of units

write.csv3 <- function(d, file) {
  c <- colnames(d)
  opts <- options(useFancyQuotes = FALSE)
  on.exit(options(opts))
  h1 <- paste(dQuote(c("", names(d))), collapse = ",")
  h2 <- paste(dQuote(c("", comment(d$a), comment(d$b))), collapse = ",")
  writeLines(paste(h1, h2, sep = "\n"), file)
  write.table(d, file, sep = ",", append = TRUE, col.names = FALSE)
}