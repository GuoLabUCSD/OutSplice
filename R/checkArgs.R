# Function to validate function arguments

checkOSArgs <- function(junction, gene_expr, rawcounts, output_file_prefix, dir) {
  stopifnot("Junction File does not exist. Check path to file." = file.exists(junction))
  stopifnot("Expression File does not exist. Check path to file." = file.exists(gene_expr))
  stopifnot("Rawcounts File does not exist. Check path to file." = file.exists(rawcounts))
  stopifnot("Output File Prefix is not a string. Check quotations." = is.character(output_file_prefix))
  stopifnot("Path to Output Directory should end with a '/'" = endsWith(dir, "/"))
  stopifnot("Output Directory does not exist" = dir.exists(dir))
}
