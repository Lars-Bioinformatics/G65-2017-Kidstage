####################################
###         Load Packages        ###

# Declare packages
packages <- c("stringr")

# Function to install and load packages
usePackage <- function( p ) {
  if ( !is.element(p, installed.packages()[,1]) )
    install.packages(p, dep = TRUE, repos = c("https://mirrors.dotsrc.org/cran/", "http://cran.us.r-project.org"))
  require(p, character.only = TRUE)
}

# Loop through packages
for ( i in 1:length(packages) ) {
  usePackage(packages[i])
}



####################################
###     Initial preparations     ###

args = commandArgs(trailingOnly=TRUE)

filepath_vcf <- args[1]
if ( !file.exists(args[2]) ) {
  file.create(args[2])
}
filepath_output <- args[2]

# Establish connection to input files
vcf = file(filepath_vcf, "r")
output = file(filepath_output, "w")

# Skip (and save) comment lines in vcf file
header_vcf <- list()
repeat {
  line_vcf <- readLines(vcf, n = 1)
  if ( str_sub(line_vcf, 1, 2) == "##" ) {
    header_vcf[[length(header_vcf) + 1]] <- line_vcf
  } else {
    header_vcf[[length(header_vcf) + 1]] <- line_vcf
    break
  }
}

# Write vcf header to output file
writeLines(unlist(header_vcf), output, "\n")



####################################
###             Main             ###

# Logical condition whether to terminate main loop
terminateMain <- FALSE

# Main Loop
repeat {
  # Read next line in vcf file
  line_vcf <- readLines(vcf, n = 1)
  # Break if 
  if ( length(line_vcf) == 0 ) {
    break
  }
  # Split line
  line_vcf <- str_split(line_vcf, "\t")
  # Advance if the variant is not flagged with PASS
  if ( !all(line_vcf[[1]][7] == "PASS") ) {
    next
  }
  # Advance if more than one alternative allele exists
  if ( length(unlist(str_split(line_vcf[[1]][5], ","))) > 1 ) {
    next
  }
  # Write variant to output file if it has passed the filtering 
  writeLines(str_c(line_vcf[[1]], collapse = "\t"), output, "\n")
}

# Terminate connection
close(vcf)
close(output)
