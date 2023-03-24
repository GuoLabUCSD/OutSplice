# Function to normalize junction data by gene expression values

normalizeJunctions <- function(junc.RPM2, junctionGenegene_expr) {
  ## Turn all the zeros in junctionGenegene_expr into 1.
  ## replace ALL zeros with 1.  Then when you divide it does not do anything.
  junctionGenegene_expr <- junctionGenegene_expr[row.names(junc.RPM2), ]
  junctionGenegene_expr2 <- junctionGenegene_expr
  junctionGenegene_expr2[junctionGenegene_expr == 0] <- 1
  junc.RPM.norm <- junc.RPM2
  junc.RPM.norm <- junc.RPM2 / junctionGenegene_expr2

  return(junc.RPM.norm)
}
