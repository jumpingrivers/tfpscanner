#' Extracts node-specific data and annotates whether they are internal or a tip in the tree
#'
#' @param   sc0   data.frame.
#' @param   tr2   phylo.
#' @param   branch_cols   Character vector defining the statistics that must be extracted from
#'   \code{sc0} for all nodes in the tree.
#' @return   data.frame. Containing a subset of the columns from \code{sc0} (including all
#'   \code{branch_cols}) and some additional columns. The additional columns are \code{lineages},
#'   \code{cocirc_summary}, \code{node}, \code{internal}, \code{lineages1}.

extract_tree_dataframe <- function(sc0,
                                   tr2,
                                   branch_cols) {
  tdvars <- unique(c(
    branch_cols,
    "logistic_growth_rate",
    "clock_outlier",
    "cluster_size",
    "date_range",
    "cluster_id",
    "region_summary",
    "cocirc_lineage_summary",
    "lineage",
    "tr2mrca"
  ))

  sc2 <- sc0[!is.na(sc0$tr2mrca), ]
  sc2$date_range <- sapply(
    seq_len(nrow(sc2)),
    function(i) glue::glue("{sc2$least_recent_tip[i]} -> {sc2$most_recent_tip[i]}")
  )

  ## tips
  td0 <- sc2[sc2$tr2mrca <= ape::Ntip(tr2), tdvars]
  td0$lineages <- td0$lineage
  td0$cocirc_summary <- td0$cocirc_lineage_summary
  td0$node <- td0$tr2mrca
  td0$internal <- "N"

  ## internal
  td1 <- sc2[sc2$tr2mrca > ape::Ntip(tr2), tdvars]
  if (nrow(td1) > 0) {
    td1$lineages <- td1$lineage
    td1$cocirc_summary <- td1$cocirc_lineage_summary
    td1$node <- td1$tr2mrca
    td1$internal <- "Y"
    td1$cluster_size <- 0
    x <- setdiff(
      (ape::Ntip(tr2) + 1):(ape::Ntip(tr2) + ape::Nnode(tr2)),
      td1$node
    ) # make sure every node represented
    td1 <- merge(td1,
      data.frame(node = x),
      all = TRUE
    )
    td <- rbind(td0, td1)
  } else {
    td <- td0
  }
  td <- td[order(td$node), ] # important

  # rescale clock ?
  td$clock_outlier <- scale(td$clock_outlier) / 2

  # interpolate missing values  &  repair cluster sizes
  td$logistic_growth_rate[(td$node <= ape::Ntip(tr2)) & (is.na(td$logistic_growth_rate))] <- 0
  td$clock_outlier[(td$node <= ape::Ntip(tr2)) & (is.na(td$clock_outlier))] <- 0
  for (ie in ape::postorder(tr2)) {
    a <- tr2$edge[ie, 1]
    u <- tr2$edge[ie, 2]
    td$cluster_size[a] <- td$cluster_size[a] + td$cluster_size[u]
    for (vn in branch_cols) {
      if (is.na(td[[vn]][a])) {
        td[[vn]][a] <- td[[vn]][u]
      }
    }
  }

  # lineages for clade labels
  td$lineages1 <- sapply(strsplit(td$lineages, split = "\\|"), "[", 1)

  td
}
