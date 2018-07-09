## Hydrographs for Pepsi 2

for (i in 1:length(pep2_qobs)) {
  namei <- names(pep2_qobs)[i]
  
  hgi <- bam_hydrograph(bamsamps[[i]], qobs = pep2_qobs[[i]]) +
    ggtitle(namei) + 
    theme_minimal(base_size = 20, base_line_size = 1)
  
  ggsave(hgi, filename = sprintf("graphs/man_hgraphs/%s.png", namei))
}


for (i in 1:length(pep2_qobs)) {
  namei <- names(pep2_qobs)[i]
  
  hgi <- bam_hydrograph(bamsamps_man_amhg[[i]], qobs = pep2_qobs[[i]]) +
    ggtitle(namei) + 
    theme_minimal(base_size = 20, base_line_size = 1)
  
  ggsave(hgi, filename = sprintf("graphs/man_amhg_hgraphs/%s.png", namei))
}
