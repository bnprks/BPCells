#' Conclusions:
#' - Across the board, the new implementation is faster than the old one, so safe to transition
#' - There are slight increases in size of packed fragments in memory (order of 1%), which I'm willing
#'   to eat
#' - For HDF5 at least, increasing the buffer size somewhat appears to be highly beneficial. 65536
#'   seems to be about where diminishing returns start to really kick in
#' - With HDF5 properly buffered, it's within a factor of two from the dir speed

devtools::load_all("~/Dropbox/greenleaf/playground/fragment_io/BPCells/")

# 220s user, 265s elapsed (very slow!)
system.time({
  open_10x_fragments("~/Downloads/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz") %>%
    write_fragments_dir2("~/Downloads/pbmc_10k_frags/")
})

# 4.6s user, 4.8s system, 12.1s elapsed, 2.7GB memory (slow SSD?)
system.time({
  f <- open_fragments_dir2("~/Downloads/pbmc_10k_frags") %>%
    write_fragments_memory2(compress = FALSE)
})

system.time({
  s <- scan_fragments(f)
})
system.time({
  n <- nucleosome_counts(f)
})

# 5.2s user, 13.7s elapsed, another 2.7GB memory (somewhat confusing to me)
system.time({
  m_u1 <- write_fragments_memory(f, compressed=FALSE)
})

# 4.7s user, 15.9s elapsed
system.time({
  write_fragments_memory2(f, compress=FALSE)
})

# 2.1s user, 3.1s elapsed, 1.0GB memory
system.time({
  m_p1 <- write_fragments_memory(f, compressed=TRUE)
})

# 2.3s user, 3.5s elapsed, 1.0GB memory
system.time({
  m_p2 <- write_fragments_memory2(f, compress=TRUE)
})

# 2.1s user, 11.98s elapsed, 2.7GB disk
system.time({
  write_fragments_dir(f, "~/Downloads/frag_test/u1", compressed = FALSE)
})
d_u1 <- open_fragments_dir("~/Downloads/frag_test/u1", compressed = FALSE)

# 6.0s user, 13.2s elapsed, 984Mb disk
system.time({
  write_fragments_dir(f, "~/Downloads/frag_test/p1", compressed = TRUE)
})
d_p1 <- open_fragments_dir("~/Downloads/frag_test/p1", compressed = TRUE)

# 1.8s user, 7.5s elapsed, 2.7GB disk
system.time({
  d_u2 <- write_fragments_dir2(f, "~/Downloads/frag_test/u2", compress = FALSE)
})

# 1.5s user, 3.4s elapsed, 997Mb disk
system.time({
  d_p2 <- write_fragments_dir2(f, "~/Downloads/frag_test/p2", compress = TRUE)
})

# 57.5s user, 71s elapsed, 2.7G on disk
system.time({
  write_fragments_h5(f, "~/Downloads/frag_test/u1.h5", "unpacked", chunk_size=1000L, compressed=FALSE)
})
h_u1 <- open_fragments_h5("~/Downloads/frag_test/u1.h5", "unpacked", compressed=FALSE)

# 101s user, 111s elapsed, 1.0G on disk
system.time({
  write_fragments_h5(f, "~/Downloads/frag_test/p1.h5", "packed", chunk_size=1000L, compressed=TRUE)
})
h_p1 <- open_fragments_h5("~/Downloads/frag_test/p1.h5", "packed", compressed=TRUE)

# 7.3s user, 21.4s elapsed, 2.7G on disk
system.time({
  h_u2 <- write_fragments_h52(f, "~/Downloads/frag_test/u2.h5", "unpacked", compress=FALSE)
})

# 4.2s user, 8.2s elapsed, 1.0G on disk
system.time({
  h_p2 <- write_fragments_h52(f, "~/Downloads/frag_test/p2.h5", "packed", compress=TRUE)
})

scan_timing <- bench::mark(
  scan_fragments(f),
  scan_fragments(m_u1),
  scan_fragments(m_p1),
  scan_fragments(m_p2),
  scan_fragments(d_u1),
  scan_fragments(d_u2),
  scan_fragments(d_p1),
  scan_fragments(d_p2),
  scan_fragments(h_u2),
  scan_fragments(h_p2)
)

#' expression                min  median `itr/sec` mem_alloc `gc/sec` n_itr  n_gc total_time result memory  
#1 scan_fragments(f)       3.87s   3.87s    0.259     4.98KB        0     1     0      3.87s <dbl … <Rprofm…
#2 scan_fragments(m_u1)    4.91s   4.91s    0.204     4.98KB        0     1     0      4.91s <dbl … <Rprofm…
#3 scan_fragments(m_p1)     1.4s    1.4s    0.717     4.98KB        0     1     0       1.4s <dbl … <Rprofm…
#4 scan_fragments(m_p2)    834ms   834ms    1.20      4.98KB        0     1     0      834ms <dbl … <Rprofm…
#5 scan_fragments(d_u1)    2.87s   2.87s    0.348     4.98KB        0     1     0      2.87s <dbl … <Rprofm…
#6 scan_fragments(d_u2)    2.43s   2.43s    0.411     4.98KB        0     1     0      2.43s <dbl … <Rprofm…
#7 scan_fragments(d_p1)    2.07s   2.07s    0.484     4.98KB        0     1     0      2.07s <dbl … <Rprofm…
#8 scan_fragments(d_p2)    1.43s   1.43s    0.698     4.98KB        0     1     0      1.43s <dbl … <Rprofm…
#9 scan_fragments(h_u2)   34.09s  34.09s    0.0293    4.98KB        0     1     0     34.09s <dbl … <Rprofm…
#10 scan_fragments(h_p2)   23.12s  23.12s    0.0433    4.98KB        0     1     0     23.12s <dbl … <Rprofm…
# … with 2 more variables: time <list>, gc <list>

h_p2@buffer_size <- 32L * 8192L # From 8192L
system.time({
  scan_fragments(h_p2)
})