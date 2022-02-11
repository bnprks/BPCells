
devtools::load_all("~/Dropbox/greenleaf/playground/fragment_io/BPCells/")

input <- open_10x_fragments("~/Dropbox/greenleaf/playground/fragment_io/BPCells/tests/data/mini_fragments.tsv.gz")

d_u <- open_fragments_dir2("~/Downloads/test_mini_fragments/")
m_u <- write_fragments_memory(input, compressed = FALSE)
d_u <- write_fragments_dir2(input, "~/Downloads/test_mini_fragments", compress=FALSE)
as(d_u, "GRanges")
all.equal(as(input, "GRanges"), as(d_u, "GRanges"))

h_u <- write_fragments_h52(d_u, "~/Downloads/test_mini_fragments.h5", compress = FALSE)
all.equal(as(input, "GRanges"), as(h_u, "GRanges"))

all.equal(as(input, "GRanges"), as(m_u, "GRanges"))



# get dir array
gda <- function(s) {
  readBin(file.path(d_u@dir, s), integer(), n=15000, endian="big")[c(-1,-2)]
}
gda("start") %>% tail()
gda("end_max") %>% tail()

m <- gda("end_max")

gr <- as(input, "GRanges")
table(gr@seqnames)

system.time({
  open_10x_fragments("~/Downloads/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz") %>%
    write_fragments_dir("~/Downloads/pbmc_10k_frags/")
})

system.time(
  open_fragments_dir("~/Downloads/pbmc_10k_frags/") %>%
    nucleosome_counts()
)
# m <- open_fragments_dir("~/Downloads/pbmc_10k_frags/") %>% write_fragments_memory()
# system.time(
#   nucleosome_counts(m)
# )


open_fragments_dir("~/Downloads/pbmc_10k_frags/") %>% 
  write_fragments_dir("~/Downloads/pbmc_10k_frags_unpacked/", compressed=FALSE)

system.time(
  x1 <- nucleosome_counts(open_fragments_dir("~/Downloads/pbmc_10k_frags/"))
)
system.time(
  x2 <- nucleosome_counts(open_fragments_dir("~/Downloads/pbmc_10k_frags_unpacked/", compressed=FALSE))
)
system.time(
  nucleosome_counts(d)
)
r <- write_fragments_memory(d)

mean(diff(r@fragments[[1]]$cell_idx)/4)

# ~23s
system.time({
raw <- open_10x_fragments("~/Downloads/pbmc_granulocyte_sorted_3k_atac_fragments.tsv.gz") %>%
  write_fragments_memory(compressed=FALSE)
})

# 28.7 seconds, 229.4MB large
# revised: 20.9s, 18s user
system.time({
  write_fragments_h5(raw, "~/Downloads/fragments_test_1K.h5", "packed", chunk_size=1000, compressed = TRUE)
})

# 30 seconds, 237MB large
# revised: 18.7s, 17.4 user
system.time({
  write_fragments_h5(raw, "~/Downloads/fragments_test_10K.h5", "packed", chunk_size=10000, compressed = TRUE)
})

# revised: 18.6s, 17.4 user, 346MB large
system.time({
  write_fragments_h5(raw, "~/Downloads/fragments_test_100K.h5", "packed", chunk_size=100000, compressed = TRUE)
})

h1 <- open_fragments_h5("~/Downloads/fragments_test_1K.h5", "packed")
h10 <- open_fragments_h5("~/Downloads/fragments_test_10K.h5", "packed")
h100 <- open_fragments_h5("~/Downloads/fragments_test_100K.h5", "packed")

# 21 seconds, but only 14 for user
# revised: 19.8s, 13 for user
system.time({
  n5 <- scan_fragments_modulo_cpp(iterate_fragments(h1))
})

# 16.4 seconds
system.time({
  n6 <- scan_fragments_modulo_cpp(iterate_fragments(h10))
})

# 11.8 seconds user
system.time({
  n7 <- scan_fragments_modulo_cpp(iterate_fragments(h100))
})

# 2.28 seconds, 1.2 user
system.time({
  write_fragments_dir(raw, "~/Downloads/fragments_test/packed")
})
d <- open_fragments_dir("~/Downloads/fragments_test/packed")

# 0.3 user, 0.7 total
system.time({
  n7 <- scan_fragments_modulo_cpp(iterate_fragments(d))
})

# 2.9 seconds, 0.57 user
system.time({
  write_fragments_dir(raw, "~/Downloads/fragments_test/unpacked", compressed = FALSE)
})
du <- open_fragments_dir("~/Downloads/fragments_test/unpacked", compressed = FALSE)

# 0.26 user, 0.836 total
system.time({
  n7 <- scan_fragments_modulo_cpp(iterate_fragments(du))
})




# ~0.9-1s
system.time({
  packed <- write_fragments_memory(raw)
})

system.time({
  raw_v1 <- write_raw_fragments(raw)
})

# 0.9-1s
system.time({
  packed_v1 <- write_packed_fragments(raw_v1)
})

# 0.2 - 0.26s
system.time({
  n1 <- scan_fragments_modulo_cpp(iterate_fragments(raw))
})

# basically identical
system.time({
  n2 <- scan_fragments_modulo_cpp(iterate_fragments(raw_v1))
})

# 0.179-0.184 (faster??)
system.time({
  n3 <- scan_fragments_modulo_cpp(iterate_packed_fragments2_cpp(packed))
})

# 0.13 - 0.14ish
system.time({
  n4 <- scan_fragments_modulo_cpp(iterate_fragments(packed_v1))
})

system.time({
  x1 <- nucleosome_counts(packed)
})

bench::mark(
  v1 = scan_fragments_modulo_cpp(iterate_fragments(packed_v1)),
  v2 = scan_fragments_modulo_cpp(iterate_fragments(packed))
)
bench::mark(
  v2 = nucleosome_counts(packed),
  v1 = nucleosome_counts(packed_v1)
)
  x1 <- nucleosome_counts(packed_v1)
})


r2 <- RawFragments2(r)

r3 <- RawFragments(r)

all.equal(as(r2, "GRanges"), r)

attributes <- c("start", "end", "cell_id")
all(purrr::map_lgl(seq_along(r2@fragments), function(i) 
  all.equal(r2@fragments[[i]][attributes], r3@fragments[[i]][attributes])))

r_return <- as(r2, "GRanges")
