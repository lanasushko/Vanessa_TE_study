# cd-hit-est run to reduce library redundancy

# concatenate individual libraries
cat Vcard_curated_sequences_NR.fa Vatal_curated_sequences_NR.fa Vtam_newcurated_sequences_NR.fa Aio_curated_sequences_NR.fa Arth_mon_curated_RMDLlib_v2.2_nowrap.fasta > full_lib_curated_concat.fa

# run cd-hit-est
cd-hit-est -i full_lib_curated_concat.fa -o clustered_full_lib.fa -T 5 -aS 0.8 -aL 0.8
