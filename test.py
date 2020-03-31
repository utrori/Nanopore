from rDNA_structure_of_long_reads import analyze_split_reads, easy_flag, analyze_direction_distance

samfile = ('temp_sam.sam')
analyze_split_reads(samfile, 200)
