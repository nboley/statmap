'''
Mappings of the enums used in libstatmap.so

Most, but not all of these, are defined in src/config.h
'''

# enum bool
(false, true) = (0, 1)

# enum CHR_SOURCE
(REFERENCE, PATERNAL, MATERNAL) = (1, 2, 3)

# enum input_file_type_t
(SANGER_FQ, ILLUMINA_v13_FQ, SOLEXA_v14_FQ, ILLUMINA_v15_FQ, SOLEXA_LOG_ODDS_FQ, TEST_SUITE_FORMAT, MARKS_SOLEXA) \
        = (1, 2, 2, 2, 3, 2, 4)

# enum assay_type_t
(CAGE, CHIP_SEQ) = (1, 2)

# enum inputfile_type (rawread.h)
(FASTQ, ) = (1, )

# UNKNOWN
(UNKNOWN, ) = (0, )
