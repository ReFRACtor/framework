list(APPEND FULL_TEST_REGULAR_LIST
    sample
    l2_full
    fts
    fts_sim
    tccon_sounding_1
    # These cases are very sensitive to certain machine architectures
    # and cause spurious differences, ie retrieval is not stable
    #tccon_sounding_2
    #tccon_sounding_3
    tccon_sounding_4
    tccon_sounding_5
    # This test uses ascii files, not really essential to test anymore
    #oco_simple_vs_26 
    oco2_2stream_fm
    oco2_sounding_1
    oco2_sounding_2
    oco2_sounding_3
    oco2_brdf
)
