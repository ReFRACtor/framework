test_data_dir = os.getenv("abs_top_srcdir") .. "/test/unit/data/"
sid_string = "2014090915251774"

l1b_fname = test_data_dir .. "in/common/l1b_example_data.h5"
l1b_hdf = HdfFile(l1b_fname)
l1b = ExampleLevel1b(l1b_hdf, sid_string)
