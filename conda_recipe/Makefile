CONDA_BUILD = conda-build --output-folder refractor-conda-channel/

message:
	@echo "======================================================"
	@echo "Often you want to first do"
	@echo "  conda clean --all"
	@echo "  conda-build purge-all"
	@echo "Without this, may be stale packages in cache"
	@echo "which will cause mysterious problems"
	@echo ""
	@echo "On the Mac, if you have fink or macports installed make"
	@echo "sure we aren't picking up those versions of system tools"
	@echo "like as. You can set your path to something like:"
	@echo ""
	@echo 'export PATH=$${CONDA_PREFIX}/bin:/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin'
	@echo "======================================================"


ENV_NAME = framework-development-env
CHANNEL = $(shell pwd)/refractor-conda-channel
UNAME := $(shell uname)
ifeq ($(UNAME), Linux)
COMPILER_LIST = gcc_linux-64 gfortran_linux-64 gxx_linux-64 
endif
ifeq ($(UNAME), Darwin)
COMPILER_LIST = clang_osx-64 clangxx_osx-64 gfortran_osx-64
endif

create-framework-development-env:
	@echo "This has been pretty much repalced with refractor-development-tools"
	@echo "We will leave this in place for a little while in case it"
	@echo "proves useful (e.g., package we forgot in refractor-development-tools)"
	-conda env remove --yes -n $(ENV_NAME)
	conda create --channel $(CHANNEL) --yes -n $(ENV_NAME) $(COMPILER_LIST) cmake ninja-build swig doxygen gsl hdf5 boost pkg-config eigen lua luabind libblitz lidort twostream first_order ptpython jupyter ipython matplotlib netcdf4 attrdict pytest pytest-xdist pytest-cov pytest-runner

Manifest_linux.yml:
	conda env export > $@

Manifest_osx.yml:
	conda env export > $@

build-development-tools:
	$(CONDA_BUILD) refractor-development-tools

build-cris: message
	$(CONDA_BUILD) refractor-cris

build-oco: message
	$(CONDA_BUILD) refractor-oco

build-omps-nm: message
	$(CONDA_BUILD) refractor-omps-nm

build-framework: message
	$(CONDA_BUILD) refractor-framework

build-ninja-build: message
	$(CONDA_BUILD) ninja-build

build-ptpython: message
	$(CONDA_BUILD) ptpython

build-luabind: message
	$(CONDA_BUILD) luabind

build-lidort: message
	$(CONDA_BUILD) lidort

build-first_order: message
	$(CONDA_BUILD) first_order

build-twostream: message
	$(CONDA_BUILD) twostream

build-libblitz: message
	$(CONDA_BUILD) libblitz


