CPP := mpicxx -O3 -std=c++11
LIBS := ARBFN/interchange.o

.PHONY:	check
check:
	@echo "Checking for libboost-json and mpicxx w/ C++11..."
	@$(CPP) tests/example_controller.cpp -c -o /dev/null

	@echo "Checking for python3..."
	@python3 --version > /dev/null

	@echo "Checking for Doxygen..."
	@which doxygen > /dev/null

	@echo "Checking for pip..."
	@pip --version > /dev/null

	@echo "Checking for mpi4py..."
	@pip list | grep mpi4py > /dev/null

	@echo "Checking for autopep8..."
	@pip list | grep autopep8 > /dev/null

	@echo "Checking for cmake..."
	@cmake --version > /dev/null

	@echo "Checking for clang-format..."
	@clang-format --version > /dev/null

	@echo "Environment is valid."

.PHONY:	docs
docs:	docs/paper/* docs/pres/*
	$(MAKE) -C docs/paper
	$(MAKE) -C docs/pres
	doxygen -q
	$(MAKE) -C latex
	cp latex/refman.pdf refman.pdf

.PHONY:	format
format:
	find . -type f \( -iname "*.cpp" -or -iname "*.hpp" \) \
		-exec clang-format -i "{}" \;
	find . -type f -iname "*.py" -exec \
		autopep8 --in-place --aggressive --aggressive "{}" \;

.PHONY:	test
test:
	@echo "Ensuring valid documentation..."
	doxygen -q

	@echo "Running package tests..."
	$(MAKE) -C tests $@

.PHONY:	clean
clean:
	find . -type f \( -iname '*.o' -or -iname '*.out' -or \
		-iname '*.so' \) -exec rm -f "{}" \;

################################################################
# Container launching stuff
################################################################

# For absolute path usage later
cwd := $(shell pwd)

.PHONY: docker
docker:
	docker build --tag 'arbfn' .
	docker run \
		--mount type=bind,source="${cwd}",target="/host" \
		-i \
		-t arbfn:latest

.PHONY: podman
podman:
	podman build --tag 'arbfn' .
	podman run \
		--mount type=bind,source="${cwd}",target="/host" \
		-i \
		-t arbfn:latest
