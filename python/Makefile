# Configuration
PACKAGE_NAME=tdavec
PYPI_DIST=dist/*
README=README.md

.PHONY: clean build testpypi pypi

# Remove build artifacts
clean:
	rm -rf build dist *.egg-info tdavec/__pycache__/ tdavec/tdavec_core.c

# Build the package (source and wheel)
build: clean
	python -m build

# Upload to Test PyPI
testpypi: build
	twine upload --repository-url https://test.pypi.org/legacy/ $(PYPI_DIST)

# Upload to real PyPI
pypi: build
	twine upload $(PYPI_DIST)
