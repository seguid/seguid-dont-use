SHELL=bash

all: test install

install:
	python -m pip install .

test:
	python -m pytest

coverage:
	python -m pytest --cov=src/

coverage-html:
	pytest --cov=src/ --cov-report=html; xdg-open htmlcov/index.html

cleanup:
	find . -type d -name "__pycache__" -exec rm -rf {} +

