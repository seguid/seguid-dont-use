SHELL=bash

all: test install

install:
	python -m pip install .

test:
	python -m pytest

coverage:
	(cd Python/; python -m pytest --cov=src/)

coverage-html:
	(cd Python/; pytest --cov=src/ --cov-report=html; xdg-open htmlcov/index.html)



