MODULE=pychopper

.PHONY: clean clean-test clean-pyc clean-build docs com help 

.DEFAULT_GOAL := help

define PRINT_HELP_PYSCRIPT
import re, sys

for line in sys.stdin:
	match = re.match(r'^([a-zA-Z_-]+):.*?## (.*)$$', line)
	if match:
		target, help = match.groups()
		print("%-20s %s" % (target, help))
endef
export PRINT_HELP_PYSCRIPT

help:
	@python -c "$$PRINT_HELP_PYSCRIPT" < $(MAKEFILE_LIST)

clean: clean-build clean-pyc clean-test ## remove all build, test, coverage and Python artifacts


clean-build: ## remove build artifacts
	rm -fr build/
	rm -fr dist/
	rm -fr .eggs/
	find . -name '*.egg-info' -exec rm -fr {} +
	find . -name '*.egg' -exec rm -f {} +

clean-pyc: ## remove Python file artifacts
	find . -name '*.pyc' -exec rm -f {} +
	find . -name '*.pyo' -exec rm -f {} +
	find . -name '*~' -exec rm -f {} +
	find . -name '__pycache__' -exec rm -fr {} +

clean-test: ## remove test and coverage artifacts
	rm -f .coverage
	rm -fr htmlcov/

lint: ## check style with flake8
	@(flake8 --max-line-length=120 $(MODULE) | grep -v "E501 line too long") || true
	@(flake8 --max-line-length=120 scripts/*.py | grep -v "E501 line too long") || true

test: ## run tests quickly with the default Python
	py.test

coverage: ## check code coverage quickly with the default Python
		coverage run --source $(MODULE) --omit="*/tests/*,*__init__.py" `which py.test`
		coverage report -m --omit="*/tests/*,*__init__.py"
		coverage html

docs: ## generate Sphinx HTML documentation, including API docs
	@cd docs; make clean html

servedocs: docs ## compile the docs watching for changes
	watchmedo shell-command -p '*.rst' -c '$(MAKE) -C docs html' -R -D .

release: clean ## package and upload a release
	python setup.py sdist upload
	python setup.py bdist_wheel upload

dist: clean ## builds source and wheel package
	python setup.py sdist
	python setup.py bdist_wheel
	ls -l dist

install: clean ## install the package to the active Python's site-packages
	python setup.py install

com: ## commit all changes to git
	git commit -a

it: ## integration test
	./scripts/cdna_classifier.py -s 95 -i fasta -b pychopper/tests/data/barcodes.fas pychopper/tests/data/ref.fas pychopper/tests/data/test_output.fas
