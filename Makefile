.PHONY: run clean-all test init

VENV_DIR ?= .env
PYTHON = $(VENV_DIR)/bin/python

run:
	clear
	$(PYTHON) samstat/samstat.py

init:
	rm -rf $(VENV_DIR)
	@$(MAKE) $(VENV_DIR)

clean:
	find . -iname "*.pyc" -delete
	find . -iname "*.pyo" -delete
	find . -iname "__pycache__" -delete
	-rm -rf SameStat.egg-info

clean-all:
	-rm -rf $(VENV_DIR)
	-rm -rf SamStat.egg-info
	-rm -rf dist
	-rm -rf build

test:
	$(PYTHON) -m unittest discover

travis-test:
	python -m unittest discover

$(VENV_DIR):
	virtualenv $(VENV_DIR)
	$(VENV_DIR)/bin/pip install --upgrade pip
	$(VENV_DIR)/bin/pip install --upgrade setuptools
	-$(VENV_DIR)/bin/pip install -r requirements.txt
