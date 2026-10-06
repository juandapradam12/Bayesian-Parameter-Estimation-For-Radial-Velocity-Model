.PHONY: install test run quick mock clean legacy-build

PYTHON ?= python3

install:
	$(PYTHON) -m pip install -r requirements.txt
	$(PYTHON) -m pip install -e .

test:
	$(PYTHON) -m pytest -q

quick:
	$(PYTHON) scripts/run_inference.py --quick --fit-ah

run:
	$(PYTHON) scripts/run_inference.py --fit-ah --fit-sigma

mock:
	$(PYTHON) scripts/mock_recovery.py

# Historical C reference (fixed build of the original homework code)
legacy-build:
	$(MAKE) -C legacy -f Tarea5.mk CurvaRotacion.x

clean:
	rm -rf results/*.png results/*.json results/*.txt results/*.npy
	rm -rf legacy/*.dat legacy/*.x legacy/*.png legacy/*.log legacy/*.aux
	find . -type d -name __pycache__ -exec rm -rf {} + 2>/dev/null || true
