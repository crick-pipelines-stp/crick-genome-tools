set shell := ["bash", "-eu", "-o", "pipefail", "-c"]

default: dev

dev:
    if [ ! -d .venv ]; then uv venv .venv; fi
    source .venv/bin/activate && uv sync --group dev && python -V && which python && exec $SHELL

test:
    source .venv/bin/activate && pytest

test-cli:
    source .venv/bin/activate && pytest -o log_cli=true -o log_cli_level=INFO

lint:
    source .venv/bin/activate && ruff check .

python-upgrade:
    asdf plugin update python
    ver="$(asdf latest python 3.13)"; \
    echo ">>> Upgrading to $ver"; \
    asdf install python "$ver"; \
    asdf set python "$ver"; \
    rm -rf .venv; \
    uv venv .venv; \
    source .venv/bin/activate && uv sync && python -V
