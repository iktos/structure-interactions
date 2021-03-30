TODO change Your-package-name
# Your-package-name

## Adapt
To use this template, find and resolve all TODO in the repo
Then remove this paragraph

## Installation

### Poetry

[After installing poetry](https://python-poetry.org/docs/#installation),
 you can install the environment:
```
poetry install
```

Don't forget to update the lock file if the dependencies have change:
```
poetry lock
```

### Pre-commit

[Pre-commit](https://pre-commit.com/) add a hook before your commit
and run flake8, black, mypy and other checks.

To install pre-commit:
```
brew install pre-commit
```

To set the hook on the repository
```bash
pre-commit install
```

## Environnement

Useful commands:
```
# Run tests
poetry run task test

# Indent your files with black
poetry run task black

# Check your code with mypy
poetry run task mypy

# Check style
poetry run task linter
```
