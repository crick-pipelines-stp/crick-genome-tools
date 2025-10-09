# Crick Genome Tools

A comprehensive Python toolkit for genomic data analysis and processing at the Francis Crick Institute.

## Table of Contents

- [Overview](#overview)
- [Development Setup](#development-setup)
- [Development Workflow](#development-workflow)
- [Architecture](#architecture)
- [Testing](#testing)
- [Coding Standards](#coding-standards)
- [Project Structure](#project-structure)
- [Contributing](#contributing)
- [Maintenance](#maintenance)

## Overview

Crick Genome Tools is a Python-based toolkit designed to streamline genomic data analysis workflows. It provides a comprehensive suite of utilities for:

- **Genomic Data Processing**: Handling FASTA, FASTQ, and other genomic file formats
- **Quality Control**: Sequence quality assessment and filtering
- **Analysis Pipelines**: Automated genomic analysis workflows
- **Data Management**: Organizing and processing large-scale genomic datasets
- **Bioinformatics Utilities**: Common genomic data manipulation tasks

**Repository**: [pipelines-tech/crick-genome-tools](https://github.com/pipelines-tech/crick-genome-tools)

---

## Development Setup

### Prerequisites

- **Python**: 3.8+
- **UV**: Package manager for fast dependency management
- **Operating System**: Linux or macOS

### Initial Setup

1. **Clone the repository**:
   ```bash
   git clone <repository-url>
   cd crick-genome-tools
   ```

2. **Activate development environment**:
   ```bash
   . .venv/bin/activate && uv sync --group dev
   ```

This command activates the virtual environment and installs all dependencies including development tools.

### Verification

```bash
# Verify installation
python -c "import asf_tools; print('Installation successful')"

# Run tests to ensure everything works
pytest
```

## Development Workflow

### Test-Driven Development (TDD)

This project follows **Test-Driven Development** principles:

1. **Write tests first** before implementing functionality
2. **Run tests** to see them fail
3. **Write minimal code** to make tests pass
4. **Refactor** while keeping tests green
5. **Repeat** for each new feature

### Running Tests

```bash
# Run all tests
pytest

# Run tests with coverage
pytest --cov=asf_tools

# Run specific test file
pytest tests/test_core.py

# Run tests in watch mode during development
pytest --watch
```

**Important**: Tests are automatically run after each code change to ensure they pass.

### Code Quality Tools

The project uses several tools to maintain code quality:

```bash
# Format code
black .
isort .

# Lint code
ruff check .

# Run all quality checks
black . && isort . && ruff check . && pytest
```

## Architecture

### Module Overview

```
asf_tools/
├── core/                   # Core genomic processing modules
├── io/                     # File I/O and data management
├── quality/                # Quality control and filtering
├── pipelines/              # Analysis pipeline components
├── utils/                  # Utility functions and helpers
└── config/                 # Configuration management
```

### Core Components

#### Data Processing (`asf_tools.core`)
- Genomic sequence parsing and manipulation
- File format handling (FASTA, FASTQ, etc.)
- Sequence transformation utilities

#### I/O Operations (`asf_tools.io`)
- File reading and writing abstractions
- Data validation and error handling
- Batch processing capabilities

#### Quality Control (`asf_tools.quality`)
- Sequence quality assessment
- Filtering and trimming operations
- Quality metrics calculation

#### Pipeline Management (`asf_tools.pipelines`)
- Workflow orchestration
- Pipeline configuration
- Progress monitoring

### Data Flow

1. **Input Processing**: Parse genomic data files
2. **Quality Assessment**: Evaluate sequence quality
3. **Filtering**: Apply quality filters and transformations
4. **Analysis**: Run genomic analysis algorithms
5. **Output Generation**: Write results in specified formats

### Configuration Management

Configuration is managed through:
- Environment variables
- Configuration files (TOML format)
- Command-line arguments
- Default settings in code

## Testing

### Test Organization

Tests are organized in a flat structure mirroring the source code:

```
tests/
├── test_core_processor.py      # Core functionality
├── test_io_file_handler.py     # I/O operations
├── test_quality_filter.py      # Quality control
└── test_pipelines_workflow.py  # Pipeline components
```

### Testing Guidelines

#### Test Structure
```python
# tests/test_core.py
import pytest
from assertpy import assert_that

class TestGenomeProcessor:
    
    def test_should_process_valid_genome_file(self):
        # Setup
        processor = GenomeProcessor()
        test_file = "tests/fixtures/sample.fasta"
        
        # Test
        result = processor.process(test_file)
        
        # Assert
        assert_that(result).is_not_none()
        assert_that(result["sequences"]).is_length(5)
```

#### Test Best Practices
- **Descriptive names**: `test_should_process_valid_genome_file` not `test_process`
- **Clear structure**: Setup, Test, Assert sections
- **Isolation**: No I/O or network calls unless mocked
- **Parametrization**: Use `pytest.mark.parametrize` for multiple scenarios
- **Fixtures**: Use pytest fixtures for common setup/teardown
- **Coverage**: Aim for 100% test coverage

#### Running Tests

```bash
# All tests
pytest

# Specific module
pytest tests/test_core_processor.py

# With coverage
pytest --cov=asf_tools --cov-report=html

# Verbose output
pytest -v

# Stop on first failure
pytest -x
```

## Coding Standards

### General Principles

- **DRY (Don't Repeat Yourself)**: Avoid code duplication
- **KISS (Keep It Simple, Stupid)**: Prefer simple, readable solutions
- **Test everything**: All new code must include tests
- **Clear naming**: Use descriptive names for functions, variables, and classes
- **Document complex logic**: Use comments for complex algorithms, docstrings for public APIs

### Python Style Guidelines

We follow [PEP 8](https://peps.python.org/pep-0008/) with these specific conventions:

#### Code Formatting
- **Indentation**: 4 spaces (no tabs)
- **Line length**: 88 characters (Black default)
- **Quotes**: Double quotes for strings, triple double quotes for docstrings
- **Imports**: Grouped (standard library, third-party, local) and only at file top

#### Type Hints and Documentation
```python
def process_genome_data(file_path: str, quality_threshold: float = 0.8) -> dict[str, Any]:
    """Process genomic data from a file.
    
    Args:
        file_path: Path to the input file
        quality_threshold: Minimum quality score for filtering
        
    Returns:
        Dictionary containing processed results
        
    Raises:
        FileNotFoundError: If the input file doesn't exist
        ValueError: If quality_threshold is not between 0 and 1
    """
    # Implementation here
    pass
```

#### Modern Python Features
- Use **f-strings** for string interpolation: `f"Processing {filename}"`
- Use **list/dict comprehensions** where appropriate
- Use **type hints** in all function signatures
- Prefer **pathlib** over os.path for file operations

### Code Organization

- **Single responsibility**: Each module should have one clear purpose
- **Dependency injection**: Prefer explicit dependencies over global state
- **Error handling**: Use specific exceptions with clear messages
- **Logging**: Use structured logging for debugging and monitoring

### When to Add New Abstractions

Only introduce new patterns or abstractions when:
1. **Clear benefit**: Reduces complexity or improves reusability
2. **Justified rationale**: Document why in code comments or PR description
3. **Not premature**: Focus on clarity and correctness first

## Project Structure

```
crick-genome-tools/
├── asf_tools/              # Main source code directory
│   ├── __init__.py         # Package initialization
│   ├── core/               # Core processing modules
│   ├── io/                 # Input/output operations
│   ├── quality/            # Quality control
│   ├── pipelines/          # Pipeline components
│   ├── utils/              # Utility functions
│   └── config/             # Configuration management
├── tests/                  # Test directory (flat structure)
│   ├── fixtures/           # Test data files
│   └── test_*.py           # Test modules
├── .github/               # GitHub workflows and configurations
├── .venv/                 # Virtual environment (local development)
├── pyproject.toml         # Project configuration and dependencies
├── pytest.ini            # Test configuration
└── README.md              # This file
```

### Key Files

- **`pyproject.toml`**: Project metadata, dependencies, and tool configuration
- **`asf_tools/`**: Main source code with modular organization
- **`tests/`**: Comprehensive test suite with pytest
- **`.github/copilot-instructions.md`**: Development guidelines and conventions

## Contributing

### Before Starting Development

1. **Discuss the task**: Always discuss with the team before writing code
2. **Start with design**: High-level architecture before implementation details
3. **Ask for clarification**: If requirements are unclear
4. **Break down complex tasks**: Into smaller, manageable pieces
5. **Get approval**: Don't start coding without a green light

### Development Process

1. **Create feature branch**: `git checkout -b feature/your-feature-name`
2. **Write tests first**: Following TDD principles
3. **Implement functionality**: Make tests pass with minimal code
4. **Run quality checks**: `black . && isort . && ruff check . && pytest`
5. **Commit changes**: With clear, descriptive messages
6. **Create Pull Request**: With comprehensive description

### Pull Request Guidelines

#### Before Submitting
- [ ] All tests pass (`pytest`)
- [ ] Code is formatted (`black . && isort .`)
- [ ] No linting errors (`ruff check .`)
- [ ] Test coverage is sufficient (aim for 100%)
- [ ] Documentation is updated if needed

#### PR Description Should Include
- **Clear summary** of changes made
- **Rationale** for new patterns or abstractions
- **API-breaking changes** explicitly highlighted
- **Testing approach** and coverage details

### Commit Message Format

```
type(scope): brief description

Detailed explanation if needed.

- API-breaking changes highlighted here
- Any special considerations
```

Example:
```
feat(genome): add quality filtering for FASTA sequences

Implement quality score filtering with configurable thresholds.

- BREAKING: GenomeProcessor.process() now requires quality_threshold parameter
- Added comprehensive test coverage for edge cases
```

## Maintenance

### Regular Tasks

- **Update dependencies**: `uv sync --upgrade`
- **Run full test suite**: `pytest`
- **Check security**: `uv audit`
- **Review coverage**: `pytest --cov=asf_tools --cov-report=html`

### Performance Considerations

- **Profile before optimizing**: Use tools like `cProfile` when needed
- **Measure impact**: Benchmark performance changes
- **Document trade-offs**: Explain performance vs. readability decisions

### Getting Help

#### Common Issues

1. **Tests failing after changes**: Run `pytest -v` for detailed output
2. **Import errors**: Ensure you're in the virtual environment
3. **Formatting issues**: Run `black . && isort .` to auto-fix
4. **Type errors**: Check type hints match actual usage

#### Resources

- **Project documentation**: Check docstrings in source code
- **Test examples**: Look at existing tests for patterns
- **GitHub Issues**: Search existing issues for similar problems
- **Team discussion**: Reach out to team members for architectural questions

---

**Happy coding!** 🧬✨

For questions or suggestions about this guide, please open an issue or reach out to the development team.
