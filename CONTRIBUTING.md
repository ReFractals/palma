# Contributing to PALMA

Thank you for your interest in contributing to PALMA! This document provides guidelines for contributors.

## Code of Conduct

Please read and follow our [Code of Conduct](CODE_OF_CONDUCT.md).

## How to Contribute

### Reporting Bugs

1. Check existing [issues](https://github.com/axiom-research/palma/issues) to avoid duplicates
2. Use the bug report template
3. Include:
   - PALMA version
   - Platform (e.g., Raspberry Pi 4, ARM64)
   - Compiler version
   - Minimal reproduction code
   - Expected vs actual behavior

### Suggesting Features

1. Open an issue with the feature request template
2. Describe the use case and motivation
3. If possible, outline a proposed implementation

### Pull Requests

1. Fork the repository
2. Create a feature branch: `git checkout -b feature/your-feature`
3. Make your changes
4. Ensure all tests pass: `make test`
5. Follow the coding style (see below)
6. Write clear commit messages
7. Submit a pull request

## Development Setup

```bash
# Clone your fork
git clone https://github.com/YOUR_USERNAME/palma.git
cd palma

# Build with debug symbols
make debug

# Run tests
make test

# Run benchmarks
make run-benchmark
```

## Coding Style

PALMA follows a consistent C style:

### Naming
- Functions: `palma_module_action()` (lowercase with underscores)
- Types: `palma_type_t` (typedef structs)
- Constants: `PALMA_CONSTANT` (uppercase)
- Local variables: `snake_case`

### Formatting
- 4-space indentation (no tabs)
- Opening braces on same line
- Maximum line length: 100 characters
- Space after keywords (`if`, `for`, `while`)

### Example
```c
palma_matrix_t* palma_matrix_create(size_t rows, size_t cols) {
    if (rows == 0 || cols == 0) {
        return NULL;
    }
    
    palma_matrix_t *mat = malloc(sizeof(palma_matrix_t));
    if (mat == NULL) {
        return NULL;
    }
    
    mat->rows = rows;
    mat->cols = cols;
    mat->data = calloc(rows * cols, sizeof(palma_val_t));
    
    return mat;
}
```

### Documentation
- All public functions must have header comments
- Use `/** ... */` for documentation comments
- Document parameters, return values, and edge cases

## Testing

- Add tests for new features
- Ensure existing tests pass
- Test on ARM platforms if modifying NEON code
- Include edge cases (empty matrices, single elements, etc.)

## Areas for Contribution

We especially welcome contributions in:

- **New examples**: Control systems, Petri nets, etc.
- **Platform support**: RISC-V, other ARM variants
- **Language bindings**: Python, Rust
- **Performance**: GPU/CUDA port, multi-threading improvements
- **Documentation**: Tutorials, mathematical background
- **Testing**: More comprehensive test coverage

## Questions?

- Open an issue for general questions
- Email the maintainer for private matters

## License

By contributing, you agree that your contributions will be licensed under the MIT License.
