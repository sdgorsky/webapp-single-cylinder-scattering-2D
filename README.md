# 2D Electromagnetic Scattering Simulator

Interactive visualization of electromagnetic plane wave scattering from an infinite dielectric cylinder. Computes and displays the scattered field using Mie theory with configurable material properties, wavelength, and polarization.

## For Developers

### Tech Stack

- **Frontend**: React + TypeScript + Vite
- **Computation**: Rust compiled to WebAssembly (wasm-pack)
- **Core Math**: Custom Bessel function implementations validated against SciPy

### Prerequisites

- Node.js (v18+)
- Rust toolchain (`rustup`)
- wasm-pack (`cargo install wasm-pack`)

### Commands

```bash
# Initial setup (install dependencies, configure git hooks)
make setup

# Run development server
make run

# Build for production
make build

# Run all tests
make test

# Run Bessel function validation tests
make test-bessel        # Full (~180K points)
make test-bessel-quick  # Quick (1000 samples)

# Lint and format
make lint
make format
```

### Project Structure

```
├── src/                  # React frontend
├── scattering-core/      # Rust WASM library
│   ├── src/
│   │   ├── bessel.rs     # Bessel/Hankel functions
│   │   ├── scattering.rs # Mie coefficients
│   │   └── field.rs      # Field computation
│   └── tests/            # Validation tests
└── Makefile              # Build commands
```
