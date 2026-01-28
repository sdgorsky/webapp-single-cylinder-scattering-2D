.PHONY: setup run build clean check-deps install-deps lint format test test-bessel test-bessel-quick help

# Colors for output
RED := \033[0;31m
GREEN := \033[0;32m
YELLOW := \033[0;33m
NC := \033[0m # No Color

help:
	@echo "Available commands:"
	@echo "  make setup    - Check and install dependencies"
	@echo "  make run      - Build and run the development server"
	@echo "  make build    - Build WASM and frontend"
	@echo "  make lint     - Run linters for Rust and TypeScript"
	@echo "  make format   - Format Rust and TypeScript code"
	@echo "  make test     - Run tests"
	@echo "  make clean    - Clean build artifacts"

check-deps:
	@echo "Checking dependencies..."
	@MISSING=""; \
	if ! command -v node &> /dev/null; then \
		echo "$(RED)✗ node is not installed$(NC)"; \
		MISSING="$$MISSING node"; \
	else \
		echo "$(GREEN)✓ node found ($$(node --version))$(NC)"; \
	fi; \
	if ! command -v npm &> /dev/null; then \
		echo "$(RED)✗ npm is not installed$(NC)"; \
		MISSING="$$MISSING npm"; \
	else \
		echo "$(GREEN)✓ npm found ($$(npm --version))$(NC)"; \
	fi; \
	if ! command -v cargo &> /dev/null; then \
		echo "$(RED)✗ cargo is not installed$(NC)"; \
		MISSING="$$MISSING cargo"; \
	else \
		echo "$(GREEN)✓ cargo found ($$(cargo --version))$(NC)"; \
	fi; \
	if ! command -v rustfmt &> /dev/null; then \
		echo "$(RED)✗ rustfmt is not installed$(NC)"; \
		MISSING="$$MISSING rustfmt"; \
	else \
		echo "$(GREEN)✓ rustfmt found$(NC)"; \
	fi; \
	if ! command -v wasm-pack &> /dev/null; then \
		echo "$(RED)✗ wasm-pack is not installed$(NC)"; \
		MISSING="$$MISSING wasm-pack"; \
	else \
		echo "$(GREEN)✓ wasm-pack found ($$(wasm-pack --version))$(NC)"; \
	fi; \
	if [ -n "$$MISSING" ]; then \
		echo ""; \
		echo "$(YELLOW)Missing dependencies:$$MISSING$(NC)"; \
		exit 1; \
	else \
		echo ""; \
		echo "$(GREEN)All dependencies are installed!$(NC)"; \
	fi

install-deps:
	@echo "Installing missing dependencies..."
	@if ! command -v rustup &> /dev/null; then \
		echo "$(YELLOW)Rust is not installed. Install it from https://rustup.rs/$(NC)"; \
		read -p "Would you like to install Rust now? [y/N] " confirm; \
		if [ "$$confirm" = "y" ] || [ "$$confirm" = "Y" ]; then \
			curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh; \
		fi; \
	fi
	@if ! command -v wasm-pack &> /dev/null; then \
		echo "$(YELLOW)wasm-pack is not installed.$(NC)"; \
		read -p "Would you like to install wasm-pack now? [y/N] " confirm; \
		if [ "$$confirm" = "y" ] || [ "$$confirm" = "Y" ]; then \
			cargo install wasm-pack; \
		fi; \
	fi
	@if ! command -v rustfmt &> /dev/null; then \
		echo "$(YELLOW)rustfmt is not installed.$(NC)"; \
		read -p "Would you like to install rustfmt now? [y/N] " confirm; \
		if [ "$$confirm" = "y" ] || [ "$$confirm" = "Y" ]; then \
			rustup component add rustfmt; \
		fi; \
	fi
	@if ! rustup component list | grep -q "clippy.*installed"; then \
		echo "$(YELLOW)clippy is not installed.$(NC)"; \
		read -p "Would you like to install clippy now? [y/N] " confirm; \
		if [ "$$confirm" = "y" ] || [ "$$confirm" = "Y" ]; then \
			rustup component add clippy; \
		fi; \
	fi
	@if [ ! -d "node_modules" ]; then \
		echo "$(YELLOW)Node modules not installed.$(NC)"; \
		read -p "Would you like to install npm dependencies now? [y/N] " confirm; \
		if [ "$$confirm" = "y" ] || [ "$$confirm" = "Y" ]; then \
			npm install; \
		fi; \
	fi

setup: check-deps install-deps
	@echo ""
	@git config core.hooksPath .githooks
	@echo "$(GREEN)Git hooks configured.$(NC)"
	@echo ""
	@echo "$(GREEN)Setup complete!$(NC)"
	@echo "Run 'make run' to start the development server."

# Build WASM module
build-wasm:
	@echo "Building WASM module..."
	cd scattering-core && wasm-pack build --target web

# Install npm dependencies if needed
npm-install:
	@if [ ! -d "node_modules" ]; then \
		echo "Installing npm dependencies..."; \
		npm install; \
	fi

# Build frontend
build-frontend: npm-install
	@echo "Building frontend..."
	npm run build

# Full build
build: build-wasm build-frontend
	@echo "$(GREEN)Build complete!$(NC)"

# Run development server
run: build-wasm npm-install
	@echo "Starting development server..."
	npm run dev

# Lint Rust code
lint-rust:
	@echo "Linting Rust code..."
	cd scattering-core && cargo clippy -- -D warnings

# Lint TypeScript code
lint-ts: npm-install
	@echo "Linting TypeScript code..."
	npm run lint

# Lint all
lint: lint-rust lint-ts
	@echo "$(GREEN)Linting complete!$(NC)"

# Format Rust code
format-rust:
	@echo "Formatting Rust code..."
	cd scattering-core && cargo fmt

# Format TypeScript code
format-ts: npm-install
	@echo "Formatting TypeScript code..."
	npm run format

# Format all
format: format-rust format-ts
	@echo "$(GREEN)Formatting complete!$(NC)"

# Check formatting without modifying files
format-check:
	@echo "Checking Rust formatting..."
	cd scattering-core && cargo fmt --check
	@echo "Checking TypeScript formatting..."
	npm run format:check

# Run Rust tests
test-rust:
	@echo "Running Rust tests..."
	cd scattering-core && cargo test

# Run all tests
test: test-rust
	@echo "$(GREEN)All tests passed!$(NC)"

# Run Bessel function validation tests against SciPy reference data
test-bessel:
	@echo "Running Bessel validation tests (full ~180K points)..."
	cd scattering-core && cargo test validate_ -- --nocapture

# Run quick Bessel validation with 1000 sample points
test-bessel-quick:
	@echo "Running Bessel validation tests (quick - 1000 samples)..."
	cd scattering-core && BESSEL_SAMPLE_SIZE=1000 cargo test validate_ -- --nocapture

# Clean build artifacts
clean:
	@echo "Cleaning build artifacts..."
	cd scattering-core && cargo clean
	rm -rf scattering-core/pkg
	rm -rf node_modules
	rm -rf dist
	@echo "$(GREEN)Clean complete!$(NC)"
