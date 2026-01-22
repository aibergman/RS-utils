# Plan for Supporting Multiple Hamiltonian Formats

## Current Architecture

The band structure module currently supports **RS-LMTO native format** only:
- **LatticeReader** (`rsutil/bands/readers.py`): Parses `lattice.nml`
- **HamiltonianReader** (`rsutil/bands/readers.py`): Parses `rs2txt.dat`

These are instantiated directly in:
- `rsutil/cli.py` (lines 67, 70 in `run_bands_subcommand`, etc.)
- All band calculation classes expect these reader objects

## Proposed Architecture

### 1. **Create Reader Interface (Base Class)**
   **Location:** `rsutil/bands/readers.py`
   
   - Define abstract base class `BaseLatticeReader` with methods:
     - `parse_lattice_file()` (abstract)
     - `get_structure_info()` (shared)
   
   - Define abstract base class `BaseHamiltonianReader` with methods:
     - `parse_hamiltonian_file()` (abstract)
     - `fourier_transform_hamiltonian(k_point)` (shared via `hamiltonian_data`)
     - Properties: `orbital_size`, `hamiltonian_data`

   - Keep existing `LatticeReader`, `HamiltonianReader` as concrete RS-LMTO implementations

### 2. **Create Format-Specific Reader Modules**
   **Location:** Create new subpackage `rsutil/bands/readers/`
   
   Structure:
   ```
   rsutil/bands/readers/
   ├── __init__.py              # Export all readers
   ├── base.py                  # Base classes
   ├── rs_lmto.py              # RS-LMTO native format (refactored from readers.py)
   ├── paoflow.py              # PAOFLOW format reader
   └── kite.py                 # KITE format reader
   ```

### 3. **PAOFLOW Support**
   **File:** `rsutil/bands/readers/paoflow.py`
   
   - Class: `PAOFLOWLatticeReader(BaseLatticeReader)`
     - Parse PAOFLOW structure output format
     - Convert to internal representation
   
   - Class: `PAOFLOWHamiltonianReader(BaseHamiltonianReader)`
     - Parse PAOFLOW tight-binding Hamiltonian
     - Handle basis information (number of orbitals, orbital types)
     - Populate `hamiltonian_data` structure

### 4. **KITE Support**
   **File:** `rsutil/bands/readers/kite.py`
   
   - Class: `KITELatticeReader(BaseLatticeReader)`
     - Parse KITE lattice definition format
   
   - Class: `KITEHamiltonianReader(BaseHamiltonianReader)`
     - Parse KITE Hamiltonian format
     - Handle orbital definitions
     - Populate `hamiltonian_data` structure

### 5. **Factory Pattern for Reader Selection**
   **Location:** `rsutil/bands/readers/__init__.py`
   
   - Function: `load_lattice(filename, format='rs-lmto')`
     - Auto-detect format if not specified
     - Return appropriate reader instance
   
   - Function: `load_hamiltonian(filename, format='rs-lmto')`
     - Auto-detect format if not specified
     - Return appropriate reader instance

### 6. **Update CLI**
   **File:** `rsutil/cli.py`
   
   Changes to `_make_common_parser()`:
   ```python
   p.add_argument('--lattice-format', choices=['rs-lmto', 'paoflow', 'kite'], 
                  default='rs-lmto', help='Lattice file format')
   p.add_argument('--hamiltonian-format', choices=['rs-lmto', 'paoflow', 'kite'],
                  default='rs-lmto', help='Hamiltonian file format')
   ```
   
   Changes to `run_bands_subcommand()` and `run_dos_subcommand()`:
   ```python
   from .bands.readers import load_lattice, load_hamiltonian
   
   lattice = load_lattice(args.lattice, format=args.lattice_format)
   hamiltonian = load_hamiltonian(args.hamiltonian, format=args.hamiltonian_format)
   ```

### 7. **Update Main Module Exports**
   **File:** `rsutil/bands/__init__.py`
   
   - Export base classes: `BaseLatticeReader`, `BaseHamiltonianReader`
   - Export format readers: `PAOFLOWLatticeReader`, `PAOFLOWHamiltonianReader`, etc.
   - Export factory: `load_lattice`, `load_hamiltonian`
   - Keep backward compatibility with `LatticeReader`, `HamiltonianReader`

## Implementation Steps

1. **Refactor readers into modular structure** (minimal breaking changes)
2. **Extract base classes**
3. **Implement PAOFLOW reader** (research PAOFLOW format first)
4. **Implement KITE reader** (research KITE format first)
5. **Add factory functions and format detection**
6. **Update CLI with format selection**
7. **Update documentation** (CLI section, examples)
8. **Add tests** for each format (using example files)

## Key Design Principles

- **Backward compatible**: Default format is 'rs-lmto'; existing code works unchanged
- **Extensible**: New formats can be added without modifying core code
- **Consistent interface**: All readers implement same base methods
- **Auto-detection**: Can attempt to detect format from file signatures
- **Error handling**: Clear error messages if format detection fails

## Files to Modify

**Create (new):**
- `rsutil/bands/readers/__init__.py`
- `rsutil/bands/readers/base.py`
- `rsutil/bands/readers/rs_lmto.py`
- `rsutil/bands/readers/paoflow.py`
- `rsutil/bands/readers/kite.py`

**Modify:**
- `rsutil/bands/readers.py` → Move to `rsutil/bands/readers/rs_lmto.py`
- `rsutil/cli.py` → Add format selection arguments
- `rsutil/bands/__init__.py` → Update exports
- Documentation files (CLI examples, module docs)

## Next Steps

1. Research PAOFLOW Hamiltonian file format
2. Research KITE Hamiltonian file format  
3. Begin with refactoring into modular structure
4. Implement readers incrementally with tests
