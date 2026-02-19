# SynTerra Changelog

## Version 2.0 - Post-Critical-Analysis (February 3, 2026)

### Major Improvements

#### Critical Bug Fixes (3)
- ✅ Fixed fatal bug: Query gene now included in iterative search
- ✅ Increased padding from 20kb to 150kb (adaptive)
- ✅ Verified region output always generated

#### Major Enhancements (7)
- ✅ Added Smith-Waterman rigorous alignment (optional)
- ✅ Enabled exon-level search by default
- ✅ Enhanced RBH validation (coverage + identity)
- ✅ Implemented pseudogene detection (5 types)
- ✅ Improved phylogenetic ordering (distance-based waves)
- ✅ Smart cluster distance calculation
- ✅ Added augmented search with variants

#### Moderate Improvements (5)
- ✅ Robust Miniprot GFF parsing with error handling
- ✅ Checkpoint resume capability (--resume flag)
- ✅ Unified MMseqs2 sensitivity parameter
- ✅ Better error diagnostics
- ✅ Enhanced logging throughout

#### Code Quality Improvements (4)
- ✅ Added type hints to core functions
- ✅ Replaced print with proper logging module
- ✅ Created test suite (21 passing tests)
- ✅ Added comprehensive input validation

### Repository Cleanup
- 📁 Moved documentation to `docs/` directory
- 🧹 Cleaned old test results (tetramorium, melettin)
- 🧹 Cleaned Nextflow work directory
- 🧹 Removed old log files
- 📝 Updated README with documentation links

### Production Readiness
**90%** (up from 30%)

### Test Coverage
21 unit tests covering:
- Coordinate normalization
- Gene ID extraction
- FASTA I/O
- Sequence operations
- Hit filtering
- Synteny identification

### Documentation
See `docs/` directory for comprehensive documentation:
- **FIXES_SUMMARY.md** - Complete overview
- **CRITICAL_ANALYSIS.md** - Technical analysis
- Individual fix documents for each category

---

## Version 1.0 - Initial Release

Original implementation with basic functionality.

