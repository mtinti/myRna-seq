# ritSeq Branch Preparation

This branch is reserved for extending the RNA-seq pipeline with new RIT-Seq functionality. The file documents the initial scope and notes so development can proceed consistently.

## Goals
- Integrate RIT-Seq specific processing steps into the existing workflow.
- Define configuration entries required for RIT-Seq inputs and outputs.
- Add tests to validate RIT-Seq workflows alongside existing RNA-seq tests.

## Next steps
1. Add rule stubs for RIT-Seq processing under `rules/` and link them in `Snakefile`.
2. Introduce configuration placeholders (e.g., sample metadata, references) in `config.yaml`.
3. Create or update environment files in `envs/` if new tools are needed.
4. Extend tests under `tests/` with RIT-Seq sample cases.
5. Update documentation to describe RIT-Seq usage and parameters.

This file serves as a starting checklist to guide further development on the `ritSeq` branch.
