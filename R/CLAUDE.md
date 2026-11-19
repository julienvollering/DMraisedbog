# R/ — code-level notes

Loaded when working with files under `R/`.

## Known Issues

**randomForestSRC + Positron crash with `do.trace`**
- **Problem**: R crashes (exit code -1073740791) when `randomForestSRC::imbalanced()` is
  called with `do.trace = TRUE` (or any numeric value) inside functions wrapped with
  `purrr::map_dfr()` in the Positron IDE.
- **Symptom**: "R exited unexpectedly" / "Error refreshing variables: The language runtime
  exited before RPC completed".
- **Solution**: always use `do.trace = FALSE` when calling `imbalanced()` inside functions
  or loops.
- **Root cause**: console trace output conflicts with Positron's variable-refresh
  mechanism in nested function contexts.
