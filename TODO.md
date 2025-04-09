
# TODO.md

## 🔧 Deployment Fixes
- [x] Fix `Tabulator` init error by assigning `columns` after widget creation
- [x] Remove unsupported `virtualization` keyword from `Tabulator`
- [ ] Ensure correct `panel serve` command with `--port` and `--address` on Render

## ⚠️ Warnings to Address
- [ ] Clean up non-critical verbose warnings
- [ ] Optimize DataFrame construction to avoid fragmentation warnings (consider using `copy()` after column insertions)