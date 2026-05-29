# Poiseuille wallVP case-selection hotfix

Fixes MATLAB struct assignment in `build_cases` when selecting a subset of Poiseuille wallVP cases.

The selected case array is now initialized as `allCases([])` so it has the same fields as the available case definitions before scalar assignment.
