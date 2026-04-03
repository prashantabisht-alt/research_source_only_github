# Source-Only Repo Organization Plan

This repo should stay simpler than the full local workspace. The idea is:

- keep the canonical folders readable
- avoid generated data and figure exports
- keep archive material available, but clearly secondary

## Best Working Rule

Use these folders first:

- `chapter1`
- `ABP in 2D`
- `jerky particle`
- `RTP`
- `TCRW`

Treat these as references:

- `Brownian Motion`
- `Brownian motion separate`
- `Brownian Motion 2`
- `new_work`

## Recommended Next Improvements

1. Add a small top-level project description on GitHub.
2. Keep future DP2 code changes inside canonical folders first.
3. If a new clean DP2 branch emerges, create a dedicated folder for it instead of extending `Brownian Motion 2`.
4. Eventually separate archive folders more explicitly if the repo grows.

## Why This Repo Exists

The full local project is better for historical completeness. This source-only repo is better for:

- GitHub sharing
- AI-assisted coding tools
- code review
- future DP2 development
