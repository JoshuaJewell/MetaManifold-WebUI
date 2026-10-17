# Copy packages into renv/library rather than linking them into renv's per-user
# cache. That cache lives under ~/.cache (~/Library/Caches on macOS), which cache
# cleaners treat as disposable; emptying it left every linked package dangling and
# the project's R library unusable until a full restore. Copies still come from
# the cache, so restores stay fast. Convert an existing linked library in place
# with renv::isolate().
options(renv.config.cache.symlinks = FALSE)
source("renv/activate.R")
