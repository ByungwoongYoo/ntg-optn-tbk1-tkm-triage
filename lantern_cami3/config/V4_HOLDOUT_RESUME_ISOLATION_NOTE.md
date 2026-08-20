# Isolated v4 holdout resume boundary

This branch is used only to execute the repaired sample-14/15 holdout workflow against the `cami3-lantern-20260819` base branch without triggering the repository-wide main-target workflow matrix. The scientific inputs, thresholds, model, and Toy pair are unchanged. The only repair is materializing immutable read artifacts as durable files rather than deleting the directories behind symbolic links before mapping.
