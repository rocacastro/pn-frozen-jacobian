# Rationale for the 15-to-31 paired repetition rule

The expensive multiprecision campaign is designed to reduce computational cost without
choosing sample size after seeing which method wins.

Every one of the 16 predeclared groups receives 15 paired timing repetitions. After all
initial groups are complete, the program applies one frozen mechanical rule:

- the paired percentile-bootstrap interval for the ratio of medians must exclude 1;
- the descriptive two-IQR separation rule must be met; and
- the direction of both diagnostics must agree with the observed ratio of medians.

Only a group satisfying all three conditions stops at 15. Every other group receives
exactly 16 additional paired repetitions, ending at 31.

The rule is intentionally conservative. It does not turn either diagnostic into a formal
universal superiority test. Its purpose is to avoid spending 31 expensive repetitions on
a group whose separation is already clear under both predeclared descriptive summaries,
while forcing extra sampling whenever either summary remains ambiguous.

The extension decisions are generated only after all 16 initial groups have completed and
are written to an immutable `extension_decisions.json` file. Manual overrides are not part
of the protocol. The methods, fields, tolerances, initial distance, precision and stopping
criterion are never changed in response to the observed winner.
