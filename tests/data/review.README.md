# Derived streaming fixtures

The three `large-*.hic` fixtures are independently assembled by
`straw/C++/tests/review_fixture.py`. Run that script with this directory as its
argument to regenerate them. They contain 70,000 source cells, exceeding the
65,536-entry aggregation buffer:

- `large-derived-count.hic`: exact count cells (2i, 2i+2)=1, derived by factor 2.
- `large-derived-score.hic`: the same geometry with float32 scores.
- `large-score-fanin.hic`: 70,000 trans score cells, all derived into one pixel.

The tests check complete output, overlapping cis windows including reflected
edge cells, and the bounded fallback for a single large score aggregation.
