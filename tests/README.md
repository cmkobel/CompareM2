# tests/

`unit/` is the test suite — pure Python, no pixi and no tools needed:

```bash
python -m pytest tests/unit -q
```

The other subdirectories hold real genomes for end-to-end runs, shipped zipped
because the unpacked `.fna` files are deliberately gitignored. `pixi run unpack`
extracts `E._faecium/`.

`E._faecium/` is the standing cross-check: `116_2.fna` and
`116_2 duplicate.fna` are the same genome twice, under a filename with a space
in it. Any tool that treats them differently is wrong — the pair must come out
at 0.00000 mash distance, 100.00% ANI, 0 SNPs and identical CDS counts.
