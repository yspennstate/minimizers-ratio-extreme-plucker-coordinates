# Minimizers of the ratio of extreme Plücker coordinates on the positive Grassmannian

The paper of the same title by Yitzchak Shmalo (`paper/main.tex`, `paper/main.pdf`) and the scripts for its numerical checks.

For a point V of the positive Grassmannian Gr>0(k,n), let E(V) be the ratio of its largest Plücker coordinate to its smallest. Ogranovich conjectured that E is minimized at the point V0 fixed by the cyclic shift, and proved this for k=2. The paper proves the conjecture for all k and n and determines all minimizers. They are the images of V0 under the positive column scalings in a compact convex polytope of dimension d − gcd(d, n/d), where d = gcd(k, n), so V0 is the unique minimizer exactly when d² divides n. For k=2 this settles the case 4 | n, which Ogranovich left open, and gives the exact range of his family of minimizers for n ≡ 2 (mod 4). The paper also shows that for (k,n) = (3,6) the version of the problem with absolute values of minors over all real matrices has a smaller minimum than the positive version.

## Running the checks

Requires Python 3 with numpy and scipy.

    sh run_checks.sh

This writes `outputs/<script>_out.txt` for each script and `outputs/versions.txt`. The outputs committed here were produced on macOS with the versions listed in `outputs/versions.txt`. The Nelder–Mead values in `ogr_check` and `abs_check` can differ slightly with other numpy and scipy versions.

| Script | What it checks (section of the paper) |
|---|---|
| `ogr_check.py` | E(V0) = P(B)/P(J) and Δ_I(V0) ∝ P(I) for (3,5), (3,6), (3,7), (3,8), (4,7), (4,8); Nelder–Mead searches for a smaller E (Sections 2 and 8) |
| `orbit_check.py` | the orbit inequality G_I/G_J ≥ P(I)/P(J) at random totally positive points (the main theorem, Section 3) |
| `col_check.py` | the dimension d − g of the set of minimizing column scalings, E under scalings in and outside the subspace, and the inequality G_I ≤ G_B (Sections 4 and 5) |
| `rank2_check.py` | the interval of minimizing alternating scalings for n ≡ 2 (mod 4) and its absence when 4 divides n (Section 6) |
| `lemma_check.py` | the three-term Plücker relation at random totally positive matrices, and that every T_r is a translate of B for 1 ≤ k < n ≤ 24 (Sections 4 and 5) |
| `ico.py` | the twenty 3×3 minors of the icosahedral 3×6 matrix (Section 7) |
| `abs_check.py` | Nelder–Mead searches for the absolute-value ratio over real 3×6 and 3×7 matrices (Section 7) |

## License

MIT.
