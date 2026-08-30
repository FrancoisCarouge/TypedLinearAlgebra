# Test Matrix {#test_matrix}

Which operations are exercised against which backends, from the `pass()` and `fail()` declarations under `test/`. Regenerated after every `ctest` run. Hover a cell for the shapes covered (`✗` = a compile-`fail` test).

@htmlonly
<style>
#tm{border-collapse:collapse;font:12px/1.3 system-ui,sans-serif}
#tm th,#tm td{border:1px solid rgba(127,127,127,.35);padding:2px 6px}
#tm thead th:not(:first-child){writing-mode:vertical-rl;transform:rotate(180deg);
  vertical-align:bottom;font-weight:600}
#tm tbody th{text-align:left;white-space:nowrap;font-weight:600}
#tm td{text-align:center}
#tm .y{color:#2e9a3e;font-weight:700}
#tm .n{color:rgba(127,127,127,.6)}
</style>
<p><span class="y">&#9679;</span> tested&nbsp;&nbsp;
<span class="n">&#8709;</span> not applicable&nbsp;&nbsp;
blank = untested</p>
<div style="overflow-x:auto"><table id="tm"><thead><tr><th>operation</th><th>au_eigen</th><th>au_std</th><th>chrono_eigen</th><th>chrono_std</th><th>eigen</th><th>eigexed</th><th>mp_units_eigen</th><th>mp_units_std</th><th>nested_typed_eigen</th><th>nholthaus_eigen</th><th>nholthaus_std</th></tr></thead><tbody>
<tr><th>addition</th><td class="y" title="1×1 1×2 ✗1×1 ✗1×2">&#9679;</td><td class="n" title="Au's lvalue-only Quantity::operator= is incompatible with the std::linalg-style add() (use operator+)">&#8709;</td><td class="y" title="1×2 ✗1×1">&#9679;</td><td></td><td class="y" title="1×2">&#9679;</td><td class="y" title="1×2">&#9679;</td><td class="y" title="1×1 1×2 ✗1×1 ✗1×2">&#9679;</td><td class="y" title="1×1 1×2">&#9679;</td><td class="y" title="1×2">&#9679;</td><td></td><td></td></tr>
<tr><th>assign</th><td></td><td></td><td></td><td></td><td class="y" title="any">&#9679;</td><td class="y" title="1×5 5×1 any">&#9679;</td><td></td><td></td><td class="y" title="1×5 5×1 any">&#9679;</td><td></td><td></td></tr>
<tr><th>at</th><td></td><td></td><td></td><td></td><td></td><td class="y" title="1×1 1×3 3×1 3×3">&#9679;</td><td></td><td></td><td class="y" title="1×1 1×3 3×1 3×3">&#9679;</td><td></td><td></td></tr>
<tr><th>common_with</th><td></td><td></td><td></td><td></td><td class="y" title="3×1">&#9679;</td><td class="y" title="3×1">&#9679;</td><td></td><td></td><td class="y" title="3×1">&#9679;</td><td></td><td></td></tr>
<tr><th>constructor</th><td></td><td></td><td class="y" title="1×1">&#9679;</td><td></td><td class="y" title="1×3 3×1 any">&#9679;</td><td class="y" title="1×1 1×3 3×1 any">&#9679;</td><td class="y" title="1×1">&#9679;</td><td class="y" title="1×1">&#9679;</td><td class="y" title="1×1 1×3 3×1 any">&#9679;</td><td></td><td></td></tr>
<tr><th>copy</th><td></td><td></td><td></td><td></td><td class="y" title="5×5">&#9679;</td><td class="y" title="5×5">&#9679;</td><td></td><td></td><td class="y" title="5×5">&#9679;</td><td></td><td></td></tr>
<tr><th>division</th><td class="y" title="1×1 2×1">&#9679;</td><td></td><td></td><td></td><td></td><td></td><td class="y" title="1×1 2×1">&#9679;</td><td></td><td></td><td></td><td></td></tr>
<tr><th>element</th><td class="y" title="1×1 1×2 2×1 2×2 ✗1×1 ✗1×2 ✗2×1 ✗2×2">&#9679;</td><td></td><td class="y" title="1×2">&#9679;</td><td></td><td></td><td></td><td class="y" title="1×1 1×2 2×1 2×2 ✗1×1 ✗1×2 ✗2×1 ✗2×2">&#9679;</td><td></td><td></td><td></td><td></td></tr>
<tr><th>equal_to</th><td class="y" title="1×1 1×2 2×3 ✗1×1 ✗1×2 ✗2×3">&#9679;</td><td class="y" title="1×1 1×2 2×1">&#9679;</td><td class="y" title="1×2">&#9679;</td><td class="y" title="1×1">&#9679;</td><td class="y" title="1×1 1×2 2×3">&#9679;</td><td class="y" title="1×1 1×2 2×3 ✗1×2/2×1 ✗any">&#9679;</td><td class="y" title="1×1 1×2 2×3 ✗1×1 ✗1×2 ✗2×3">&#9679;</td><td class="y" title="1×1 1×2 2×1">&#9679;</td><td class="y" title="1×1 1×2 2×3 ✗1×2/2×1 ✗any">&#9679;</td><td class="y" title="1×1 1×2 ✗1×1">&#9679;</td><td class="y" title="1×1">&#9679;</td></tr>
<tr><th>format</th><td></td><td></td><td></td><td></td><td class="y" title="1×1 1×3 3×1 3×3 any">&#9679;</td><td class="y" title="1×1 1×3 3×1 3×3 any">&#9679;</td><td></td><td></td><td class="y" title="1×1 1×3 3×1 3×3 any">&#9679;</td><td></td><td></td></tr>
<tr><th>magnitude</th><td class="y" title="1×2 3×1">&#9679;</td><td></td><td class="y" title="3×1">&#9679;</td><td></td><td></td><td class="y" title="1×2 2×1 ✗2×2">&#9679;</td><td class="y" title="1×2 3×1">&#9679;</td><td></td><td class="y" title="1×2 2×1">&#9679;</td><td></td><td></td></tr>
<tr><th>matrix_product</th><td></td><td class="y" title="1×1 1×2/2×1 2×1/1×2 2×2">&#9679;</td><td></td><td></td><td></td><td></td><td></td><td class="y" title="1×1 1×2/2×1 2×1/1×2 2×2">&#9679;</td><td></td><td></td><td></td></tr>
<tr><th>matrix_vector_product</th><td></td><td class="y" title="2×2 3×2 any ✗any">&#9679;</td><td></td><td></td><td></td><td></td><td></td><td class="y" title="2×2 3×2 any ✗any">&#9679;</td><td></td><td></td><td></td></tr>
<tr><th>minus</th><td class="y" title="1×1">&#9679;</td><td></td><td></td><td></td><td></td><td class="y" title="1×1 1×2 2×3 3×1">&#9679;</td><td class="y" title="1×1">&#9679;</td><td></td><td class="y" title="1×1 1×2 2×3 3×1">&#9679;</td><td></td><td></td></tr>
<tr><th>mp_units</th><td></td><td></td><td></td><td></td><td></td><td></td><td class="y" title="any">&#9679;</td><td></td><td></td><td></td><td></td></tr>
<tr><th>multiplication</th><td class="y" title="1×1 2×2 any">&#9679;</td><td class="y" title="1×1">&#9679;</td><td></td><td></td><td class="y" title="2×2/2×1 any">&#9679;</td><td class="y" title="1×1 1×2/2×1 2×2/2×1 any">&#9679;</td><td class="y" title="1×1 2×2 any">&#9679;</td><td class="y" title="1×1">&#9679;</td><td class="y" title="1×1 1×2/2×1 2×2/2×1 any">&#9679;</td><td></td><td></td></tr>
<tr><th>nested</th><td></td><td></td><td></td><td></td><td></td><td></td><td></td><td></td><td class="y" title="3×3">&#9679;</td><td></td><td></td></tr>
<tr><th>operator</th><td></td><td></td><td></td><td></td><td class="y" title="any">&#9679;</td><td class="y" title="1×1 1×3 2×2 3×1 3×3 any">&#9679;</td><td></td><td></td><td class="y" title="1×1 1×3 2×2 3×1 3×3 any">&#9679;</td><td></td><td></td></tr>
<tr><th>scale</th><td></td><td class="y" title="1×1 1×2 2×2 any">&#9679;</td><td></td><td class="y" title="any">&#9679;</td><td></td><td></td><td></td><td class="y" title="1×1 1×2 2×2 any">&#9679;</td><td></td><td></td><td class="y" title="1×2">&#9679;</td></tr>
<tr><th>structured_bindings</th><td></td><td></td><td></td><td></td><td class="y" title="1×1 1×3 3×1">&#9679;</td><td class="y" title="1×1 1×3 3×1">&#9679;</td><td></td><td></td><td></td><td></td><td></td></tr>
<tr><th>substraction</th><td class="y" title="1×1 1×2 ✗1×1 ✗1×2">&#9679;</td><td class="n" title="Au's lvalue-only Quantity::operator= is incompatible with the std::linalg-style substract() (use operator-)">&#8709;</td><td class="y" title="1×2">&#9679;</td><td></td><td class="y" title="1×2">&#9679;</td><td class="y" title="1×2">&#9679;</td><td class="y" title="1×1 1×2 ✗1×1 ✗1×2">&#9679;</td><td></td><td class="y" title="1×2">&#9679;</td><td></td><td></td></tr>
<tr><th>transposed</th><td class="y" title="2×3">&#9679;</td><td class="y" title="1×1 2×1">&#9679;</td><td></td><td></td><td></td><td class="y" title="2×3">&#9679;</td><td class="y" title="2×3">&#9679;</td><td class="y" title="1×1 2×1">&#9679;</td><td></td><td></td><td></td></tr>
<tr><th>underlying</th><td></td><td></td><td></td><td></td><td></td><td class="y" title="3×3">&#9679;</td><td></td><td></td><td class="y" title="3×3">&#9679;</td><td></td><td></td></tr>
</tbody></table></div>
@endhtmlonly

## Not applicable {#test_matrix_na}

- `addition` × `au_std` — Au's lvalue-only Quantity::operator= is incompatible with the std::linalg-style add() (use operator+)
- `substraction` × `au_std` — Au's lvalue-only Quantity::operator= is incompatible with the std::linalg-style substract() (use operator-)

<em>191 test files, 280 compiled cases, 23 operations × 11 backends.</em>
