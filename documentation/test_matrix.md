# Test Matrix {#test_matrix}

Which algorithms are exercised against which backends, from the `pass()` declarations under `test/`. A backend column is split into the strong type it integrates (header) and the linear algebra library it composes (footer, below the data, grouped across its types). Regenerated after every `ctest` run. Hover a cell for the shapes covered; a number means that many `pass()` tests cover the combination.

@htmlonly
<p><span class="y">&nbsp;&nbsp;&nbsp;</span> tested&nbsp;&nbsp;
blank = untested</p>
@endhtmlonly

@htmlinclude test_matrix.html
