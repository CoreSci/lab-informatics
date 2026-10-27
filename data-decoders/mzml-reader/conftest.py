"""Shared pytest configuration for the mzML tools.

Every test carries a ``@pytest.mark.req(...)`` marker naming the requirements it
verifies (``traceability/requirements.tsv``). This hook copies those ids into
the test's JUnit XML ``<properties>`` so that ``traceability/build_matrix.py``
can build the requirements -> tests -> results matrix from a real test run.
"""


def pytest_runtest_setup(item):
    for marker in item.iter_markers(name="req"):
        for req_id in marker.args:
            item.user_properties.append(("requirement", req_id))
