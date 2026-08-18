import lcurve
import numpy as np
import pytest

def compare_cpp_rust_solution(
    modfile: str,
    cpp_solution_file: str
) -> None:

    # Load the mod file and the solution
    binary_model = lcurve.BinaryModel.from_file(modfile)
    cpp_solution = np.loadtxt(cpp_solution_file)

    # Prepare inputs and cpp flux measurements
    time = np.ascontiguousarray(cpp_solution[:,0])
    flux = np.ascontiguousarray(cpp_solution[:,1])
    t_exp = np.ones_like(time) * 1/len(time)
    n_div = np.ones_like(time)

    # Calculate _Rust_ light curve
    lc = binary_model.compute_light_curve(
        time,
        t_exp,
        n_div
    )

    # Compare the Rust and C++ solution
    diff = lc.total - flux
    # Use approximate values of the difference to avoid floating point problems
    assert pytest.approx(diff) == np.zeros_like(diff)
    
    return None

def test_PTF1J0823():
    compare_cpp_rust_solution(
        "tests/modfiles/PTF1J0823+0819.mod",
        "tests/cpp_solutions/PTF1J0823+0819.dat"
    )
    return None

def test_FBS_2303_tess():
    compare_cpp_rust_solution(
        "tests/modfiles/FBS_2303+344_TESS.mod",
        "tests/cpp_solutions/FBS_2303+344_TESS.dat"
    )
    return None

def test_FBS_2303_zr():
    compare_cpp_rust_solution(
        "tests/modfiles/FBS_2303+344_zr.mod",
        "tests/cpp_solutions/FBS_2303+344_zr.dat"
    )
    return None
def test_PG0941():
    compare_cpp_rust_solution(
        "tests/modfiles/PG0941+280.mod",
        "tests/cpp_solutions/PG0941+280.dat"
    )
    return None
