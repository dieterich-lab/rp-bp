
import pytest


def check_hash(h1, h2, hash_func):
    return hash_func(h1) == hash_func(h2)
    
#@pytest.mark.parametrize('files, ref_files, checker_function', 
                         #[(files, ref_files, checker_function) for files, ref_files, checker_function in get_files])
#def test_output(files, ref_files, checker_function):
    #assert check_hash(files, ref_files, checker_function)
    
def test_output(get_files):
    for files, ref_files, checker_function in get_files:
        assert check_hash(files, ref_files, checker_function)
        
        
        
