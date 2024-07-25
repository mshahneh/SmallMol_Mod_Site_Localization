

def is_shifted(val1, val2, ppm = None, mz_tol = None):
    """
    Check if two values are shifted by a given ppm value or mz_tol value
    """
    diff = abs(val1 - val2)
    if ppm is None and mz_tol is None:
        raise ValueError("Either ppm or mz_tol must be provided")
    if mz_tol is not None and diff > mz_tol:
        return True
    if ppm is not None and diff > max(val1, val2) * ppm / 1e6:
        return True
    return False