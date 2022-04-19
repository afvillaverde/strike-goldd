import sympy as sp
def rationalize_all_numbers(expr):
    """
    converts all numbers in expr to sp.Rational-objects. This does not change Integers
    :param expr:
    :return:
    """
    numbers_atoms = list(expr.atoms(sp.Number))
    rationalized_number_tpls = [(n, sp.Rational(n)) for n in numbers_atoms]
    return expr.subs(rationalized_number_tpls)