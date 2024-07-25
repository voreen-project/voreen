# This file is part of the OpenLB library
#
# Copyright (C) 2021 Adrian Kummerlaender
# E-mail contact: info@openlb.net
# The most recent release of OpenLB can be downloaded at
# <http://www.openlb.net/>
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public
# License along with this program; if not, write to the Free
# Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
# Boston, MA  02110-1301, USA.

from sympy import *

from sympy.codegen.ast import CodeBlock, Assignment, Return
from sympy.printing.c import C11CodePrinter
from sympy.codegen.rewriting import ReplaceOptim
from sympy.simplify import cse_main
from sympy.utilities.iterables import numbered_symbols

from bindings import olb

import itertools

import data as symbolic

class SymbolGenerator:
    def __init__(self, prefix):
        self.prefix = prefix
        self.symbols = set()

    def __iter__(self):
        self.generator = numbered_symbols(cls=Symbol, prefix=self.prefix)
        return self

    def __next__(self):
        symbol = next(self.generator)
        self.symbols.add(symbol)
        return symbol

class CodeBlockPrinter(C11CodePrinter):
    def __init__(self, subexprs):
        super(CodeBlockPrinter, self).__init__()
        self.known_functions['Abs'] = 'util::fabs'
        self.subexprs = subexprs

    def _print_Indexed(self, expr):
        if len(expr.indices) == 1:
            return f"{expr.base.name}[{expr.indices[0]}]"
        elif len(expr.indices) == 2:
            return f"{expr.base.name}[{expr.indices[0]}][{expr.indices[1]}]"

    def _print_Float(self, flt):
        return "V{%.15g}" % flt.evalf()

    def _print_Pow(self, expr):
        if expr.exp == -0.5:
            return "V{1} / util::sqrt(%s)" % self.doprint(expr.base)
        elif expr.exp == -1:
            return "V{1} / (%s)" % self.doprint(expr.base)
        elif expr.exp == -2:
            x = self.doprint(expr.base)
            return "V{1} / ((%s)*(%s))" % (x, x)
        elif expr.exp == 0.5:
            return "util::sqrt(%s)" % self.doprint(expr.base)
        elif expr.exp == 2:
            x = self.doprint(expr.base)
            return "((%s)*(%s))" % (x, x)
        else:
            return "util::pow(%s, %s)" % (self.doprint(expr.base), self.doprint(expr.exp))

    def _print_Assignment(self, expr):
        if expr.lhs in self.subexprs:
            return f"auto {self.doprint(expr.lhs)} = {self.doprint(expr.rhs.evalf())};"
        else:
            return f"{self.doprint(expr.lhs)} = {self.doprint(expr.rhs.evalf())};"

def code(expr, symbols={ }):
    return CodeBlockPrinter(symbols).doprint(expr)

def cse(block, generator):
    expand_pos_square = ReplaceOptim(
        lambda e: e.is_Pow and e.exp.is_integer and e.exp == 2,
        lambda p: UnevaluatedExpr(Mul(p.base, p.base, evaluate = False)),
    )
    custom_opti = cse_main.basic_optimizations + [
        (expand_pos_square, expand_pos_square)
    ]
    return block.cse(symbols=generator, optimizations=custom_opti, order='none')

def import_expr(expr, symbols={ }):
    serialized = expr.describe()
    if isinstance(serialized, str):
        return eval(serialized, globals(), symbols)
    else:
        return eval(serialized.decode('utf-8'), globals(), symbols)

def collision_cse(collision, descriptor, momenta, equilibrium, simple_post_collision = True, concrete = None):
    if concrete is None:
        concrete = collision.type[descriptor,momenta,equilibrium]
    print(f"Generating { concrete.__cpp_name__ }")
    # Distinct symbols for all extracted subexpressions and helper variables
    generator = iter(SymbolGenerator('x'))
    # Setup symbolic cell and parameter for collision
    cell = symbolic.Cell(descriptor, 'cell', generator)
    parameters = symbolic.Parameters(descriptor, collision.parameters, 'parameters', generator)
    # Collect symbols used by cell and parameters as evaluation context
    symbols = cell.symbols | parameters.symbols
    # Apply collision operator
    result = concrete().apply(cell.expr, parameters.expr)
    # Import expression tree from C++ into Python context
    result_cell = import_expr(cell, symbols)
    result_rho  = import_expr(result.rho, symbols)
    result_uSqr = import_expr(result.uSqr, symbols)
    # Collect population and optional field / parameter assignments
    assignments = result_cell
    optional_assignments = parameters.optional_assignments + cell.optional_assignments
    # Substitions between internal placeholder variables and symbolic arrays
    # These are required due to current CSE limititations
    post_collision = [ ]
    if simple_post_collision:
        # Simple dependency between pre and post collision values
        substitutions = cell.substitutions + parameters.substitutions
    else:
        substitutions = parameters.substitutions
        # Use additional intermediary variables if post collision values
        # depend on different pre collision populations
        for subs in cell.substitutions:
            post_collision.append(Assignment(subs[1], subs[0]))
    # Assign cell statistics placeholders
    assignments.append(Assignment(Symbol("x_rho"), result_rho))
    assignments.append(Assignment(Symbol("x_uSqr"), result_uSqr))
    # Apply common subexpression elimination and back-substitution of placeholders
    block = cse(CodeBlock(*assignments), generator).subs(substitutions)
    # Detect which optional assignments are required to close the symbolic collision
    optional = [ ]
    for symbol in block.free_symbols:
        if symbol in { assgn.lhs for assgn in optional_assignments }:
            for assgn in optional_assignments:
                if assgn.lhs is symbol:
                    optional.append(assgn)
    # Generate C++ template for CSE-ified collision operator
    return '\n'.join([
        f"""template <typename CELL, typename PARAMETERS, typename V=typename CELL::value_t>""",
        f"""CellStatistic<V> apply(CELL& cell, PARAMETERS& parameters) any_platform""",
        "{",
        code(CodeBlock(*optional, *block.args[:-2], *post_collision), symbols=generator.symbols),
        "return { %s, %s };" % (code(block.args[-2].rhs), code(block.args[-1].rhs)),
        "}"
    ])

def cell_operator_cse(descriptor, struct, method, **args):
    print(f"Generating { struct.__cpp_name__ }::{ method }")
    # Distinct symbols for all extracted subexpressions and helper variables
    generator = iter(SymbolGenerator('x'))
    # Resolve operator callable
    operator = getattr(struct, method)
    # Generate symbols for arguments
    symbols = { }
    for name, arg in args.items():
        arg.realize(name, generator)
        symbols = symbols | arg.symbols
    # Apply operator to arguments
    result = operator(*[ arg.expr for arg in args.values() ])
    # Import returned expression tree from C++ into Python context
    if result != None:
        result = import_expr(result, symbols)
    # Import expression trees of modified arguments from C++ into Python context
    assignments = [ ]
    substitutions = [ ]
    for arg in filter(lambda arg: arg.isChanged(), args.values()):
        description = import_expr(arg, symbols)
        assignments = assignments + description
        substitutions = substitutions + arg.substitutions
    # Generate C++ result type depending on operator results
    result_type = "void"
    if result != None:
        result_type = "auto"
        assignments.append(Assignment(Symbol("result"), result))
    # Apply common subexpression elimination and back-substitutions
    block = cse(CodeBlock(*assignments), generator).subs(substitutions)
    # Generate C++ return statement if required
    optional = ""
    if result != None:
        optional = f"return { code(block.args[-1].rhs) };"
        block = CodeBlock(*block.args[:-1])
    # Generate C++ template for CSE-ified operator
    return '\n'.join([
        f"""template <{ ', '.join([ "typename " + arg.name.upper() for arg in args.values() ]) }, typename V=typename CELL::value_t>""",
        f"""static { result_type } { method }({ ', '.join([ arg.signature() for arg in args.values() ]) }) any_platform""",
        "{",
        code(block, symbols=generator.symbols),
        optional,
        "}"
    ])

def cell_method_cse(descriptor, struct, method, **args):
    print(f"Generating { struct.__cpp_name__ }::{ method }")
    # Distinct symbols for all extracted subexpressions and helper variables
    generator = iter(SymbolGenerator('x'))
    # Resolve operator callable
    operator = getattr(struct(), method)
    # Generate symbols for arguments
    symbols = { }
    optional_assignments = [ ]
    for name, arg in args.items():
        arg.realize(name, generator)
        symbols = symbols | arg.symbols
        if arg.optional_assignments:
            optional_assignments = optional_assignments + arg.optional_assignments
    # Apply operator to arguments
    result = operator(*[ arg.expr for arg in args.values() ])
    # Import returned expression tree from C++ into Python context
    if result != None:
        result = import_expr(result, symbols)
    # Import expression trees of modified arguments from C++ into Python context
    assignments = [ ]
    substitutions = [ ]
    for arg in filter(lambda arg: arg.isChanged(), args.values()):
        description = import_expr(arg, symbols)
        assignments = assignments + description
        substitutions = substitutions + arg.substitutions
    # Generate C++ result type depending on operator results
    result_type = "void"
    if result != None:
        result_type = "auto"
        assignments.append(Assignment(Symbol("result"), result))
    # Apply common subexpression elimination and back-substitutions
    block = cse(CodeBlock(*assignments), generator).subs(substitutions)
    # Dependency assignments
    dependencies = [ ]
    for symbol in block.free_symbols:
        if symbol in { assgn.lhs for assgn in optional_assignments }:
            for assgn in optional_assignments:
                if assgn.lhs is symbol:
                    dependencies.append(assgn)
    # Generate C++ return statement if required
    optional = ""
    if result != None:
        optional = f"return { code(block.args[-1].rhs) };"
        block = CodeBlock(*block.args[:-1])
    # Generate C++ template for CSE-ified operator
    return '\n'.join([
        f"""template <{ ', '.join([ "typename " + arg.name.upper() for arg in args.values() ]) }, typename V=typename CELL::value_t>""",
        f"""{ result_type } { method }({ ', '.join([ arg.signature() for arg in args.values() ]) }) any_platform""",
        "{",
        code(CodeBlock(*dependencies, *block.args), symbols=generator.symbols),
        optional,
        "}"
    ])
