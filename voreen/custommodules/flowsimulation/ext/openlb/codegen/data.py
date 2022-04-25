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
from sympy.codegen.ast import Assignment
from bindings import olb

import itertools

class Cell:
    def __init__(self, descriptor, name=None, generator=None):
        self.descriptor = descriptor
        self.data_fields = [ descriptor.fields_t.get[i] for i in range(descriptor.fields_t.size) ]
        if name and generator:
            self.realize(name, generator)

    def realize(self, name, generator):
        self.name = name
        self.symbol = IndexedBase(name, self.descriptor.d)
        self.expr = olb.CellD[olb.Expr,self.descriptor]()
        self.symbols = { name: self.symbol }
        self.aliases = [ ]
        self.substitutions = [ ]
        self.optional_assignments = [ ]
        self.optional_writebacks = [ ]
        for iPop in range(self.descriptor.q):
            self.expr[iPop] = olb.Expr(str(self.symbol[iPop]))
            alias = next(generator)
            self.symbols[alias.name] = alias
            self.aliases.append(alias)
            self.substitutions.append((alias, self.symbol[iPop]))
        for field in filter(lambda f: f != olb.descriptors.POPULATION, self.data_fields):
            for iD in range(olb.FieldD[olb.Expr,self.descriptor,field].d):
                alias = next(generator)
                self.symbols[alias.name] = alias
                self.expr.getFieldPointer[field]()[iD] = olb.Expr(alias.name)
                self.optional_assignments.append(Assignment(alias, Symbol(f"{ name }.template getFieldComponent<{ field.__cpp_name__ }>({ iD })")))
                self.optional_writebacks.append((field, iD, olb.Expr(alias.name).describe()))
        self.original = self.describe()

    def describe(self):
        assignments = [ ]
        for iPop in range(self.descriptor.q):
            assignments.append(f"Assignment({ self.aliases[iPop] }, { self.expr[iPop].describe().decode('utf-8') })")
        for field, iD, original in self.optional_writebacks:
            current = self.expr.getFieldComponent[field](iD).describe()
            if current != original:
                assignments.append(f"Assignment(Symbol(\"{ self.name }.template getFieldPointer<{ field.__cpp_name__ }>()[{ iD }]\"), { current })")
        return '[' + ','.join(assignments) + ']'

    def isChanged(self):
        return self.original != self.describe()

    def signature(self):
        return f"CELL& { self.name }"

class Scalar:
    def realize(self, name, generator):
        self.name = name
        self.alias = next(generator)
        self.symbol = Symbol(name)
        self.expr = olb.Expr(name)
        self.symbols = { name: self.symbol, self.alias.name: self.alias }
        self.substitutions = [ (self.alias, self.symbol) ]
        self.original = self.describe()

    def describe(self):
        return f"[ Assignment({ self.alias.name }, { self.expr.describe().decode('utf-8') }) ]"

    def isChanged(self):
        return self.original != self.describe()

    def signature(self):
        return f"{ self.name.upper() }& { self.name }"

class Vector:
    def __init__(self, size):
        self.size = size

    def realize(self, name, generator):
        self.name = name
        self.symbol = IndexedBase(name, self.size)
        self.expr = olb.Vector[olb.Expr,self.size]()
        self.symbols = { name: self.symbol }
        self.aliases = [ ]
        self.substitutions = [ ]
        for iD in range(self.size):
            self.expr[iD] = olb.Expr(str(self.symbol[iD]))
            alias = next(generator)
            self.symbols[alias.name] = alias
            self.aliases.append(alias)
            self.substitutions.append((alias, self.symbol[iD]))
        self.original = self.describe()

    def describe(self):
        return '[' + ','.join([ f"Assignment({ self.aliases[iD] }, { self.expr[iD].describe().decode('utf-8') })" for iD in range(self.size) ]) + ']'

    def isChanged(self):
        return self.original != self.describe()

    def signature(self):
        return f"{ self.name.upper() }& { self.name }"

class SquareMatrix:
    def __init__(self, size):
        self.size = size

    def realize(self, name):
        self.name = name
        self.symbol = IndexedBase(name, shape=(self.size, self.size))
        self.expr = olb.Vector[olb.Vector[olb.Expr,self.size],self.size]()
        self.symbols = { name: self.symbol }
        self.aliases = [ ]
        self.substitutions = [ ]
        for i in range(self.size):
            for j in range(self.size):
                self.expr[i][j] = olb.Expr(str(self.symbol[i,j]))
                alias = next(generator)
                self.symbols[alias.name] = alias
                self.aliases.append(alias)
                self.substitutions.append((alias, self.symbol[i,j]))
        self.original = self.describe()

    def describe(self):
        return '[' + ','.join([ f"Assignment({ self.aliases[self.size*i+j] }, { self.expr[i][j].describe().decode('utf-8') })" for i, j in itertools.product(range(self.size), repeat=2) ]) + ']'

    def isChanged(self):
        return self.original != self.describe()

    def signature(self):
        return f"{ self.name.upper() }& { self.name }"

class Parameters:
    def __init__(self, descriptor, fields, name=None, generator=None):
        self.descriptor = descriptor
        self.fields = fields
        if name and generator:
            self.realize(name, generator)

    def realize(self, name, generator):
        self.name = name
        self.expr = self.fields.decompose_into[olb.ParametersD[olb.Expr,self.descriptor].include_fields]()
        self.symbols = { }
        self.optional_assignments = [ ]
        self.substitutions = [ ]
        for i in range(self.fields.size):
            field = self.fields.get[i]
            size = olb.FieldD[olb.Expr,self.descriptor,field].d
            if size == 1:
                alias = next(generator)
                self.expr.set[field](olb.Expr(alias.name))
                self.symbols[alias.name] = alias
                self.optional_assignments.append(
                    Assignment(alias, Symbol(f"{ self.name }.template get<{ field.__cpp_name__ }>()")))
            else:
                data = olb.FieldD[olb.Expr,self.descriptor,field]()
                for iD in range(size):
                    alias = next(generator)
                    data[iD] = olb.Expr(alias.name)
                    self.symbols[alias.name] = alias
                    self.optional_assignments.append(
                        Assignment(alias, Symbol(f"{ self.name }.template get<{ field.__cpp_name__ }>({ iD })")))
                self.expr.set[field](data)

    def isChanged(self):
        return False

    def signature(self):
        return f"PARAMETERS& { self.name }"
