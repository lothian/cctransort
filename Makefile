#
#@BEGIN LICENSE
#
# . by Psi4 Developer, a plugin to:
#
# PSI4: an ab initio quantum chemistry software package
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#
#@END LICENSE
#

#
# Plugin Makefile generated by Psi4.
#
# You shouldn't need to modify anything in this file.
#

# The name of your plugin. Taken from the directory name.
NAME = $(shell basename `pwd`)

# C++ source files for your plugin. By default we grab all *.cc files.
CXXSRC = $(notdir $(wildcard *.cc))

# Flags that were used to compile Psi4.
CXX = /usr/bin/g++
CXXDEFS = -DHAVE_PCMSOLVER -DHAVE_DKH -DHAVE_MM_MALLOC_H -DHAVE_SYSTEM_NATIVE_LAPACK -DHAVE_SYSTEM_NATIVE_BLAS -DHAS_CXX11_VARIADIC_TEMPLATES -DHAS_CXX11_STATIC_ASSERT -DHAS_CXX11_SIZEOF_MEMBER -DHAS_CXX11_RVALUE_REFERENCES -DHAS_CXX11_LIB_REGEX -DHAS_CXX11_NULLPTR -DHAS_CXX11_LONG_LONG -DHAS_CXX11_LAMBDA -DHAS_CXX11_INITIALIZER_LIST -DHAS_CXX11_DECLTYPE -DHAS_CXX11_CSTDINT_H -DHAS_CXX11_CONSTEXPR -DHAS_CXX11_AUTO_RET_TYPE -DHAS_CXX11_AUTO -DHAS_CXX11_FUNC -DHAS_CXX11 -DVAR_MFDS -DSYS_DARWIN -DUSE_FCMANGLE_H
CXXFLAGS = -DRESTRICT=__restrict__ -fPIC -std=c++11 -O0 -g -DDEBUG -Wall -Wextra -Winit-self -Woverloaded-virtual -Wuninitialized -Wmissing-declarations -Wwrite-strings -Weffc++ -Wdocumentation -Wno-unknown-pragmas
LDFLAGS = -Wl,-search_paths_first -Wl,-headerpad_max_install_names
INCLUDES = -I/Users/crawdad/src/psi4/objdir/interfaces/include -I/Users/crawdad/src/psi4/objdir/src/lib -I/Users/crawdad/src/psi4/src/lib -I/Users/crawdad/src/psi4/include -I/Users/crawdad/src/psi4/objdir/include -I/usr/local/include -I/usr/local/Cellar/python/2.7.10_2/Frameworks/Python.framework/Versions/2.7/include/python2.7 -I/usr/include

# Used to determine linking flags.
UNAME = $(shell uname)

DEPENDINCLUDE = $(notdir $(wildcard *.h*))

PSITARGET = $(NAME).so

# Start the compilation rules
default:: $(PSITARGET)

# Add the flags needed for shared library creation
ifeq ($(UNAME), Linux)
    LDFLAGS += -shared
endif
ifeq ($(UNAME), Darwin)
    LDFLAGS += -shared -undefined dynamic_lookup
    CXXFLAGS += -fno-common
endif

# The object files
BINOBJ = $(CXXSRC:%.cc=%.o)

%.o: %.cc
	$(CXX) $(CXXDEFS) $(CXXFLAGS) $(INCLUDES) -c $<

$(PSITARGET): $(BINOBJ)
	$(CXX) $(LDFLAGS) -o $@ $^ $(CXXDEFS)

# Erase all compiled intermediate files
clean:
	rm -f $(BINOBJ) $(PSITARGET) *.d *.pyc *.test output.dat psi.timer.dat

