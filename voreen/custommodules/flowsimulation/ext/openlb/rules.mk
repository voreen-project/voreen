#  This file is part of the OpenLB library
#
#  Copyright (C) 2017 Markus Mohrhard
#  E-mail contact: info@openlb.net
#  The most recent release of OpenLB can be downloaded at
#  <http://www.openlb.net/>
#
#  This program is free software; you can redistribute it and/or
#  modify it under the terms of the GNU General Public License
#  as published by the Free Software Foundation; either version 2
#  of the License, or (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public
#  License along with this program; if not, write to the Free
#  Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
#  Boston, MA  02110-1301, USA.

###########################################################################

include global.mk

###############################################################
## Handling of our samples

# documentation of the arguments:
# $(1) the source file name
# $(2) the object file name
define sample_object
$(2) : $(1)
	$(CXX) $(CXXFLAGS) $(INCLUDEDIR) -c $$< -o $$@

$(DEPENDDIR)/$(2:.o=.d) : $(1)
	@mkdir -p $$(dir $$@)
	@$(SHELL) -ec '$(CXX) -M $(CXXFLAGS) $(INCLUDEDIR) $(EXTRA_IDIR) $$< \
	    | sed -e "s!$(notdir $(2))!$(2)!1" > $$@;'

include $(wildcard $(DEPENDDIR)/$(2:.o=.d))

endef

# documentation of the arguments:
# $(1) the sample name
# $(2) the list of source files
# $(3) the directory name
define sample

$(foreach source,$(2),$(eval $(call sample_object,$(source),$(source:.cpp=.o))))

$(1) : $(2:.cpp=.o) $(LIBDIR)/lib$(LIB).a $(LIBDIR)/libz.a | $(DEPENDDIR)/$(2:.cpp=.d)
	$(CXX) $(LDFLAGS) $(2:.cpp=.o) -L./$(LIBDIR) $(LIBS) -o $$@

$(1)_clean :
	@rm -f $(2:.cpp=.o) &> /dev/null || true
	@rm -f $(1) &> /dev/null || true
	@rm -f $(DIR)*.d

$(1)_samples_clean :
	$(eval DIR=$(dir $(2)))
	@rm -f $(DIR)core $(DIR).tmpfile $(DIR)tmp/*.*
	@rm -f $(DIR)tmp/vtkData/*.* $(DIR)tmp/vtkData/data/*.* $(DIR)tmp/imageData/*.* $(DIR)tmp/gnuplotData/*.*
	@rm -f $(DIR)*.d

CLEANTARGETS += $(1)_clean

SAMPLESCLEAN += $(1)_clean \
								$(1)_samples_clean

SAMPLES += $(1)

endef
