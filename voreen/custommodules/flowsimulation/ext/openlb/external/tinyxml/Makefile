AR := ar
ARFLAGS := rc

TIXML_OBJDIR:=./build

TIXML_SRC := $(wildcard *.cpp)
TIXML_OBJ := $(TIXML_SRC:.cpp=.o)
TIXML_OBJ := $(foreach file,$(TIXML_OBJ), $(addprefix $(TIXML_OBJDIR)/,$(file)))

all: $(TIXML_OBJDIR) $(TIXML_OBJDIR)/libtinyxml.a

$(TIXML_OBJDIR):
	mkdir -p $(TIXML_OBJDIR)

$(TIXML_OBJDIR)/%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(TIXML_OBJDIR)/libtinyxml.a: $(TIXML_OBJ)
	ar rc $@ $(TIXML_OBJ)

clean:
	@rm $(TIXML_OBJDIR)/* &> /dev/null || true
